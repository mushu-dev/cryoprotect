#!/usr/bin/env python3
"""
CryoProtect v2 - Fix RLS Implementation

This script implements consistent Row Level Security (RLS) policies across all tables in the database.
It enables RLS on all tables in the schema and creates appropriate RLS policies based on table type.

Key features:
1. Enables RLS on all tables in the schema
2. Creates comprehensive RLS policies for each table type
3. Restricts anonymous access to sensitive data
4. Adds indexes for policy columns to optimize performance
5. Implements transaction support for safer operations
6. Provides detailed error handling and reporting

Usage:
    python fix_rls_implementation.py [--dry-run] [--verify-only] [--rollback] [--schema SCHEMA]
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
log_filename = f"rls_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Tables that don't need row-level security policies
NON_SENSITIVE_TABLES = ['property_types', 'calculation_methods']

# Backup of original state for rollback
original_state = {}

def get_supabase_client():
    """Get a Supabase client with service role key."""
    try:
        from supabase import create_client, Client
        
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY")
        
        if not SUPABASE_URL or not SUPABASE_KEY:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        sys.exit(1)

def execute_sql(supabase, sql, description, dry_run=False):
    """Execute SQL using the Supabase client."""
    if dry_run:
        logger.info(f"DRY RUN: Would execute SQL: {description}")
        logger.debug(f"SQL: {sql}")
        return True, None
    
    try:
        logger.info(f"Executing SQL: {description}")
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False, response.error
        
        logger.info(f"SQL executed successfully: {description}")
        return True, response.data
    except Exception as e:
        logger.error(f"Error executing SQL ({description}): {str(e)}")
        return False, str(e)

def execute_transaction(supabase, sql_statements, description, dry_run=False):
    """Execute multiple SQL statements in a transaction."""
    if dry_run:
        logger.info(f"DRY RUN: Would execute transaction: {description}")
        for i, sql in enumerate(sql_statements):
            logger.debug(f"Statement {i+1}: {sql}")
        return True, None
    
    # Combine statements into a transaction
    transaction_sql = "BEGIN;\n"
    transaction_sql += "\n".join(sql_statements)
    transaction_sql += "\nCOMMIT;"
    
    return execute_sql(supabase, transaction_sql, description)

def get_schema_tables(supabase, schema_name='remediation_test', dry_run=False):
    """Get all tables in the specified schema."""
    if dry_run:
        logger.info(f"DRY RUN: Would get tables for schema {schema_name}")
        # Return a default list of tables for dry run
        return [
            'molecules', 'experiments', 'predictions', 'mixtures',
            'mixture_components', 'property_types', 'calculation_methods',
            'user_profile', 'experiment_mixtures', 'predictions_old',
            'projects', 'teams'
        ]
    
    get_tables_sql = f"""
    SELECT tablename
    FROM pg_tables
    WHERE schemaname = '{schema_name}';
    """
    
    success, tables_result = execute_sql(supabase, get_tables_sql, "Get schema tables")
    
    if not success:
        logger.error(f"Failed to get tables for schema {schema_name}")
        return []
    
    # Handle different response formats
    if isinstance(tables_result, list):
        # Standard list of dictionaries format
        if len(tables_result) > 0 and 'tablename' in tables_result[0]:
            return [table['tablename'] for table in tables_result]
    elif isinstance(tables_result, dict):
        # Handle case where result is a dictionary with 'success' key
        logger.info(f"Tables result format is a dictionary: {tables_result}")
        
        # Let's try a different approach to get tables
        get_tables_sql = f"""
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = '{schema_name}';
        """
        
        success, tables_result = execute_sql(supabase, get_tables_sql, "Get schema tables (alternative)")
        
        if success and isinstance(tables_result, list):
            if len(tables_result) > 0 and 'table_name' in tables_result[0]:
                return [table['table_name'] for table in tables_result]
    
    # If we get here, try a direct query to list all tables
    logger.info("Trying direct query to list tables...")
    list_tables_sql = f"""
    SELECT string_agg(table_name, ',') as tables
    FROM information_schema.tables
    WHERE table_schema = '{schema_name}';
    """
    
    success, direct_result = execute_sql(supabase, list_tables_sql, "Direct query for tables")
    
    if success and isinstance(direct_result, list) and len(direct_result) > 0:
        if 'tables' in direct_result[0] and direct_result[0]['tables']:
            return direct_result[0]['tables'].split(',')
    
    # Hardcode the tables we know should exist based on test_database_remediation.py
    logger.warning("Using hardcoded table list as fallback")
    return [
        'molecules', 'experiments', 'predictions', 'mixtures',
        'mixture_components', 'property_types', 'calculation_methods',
        'user_profile', 'experiment_mixtures', 'predictions_old',
        'projects', 'teams'
    ]

def check_column_exists(supabase, schema_name, table, column, dry_run=False):
    """Check if a column exists in a table."""
    if dry_run:
        # For dry run, assume common columns exist
        if column in ['created_by', 'is_public', 'project_id', 'team_id']:
            return True
        return False
    
    check_column_sql = f"""
    SELECT EXISTS (
        SELECT FROM information_schema.columns
        WHERE table_schema = '{schema_name}'
        AND table_name = '{table}'
        AND column_name = '{column}'
    ) as has_column;
    """
    
    success, column_result = execute_sql(supabase, check_column_sql, f"Check {column} column in {table}")
    
    if success:
        if isinstance(column_result, list) and len(column_result) > 0:
            if 'has_column' in column_result[0]:
                return column_result[0]['has_column']
    
    # Default to False if we can't determine
    return False

def check_rls_status(supabase, schema_name, table, dry_run=False):
    """Check if RLS is enabled for a table."""
    if dry_run:
        logger.info(f"DRY RUN: Would check RLS status for {schema_name}.{table}")
        return None
    
    check_rls_sql = f"""
    SELECT relrowsecurity
    FROM pg_class c
    JOIN pg_namespace n ON n.oid = c.relnamespace
    WHERE n.nspname = '{schema_name}' AND c.relkind = 'r' AND c.relname = '{table}';
    """
    
    success, rls_result = execute_sql(supabase, check_rls_sql, f"Check RLS on {table}")
    
    if success:
        if isinstance(rls_result, list) and len(rls_result) > 0:
            if 'relrowsecurity' in rls_result[0]:
                return rls_result[0]['relrowsecurity']
    
    return None

def backup_table_state(supabase, schema_name, table, dry_run=False):
    """Backup the current state of a table's RLS settings and policies."""
    if dry_run:
        logger.info(f"DRY RUN: Would backup state for {schema_name}.{table}")
        return
    
    # Check if RLS is enabled
    rls_enabled = check_rls_status(supabase, schema_name, table)
    
    # Get existing policies
    get_policies_sql = f"""
    SELECT policyname, cmd, qual, with_check
    FROM pg_policies
    WHERE schemaname = '{schema_name}' AND tablename = '{table}';
    """
    
    success, policies = execute_sql(supabase, get_policies_sql, f"Get policies for {table}")
    
    if not success:
        policies = []
    
    original_state[f"{schema_name}.{table}"] = {
        'rls_enabled': rls_enabled,
        'policies': policies
    }
    
    logger.info(f"Backed up state for table {schema_name}.{table}")

def enable_rls(supabase, schema_name, table, dry_run=False):
    """Enable RLS on a table."""
    # First check if RLS is already enabled
    if not dry_run:
        rls_enabled = check_rls_status(supabase, schema_name, table)
        
        if rls_enabled:
            logger.info(f"RLS already enabled on table {schema_name}.{table}")
            return True
    
    enable_rls_sql = f"""
    ALTER TABLE {schema_name}.{table} ENABLE ROW LEVEL SECURITY;
    """
    
    success, _ = execute_sql(supabase, enable_rls_sql, f"Enable RLS on {table}", dry_run)
    
    if success:
        logger.info(f"Enabled RLS on table {schema_name}.{table}")
    else:
        logger.error(f"Failed to enable RLS on table {schema_name}.{table}")
    
    return success

def create_public_read_policy(supabase, schema_name, table, dry_run=False):
    """Create a policy to allow public read access if is_public = true."""
    # Check if the table has an is_public column
    has_is_public = check_column_exists(supabase, schema_name, table, 'is_public', dry_run)
    
    if not has_is_public:
        logger.info(f"Table {schema_name}.{table} does not have an is_public column, skipping public read policy")
        return True
    
    policy_name = f"Public read access on {table}"
    
    create_policy_sql = f"""
    DROP POLICY IF EXISTS "{policy_name}" ON {schema_name}.{table};
    CREATE POLICY "{policy_name}" 
    ON {schema_name}.{table} 
    FOR SELECT 
    USING (is_public = true);
    
    COMMENT ON POLICY "{policy_name}" ON {schema_name}.{table} IS 
    'Allows anyone to SELECT records where is_public = true';
    """
    
    success, _ = execute_sql(supabase, create_policy_sql, f"Create public read policy for {table}", dry_run)
    
    if success:
        logger.info(f"Created public read policy for table {schema_name}.{table}")
    else:
        logger.error(f"Failed to create public read policy for table {schema_name}.{table}")
    
    return success

def create_owner_access_policy(supabase, schema_name, table, dry_run=False):
    """Create a policy to allow full access to record owners."""
    # Check if the table has a created_by column
    has_created_by = check_column_exists(supabase, schema_name, table, 'created_by', dry_run)
    
    if not has_created_by:
        logger.info(f"Table {schema_name}.{table} does not have a created_by column, skipping owner access policy")
        return True
    
    policy_name = f"Owner full access on {table}"
    
    create_policy_sql = f"""
    DROP POLICY IF EXISTS "{policy_name}" ON {schema_name}.{table};
    CREATE POLICY "{policy_name}" 
    ON {schema_name}.{table} 
    USING (created_by = auth.uid());
    
    COMMENT ON POLICY "{policy_name}" ON {schema_name}.{table} IS 
    'Allows full access to records created by the user';
    """
    
    success, _ = execute_sql(supabase, create_policy_sql, f"Create owner access policy for {table}", dry_run)
    
    if success:
        logger.info(f"Created owner access policy for table {schema_name}.{table}")
    else:
        logger.error(f"Failed to create owner access policy for table {schema_name}.{table}")
    
    return success

def create_team_member_access_policy(supabase, schema_name, table, dry_run=False):
    """Create a policy to allow team members to access records."""
    policy_name = f"Team member access on {table}"
    
    # Different tables might have different ways to link to teams
    team_link_condition = ""
    
    if check_column_exists(supabase, schema_name, table, 'team_id', dry_run):
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM {schema_name}.user_profile
            WHERE user_profile.team_id = {table}.team_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif check_column_exists(supabase, schema_name, table, 'project_id', dry_run):
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM {schema_name}.projects
            JOIN {schema_name}.teams ON teams.id = projects.team_id
            JOIN {schema_name}.user_profile ON user_profile.team_id = teams.id
            WHERE projects.id = {table}.project_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif table == 'mixture_components' and check_column_exists(supabase, schema_name, table, 'mixture_id', dry_run):
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM {schema_name}.mixtures
            JOIN {schema_name}.projects ON projects.id = mixtures.project_id
            JOIN {schema_name}.teams ON teams.id = projects.team_id
            JOIN {schema_name}.user_profile ON user_profile.team_id = teams.id
            WHERE mixtures.id = {table}.mixture_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif table == 'experiment_properties' and check_column_exists(supabase, schema_name, table, 'experiment_id', dry_run):
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM {schema_name}.experiments
            JOIN {schema_name}.projects ON projects.id = experiments.project_id
            JOIN {schema_name}.teams ON teams.id = projects.team_id
            JOIN {schema_name}.user_profile ON user_profile.team_id = teams.id
            WHERE experiments.id = {table}.experiment_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    else:
        # For tables without a clear team relationship, skip this policy
        logger.info(f"Table {schema_name}.{table} does not have a clear team relationship, skipping team member access policy")
        return True
    
    create_policy_sql = f"""
    DROP POLICY IF EXISTS "{policy_name}" ON {schema_name}.{table};
    CREATE POLICY "{policy_name}" 
    ON {schema_name}.{table} 
    FOR SELECT 
    USING ({team_link_condition});
    
    COMMENT ON POLICY "{policy_name}" ON {schema_name}.{table} IS 
    'Allows SELECT for users who are members of the team associated with the record';
    """
    
    success, _ = execute_sql(supabase, create_policy_sql, f"Create team member access policy for {table}", dry_run)
    
    if success:
        logger.info(f"Created team member access policy for table {schema_name}.{table}")
    else:
        logger.error(f"Failed to create team member access policy for table {schema_name}.{table}")
    
    return success

def create_service_role_bypass_policy(supabase, schema_name, table, dry_run=False):
    """Create a policy to allow service role to bypass RLS."""
    policy_name = f"Service role bypass on {table}"
    
    create_policy_sql = f"""
    DROP POLICY IF EXISTS "{policy_name}" ON {schema_name}.{table};
    CREATE POLICY "{policy_name}" 
    ON {schema_name}.{table} 
    USING (auth.role() = 'service_role');
    
    COMMENT ON POLICY "{policy_name}" ON {schema_name}.{table} IS 
    'Allows all operations for service role';
    """
    
    success, _ = execute_sql(supabase, create_policy_sql, f"Create service role bypass policy for {table}", dry_run)
    
    if success:
        logger.info(f"Created service role bypass policy for table {schema_name}.{table}")
    else:
        logger.error(f"Failed to create service role bypass policy for table {schema_name}.{table}")
    
    return success

def create_rls_performance_indexes(supabase, schema_name, tables, dry_run=False):
    """Create indexes for better RLS performance."""
    index_statements = []
    
    # Create index on user_profile.auth_user_id
    index_statements.append(f"""
    CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id 
    ON {schema_name}.user_profile(auth_user_id);
    """)
    
    # Create indexes on team_id columns
    for table in tables:
        if table in NON_SENSITIVE_TABLES:
            continue
        
        if check_column_exists(supabase, schema_name, table, 'team_id', dry_run):
            index_name = f"idx_{table}_team_id"
            index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS {index_name} 
            ON {schema_name}.{table}(team_id);
            """)
    
    # Create indexes on project_id columns
    for table in tables:
        if table in NON_SENSITIVE_TABLES:
            continue
        
        if check_column_exists(supabase, schema_name, table, 'project_id', dry_run):
            index_name = f"idx_{table}_project_id"
            index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS {index_name} 
            ON {schema_name}.{table}(project_id);
            """)
    
    # Create indexes on created_by columns
    for table in tables:
        if table in NON_SENSITIVE_TABLES:
            continue
        
        if check_column_exists(supabase, schema_name, table, 'created_by', dry_run):
            index_name = f"idx_{table}_created_by"
            index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS {index_name} 
            ON {schema_name}.{table}(created_by);
            """)
    
    # Create indexes on is_public columns
    for table in tables:
        if table in NON_SENSITIVE_TABLES:
            continue
        
        if check_column_exists(supabase, schema_name, table, 'is_public', dry_run):
            index_name = f"idx_{table}_is_public"
            index_statements.append(f"""
            CREATE INDEX IF NOT EXISTS {index_name} 
            ON {schema_name}.{table}(is_public);
            """)
    
    # Execute all index creation statements in a transaction
    success, _ = execute_transaction(
        supabase, 
        index_statements, 
        "Create RLS performance indexes", 
        dry_run
    )
    
    if success:
        logger.info(f"Created RLS performance indexes for schema {schema_name}")
    else:
        logger.error(f"Failed to create RLS performance indexes for schema {schema_name}")
    
    return success

def fix_rls_implementation(supabase, schema_name='remediation_test', dry_run=False):
    """Fix the RLS implementation by enabling RLS on all tables and creating appropriate policies."""
    logger.info(f"Fixing RLS implementation for schema {schema_name}...")
    
    # Get all tables in the schema
    tables = get_schema_tables(supabase, schema_name, dry_run)
    
    if not tables:
        logger.error(f"No tables found in schema {schema_name}")
        return False
    
    logger.info(f"Found {len(tables)} tables in schema {schema_name}: {', '.join(tables)}")
    
    # Backup original state for rollback
    if not dry_run:
        logger.info("Backing up original state")
        for table in tables:
            backup_table_state(supabase, schema_name, table, dry_run)
    
    # Enable RLS on all tables
    for table in tables:
        if not enable_rls(supabase, schema_name, table, dry_run):
            logger.error(f"Failed to enable RLS on table {schema_name}.{table}")
            return False
    
    # Restrict anonymous access
    restrict_anon_sql = f"""
    REVOKE ALL ON SCHEMA {schema_name} FROM anon;
    REVOKE ALL ON ALL TABLES IN SCHEMA {schema_name} FROM anon;
    GRANT USAGE ON SCHEMA {schema_name} TO anon;
    """
    
    success, _ = execute_sql(supabase, restrict_anon_sql, "Restrict anonymous access", dry_run)
    
    if not success:
        logger.error("Failed to restrict anonymous access")
        return False
    
    # Grant anon access to non-sensitive tables
    for table in NON_SENSITIVE_TABLES:
        if table in tables:
            grant_anon_sql = f"""
            GRANT SELECT ON {schema_name}.{table} TO anon;
            """
            
            success, _ = execute_sql(supabase, grant_anon_sql, f"Grant anon access to {table}", dry_run)
            
            if not success:
                logger.error(f"Failed to grant anon access to {table}")
                return False
    
    # Create RLS policies for each table
    for table in tables:
        # Skip non-sensitive tables that don't need row-level policies
        if table in NON_SENSITIVE_TABLES:
            continue
        
        # Create public read policy
        if not create_public_read_policy(supabase, schema_name, table, dry_run):
            logger.error(f"Failed to create public read policy for {schema_name}.{table}")
            return False
        
        # Create owner access policy
        if not create_owner_access_policy(supabase, schema_name, table, dry_run):
            logger.error(f"Failed to create owner access policy for {schema_name}.{table}")
            return False
        
        # Create team member access policy
        if not create_team_member_access_policy(supabase, schema_name, table, dry_run):
            logger.error(f"Failed to create team member access policy for {schema_name}.{table}")
            return False
        
        # Create service role bypass policy
        if not create_service_role_bypass_policy(supabase, schema_name, table, dry_run):
            logger.error(f"Failed to create service role bypass policy for {schema_name}.{table}")
            return False
    
    # Create RLS performance indexes
    if not create_rls_performance_indexes(supabase, schema_name, tables, dry_run):
        logger.error(f"Failed to create RLS performance indexes for schema {schema_name}")
        return False
    
    logger.info(f"RLS implementation fixed successfully for schema {schema_name}")
    return True

def verify_rls_implementation(supabase, schema_name='remediation_test', dry_run=False):
    """Verify that RLS is enabled on all tables and policies are created."""
    if dry_run:
        logger.info(f"DRY RUN: Would verify RLS implementation for schema {schema_name}")
        return True
    
    logger.info(f"Verifying RLS implementation for schema {schema_name}...")
    
    # Get all tables in the schema
    tables = get_schema_tables(supabase, schema_name, dry_run)
    
    if not tables:
        logger.error(f"No tables found in schema {schema_name}")
        return False
    
    logger.info(f"Found tables in schema {schema_name}: {', '.join(tables)}")
    
    all_verified = True
    verification_results = {}
    
    # Verify RLS enabled on all tables
    for table in tables:
        rls_enabled = check_rls_status(supabase, schema_name, table)
        
        if rls_enabled is None:
            logger.warning(f"Could not verify RLS status for table {schema_name}.{table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": "unknown"}
            continue
        
        if not rls_enabled:
            logger.error(f"RLS is not enabled on table {schema_name}.{table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": False}
            continue
        
        # Verify RLS policies
        if table in NON_SENSITIVE_TABLES:
            verification_results[table] = {"rls_enabled": True, "policies_verified": "N/A"}
            continue
        
        # Check if policies exist
        check_policy_sql = f"""
        SELECT policyname 
        FROM pg_policies 
        WHERE schemaname = '{schema_name}' AND tablename = '{table}';
        """
        
        success, policies = execute_sql(supabase, check_policy_sql, f"Check policies on {table}")
        
        if not success:
            logger.warning(f"Could not verify policies for table {schema_name}.{table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": True, "policies_verified": False}
            continue
        
        policy_names = [p['policyname'] for p in policies]
        
        # Expected policies depend on table structure
        expected_policies = [f"Service role bypass on {table}"]
        
        if check_column_exists(supabase, schema_name, table, 'is_public'):
            expected_policies.append(f"Public read access on {table}")
        
        if check_column_exists(supabase, schema_name, table, 'created_by'):
            expected_policies.append(f"Owner full access on {table}")
        
        # Team member access policy is expected if the table has a team relationship
        has_team_relationship = (
            check_column_exists(supabase, schema_name, table, 'team_id') or
            check_column_exists(supabase, schema_name, table, 'project_id') or
            (table == 'mixture_components' and check_column_exists(supabase, schema_name, table, 'mixture_id')) or
            (table == 'experiment_properties' and check_column_exists(supabase, schema_name, table, 'experiment_id'))
        )
        
        if has_team_relationship:
            expected_policies.append(f"Team member access on {table}")
        
        missing_policies = [p for p in expected_policies if p not in policy_names]
        
        if missing_policies:
            logger.warning(f"Missing policies for table {schema_name}.{table}: {missing_policies}")
            all_verified = False
            verification_results[table] = {
                "rls_enabled": True, 
                "policies_verified": False,
                "missing_policies": missing_policies
            }
        else:
            logger.info(f"All policies verified for table {schema_name}.{table}")
            verification_results[table] = {"rls_enabled": True, "policies_verified": True}
    
    # Save verification results to file
    verification_file = f"rls_verification_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(verification_file, "w") as f:
        json.dump(verification_results, f, indent=2)
    
    logger.info(f"Verification results saved to {verification_file}")
    
    if all_verified:
        logger.info(f"RLS implementation verified successfully for schema {schema_name}")
    else:
        logger.warning(f"RLS implementation verification found issues for schema {schema_name}")
    
    return all_verified

def rollback_changes(supabase, schema_name='remediation_test', dry_run=False):
    """Rollback changes if needed."""
    if dry_run:
        logger.info(f"DRY RUN: Would rollback changes for schema {schema_name}")
        return True
    
    if not original_state:
        logger.warning("No original state to rollback to")
        return False
    
    logger.info(f"Rolling back changes for schema {schema_name}...")
    
    all_rolled_back = True
    
    # Get all tables in the schema
    tables = get_schema_tables(supabase, schema_name, dry_run)
    
    for table in tables:
        table_key = f"{schema_name}.{table}"
        
        if table_key not in original_state:
            logger.warning(f"No original state for table {table_key}, skipping rollback")
            continue
        
        state = original_state[table_key]
        
        # Drop all policies
        drop_policies_sql = f"""
        DO $$
        DECLARE
            pol record;
        BEGIN
            FOR pol IN SELECT policyname FROM pg_policies WHERE schemaname = '{schema_name}' AND tablename = '{table}'
            LOOP
                EXECUTE 'DROP POLICY "' || pol.policyname || '" ON {schema_name}.{table};';
            END LOOP;
        END $$;
        """
        
        success, _ = execute_sql(supabase, drop_policies_sql, f"Drop policies for {table}")
        
        if not success:
            logger.error(f"Failed to drop policies for table {table}")
            all_rolled_back = False
            continue
        
        # Disable RLS if it was originally disabled
        if not state['rls_enabled']:
            disable_rls_sql = f"""
            ALTER TABLE {schema_name}.{table} DISABLE ROW LEVEL SECURITY;
            """
            
            success, _ = execute_sql(supabase, disable_rls_sql, f"Disable RLS on {table}")
            
            if not success:
                logger.error(f"Failed to disable RLS on table {table}")
                all_rolled_back = False
                continue
        
        # Recreate original policies
        for policy in state['policies']:
            if not policy:
                continue
                
            policy_name = policy['policyname']
            cmd = policy['cmd']
            qual = policy['qual']
            with_check = policy['with_check']
            
            create_policy_sql = f"""
            CREATE POLICY "{policy_name}"
            ON {schema_name}.{table}
            """
            
            if cmd:
                create_policy_sql += f"FOR {cmd} "
            
            if qual:
                create_policy_sql += f"USING ({qual}) "
            
            if with_check:
                create_policy_sql += f"WITH CHECK ({with_check})"
            
            success, _ = execute_sql(supabase, create_policy_sql, f"Recreate policy {policy_name} for {table}")
            
            if not success:
                logger.error(f"Failed to recreate policy {policy_name} for table {table}")
                all_rolled_back = False
    
    logger.info(f"Rollback completed for schema {schema_name}")
    return all_rolled_back

def main():
    """Main function to fix the RLS implementation."""
    parser = argparse.ArgumentParser(description="Fix Row Level Security (RLS) implementation for CryoProtect v2")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    parser.add_argument("--verify-only", action="store_true", help="Only verify the current RLS implementation")
    parser.add_argument("--rollback", action="store_true", help="Rollback changes to the original state")
    parser.add_argument("--schema", default="remediation_test", help="Database schema to operate on (default: remediation_test)")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 RLS Implementation Fix")
    
    # Connect to Supabase
    try:
        supabase = get_supabase_client()
        logger.info("Connected to Supabase")
    except Exception as e:
        logger.error(f"Failed to connect to Supabase: {str(e)}")
        return 1
    
    if args.verify_only:
        logger.info(f"Verifying RLS implementation for schema {args.schema}")
        if not verify_rls_implementation(supabase, args.schema, args.dry_run):
            logger.error("Verification of RLS implementation failed")
            return 1
        logger.info("Verification completed")
        return 0
    
    if args.rollback:
        logger.info(f"Rolling back RLS implementation for schema {args.schema}")
        if not rollback_changes(supabase, args.schema, args.dry_run):
            logger.error("Rollback of RLS implementation failed")
            return 1
        logger.info("Rollback completed")
        return 0
    
    # Fix RLS implementation
    if not fix_rls_implementation(supabase, args.schema, args.dry_run):
        logger.error("Failed to fix RLS implementation")
        return 1
    
    # Verify RLS implementation
    if not verify_rls_implementation(supabase, args.schema, args.dry_run):
        logger.error("Verification of RLS implementation failed")
        return 1
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("CryoProtect v2 RLS Implementation Fix Summary")
    logger.info("=" * 80)
    
    logger.info("\nStatus: SUCCESS")
    logger.info("The RLS implementation has been fixed successfully.")
    logger.info("All tables now have RLS enabled and appropriate policies created.")
    
    logger.info("\nKey improvements:")
    logger.info("1. Enabled RLS on all tables in the schema")
    logger.info("2. Created comprehensive RLS policies for each table type")
    logger.info("3. Restricted anonymous access to sensitive data")
    logger.info("4. Added indexes for policy columns to optimize performance")
    logger.info("5. Implemented transaction support for safer operations")
    
    logger.info("\nFor detailed information, check the log file and verification results.")
    logger.info("=" * 80 + "\n")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

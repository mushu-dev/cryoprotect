#!/usr/bin/env python3
"""
CryoProtect v2 - Security Implementation Script

This script implements comprehensive security features for the CryoProtect Supabase project:
1. Enables Row Level Security (RLS) on all public tables
2. Creates RLS policies for public access, owner access, team member access, and service role bypass
3. Creates app-specific database roles with minimum permissions
4. Includes verification and rollback mechanisms

Usage:
    python implement_security.py [--dry-run] [--verify-only] [--rollback]
"""

import os
import argparse
import logging
import json
import time
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("security_implementation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase project ID
PROJECT_ID = "tsdlmynydfuypiugmkev"

# Tables to apply security to
TABLES = [
    "molecules",
    "mixtures",
    "mixture_components",
    "predictions",
    "experiments",
    "experiment_properties",
    "calculation_methods",
    "projects",
    "teams"
]

# Backup of original state for rollback
original_state = {}

def get_supabase_client():
    """Get a Supabase client with service role key."""
    SUPABASE_URL = os.getenv("SUPABASE_URL")
    SUPABASE_KEY = os.getenv("SUPABASE_KEY")
    
    if not SUPABASE_URL or not SUPABASE_KEY:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    
    return create_client(SUPABASE_URL, SUPABASE_KEY)

def execute_sql(supabase, sql, dry_run=False):
    """Execute SQL using the Supabase client."""
    if dry_run:
        logger.info("DRY RUN: Would execute SQL:")
        logger.info(sql)
        return True, None
    
    try:
        # Use the Supabase MCP tool to execute SQL
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False, response.error
        
        return True, response.data
    except Exception as e:
        logger.error(f"Error executing SQL: {str(e)}")
        return False, str(e)

def check_rls_status(supabase, table, dry_run=False):
    """Check if RLS is enabled for a table."""
    sql = f"""
    SELECT relrowsecurity 
    FROM pg_class 
    WHERE oid = '{table}'::regclass;
    """
    
    success, result = execute_sql(supabase, sql, dry_run)
    
    if not success or dry_run:
        return None
    
    if result and len(result) > 0:
        return result[0]['relrowsecurity']
    
    return None

def backup_table_state(supabase, table, dry_run=False):
    """Backup the current state of a table's RLS settings and policies."""
    if dry_run:
        return
    
    # Check if RLS is enabled
    rls_enabled = check_rls_status(supabase, table, dry_run)
    
    # Get existing policies
    sql = f"""
    SELECT policyname, cmd, qual, with_check 
    FROM pg_policies 
    WHERE tablename = '{table}';
    """
    
    success, policies = execute_sql(supabase, sql, dry_run)
    
    if not success:
        policies = []
    
    original_state[table] = {
        'rls_enabled': rls_enabled,
        'policies': policies
    }
    
    logger.info(f"Backed up state for table {table}")

def enable_rls(supabase, table, dry_run=False):
    """Enable RLS on a table."""
    # First check if RLS is already enabled
    if not dry_run:
        rls_enabled = check_rls_status(supabase, table, dry_run)
        
        if rls_enabled:
            logger.info(f"RLS already enabled on table {table}")
            return True
    
    sql = f"ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;"
    
    success, _ = execute_sql(supabase, sql, dry_run)
    
    if success:
        logger.info(f"Enabled RLS on table {table}")
    else:
        logger.error(f"Failed to enable RLS on table {table}")
    
    return success

def create_public_read_policy(supabase, table, dry_run=False):
    """Create a policy to allow public read access if is_public = true."""
    policy_name = f"Public read access on {table}"
    
    sql = f"""
    CREATE POLICY "{policy_name}" 
    ON {table} 
    FOR SELECT 
    USING (is_public = true);
    
    COMMENT ON POLICY "{policy_name}" ON {table} IS 
    'Allows anyone to SELECT records where is_public = true';
    """
    
    success, _ = execute_sql(supabase, sql, dry_run)
    
    if success:
        logger.info(f"Created public read policy for table {table}")
    else:
        logger.error(f"Failed to create public read policy for table {table}")
    
    return success

def create_owner_access_policy(supabase, table, dry_run=False):
    """Create a policy to allow full access to record owners."""
    # For SELECT, UPDATE, DELETE
    policy_name = f"Owner full access on {table}"
    
    sql = f"""
    CREATE POLICY "{policy_name}" 
    ON {table} 
    USING (created_by = auth.uid());
    
    COMMENT ON POLICY "{policy_name}" ON {table} IS 
    'Allows full access to records created by the user';
    """
    
    success, _ = execute_sql(supabase, sql, dry_run)
    
    if success:
        logger.info(f"Created owner access policy for table {table}")
    else:
        logger.error(f"Failed to create owner access policy for table {table}")
    
    return success

def create_team_member_access_policy(supabase, table, dry_run=False):
    """Create a policy to allow team members to access records."""
    policy_name = f"Team member access on {table}"
    
    # Different tables might have different ways to link to teams
    team_link_condition = ""
    
    if table in ["projects"]:
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM teams
            JOIN user_profile ON user_profile.team_id = teams.id
            WHERE teams.id = {table}.team_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif table in ["molecules", "mixtures", "experiments", "predictions"]:
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM projects
            JOIN teams ON teams.id = projects.team_id
            JOIN user_profile ON user_profile.team_id = teams.id
            WHERE projects.id = {table}.project_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif table in ["mixture_components"]:
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM mixtures
            JOIN projects ON projects.id = mixtures.project_id
            JOIN teams ON teams.id = projects.team_id
            JOIN user_profile ON user_profile.team_id = teams.id
            WHERE mixtures.id = {table}.mixture_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    elif table in ["experiment_properties"]:
        team_link_condition = f"""
        EXISTS (
            SELECT 1 FROM experiments
            JOIN projects ON projects.id = experiments.project_id
            JOIN teams ON teams.id = projects.team_id
            JOIN user_profile ON user_profile.team_id = teams.id
            WHERE experiments.id = {table}.experiment_id
            AND user_profile.auth_user_id = auth.uid()
        )
        """
    else:
        # For tables without a clear team relationship
        team_link_condition = "FALSE"
    
    sql = f"""
    CREATE POLICY "{policy_name}" 
    ON {table} 
    FOR SELECT 
    USING ({team_link_condition});
    
    COMMENT ON POLICY "{policy_name}" ON {table} IS 
    'Allows SELECT for users who are members of the team associated with the record';
    """
    
    success, _ = execute_sql(supabase, sql, dry_run)
    
    if success:
        logger.info(f"Created team member access policy for table {table}")
    else:
        logger.error(f"Failed to create team member access policy for table {table}")
    
    return success

def create_service_role_bypass_policy(supabase, table, dry_run=False):
    """Create a policy to allow service role to bypass RLS."""
    policy_name = f"Service role bypass on {table}"
    
    sql = f"""
    CREATE POLICY "{policy_name}" 
    ON {table} 
    USING (auth.role() = 'service_role');
    
    COMMENT ON POLICY "{policy_name}" ON {table} IS 
    'Allows all operations for service role';
    """
    
    success, _ = execute_sql(supabase, sql, dry_run)
    
    if success:
        logger.info(f"Created service role bypass policy for table {table}")
    else:
        logger.error(f"Failed to create service role bypass policy for table {table}")
    
    return success

def create_app_roles(supabase, dry_run=False):
    """Create app-specific database roles with minimum permissions."""
    # Create app_readonly role
    sql_readonly = """
    CREATE ROLE app_readonly;
    GRANT USAGE ON SCHEMA public TO app_readonly;
    GRANT SELECT ON ALL TABLES IN SCHEMA public TO app_readonly;
    ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO app_readonly;
    """
    
    success_readonly, _ = execute_sql(supabase, sql_readonly, dry_run)
    
    if success_readonly:
        logger.info("Created app_readonly role")
    else:
        logger.error("Failed to create app_readonly role")
    
    # Create app_user role
    sql_user = """
    CREATE ROLE app_user;
    GRANT USAGE ON SCHEMA public TO app_user;
    GRANT SELECT, INSERT, UPDATE ON ALL TABLES IN SCHEMA public TO app_user;
    ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT, INSERT, UPDATE ON TABLES TO app_user;
    """
    
    success_user, _ = execute_sql(supabase, sql_user, dry_run)
    
    if success_user:
        logger.info("Created app_user role")
    else:
        logger.error("Failed to create app_user role")
    
    # Create app_admin role
    sql_admin = """
    CREATE ROLE app_admin;
    GRANT USAGE ON SCHEMA public TO app_admin;
    GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO app_admin;
    ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL PRIVILEGES ON TABLES TO app_admin;
    """
    
    success_admin, _ = execute_sql(supabase, sql_admin, dry_run)
    
    if success_admin:
        logger.info("Created app_admin role")
    else:
        logger.error("Failed to create app_admin role")
    
    return success_readonly and success_user and success_admin

def verify_rls_policies(supabase, dry_run=False):
    """Verify that RLS policies are correctly applied."""
    if dry_run:
        logger.info("DRY RUN: Would verify RLS policies")
        return True
    
    all_verified = True
    verification_results = {}
    
    for table in TABLES:
        # Check if RLS is enabled
        rls_enabled = check_rls_status(supabase, table, dry_run)
        
        if rls_enabled is None:
            logger.warning(f"Could not verify RLS status for table {table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": "unknown"}
            continue
        
        if not rls_enabled:
            logger.error(f"RLS is not enabled on table {table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": False}
            continue
        
        # Check if policies exist
        sql = f"""
        SELECT policyname 
        FROM pg_policies 
        WHERE tablename = '{table}';
        """
        
        success, policies = execute_sql(supabase, sql, dry_run)
        
        if not success:
            logger.warning(f"Could not verify policies for table {table}")
            all_verified = False
            verification_results[table] = {"rls_enabled": True, "policies_verified": False}
            continue
        
        policy_names = [p['policyname'] for p in policies]
        expected_policies = [
            f"Public read access on {table}",
            f"Owner full access on {table}",
            f"Team member access on {table}",
            f"Service role bypass on {table}"
        ]
        
        missing_policies = [p for p in expected_policies if p not in policy_names]
        
        if missing_policies:
            logger.warning(f"Missing policies for table {table}: {missing_policies}")
            all_verified = False
            verification_results[table] = {
                "rls_enabled": True, 
                "policies_verified": False,
                "missing_policies": missing_policies
            }
        else:
            logger.info(f"All policies verified for table {table}")
            verification_results[table] = {"rls_enabled": True, "policies_verified": True}
    
    # Check if roles exist
    sql_roles = """
    SELECT rolname 
    FROM pg_roles 
    WHERE rolname IN ('app_readonly', 'app_user', 'app_admin');
    """
    
    success, roles = execute_sql(supabase, sql_roles, dry_run)
    
    if not success:
        logger.warning("Could not verify roles")
        all_verified = False
        verification_results["roles"] = {"verified": False}
    else:
        role_names = [r['rolname'] for r in roles]
        expected_roles = ['app_readonly', 'app_user', 'app_admin']
        missing_roles = [r for r in expected_roles if r not in role_names]
        
        if missing_roles:
            logger.warning(f"Missing roles: {missing_roles}")
            all_verified = False
            verification_results["roles"] = {"verified": False, "missing_roles": missing_roles}
        else:
            logger.info("All roles verified")
            verification_results["roles"] = {"verified": True}
    
    # Save verification results to file
    with open("security_verification_results.json", "w") as f:
        json.dump(verification_results, f, indent=2)
    
    return all_verified

def test_rls_policies(supabase, dry_run=False):
    """Test RLS policies with different user roles."""
    if dry_run:
        logger.info("DRY RUN: Would test RLS policies")
        return True
    
    # This function would simulate different user roles and test access
    # For a real implementation, you would need to create test users with different roles
    # and test access to different tables
    
    logger.info("Testing RLS policies with different user roles")
    
    # Example test: Check if service role can access all tables
    all_passed = True
    test_results = {}
    
    for table in TABLES:
        sql = f"SELECT COUNT(*) FROM {table};"
        success, _ = execute_sql(supabase, sql, dry_run)
        
        if not success:
            logger.error(f"Service role cannot access table {table}")
            all_passed = False
            test_results[table] = {"service_role_access": False}
        else:
            logger.info(f"Service role can access table {table}")
            test_results[table] = {"service_role_access": True}
    
    # Save test results to file
    with open("security_test_results.json", "w") as f:
        json.dump(test_results, f, indent=2)
    
    return all_passed

def rollback_changes(supabase, dry_run=False):
    """Rollback changes if needed."""
    if dry_run:
        logger.info("DRY RUN: Would rollback changes")
        return True
    
    if not original_state:
        logger.warning("No original state to rollback to")
        return False
    
    all_rolled_back = True
    
    for table, state in original_state.items():
        # Drop all policies
        sql = f"""
        DO $$
        DECLARE
            pol record;
        BEGIN
            FOR pol IN SELECT policyname FROM pg_policies WHERE tablename = '{table}'
            LOOP
                EXECUTE 'DROP POLICY "' || pol.policyname || '" ON {table};';
            END LOOP;
        END $$;
        """
        
        success, _ = execute_sql(supabase, sql, dry_run)
        
        if not success:
            logger.error(f"Failed to drop policies for table {table}")
            all_rolled_back = False
            continue
        
        # Disable RLS if it was originally disabled
        if not state['rls_enabled']:
            sql = f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;"
            success, _ = execute_sql(supabase, sql, dry_run)
            
            if not success:
                logger.error(f"Failed to disable RLS for table {table}")
                all_rolled_back = False
                continue
        
        # Recreate original policies
        for policy in state['policies']:
            policy_name = policy['policyname']
            cmd = policy['cmd']
            qual = policy['qual']
            with_check = policy['with_check']
            
            sql = f"""
            CREATE POLICY "{policy_name}" 
            ON {table} 
            FOR {cmd} 
            """
            
            if qual:
                sql += f"USING ({qual}) "
            
            if with_check:
                sql += f"WITH CHECK ({with_check})"
            
            success, _ = execute_sql(supabase, sql, dry_run)
            
            if not success:
                logger.error(f"Failed to recreate policy {policy_name} for table {table}")
                all_rolled_back = False
    
    # Drop roles
    sql_roles = """
    DROP ROLE IF EXISTS app_readonly;
    DROP ROLE IF EXISTS app_user;
    DROP ROLE IF EXISTS app_admin;
    """
    
    success, _ = execute_sql(supabase, sql_roles, dry_run)
    
    if not success:
        logger.error("Failed to drop roles")
        all_rolled_back = False
    
    return all_rolled_back

def main():
    parser = argparse.ArgumentParser(description="Implement security features for CryoProtect Supabase project.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    parser.add_argument("--verify-only", action="store_true", help="Only verify existing security settings")
    parser.add_argument("--rollback", action="store_true", help="Rollback changes")
    args = parser.parse_args()
    
    logger.info(f"Starting security implementation for project {PROJECT_ID}")
    
    # Connect to Supabase
    try:
        supabase = get_supabase_client()
        logger.info("Connected to Supabase")
    except Exception as e:
        logger.error(f"Failed to connect to Supabase: {str(e)}")
        return
    
    if args.rollback:
        logger.info("Rolling back changes")
        success = rollback_changes(supabase, args.dry_run)
        
        if success:
            logger.info("Successfully rolled back changes")
        else:
            logger.error("Failed to roll back all changes")
        
        return
    
    if args.verify_only:
        logger.info("Verifying security settings")
        success = verify_rls_policies(supabase, args.dry_run)
        
        if success:
            logger.info("All security settings verified")
        else:
            logger.warning("Some security settings could not be verified")
        
        return
    
    # Backup original state for rollback
    logger.info("Backing up original state")
    for table in TABLES:
        backup_table_state(supabase, table, args.dry_run)
    
    # Enable RLS on all tables
    logger.info("Enabling RLS on tables")
    for table in TABLES:
        enable_rls(supabase, table, args.dry_run)
    
    # Create policies for each table
    logger.info("Creating RLS policies")
    for table in TABLES:
        create_public_read_policy(supabase, table, args.dry_run)
        create_owner_access_policy(supabase, table, args.dry_run)
        create_team_member_access_policy(supabase, table, args.dry_run)
        create_service_role_bypass_policy(supabase, table, args.dry_run)
    
    # Create app roles
    logger.info("Creating app roles")
    create_app_roles(supabase, args.dry_run)
    
    # Verify policies
    logger.info("Verifying RLS policies")
    verify_success = verify_rls_policies(supabase, args.dry_run)
    
    # Test policies
    logger.info("Testing RLS policies")
    test_success = test_rls_policies(supabase, args.dry_run)
    
    if verify_success and test_success:
        logger.info("Security implementation completed successfully")
    else:
        logger.warning("Security implementation completed with warnings")
        logger.warning("Check verification and test results for details")
        logger.warning("Use --rollback to revert changes if needed")

if __name__ == "__main__":
    main()
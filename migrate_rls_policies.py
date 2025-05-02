#!/usr/bin/env python3
"""
CryoProtect v2 - RLS Policy Migration Script

This script migrates Row Level Security (RLS) policies from singular-named tables to plural-named tables.
It ensures that security is maintained during and after the database migration process.

Tables to migrate:
1. molecule → molecules
2. mixture → mixtures
3. experiment → experiments
4. prediction → predictions
5. project → projects

Key features:
1. Identifies existing RLS policies on singular-named tables
2. Creates equivalent policies on plural-named tables
3. Updates policy expressions that reference table names
4. Verifies that policies are correctly applied
5. Supports rollback in case of failure
6. Provides detailed logging
7. Supports dry-run mode for testing

Author: CryoProtect Team
Date: 2025-04-18
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime
import psycopg2
from psycopg2 import sql
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Configure logging
log_filename = f"rls_migration_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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

# Database connection parameters
DB_URL = os.getenv("DATABASE_URL") or os.getenv("SUPABASE_DB_URL")
if not DB_URL:
    logger.error("Database URL not found in environment variables")
    sys.exit(1)

# Migration configuration
TABLE_MAPPINGS = [
    {"singular": "molecule", "plural": "molecules"},
    {"singular": "mixture", "plural": "mixtures"},
    {"singular": "experiment", "plural": "experiments"},
    {"singular": "prediction", "plural": "predictions"},
    {"singular": "project", "plural": "projects"}
]

def connect_to_db():
    """
    Establish a connection to the database.
    
    Returns:
        tuple: (connection, cursor) or (None, None) if connection fails
    """
    try:
        conn = psycopg2.connect(DB_URL)
        conn.autocommit = False  # We'll manage transactions explicitly
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        logger.info("Successfully connected to the database")
        return conn, cursor
    except Exception as e:
        logger.error(f"Error connecting to the database: {str(e)}")
        return None, None

def get_rls_policies(cursor, schema_name='public'):
    """
    Get all RLS policies in the specified schema.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Dictionary mapping table names to their policies
    """
    try:
        logger.info(f"Getting RLS policies for schema {schema_name}")
        
        # Query to get all policies
        query = """
        SELECT 
            schemaname, 
            tablename, 
            policyname, 
            permissive,
            roles,
            cmd, 
            qual, 
            with_check
        FROM pg_policies
        WHERE schemaname = %s
        ORDER BY tablename, policyname
        """
        
        cursor.execute(query, (schema_name,))
        policies = cursor.fetchall()
        
        # Organize policies by table
        policies_by_table = {}
        for policy in policies:
            table_name = policy['tablename']
            if table_name not in policies_by_table:
                policies_by_table[table_name] = []
            policies_by_table[table_name].append(policy)
        
        logger.info(f"Found policies for {len(policies_by_table)} tables")
        return policies_by_table
    
    except Exception as e:
        logger.error(f"Error getting RLS policies: {str(e)}")
        return {}

def check_rls_enabled(cursor, schema_name, table_name):
    """
    Check if RLS is enabled for a table.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        table_name (str): Table name
        
    Returns:
        bool: True if RLS is enabled, False otherwise
    """
    try:
        query = """
        SELECT relrowsecurity
        FROM pg_class c
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE n.nspname = %s AND c.relkind = 'r' AND c.relname = %s
        """
        
        cursor.execute(query, (schema_name, table_name))
        result = cursor.fetchone()
        
        if result:
            return result['relrowsecurity']
        return False
    
    except Exception as e:
        logger.error(f"Error checking if RLS is enabled for {table_name}: {str(e)}")
        return False

def enable_rls(cursor, schema_name, table_name):
    """
    Enable RLS on a table.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        table_name (str): Table name
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Enabling RLS on {schema_name}.{table_name}")
        
        # Check if RLS is already enabled
        if check_rls_enabled(cursor, schema_name, table_name):
            logger.info(f"RLS already enabled on {schema_name}.{table_name}")
            return True
        
        # Enable RLS
        cursor.execute(
            sql.SQL("ALTER TABLE {}.{} ENABLE ROW LEVEL SECURITY").format(
                sql.Identifier(schema_name),
                sql.Identifier(table_name)
            )
        )
        
        logger.info(f"Enabled RLS on {schema_name}.{table_name}")
        return True
    
    except Exception as e:
        logger.error(f"Error enabling RLS on {schema_name}.{table_name}: {str(e)}")
        return False

def update_policy_expression(expression, table_mappings):
    """
    Update policy expression to reference plural table names.
    
    Args:
        expression (str): SQL expression used in policy
        table_mappings (list): List of dictionaries mapping singular to plural table names
        
    Returns:
        str: Updated expression
    """
    if not expression:
        return expression
    
    updated_expression = expression
    
    for mapping in table_mappings:
        singular = mapping['singular']
        plural = mapping['plural']
        
        # Replace table references in expressions
        # Be careful to only replace full table names, not partial matches
        updated_expression = updated_expression.replace(f" {singular} ", f" {plural} ")
        updated_expression = updated_expression.replace(f" {singular}.", f" {plural}.")
        updated_expression = updated_expression.replace(f"({singular}.", f"({plural}.")
        updated_expression = updated_expression.replace(f"FROM {singular}", f"FROM {plural}")
        updated_expression = updated_expression.replace(f"JOIN {singular}", f"JOIN {plural}")
        updated_expression = updated_expression.replace(f"EXISTS (SELECT 1 FROM {singular}", f"EXISTS (SELECT 1 FROM {plural}")
    
    return updated_expression

def create_policy(cursor, schema_name, table_name, policy_data):
    """
    Create a policy on a table.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        table_name (str): Table name
        policy_data (dict): Policy data
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        policy_name = policy_data['policyname']
        logger.info(f"Creating policy '{policy_name}' on {schema_name}.{table_name}")
        
        # Build the CREATE POLICY statement
        create_policy_sql = sql.SQL("CREATE POLICY {} ON {}.{}").format(
            sql.Identifier(policy_name),
            sql.Identifier(schema_name),
            sql.Identifier(table_name)
        )
        
        # Add AS clause if roles are specified
        if policy_data['roles'] and policy_data['roles'] != '{public}':
            roles_list = policy_data['roles'][1:-1].split(',')  # Remove {} and split
            create_policy_sql = sql.SQL("{} AS {}").format(
                create_policy_sql,
                sql.SQL(', '.join(roles_list))
            )
        
        # Add FOR clause if command is specified
        if policy_data['cmd'] and policy_data['cmd'] != 'ALL':
            create_policy_sql = sql.SQL("{} FOR {}").format(
                create_policy_sql,
                sql.SQL(policy_data['cmd'])
            )
        
        # Add USING clause if qual is specified
        if policy_data['qual']:
            create_policy_sql = sql.SQL("{} USING ({})").format(
                create_policy_sql,
                sql.SQL(policy_data['qual'])
            )
        
        # Add WITH CHECK clause if with_check is specified
        if policy_data['with_check']:
            create_policy_sql = sql.SQL("{} WITH CHECK ({})").format(
                create_policy_sql,
                sql.SQL(policy_data['with_check'])
            )
        
        # Execute the CREATE POLICY statement
        cursor.execute(create_policy_sql)
        
        logger.info(f"Created policy '{policy_name}' on {schema_name}.{table_name}")
        return True
    
    except Exception as e:
        logger.error(f"Error creating policy on {schema_name}.{table_name}: {str(e)}")
        return False

def migrate_policies(cursor, schema_name, policies_by_table, table_mappings, dry_run=False):
    """
    Migrate policies from singular to plural tables.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        policies_by_table (dict): Dictionary mapping table names to their policies
        table_mappings (list): List of dictionaries mapping singular to plural table names
        dry_run (bool): If True, don't execute any changes
        
    Returns:
        dict: Migration report
    """
    migration_report = {
        "timestamp": datetime.now().isoformat(),
        "migrated_policies": [],
        "failed_policies": [],
        "skipped_tables": []
    }
    
    # Create a mapping from singular to plural table names for easier lookup
    table_name_map = {mapping['singular']: mapping['plural'] for mapping in table_mappings}
    
    # Process each singular table
    for singular_table, policies in policies_by_table.items():
        # Skip if this is not a table we need to migrate
        if singular_table not in table_name_map:
            continue
        
        plural_table = table_name_map[singular_table]
        logger.info(f"Migrating policies from {singular_table} to {plural_table}")
        
        # Check if plural table exists
        cursor.execute(
            """
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
            """,
            (schema_name, plural_table)
        )
        table_exists = cursor.fetchone()['exists']
        
        if not table_exists:
            logger.warning(f"Plural table {plural_table} does not exist, skipping")
            migration_report["skipped_tables"].append({
                "singular": singular_table,
                "plural": plural_table,
                "reason": "Plural table does not exist"
            })
            continue
        
        # Enable RLS on plural table if not already enabled
        if not dry_run:
            if not enable_rls(cursor, schema_name, plural_table):
                logger.error(f"Failed to enable RLS on {plural_table}, skipping")
                migration_report["skipped_tables"].append({
                    "singular": singular_table,
                    "plural": plural_table,
                    "reason": "Failed to enable RLS"
                })
                continue
        else:
            logger.info(f"DRY RUN: Would enable RLS on {schema_name}.{plural_table}")
        
        # Migrate each policy
        for policy in policies:
            policy_name = policy['policyname']
            logger.info(f"Migrating policy '{policy_name}' from {singular_table} to {plural_table}")
            
            # Update policy expressions to reference plural tables
            updated_qual = update_policy_expression(policy['qual'], table_mappings)
            updated_with_check = update_policy_expression(policy['with_check'], table_mappings)
            
            # Create a copy of the policy with updated expressions
            updated_policy = policy.copy()
            updated_policy['qual'] = updated_qual
            updated_policy['with_check'] = updated_with_check
            
            # Create the policy on the plural table
            if not dry_run:
                # First drop the policy if it already exists
                try:
                    cursor.execute(
                        sql.SQL("DROP POLICY IF EXISTS {} ON {}.{}").format(
                            sql.Identifier(policy_name),
                            sql.Identifier(schema_name),
                            sql.Identifier(plural_table)
                        )
                    )
                except Exception as e:
                    logger.warning(f"Error dropping existing policy '{policy_name}' on {plural_table}: {str(e)}")
                
                # Create the policy
                if create_policy(cursor, schema_name, plural_table, updated_policy):
                    migration_report["migrated_policies"].append({
                        "table": singular_table,
                        "plural_table": plural_table,
                        "policy_name": policy_name,
                        "original_qual": policy['qual'],
                        "updated_qual": updated_qual,
                        "original_with_check": policy['with_check'],
                        "updated_with_check": updated_with_check
                    })
                else:
                    migration_report["failed_policies"].append({
                        "table": singular_table,
                        "plural_table": plural_table,
                        "policy_name": policy_name,
                        "reason": "Failed to create policy"
                    })
            else:
                logger.info(f"DRY RUN: Would create policy '{policy_name}' on {schema_name}.{plural_table}")
                logger.info(f"  Original USING: {policy['qual']}")
                logger.info(f"  Updated USING: {updated_qual}")
                logger.info(f"  Original WITH CHECK: {policy['with_check']}")
                logger.info(f"  Updated WITH CHECK: {updated_with_check}")
                
                migration_report["migrated_policies"].append({
                    "table": singular_table,
                    "plural_table": plural_table,
                    "policy_name": policy_name,
                    "original_qual": policy['qual'],
                    "updated_qual": updated_qual,
                    "original_with_check": policy['with_check'],
                    "updated_with_check": updated_with_check
                })
    
    return migration_report

def verify_policies(cursor, schema_name, migration_report):
    """
    Verify that policies were correctly migrated.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        migration_report (dict): Migration report
        
    Returns:
        dict: Verification report
    """
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "verified_policies": [],
        "failed_verifications": []
    }
    
    # Get current policies
    policies_by_table = get_rls_policies(cursor, schema_name)
    
    # Verify each migrated policy
    for policy_info in migration_report["migrated_policies"]:
        plural_table = policy_info["plural_table"]
        policy_name = policy_info["policy_name"]
        
        # Check if the policy exists on the plural table
        if plural_table in policies_by_table:
            table_policies = policies_by_table[plural_table]
            matching_policies = [p for p in table_policies if p['policyname'] == policy_name]
            
            if matching_policies:
                policy = matching_policies[0]
                
                # Verify policy expressions
                qual_matches = policy['qual'] == policy_info['updated_qual']
                with_check_matches = (
                    (policy['with_check'] is None and policy_info['updated_with_check'] is None) or
                    policy['with_check'] == policy_info['updated_with_check']
                )
                
                if qual_matches and with_check_matches:
                    verification_report["verified_policies"].append({
                        "table": plural_table,
                        "policy_name": policy_name,
                        "status": "verified"
                    })
                else:
                    verification_report["failed_verifications"].append({
                        "table": plural_table,
                        "policy_name": policy_name,
                        "status": "expression_mismatch",
                        "expected_qual": policy_info['updated_qual'],
                        "actual_qual": policy['qual'],
                        "expected_with_check": policy_info['updated_with_check'],
                        "actual_with_check": policy['with_check']
                    })
            else:
                verification_report["failed_verifications"].append({
                    "table": plural_table,
                    "policy_name": policy_name,
                    "status": "policy_not_found"
                })
        else:
            verification_report["failed_verifications"].append({
                "table": plural_table,
                "policy_name": policy_name,
                "status": "table_has_no_policies"
            })
    
    return verification_report

def rollback_migration(cursor, schema_name, migration_report):
    """
    Rollback the policy migration.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        migration_report (dict): Migration report
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info("Rolling back policy migration")
    
    success = True
    
    # Drop migrated policies
    for policy_info in migration_report["migrated_policies"]:
        plural_table = policy_info["plural_table"]
        policy_name = policy_info["policy_name"]
        
        logger.info(f"Dropping policy '{policy_name}' from {schema_name}.{plural_table}")
        
        try:
            cursor.execute(
                sql.SQL("DROP POLICY IF EXISTS {} ON {}.{}").format(
                    sql.Identifier(policy_name),
                    sql.Identifier(schema_name),
                    sql.Identifier(plural_table)
                )
            )
        except Exception as e:
            logger.error(f"Error dropping policy '{policy_name}' from {schema_name}.{plural_table}: {str(e)}")
            success = False
    
    return success

def main():
    """Main function to run the RLS policy migration script."""
    parser = argparse.ArgumentParser(description='CryoProtect v2 RLS Policy Migration Script')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without committing changes')
    parser.add_argument('--verify-only', action='store_true', help='Only verify existing policies without migration')
    parser.add_argument('--rollback', action='store_true', help='Rollback migration using report file')
    parser.add_argument('--report', type=str, help='Path to migration report for verification or rollback')
    parser.add_argument('--schema', type=str, default='public', help='Database schema (default: public)')
    args = parser.parse_args()
    
    # Connect to the database
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        sys.exit(1)
    
    try:
        schema_name = args.schema
        
        if args.rollback:
            if not args.report:
                logger.error("Migration report is required for rollback")
                sys.exit(1)
            
            # Load migration report
            with open(args.report, 'r') as f:
                migration_report = json.load(f)
            
            # Perform rollback
            if rollback_migration(cursor, schema_name, migration_report):
                logger.info("Migration rollback completed successfully")
                conn.commit()
            else:
                logger.error("Migration rollback failed")
                conn.rollback()
                sys.exit(1)
        
        elif args.verify_only:
            if not args.report:
                logger.error("Migration report is required for verification")
                sys.exit(1)
            
            # Load migration report
            with open(args.report, 'r') as f:
                migration_report = json.load(f)
            
            # Verify policies
            verification_report = verify_policies(cursor, schema_name, migration_report)
            
            # Save verification report
            verification_file = f"rls_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(verification_file, 'w') as f:
                json.dump(verification_report, f, indent=2)
            
            logger.info(f"Verification report saved to {verification_file}")
            
            # Print summary
            verified_count = len(verification_report["verified_policies"])
            failed_count = len(verification_report["failed_verifications"])
            total_count = verified_count + failed_count
            
            logger.info(f"Verification summary: {verified_count}/{total_count} policies verified")
            
            if failed_count > 0:
                logger.warning(f"{failed_count} policies failed verification")
                for failed in verification_report["failed_verifications"]:
                    logger.warning(f"  {failed['table']}.{failed['policy_name']}: {failed['status']}")
        
        else:
            # Get existing policies
            policies_by_table = get_rls_policies(cursor, schema_name)
            
            # Migrate policies
            migration_report = migrate_policies(
                cursor, 
                schema_name, 
                policies_by_table, 
                TABLE_MAPPINGS, 
                dry_run=args.dry_run
            )
            
            # Save migration report
            report_file = f"rls_migration_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(report_file, 'w') as f:
                json.dump(migration_report, f, indent=2)
            
            logger.info(f"Migration report saved to {report_file}")
            
            # Print summary
            migrated_count = len(migration_report["migrated_policies"])
            failed_count = len(migration_report["failed_policies"])
            skipped_count = len(migration_report["skipped_tables"])
            
            logger.info(f"Migration summary:")
            logger.info(f"  {migrated_count} policies migrated")
            logger.info(f"  {failed_count} policies failed to migrate")
            logger.info(f"  {skipped_count} tables skipped")
            
            if args.dry_run:
                logger.info("Dry run completed, rolling back changes")
                conn.rollback()
            else:
                # Verify policies
                verification_report = verify_policies(cursor, schema_name, migration_report)
                
                # Save verification report
                verification_file = f"rls_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
                with open(verification_file, 'w') as f:
                    json.dump(verification_report, f, indent=2)
                
                logger.info(f"Verification report saved to {verification_file}")
                
                # Print verification summary
                verified_count = len(verification_report["verified_policies"])
                failed_count = len(verification_report["failed_verifications"])
                total_count = verified_count + failed_count
                
                logger.info(f"Verification summary: {verified_count}/{total_count} policies verified")
                
                if failed_count > 0:
                    logger.warning(f"{failed_count} policies failed verification")
                    for failed in verification_report["failed_verifications"]:
                        logger.warning(f"  {failed['table']}.{failed['policy_name']}: {failed['status']}")
                    
                    # Ask for confirmation before committing
                    confirm = input("Some policies failed verification. Commit anyway? (y/n): ")
                    if confirm.lower() != 'y':
                        logger.info("Changes rolled back")
                        conn.rollback()
                        sys.exit(1)
                
                logger.info("Committing changes")
                conn.commit()
    
    except Exception as e:
        logger.error(f"Error during migration: {str(e)}")
        conn.rollback()
        sys.exit(1)
    
    finally:
        # Close database connection
        if conn:
            conn.close()

if __name__ == "__main__":
    main()
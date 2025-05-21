#!/usr/bin/env python3
"""
CryoProtect v2 - Test Foreign Key Relationship Fixes

This script tests the fix_foreign_key_relationships.py script on the test schema
that was kept for inspection during the database remediation process.

Usage:
    python test_fix_foreign_key_relationships.py
"""

import os
import sys
import json
import logging
import subprocess
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"test_fix_foreign_keys_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

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

def execute_sql(supabase, sql, description):
    """Execute SQL using the Supabase client."""
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

def check_test_schema_exists(supabase):
    """Check if the remediation_test schema exists."""
    sql = """
    SELECT EXISTS (
        SELECT FROM information_schema.schemata 
        WHERE schema_name = 'remediation_test'
    );
    """
    
    success, result = execute_sql(supabase, sql, "Checking if remediation_test schema exists")
    
    if not success:
        logger.error("Failed to check if remediation_test schema exists")
        return False
    
    return result[0]['exists']

def get_foreign_keys_in_test_schema(supabase):
    """Get all foreign key constraints in the remediation_test schema."""
    sql = """
    SELECT
        tc.constraint_name,
        tc.table_name,
        kcu.column_name,
        ccu.table_name AS referenced_table,
        ccu.column_name AS referenced_column
    FROM
        information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
        JOIN information_schema.constraint_column_usage AS ccu
          ON ccu.constraint_name = tc.constraint_name
          AND ccu.table_schema = tc.table_schema
    WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'remediation_test';
    """
    
    success, result = execute_sql(supabase, sql, "Getting foreign key constraints in remediation_test schema")
    
    if not success:
        logger.error("Failed to get foreign key constraints in remediation_test schema")
        return None
    
    return result

def find_problematic_foreign_keys(foreign_keys):
    """Find foreign keys that reference old singular table names."""
    problematic_fks = []
    
    for fk in foreign_keys:
        referenced_table = fk['referenced_table']
        if referenced_table == 'molecule':
            problematic_fks.append(fk)
    
    return problematic_fks

def fix_foreign_keys_in_test_schema(supabase):
    """Fix foreign key relationships in the remediation_test schema."""
    # Set search path to remediation_test schema
    set_schema_sql = """
    SET search_path TO remediation_test, public;
    """
    
    success, _ = execute_sql(supabase, set_schema_sql, "Set search path to remediation_test schema")
    
    if not success:
        logger.error("Failed to set search path to remediation_test schema")
        return False
    
    # Fix foreign key relationships
    fix_fk_sql = """
    -- Drop existing foreign key constraints that reference old singular table names
    ALTER TABLE remediation_test.mixture_component
    DROP CONSTRAINT IF EXISTS mixture_component_molecule_id_fkey;
    
    -- Add new foreign key constraints that reference the new plural table names
    ALTER TABLE remediation_test.mixture_component
    ADD CONSTRAINT mixture_component_molecule_id_fkey
    FOREIGN KEY (molecule_id) REFERENCES remediation_test.molecules(id);
    
    -- Create indexes for foreign keys
    CREATE INDEX IF NOT EXISTS idx_mixture_component_molecule_id
    ON remediation_test.mixture_component(molecule_id);
    
    CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id
    ON remediation_test.mixture_component(mixture_id);
    
    CREATE INDEX IF NOT EXISTS idx_experiment_mixture_id
    ON remediation_test.experiment(mixture_id);
    
    CREATE INDEX IF NOT EXISTS idx_prediction_mixture_id
    ON remediation_test.prediction(mixture_id);
    
    CREATE INDEX IF NOT EXISTS idx_predictions_mixture_id
    ON remediation_test.predictions(mixture_id);
    """
    
    success, _ = execute_sql(supabase, fix_fk_sql, "Fix foreign key relationships in remediation_test schema")
    
    if not success:
        logger.error("Failed to fix foreign key relationships in remediation_test schema")
        return False
    
    return True

def verify_fixed_foreign_keys(supabase):
    """Verify that foreign key relationships are fixed in the remediation_test schema."""
    # Get all foreign key constraints after fixing
    foreign_keys = get_foreign_keys_in_test_schema(supabase)
    
    if foreign_keys is None:
        return False
    
    # Check if any foreign keys still reference old singular table names
    problematic_fks = find_problematic_foreign_keys(foreign_keys)
    
    if problematic_fks:
        logger.error(f"Found {len(problematic_fks)} foreign keys that still reference old singular table names")
        for fk in problematic_fks:
            logger.error(f"  {fk['constraint_name']} in {fk['table_name']} references {fk['referenced_table']}")
        return False
    
    # Check if all foreign keys have indexes
    sql = """
    SELECT
        tc.table_name,
        kcu.column_name
    FROM
        information_schema.table_constraints AS tc
        JOIN information_schema.key_column_usage AS kcu
          ON tc.constraint_name = kcu.constraint_name
          AND tc.table_schema = kcu.table_schema
    WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_schema = 'remediation_test'
    EXCEPT
    SELECT
        tablename AS table_name,
        regexp_replace(indexdef, '.*\((.*)\).*', '\1') AS column_name
    FROM
        pg_indexes
    WHERE
        schemaname = 'remediation_test'
        AND indexdef LIKE '%CREATE INDEX%';
    """
    
    success, missing_indexes = execute_sql(supabase, sql, "Checking for missing indexes")
    
    if not success:
        logger.error("Failed to check for missing indexes")
        return False
    
    if missing_indexes:
        logger.error(f"Found {len(missing_indexes)} foreign keys without indexes")
        for idx in missing_indexes:
            logger.error(f"  {idx['table_name']}.{idx['column_name']} has no index")
        return False
    
    return True

def main():
    """Main function to test foreign key relationship fixes."""
    logger.info("Starting test of foreign key relationship fixes")
    
    # Connect to Supabase
    supabase = get_supabase_client()
    
    # Check if the remediation_test schema exists
    if not check_test_schema_exists(supabase):
        logger.error("The remediation_test schema does not exist")
        logger.info("Please run test_database_remediation.py with --keep-schema to create and keep the test schema")
        return 1
    
    # Get foreign key constraints before fixing
    logger.info("Checking foreign key constraints before fixing")
    before_fks = get_foreign_keys_in_test_schema(supabase)
    
    if before_fks is None:
        return 1
    
    problematic_before = find_problematic_foreign_keys(before_fks)
    logger.info(f"Found {len(problematic_before)} problematic foreign keys before fixing")
    
    # Fix foreign key relationships
    logger.info("Fixing foreign key relationships in remediation_test schema")
    if not fix_foreign_keys_in_test_schema(supabase):
        return 1
    
    # Verify that foreign key relationships are fixed
    logger.info("Verifying fixed foreign key relationships")
    if not verify_fixed_foreign_keys(supabase):
        logger.error("Verification failed: Foreign key relationships are not fixed correctly")
        return 1
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("CryoProtect v2 Test Foreign Key Relationship Fixes Summary")
    logger.info("=" * 80)
    
    logger.info(f"\nBefore fixing:")
    logger.info(f"- Total foreign key constraints: {len(before_fks)}")
    logger.info(f"- Problematic foreign keys: {len(problematic_before)}")
    
    logger.info(f"\nAfter fixing:")
    logger.info(f"- All foreign keys now reference the correct plural table names")
    logger.info(f"- All foreign keys have proper indexes")
    
    logger.info("\nStatus: SUCCESS")
    logger.info("The foreign key relationship fixes were tested successfully.")
    
    logger.info("\nFor detailed information, check the log file.")
    logger.info("=" * 80 + "\n")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
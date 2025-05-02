#!/usr/bin/env python3
"""
CryoProtect v2 - Apply Service Role RLS Migration

This script applies the service role RLS migration to allow the service role
to insert data into tables with RLS enabled.

Usage:
    python apply_service_role_rls.py [--dry-run]
"""

import os
import argparse
import logging
from dotenv import load_dotenv
from service_role_helper import get_supabase_client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("service_role_rls_migration.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def read_migration_file(file_path):
    """Read the migration SQL file."""
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading migration file: {str(e)}")
        return None

def apply_migration(supabase, sql, dry_run=False):
    """Apply the migration SQL to the database."""
    if dry_run:
        logger.info("DRY RUN: Would execute the following SQL:")
        logger.info(sql)
        return True
    
    try:
        # Execute the SQL directly using the Supabase client
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error applying migration: {response.error}")
            return False
        
        logger.info("Migration applied successfully")
        return True
    except Exception as e:
        logger.error(f"Error applying migration: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Apply service role RLS migration.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
    args = parser.parse_args()
    
    logger.info("Starting Service Role RLS Migration")
    
    # Connect to Supabase using service role
    supabase = get_supabase_client()
    
    # Read migration file
    migration_path = os.path.join('migrations', '007_service_role_rls.sql')
    sql = read_migration_file(migration_path)
    
    if not sql:
        logger.error(f"Failed to read migration file: {migration_path}")
        return
    
    # Apply migration
    success = apply_migration(supabase, sql, args.dry_run)
    
    if success:
        logger.info("Service Role RLS Migration complete")
        
        # Verify the policies were created
        if not args.dry_run:
            try:
                # Check if the policy exists for the molecule table
                response = supabase.rpc('exec_sql', {
                    'query': "SELECT * FROM pg_policies WHERE tablename = 'molecule' AND policyname = 'Allow service role inserts on molecule'"
                }).execute()
                
                if hasattr(response, 'data') and response.data:
                    logger.info("Verified policy creation: Allow service role inserts on molecule")
                else:
                    logger.warning("Could not verify policy creation. Please check manually.")
            except Exception as e:
                logger.warning(f"Error verifying policy creation: {str(e)}")
    else:
        logger.error("Service Role RLS Migration failed")

if __name__ == "__main__":
    main()
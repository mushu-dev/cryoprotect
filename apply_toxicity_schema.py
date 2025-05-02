#!/usr/bin/env python3
"""
Apply the toxicity schema migration to the database.

This script applies the 012_toxicity_schema.sql migration to the database,
which sets up the schema for Tox21 toxicity data integration.
"""

import os
import sys
import logging
import argparse
from supabase import create_client, Client
from config import Config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def read_migration_file(file_path):
    """Read the migration SQL file."""
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading migration file: {str(e)}")
        raise

def apply_migration(supabase, sql, dry_run=False):
    """Apply the migration SQL to the database."""
    if dry_run:
        logger.info("Dry run - would execute the following SQL:")
        logger.info(sql)
        return True
    
    # Split the SQL into individual statements
    # This is a simple split by semicolon, which may not work for all SQL
    # but should be sufficient for our migration script
    statements = []
    current_statement = ""
    
    for line in sql.split('\n'):
        # Skip comments and empty lines when determining statement boundaries
        stripped_line = line.strip()
        if not stripped_line or stripped_line.startswith('--'):
            current_statement += line + '\n'
            continue
            
        current_statement += line + '\n'
        
        if stripped_line.endswith(';'):
            statements.append(current_statement)
            current_statement = ""
    
    # Add any remaining statement
    if current_statement.strip():
        statements.append(current_statement)
    
    logger.info(f"Split migration into {len(statements)} statements")
    
    # Execute each statement
    success_count = 0
    for i, statement in enumerate(statements):
        if not statement.strip():
            continue
            
        try:
            # For statements that don't return results, we need to use a different approach
            # We'll use the REST API directly to execute the SQL
            response = supabase.postgrest.rpc('exec_sql', {'query': f"SELECT 1 AS result FROM (SELECT 1) AS t WHERE EXISTS (DO $$ BEGIN {statement} END $$)"}).execute()
            
            if hasattr(response, 'error') and response.error:
                logger.error(f"Error executing statement {i+1}: {response.error}")
                logger.error(f"Statement: {statement}")
                return False
            
            success_count += 1
            
        except Exception as e:
            logger.error(f"Error executing statement {i+1}: {str(e)}")
            logger.error(f"Statement: {statement}")
            return False
    
    logger.info(f"Successfully executed {success_count} statements")
    return True

def main():
    """Main function to apply the toxicity schema migration."""
    parser = argparse.ArgumentParser(description='Apply toxicity schema migration')
    parser.add_argument('--dry-run', action='store_true', help='Print SQL without executing')
    args = parser.parse_args()
    
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Read migration SQL
    migration_path = os.path.join('migrations', '012_toxicity_schema.sql')
    try:
        sql = read_migration_file(migration_path)
    except Exception as e:
        logger.error(f"Failed to read migration file: {str(e)}")
        return False
    
    # Apply migration
    success = apply_migration(supabase, sql, args.dry_run)
    
    if success:
        logger.info("Toxicity schema migration completed successfully")
        return True
    else:
        logger.error("Toxicity schema migration failed")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
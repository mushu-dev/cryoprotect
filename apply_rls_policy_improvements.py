#!/usr/bin/env python
"""
apply_rls_policy_improvements.py: Applies the RLS policy improvements migration
to improve security and performance.

This script can be run directly or via the Supabase MCP to apply the 
improvements to the database.
"""

import os
import sys
import argparse
import psycopg2
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f"logs/rls_improvements_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
os.makedirs("logs", exist_ok=True)

def get_db_connection():
    """Create a database connection using environment variables."""
    try:
        # Get database connection parameters from environment variables
        db_params = {
            "host": os.environ.get("SUPABASE_DB_HOST"),
            "port": os.environ.get("SUPABASE_DB_PORT", "5432"),
            "database": os.environ.get("SUPABASE_DB_NAME", "postgres"),
            "user": os.environ.get("SUPABASE_DB_USER"),
            "password": os.environ.get("SUPABASE_DB_PASSWORD"),
            "sslmode": os.environ.get("SUPABASE_DB_SSLMODE", "require")
        }
        
        # Check if all required parameters are available
        for key, value in db_params.items():
            if not value and key != "sslmode":
                logger.error(f"Missing required database parameter: {key}")
                print(f"ERROR: Missing required database parameter: {key}")
                print("Please set all required environment variables:")
                print("  SUPABASE_DB_HOST")
                print("  SUPABASE_DB_PORT (defaults to 5432)")
                print("  SUPABASE_DB_NAME (defaults to postgres)")
                print("  SUPABASE_DB_USER")
                print("  SUPABASE_DB_PASSWORD")
                print("  SUPABASE_DB_SSLMODE (defaults to require)")
                sys.exit(1)
        
        # Create a connection to the database
        conn = psycopg2.connect(**db_params)
        return conn
    
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        print(f"ERROR: Failed to connect to the database: {e}")
        sys.exit(1)

def apply_migration_file(conn, file_path):
    """Apply the migration file to the database."""
    try:
        # Read the migration file
        with open(file_path, 'r') as file:
            sql = file.read()
        
        # Execute the SQL
        cursor = conn.cursor()
        
        # Split on semicolons to execute multiple statements
        # This approach is simple but doesn't handle complex SQL correctly
        # For a more robust solution, consider using a proper SQL parser
        for statement in sql.split(';'):
            # Skip empty statements
            if statement.strip():
                try:
                    cursor.execute(statement)
                    logger.info(f"Executed SQL statement: {statement[:100]}...")
                except Exception as e:
                    logger.error(f"Error executing SQL statement: {statement[:100]}...\nError: {e}")
                    raise
        
        # Commit the transaction
        conn.commit()
        cursor.close()
        
        return True
    
    except Exception as e:
        logger.error(f"Migration error: {e}")
        conn.rollback()
        return False

def record_migration(conn, migration_name):
    """Record the migration in the migrations table."""
    try:
        cursor = conn.cursor()
        
        # Check if migrations table exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' 
                AND table_name = 'migrations'
            );
        """)
        
        if not cursor.fetchone()[0]:
            # Create migrations table if it doesn't exist
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS migrations (
                    id SERIAL PRIMARY KEY,
                    name VARCHAR(255) NOT NULL,
                    applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
                );
            """)
        
        # Record the migration
        cursor.execute("""
            INSERT INTO migrations (name, applied_at)
            VALUES (%s, NOW());
        """, (migration_name,))
        
        conn.commit()
        cursor.close()
        
        return True
    
    except Exception as e:
        logger.error(f"Error recording migration: {e}")
        conn.rollback()
        return False

def main():
    """Main function to apply the RLS policy improvements."""
    parser = argparse.ArgumentParser(description='Apply RLS policy improvements to the database.')
    parser.add_argument('--migration-file', default='migrations/020_rls_policy_improvements.sql',
                        help='Path to the migration file (default: migrations/020_rls_policy_improvements.sql)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Dry run mode - prints statements but does not execute them')
    args = parser.parse_args()
    
    # Validate the migration file exists
    if not os.path.exists(args.migration_file):
        logger.error(f"Migration file not found: {args.migration_file}")
        print(f"ERROR: Migration file not found: {args.migration_file}")
        sys.exit(1)
    
    # Get the migration filename without path for recording
    migration_name = os.path.basename(args.migration_file)
    
    # If dry run mode, just print the file content
    if args.dry_run:
        with open(args.migration_file, 'r') as file:
            print(f"DRY RUN: Would apply migration {migration_name}")
            print(f"\n{file.read()}")
        return
    
    # Get a database connection
    conn = get_db_connection()
    
    try:
        # Apply the migration
        logger.info(f"Applying migration {migration_name}...")
        if apply_migration_file(conn, args.migration_file):
            logger.info(f"Successfully applied migration {migration_name}")
            
            # Record the migration
            if record_migration(conn, migration_name):
                logger.info(f"Successfully recorded migration {migration_name}")
                print(f"Successfully applied and recorded migration {migration_name}")
            else:
                logger.warning(f"Migration {migration_name} was applied but could not be recorded")
                print(f"WARNING: Migration {migration_name} was applied but could not be recorded")
        else:
            logger.error(f"Failed to apply migration {migration_name}")
            print(f"ERROR: Failed to apply migration {migration_name}")
            sys.exit(1)
    
    finally:
        # Close the database connection
        conn.close()

if __name__ == "__main__":
    main()
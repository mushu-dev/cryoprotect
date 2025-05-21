#!/usr/bin/env python3
"""
fix_consolidated_migration.py: Apply consolidated molecule migrations directly

This script directly applies the SQL migrations needed for the consolidated molecules feature,
including fixing invalid molecule_status values before applying the constraints.
"""

import os
import sys
import getpass
import argparse
import psycopg2
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_db_connection(host=None, port=None, dbname=None, user=None, password=None):
    """Create a database connection."""
    try:
        # Use provided values or try to find them
        db_host = host or os.environ.get("SUPABASE_DB_HOST") or os.environ.get("DB_HOST")
        db_port = port or os.environ.get("SUPABASE_DB_PORT") or os.environ.get("DB_PORT", "5432")
        db_name = dbname or os.environ.get("SUPABASE_DB_NAME") or os.environ.get("DB_NAME", "postgres")
        db_user = user or os.environ.get("SUPABASE_DB_USER") or os.environ.get("DB_USER")
        db_password = password or os.environ.get("SUPABASE_DB_PASSWORD") or os.environ.get("DB_PASSWORD")
        
        # If any required parameter is missing, prompt for it
        if not db_host:
            db_host = input("Database host: ")
        if not db_port:
            db_port = input("Database port (5432): ") or "5432"
        if not db_name:
            db_name = input("Database name (postgres): ") or "postgres"
        if not db_user:
            db_user = input("Database user: ")
        if not db_password:
            db_password = getpass.getpass("Database password: ")
        
        # Create a connection to the database
        logger.info(f"Connecting to database {db_name} at {db_host}:{db_port} as {db_user}")
        print(f"Connecting to database {db_name} at {db_host}:{db_port} as {db_user}")
        
        conn = psycopg2.connect(
            host=db_host,
            port=db_port,
            database=db_name,
            user=db_user,
            password=db_password
        )
        return conn
    
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        print(f"ERROR: Failed to connect to the database: {e}")
        sys.exit(1)

def fix_molecule_status(conn):
    """Fix invalid molecule_status values in the consolidated_molecules table."""
    try:
        conn.autocommit = True
        cursor = conn.cursor()
        
        try:
            # Check if table exists
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_name = 'consolidated_molecules'
                )
            """)
            table_exists = cursor.fetchone()[0]
            
            if not table_exists:
                logger.error("consolidated_molecules table does not exist")
                print("ERROR: consolidated_molecules table does not exist")
                return False
            
            # Find invalid molecule_status values
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            invalid_count = cursor.fetchone()[0]
            
            if invalid_count == 0:
                logger.info("No invalid molecule_status values found")
                print("No invalid molecule_status values found. Table is already valid.")
                return True
            
            # Print the invalid values for reporting
            cursor.execute("""
                SELECT molecule_id, molecule_status 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
                LIMIT 10
            """)
            invalid_examples = cursor.fetchall()
            
            print(f"Found {invalid_count} invalid molecule_status values.")
            print("Examples of invalid values:")
            for example in invalid_examples:
                print(f"  Molecule ID: {example[0]}, Status: {example[1]}")
            
            # Update invalid values to 'original'
            cursor.execute("""
                UPDATE consolidated_molecules 
                SET molecule_status = 'original' 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            
            print(f"Updated {invalid_count} records to have molecule_status = 'original'")
            logger.info(f"Updated {invalid_count} records to have molecule_status = 'original'")
            
            # Verify fix
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            remaining_invalid = cursor.fetchone()[0]
            
            if remaining_invalid == 0:
                print("Verification successful. All molecule_status values are now valid.")
                logger.info("Verification successful. All molecule_status values are now valid.")
                return True
            else:
                print(f"ERROR: {remaining_invalid} invalid values still remain after update.")
                logger.error(f"{remaining_invalid} invalid values still remain after update.")
                return False
                
        except Exception as e:
            logger.error(f"Error fixing molecule_status values: {e}")
            print(f"ERROR: Failed to fix molecule_status values: {e}")
            return False
        finally:
            cursor.close()
    
    except Exception as e:
        logger.error(f"Error: {e}")
        print(f"ERROR: {e}")
        return False

def apply_migration(conn, migration_file_path):
    """Apply a migration file."""
    try:
        # Read the migration file
        with open(migration_file_path, 'r') as f:
            sql = f.read()
        
        conn.autocommit = True
        cursor = conn.cursor()
        
        try:
            # Execute the SQL
            logger.info(f"Executing SQL from {migration_file_path}")
            print(f"Applying migration: {migration_file_path}")
            cursor.execute(sql)
            print(f"Successfully applied migration: {migration_file_path}")
            return True
        except Exception as e:
            logger.error(f"Error executing SQL: {e}")
            print(f"ERROR: Failed to apply migration: {e}")
            return False
        finally:
            cursor.close()
    
    except Exception as e:
        logger.error(f"Error reading migration file: {e}")
        print(f"ERROR: Failed to read migration file: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Fix and apply consolidated molecule migrations')
    parser.add_argument('--host', help='Database host')
    parser.add_argument('--port', help='Database port')
    parser.add_argument('--dbname', help='Database name')
    parser.add_argument('--user', help='Database user')
    parser.add_argument('--password', help='Database password')
    parser.add_argument('--skip-fix', action='store_true', help='Skip fixing invalid molecule_status values')
    args = parser.parse_args()
    
    # Get password securely if not provided
    if args.password is None and all([args.host, args.port, args.dbname, args.user]):
        try:
            args.password = getpass.getpass("Database password: ")
        except (EOFError, Exception) as e:
            logger.warning(f"Error getting password: {e}")
            # Prompt for password manually
            print("Error reading password input. Please provide password as command-line argument.")
            print("Example: python fix_consolidated_migration.py --host localhost --password yourpassword")
            sys.exit(1)
    
    # Get database connection
    conn = get_db_connection(
        host=args.host,
        port=args.port,
        dbname=args.dbname,
        user=args.user,
        password=args.password
    )
    
    try:
        # Step 1: Fix invalid molecule_status values
        if not args.skip_fix:
            print("\nStep 1: Fixing invalid molecule_status values...")
            if not fix_molecule_status(conn):
                print("Failed to fix invalid molecule_status values. Cannot continue.")
                sys.exit(1)
        else:
            print("\nSkipping fix of invalid molecule_status values as requested.")
        
        # Step 2: Apply consolidated molecule constraints and indexes
        print("\nStep 2: Applying consolidated molecule constraints and indexes...")
        if not apply_migration(conn, "migrations/026_consolidated_molecule_constraints_indexes.sql"):
            print("Failed to apply consolidated molecule constraints and indexes. Cannot continue.")
            sys.exit(1)
        
        # Step 3: Apply scientific data audit table
        print("\nStep 3: Applying scientific data audit table...")
        if not apply_migration(conn, "migrations/027_create_scientific_data_audit.sql"):
            print("Failed to apply scientific data audit table. Cannot continue.")
            sys.exit(1)
        
        print("\nAll migrations successfully applied!")
        print("The consolidated molecules implementation is ready for the next steps.")
        
    finally:
        conn.close()
        print("\nDatabase connection closed.")

if __name__ == "__main__":
    main()
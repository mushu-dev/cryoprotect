#!/usr/bin/env python3
"""
Run table consolidation migration with data integrity tests.

This script executes the SQL migration to consolidate redundant tables
and verifies data integrity before and after the migration.
"""

import os
import sys
import argparse
import psycopg2
import subprocess
from datetime import datetime
import unittest
import importlib.util
import uuid

# Add the tests directory to the path to import db_integrity_test
sys.path.append(os.path.join(os.path.dirname(__file__), 'tests'))

def load_test_module(module_path):
    """Dynamically load a test module given its path."""
    spec = importlib.util.spec_from_file_location("db_integrity_test", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def run_pre_migration_tests():
    """Run integrity tests before migration."""
    print(f"Running pre-migration integrity tests at {datetime.now().isoformat()}")
    
    # Load test module
    db_test_module = load_test_module("tests/db_integrity_test.py")
    
    # Run the tests
    test_suite = unittest.TestLoader().loadTestsFromModule(db_test_module)
    test_runner = unittest.TextTestRunner(verbosity=2)
    result = test_runner.run(test_suite)
    
    if not result.wasSuccessful():
        print("Pre-migration tests failed. Aborting migration.")
        sys.exit(1)
    
    print("Pre-migration tests passed.")

def run_migration(db_params, migration_path):
    """Run the SQL migration script."""
    print(f"Running migration from {migration_path} at {datetime.now().isoformat()}")
    
    # Connect to the database
    conn = psycopg2.connect(**db_params)
    conn.autocommit = False
    
    try:
        # Create a cursor
        cursor = conn.cursor()
        
        # Read the migration script
        with open(migration_path, 'r') as f:
            migration_sql = f.read()
        
        # Execute the migration
        cursor.execute(migration_sql)
        
        # Check if the operation succeeded
        cursor.execute("SELECT EXISTS(SELECT 1 FROM migrations WHERE name = '023_consolidate_calculation_tables')")
        migration_exists = cursor.fetchone()[0]
        
        if migration_exists:
            print("Migration successfully applied and recorded in the migrations table.")
        else:
            print("Warning: Migration applied but not recorded in migrations table.")
        
        # Commit the transaction
        conn.commit()
        print("Migration committed successfully.")
    
    except psycopg2.Error as e:
        # Roll back the transaction
        conn.rollback()
        print(f"Error during migration: {e}")
        sys.exit(1)
    
    finally:
        # Close the connection
        if conn:
            conn.close()

def run_post_migration_tests():
    """Run integrity tests after migration."""
    print(f"Running post-migration integrity tests at {datetime.now().isoformat()}")
    
    # Load test module
    db_test_module = load_test_module("tests/db_integrity_test.py")
    
    # Run the tests
    test_suite = unittest.TestLoader().loadTestsFromModule(db_test_module)
    test_runner = unittest.TextTestRunner(verbosity=2)
    result = test_runner.run(test_suite)
    
    if not result.wasSuccessful():
        print("Post-migration tests failed. Migration may have introduced issues.")
        sys.exit(1)
    
    print("Post-migration tests passed. Migration completed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Run table consolidation migration with integrity tests")
    parser.add_argument("--host", default=os.environ.get("DB_HOST", "aws-0-us-east-1.pooler.supabase.com"), 
                        help="Database host")
    parser.add_argument("--port", default=os.environ.get("DB_PORT", "5432"), 
                        help="Database port")
    parser.add_argument("--dbname", default=os.environ.get("DB_NAME", "postgres"), 
                        help="Database name")
    parser.add_argument("--user", default=os.environ.get("DB_USER", "postgres.tsdlmynydfuypiugmkev"), 
                        help="Database user")
    parser.add_argument("--password", default=os.environ.get("DB_PASSWORD"), 
                        help="Database password")
    parser.add_argument("--migration", default="migrations/023_consolidate_calculation_tables.sql", 
                        help="Path to migration SQL file")
    parser.add_argument("--skip-tests", action="store_true", 
                        help="Skip integrity tests")
    
    args = parser.parse_args()
    
    # Prepare database connection parameters
    db_params = {
        'host': args.host,
        'port': args.port,
        'dbname': args.dbname,
        'user': args.user,
        'password': args.password,
        'sslmode': 'require'
    }
    
    # Ensure password is provided
    if not db_params['password']:
        print("Error: Database password is required. Set DB_PASSWORD environment variable or use --password.")
        sys.exit(1)
    
    # Export environment variables for test module
    os.environ["DB_HOST"] = args.host
    os.environ["DB_PORT"] = args.port
    os.environ["DB_NAME"] = args.dbname
    os.environ["DB_USER"] = args.user
    os.environ["DB_PASSWORD"] = args.password
    
    # Create a timestamped backup file name
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file = f"database_backup_{timestamp}.sql"
    
    print(f"Starting table consolidation process at {datetime.now().isoformat()}")
    
    # Step 1: Backup the database (optional)
    print(f"Creating database backup to {backup_file}")
    try:
        # Use pg_dump to create a backup
        subprocess.run([
            "pg_dump", 
            "-h", args.host, 
            "-p", args.port, 
            "-U", args.user, 
            "-d", args.dbname, 
            "-n", "public", 
            "--no-owner", 
            "--no-acl", 
            "-f", backup_file
        ], env={"PGPASSWORD": args.password}, check=True)
        print("Backup created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Warning: Failed to create backup: {e}")
        response = input("Continue without backup? (y/n): ")
        if response.lower() != 'y':
            sys.exit(1)
    
    # Step 2: Run pre-migration tests
    if not args.skip_tests:
        run_pre_migration_tests()
    
    # Step 3: Run the migration
    run_migration(db_params, args.migration)
    
    # Step 4: Run post-migration tests
    if not args.skip_tests:
        run_post_migration_tests()
    
    print(f"Table consolidation process completed successfully at {datetime.now().isoformat()}")

if __name__ == "__main__":
    main()
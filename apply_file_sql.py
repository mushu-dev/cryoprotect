#!/usr/bin/env python3
"""
apply_file_sql.py: Apply SQL from a file to the Supabase database
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
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

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

def apply_sql_file(file_path):
    """Apply SQL from a file to the database."""
    try:
        # Read the SQL file
        with open(file_path, 'r') as file:
            sql = file.read()
        
        # Get a database connection
        conn = get_db_connection()
        conn.autocommit = True
        cursor = conn.cursor()
        
        try:
            # Execute the SQL
            logger.info(f"Executing SQL from {file_path}")
            cursor.execute(sql)
            print(f"Successfully executed SQL from {file_path}")
        except Exception as e:
            logger.error(f"Error executing SQL: {e}")
            print(f"ERROR: Failed to execute SQL: {e}")
            return False
        finally:
            cursor.close()
            conn.close()
        
        return True
    
    except Exception as e:
        logger.error(f"Error: {e}")
        print(f"ERROR: {e}")
        return False

def main():
    """Main function to apply SQL from a file."""
    parser = argparse.ArgumentParser(description='Apply SQL from a file to the database.')
    parser.add_argument('file_path', help='Path to the SQL file')
    args = parser.parse_args()
    
    if not os.path.exists(args.file_path):
        logger.error(f"File not found: {args.file_path}")
        print(f"ERROR: File not found: {args.file_path}")
        sys.exit(1)
    
    success = apply_sql_file(args.file_path)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
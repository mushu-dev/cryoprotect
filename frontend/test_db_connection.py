#!/usr/bin/env python3
"""
Simple script to test database connection.
This will attempt to connect using the configured database parameters.
"""

import sys
import os
import logging
import db_utils
import db_config
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def print_connection_info():
    """Print the current connection configuration."""
    # Get connection parameters
    params = db_config.get_db_connection_params()
    
    # Print connection info (excluding password)
    print("\nCurrent Database Connection Configuration:")
    print(f"Host: {params.get('host', 'None')}")
    print(f"Port: {params.get('port', 'None')}")
    print(f"Database: {params.get('dbname', 'None')}")
    print(f"User: {params.get('user', 'None')}")
    print(f"SSL Mode: {params.get('sslmode', 'None')}")
    print(f"Application Name: {params.get('application_name', 'None')}")
    
    # Check direct DB config
    print(f"\nDirect DB Configured: {db_config.is_direct_db_configured()}")
    
    # Check Supabase config
    print(f"Supabase Configured: {db_config.is_supabase_configured()}")
    if db_config.is_supabase_configured():
        supabase_config = db_config.get_supabase_config()
        print(f"Supabase URL: {supabase_config.get('url', 'None')}")
        print(f"Supabase Project ID: {supabase_config.get('project_id', 'None')}")
    
    # Connection pool config
    pool_config = db_config.get_connection_pool_config()
    print(f"\nConnection Pool Min: {pool_config.get('minconn', 'None')}")
    print(f"Connection Pool Max: {pool_config.get('maxconn', 'None')}")
    print(f"Connection Timeout: {pool_config.get('connection_timeout', 'None')}")

def main():
    """Main function to test database connection."""
    print("Testing database connection...")
    print_connection_info()
    
    # Test the connection
    print("\nAttempting to connect to the database...")
    try:
        if db_utils.test_connection():
            print("\nSUCCESS: Database connection successful")
            
            # Try a simple query
            print("\nRunning test query...")
            result = db_utils.execute_query(
                "SELECT current_database() as db, current_user as user, version() as version",
                cursor_factory=db_utils.RealDictCursor
            )
            
            if result:
                print(f"Database: {result[0]['db']}")
                print(f"User: {result[0]['user']}")
                print(f"Version: {result[0]['version']}")
            
            return 0
        else:
            print("\nERROR: Database connection test failed")
            return 1
    except Exception as e:
        print(f"\nERROR: Database connection failed: {e}")
        print("\nTroubleshooting Tips:")
        print("1. Check if PostgreSQL is running")
        print("2. Verify username and password are correct")
        print("3. Confirm database exists")
        print("4. Check firewall settings")
        print("5. Verify connection parameters in .env file")
        return 1

if __name__ == "__main__":
    sys.exit(main())
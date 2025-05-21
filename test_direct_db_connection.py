#\!/usr/bin/env python3
"""
Direct Database Connection Test

This script tests a direct connection to the Supabase PostgreSQL database
without using any of the CryoProtect adapter or configuration classes.
This helps identify whether any issues are with the connection credentials
themselves or with the CryoProtect database abstraction layer.
"""

import os
import sys
import json
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor

# Default Supabase connection parameters from MCP
DEFAULT_DATABASE_HOST = "db.tsdlmynydfuypiugmkev.supabase.co"
DEFAULT_DATABASE_PORT = 5432
DEFAULT_DATABASE_NAME = "postgres"
DEFAULT_DATABASE_USER = "postgres"
DEFAULT_SUPABASE_URL = "https://tsdlmynydfuypiugmkev.supabase.co"
DEFAULT_SUPABASE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRzZGxteW55ZGZ1eXBpdWdta2V2Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDQ3ODAxODMsImV4cCI6MjA2MDM1NjE4M30.T6udmD-3lhlTy4bY2p0y2lX5-11yvyn425PWrlnIPLU"

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Test direct database connection')
    parser.add_argument('--host', default=DEFAULT_DATABASE_HOST,
                        help=f'Database host (default: {DEFAULT_DATABASE_HOST})')
    parser.add_argument('--port', type=int, default=DEFAULT_DATABASE_PORT,
                        help=f'Database port (default: {DEFAULT_DATABASE_PORT})')
    parser.add_argument('--database', default=DEFAULT_DATABASE_NAME,
                        help=f'Database name (default: {DEFAULT_DATABASE_NAME})')
    parser.add_argument('--user', default=DEFAULT_DATABASE_USER,
                        help=f'Database user (default: {DEFAULT_DATABASE_USER})')
    parser.add_argument('--password', help='Database password (will prompt if not provided)')
    parser.add_argument('--config', help='Path to a JSON config file with connection details')
    return parser.parse_args()

def load_config_from_file(config_path):
    """Load connection parameters from a config file."""
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        # Try to extract Supabase connection details
        db_config = config.get('database', {}).get('connection', {}).get('supabase', {})
        if not db_config:
            print("Warning: Could not find Supabase configuration in the config file.")
            return {}
        
        return {
            'host': db_config.get('host', DEFAULT_DATABASE_HOST),
            'port': db_config.get('port', DEFAULT_DATABASE_PORT),
            'database': db_config.get('database', DEFAULT_DATABASE_NAME),
            'user': db_config.get('user', DEFAULT_DATABASE_USER),
            'password': db_config.get('password', '')
        }
    except Exception as e:
        print(f"Error loading config file: {str(e)}")
        return {}

def get_connection_params():
    """Get connection parameters from command line or config file."""
    args = parse_arguments()
    
    # Start with default values
    params = {
        'host': DEFAULT_DATABASE_HOST,
        'port': DEFAULT_DATABASE_PORT,
        'database': DEFAULT_DATABASE_NAME,
        'user': DEFAULT_DATABASE_USER,
        'password': '',
        'application_name': 'CryoProtect-DirectTest'
    }
    
    # Override with config file if provided
    if args.config:
        config_params = load_config_from_file(args.config)
        params.update(config_params)
    
    # Override with command-line arguments
    if args.host:
        params['host'] = args.host
    if args.port:
        params['port'] = args.port
    if args.database:
        params['database'] = args.database
    if args.user:
        params['user'] = args.user
    if args.password:
        params['password'] = args.password
    
    # Prompt for password if not provided
    if not params['password']:
        import getpass
        params['password'] = getpass.getpass(f"Enter password for {params['user']}@{params['host']}: ")
    
    return params

def test_connection(params):
    """Test direct connection to the database."""
    print(f"Connecting to {params['host']}:{params['port']} as {params['user']}...")
    
    try:
        # Connect directly using psycopg2
        connection = psycopg2.connect(
            host=params['host'],
            port=params['port'],
            database=params['database'],
            user=params['user'],
            password=params['password'],
            application_name=params['application_name']
        )
        
        print("✓ Connection successful\!")
        
        # Test a simple query
        with connection.cursor(cursor_factory=RealDictCursor) as cursor:
            print("\nTesting basic query: SELECT current_database()")
            cursor.execute("SELECT current_database();")
            result = cursor.fetchone()
            print(f"Current database: {result['current_database']}")
            
            print("\nTesting molecules count query")
            cursor.execute("SELECT COUNT(*) FROM molecules;")
            result = cursor.fetchone()
            print(f"Number of molecules in database: {result['count']}")
            
            print("\nChecking database version")
            cursor.execute("SELECT version();")
            result = cursor.fetchone()
            version = result['version']
            print(f"Database version: {version}")
            
            print("\nListing some tables from public schema")
            cursor.execute("SELECT tablename FROM pg_tables WHERE schemaname = 'public' LIMIT 5;")
            tables = cursor.fetchall()
            print("Sample tables in database:")
            for table in tables:
                print(f"  - {table['tablename']}")
            
            print("\nChecking server timezone")
            cursor.execute("SHOW timezone;")
            result = cursor.fetchone()
            print(f"Server timezone: {result['timezone']}")
            
            print("\nTesting transaction")
            connection.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_READ_COMMITTED)
            cursor.execute("BEGIN;")
            cursor.execute("SELECT COUNT(*) FROM molecules;")
            result = cursor.fetchone()
            print(f"Transaction query result: {result['count']} molecules")
            cursor.execute("ROLLBACK;")
            print("Transaction rolled back successfully")
        
        # Close connection
        connection.close()
        print("\n✓ All tests completed successfully\!")
        return True
    
    except Exception as e:
        print(f"\n✗ Error connecting to database: {str(e)}")
        
        # Provide more detailed error information based on the exception type
        if isinstance(e, psycopg2.OperationalError):
            print("\nThis is an operational error, which could indicate:")
            print("- Incorrect connection parameters (host, port, user, password)")
            print("- Network connectivity issues")
            print("- Firewall blocking the connection")
            print("- Database server not running")
            
            # Check if it's an authentication error
            if "password authentication failed" in str(e):
                print("\nAuthentication failed. Possible reasons:")
                print("- Incorrect password")
                print("- User does not have permission to connect")
                print("- User does not exist")
                
            # Check if it's a connection refused error
            elif "could not connect to server" in str(e):
                print("\nCould not connect to server. Possible reasons:")
                print("- Database server is not running")
                print("- Network connectivity issues")
                print("- Firewall blocking the connection")
                print("- Incorrect host or port")
                
            # Check if it's a database not found error
            elif "database" in str(e) and "does not exist" in str(e):
                print("\nDatabase does not exist. Possible reasons:")
                print("- Incorrect database name")
                print("- Database has not been created")
        
        return False

def main():
    """Main function."""
    print("=" * 60)
    print("Direct Database Connection Test")
    print("=" * 60)
    
    # Get connection parameters
    params = get_connection_params()
    
    # Display parameters (except password)
    print("Using connection parameters:")
    for key, value in params.items():
        if key \!= 'password':
            print(f"  - {key}: {value}")
    print("-" * 60)
    
    # Test connection
    test_connection(params)
    
    print("=" * 60)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Test script to check PostgreSQL connection with different credentials.
"""

import os
import sys
import psycopg2
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

def test_connection(host, port, dbname, user, password):
    """Test a PostgreSQL connection with the given parameters."""
    try:
        print(f"Trying to connect to PostgreSQL: host={host}, port={port}, dbname={dbname}, user={user}, password={'*' * len(password)}")
        conn = psycopg2.connect(
            host=host,
            port=port,
            dbname=dbname,
            user=user,
            password=password
        )
        print("Connection successful!")
        
        # Test a simple query
        with conn.cursor() as cursor:
            cursor.execute("SELECT version();")
            version = cursor.fetchone()
            print(f"PostgreSQL version: {version[0]}")
        
        conn.close()
        return True
    except Exception as e:
        print(f"Connection failed: {str(e)}")
        return False

def main():
    # Test with environment variables from .env
    print("\n=== Testing with SUPABASE_DB_* environment variables ===")
    supabase_host = os.getenv('SUPABASE_DB_HOST', 'localhost')
    supabase_port = int(os.getenv('SUPABASE_DB_PORT', 5432))
    supabase_dbname = os.getenv('SUPABASE_DB_NAME', 'postgres')
    supabase_user = os.getenv('SUPABASE_DB_USER', 'postgres')
    supabase_password = os.getenv('SUPABASE_DB_PASSWORD', 'postgres')
    
    test_connection(
        host=supabase_host,
        port=supabase_port,
        dbname=supabase_dbname,
        user=supabase_user,
        password=supabase_password
    )
    
    # Test with localhost and environment variables
    print("\n=== Testing with localhost and SUPABASE_DB_* credentials ===")
    test_connection(
        host='localhost',
        port=supabase_port,
        dbname=supabase_dbname,
        user=supabase_user,
        password=supabase_password
    )
    
    # Test with localhost and default credentials
    print("\n=== Testing with localhost and default credentials ===")
    test_connection(
        host='localhost',
        port=5432,
        dbname='postgres',
        user='postgres',
        password='postgres'
    )
    
    # Test with a different password
    print("\n=== Testing with localhost and alternative password ===")
    test_connection(
        host='localhost',
        port=5432,
        dbname='postgres',
        user='postgres',
        password='cryoprotect'
    )
    
    # Test with empty password
    print("\n=== Testing with localhost and empty password ===")
    test_connection(
        host='localhost',
        port=5432,
        dbname='postgres',
        user='postgres',
        password=''
    )

if __name__ == "__main__":
    main()
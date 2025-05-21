#!/usr/bin/env python3
"""
Direct connection test script bypassing the configuration system.
"""

import sys
import psycopg2
from psycopg2.extras import RealDictCursor

# Direct connection parameters - update these values
DB_HOST = "db.tsdlmynydfuypiugmkev.supabase.co"
DB_PORT = "5432"
DB_NAME = "postgres"
DB_USER = "postgres"
DB_PASSWORD = "LDHt$rkaM&Gmf3X@LQ37"

def test_direct_connection():
    """Test a direct database connection."""
    print("Testing direct database connection...")
    
    connection_params = {
        "host": DB_HOST,
        "port": DB_PORT, 
        "dbname": DB_NAME,
        "user": DB_USER,
        "password": DB_PASSWORD,
        "sslmode": "require"
    }
    
    print(f"Connecting to: {DB_HOST}:{DB_PORT}/{DB_NAME} as {DB_USER}")
    print("SSL Mode: require")
    
    try:
        # Connect to the database
        conn = psycopg2.connect(**connection_params)
        
        # Create a cursor
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # Execute a test query
        cursor.execute("SELECT current_database() as db, current_user as user, version() as version")
        
        # Fetch the results
        result = cursor.fetchone()
        
        print("\nSUCCESS: Database connection successful")
        print(f"Database: {result['db']}")
        print(f"User: {result['user']}")
        print(f"Version: {result['version']}")
        
        # Close the cursor and connection
        cursor.close()
        conn.close()
        
        return 0
    except Exception as e:
        print(f"\nERROR: Database connection failed: {e}")
        print("\nTroubleshooting Tips:")
        print("1. Verify hostname is correct")
        print("2. Check if username and password are correct")
        print("3. Confirm database name exists")
        print("4. Ensure SSL/TLS is properly configured")
        print("5. Check network connectivity to the database server")
        return 1

if __name__ == "__main__":
    sys.exit(test_direct_connection())
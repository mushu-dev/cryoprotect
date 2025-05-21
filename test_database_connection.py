#!/usr/bin/env python3
"""
Test database connection directly.

This script tests the database connection directly, bypassing the Flask app.
"""

import sys
import json
from pathlib import Path
from psycopg2.extras import RealDictCursor
from database import db, db_service_role, db_public

def load_config():
    """Load database configuration."""
    config_path = Path('config/config.json')
    if not config_path.exists():
        print(f"Error: Configuration file not found at {config_path}")
        return None
    
    with open(config_path, 'r') as f:
        config_data = json.load(f)
    
    if 'database' in config_data and 'connection' in config_data['database']:
        connection_config = config_data['database']['connection']
        if 'supabase' in connection_config:
            config = connection_config['supabase'].copy()
            # Add pooling settings if available
            if 'pooling' in config_data['database']:
                pooling = config_data['database']['pooling']
                if 'min_connections' in pooling:
                    config['min_connections'] = pooling['min_connections']
                if 'max_connections' in pooling:
                    config['max_connections'] = pooling['max_connections']
            return config
    
    print("Error: No valid database configuration found")
    return None

def test_database_connections():
    """Test database connections."""
    print("Testing database connections...")
    
    # Load configuration
    config = load_config()
    if not config:
        print("Failed to load configuration")
        return False
    
    # Initialize database modules
    print("Initializing database modules...")
    db.init_connection_pool(config=config)
    
    # Create service role config
    service_role_config = config.copy()
    service_role_config['options'] = "-c role=service_role"
    db_service_role.init_connection_pool(config=service_role_config)
    
    # Create public config
    public_config = config.copy()
    db_public.init_connection_pool(config=public_config)
    
    # Test main db module
    print("\nTesting main db module...")
    try:
        conn = db.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db, current_user as user")
            result = cursor.fetchone()
            print(f"✓ Main db connection successful")
            print(f"  Database: {result['db']}")
            print(f"  User: {result['user']}")
            
            cursor.execute("SELECT count(*) FROM molecules")
            result = cursor.fetchone()
            print(f"✓ Found {result['count']} molecules in the database")
        db.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with main db connection: {e}")
        return False
    
    # Test service role module
    print("\nTesting service role module...")
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db, current_user as user")
            result = cursor.fetchone()
            print(f"✓ Service role connection successful")
            print(f"  Database: {result['db']}")
            print(f"  User: {result['user']}")
            
            cursor.execute("SELECT count(*) FROM molecules")
            result = cursor.fetchone()
            print(f"✓ Found {result['count']} molecules in the database (service role)")
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with service role connection: {e}")
        return False
    
    # Test public module
    print("\nTesting public module...")
    try:
        conn = db_public.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db, current_user as user")
            result = cursor.fetchone()
            print(f"✓ Public connection successful")
            print(f"  Database: {result['db']}")
            print(f"  User: {result['user']}")
            
            cursor.execute("SELECT count(*) FROM molecules")
            result = cursor.fetchone()
            print(f"✓ Found {result['count']} molecules in the database (public)")
        db_public.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with public connection: {e}")
        return False
    
    print("\nAll database connections successful!")
    return True

if __name__ == "__main__":
    success = test_database_connections()
    sys.exit(0 if success else 1)
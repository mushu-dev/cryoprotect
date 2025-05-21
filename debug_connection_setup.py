#!/usr/bin/env python3
"""
Debug the database connection setup specifically.

This script focuses exclusively on fixing the database.db module connection issue.
"""

import os
import sys
import logging
import traceback
import json
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def fix_connection_config():
    """Fix the connection configuration."""
    print("Fixing database connection configuration...")
    
    # Ensure environment variables are loaded from .env
    load_dotenv()
    
    # Create connection configuration from the environment
    config = {
        'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
        'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
        'application_name': 'CryoProtect',
        'min_connections': 1,
        'max_connections': 10
    }
    
    # Print out the configuration (without password)
    safe_config = config.copy()
    safe_config['password'] = '******'
    print(f"Connection configuration: {json.dumps(safe_config, indent=2)}")
    
    print("Importing database modules...")
    try:
        # Import necessary database modules
        from database import db, db_service_role, db_public
        
        # Try to initialize each module's connection pool
        print("\nInitializing db module connection pool...")
        db_success = db.init_connection_pool(config=config)
        print(f"Success: {db_success}")
        
        print("\nInitializing db_service_role module connection pool...")
        service_role_success = db_service_role.init_connection_pool(config=config)
        print(f"Success: {service_role_success}")
        
        print("\nInitializing db_public module connection pool...")
        public_success = db_public.init_connection_pool(config=config)
        print(f"Success: {public_success}")
        
        # Test connection from each module
        print("\nTesting db module connection...")
        try:
            conn = db.get_connection()
            if conn:
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1 as test")
                    print("✓ db connection successful!")
                db.release_connection(conn)
            else:
                print("✗ db.get_connection() returned None")
        except Exception as e:
            print(f"✗ Error testing db connection: {e}")
            traceback.print_exc()
        
        print("\nTesting db_service_role module connection...")
        try:
            conn = db_service_role.get_connection()
            if conn:
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1 as test")
                    print("✓ db_service_role connection successful!")
                db_service_role.release_connection(conn)
            else:
                print("✗ db_service_role.get_connection() returned None")
        except Exception as e:
            print(f"✗ Error testing db_service_role connection: {e}")
            traceback.print_exc()
        
        print("\nTesting db_public module connection...")
        try:
            conn = db_public.get_connection()
            if conn:
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1 as test")
                    print("✓ db_public connection successful!")
                db_public.release_connection(conn)
            else:
                print("✗ db_public.get_connection() returned None")
        except Exception as e:
            print(f"✗ Error testing db_public connection: {e}")
            traceback.print_exc()
        
    except Exception as e:
        print(f"Error during database module initialization: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    fix_connection_config()
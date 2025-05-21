#!/usr/bin/env python3
"""
Sync database modules configuration.

This script fixes the service_role and public database modules to use the same
configuration as the main db module. It also creates a flattened db_config.json
file for backward compatibility.
"""

import os
import sys
import json
import logging
import traceback
import shutil
from pathlib import Path
from psycopg2.extras import RealDictCursor

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_db_config_json(config_data):
    """Create a flattened db_config.json file for backward compatibility."""
    print("\nCreating db_config.json for backward compatibility...")
    
    # Extract the database connection parameters from the config.json format
    if 'database' in config_data and 'connection' in config_data['database']:
        connection_config = config_data['database']['connection']
        if 'supabase' in connection_config:
            # Create a flattened config
            db_config = connection_config['supabase'].copy()
            
            # Add pooling settings if available
            if 'pooling' in config_data['database']:
                pooling = config_data['database']['pooling']
                if 'min_connections' in pooling:
                    db_config['min_connections'] = pooling['min_connections']
                if 'max_connections' in pooling:
                    db_config['max_connections'] = pooling['max_connections']
            
            # Ensure config directory exists
            config_dir = Path('config')
            config_dir.mkdir(exist_ok=True)
            
            # Backup existing db_config.json if it exists
            db_config_path = config_dir / 'db_config.json'
            if db_config_path.exists():
                backup_path = db_config_path.with_suffix('.json.bak')
                shutil.copy2(db_config_path, backup_path)
                print(f"✓ Backed up existing db_config.json to {backup_path}")
            
            # Write the flattened config to db_config.json
            with open(db_config_path, 'w') as f:
                json.dump(db_config, f, indent=2)
            
            print(f"✓ Created db_config.json at {db_config_path}")
            return True
        else:
            print("Error: No supabase configuration found in config.json")
    else:
        print("Error: No database connection configuration found in config.json")
    
    return False

def sync_modules():
    """Synchronize database modules to use the same configuration."""
    print("Syncing database modules to use the same configuration...")
    
    # Load the configuration from the main config.json file
    config_path = Path('config/config.json')
    if not config_path.exists():
        print(f"Error: Configuration file not found at {config_path}")
        return False
    
    # Import all database modules
    try:
        from database import db, db_service_role, db_public
        print("✓ Successfully imported all database modules")
    except ImportError as e:
        print(f"Error importing database modules: {e}")
        return False
    
    # Load the configuration
    with open(config_path, 'r') as f:
        config_data = json.load(f)
    
    # Create db_config.json for backward compatibility
    create_db_config_json(config_data)
    
    # Extract the database connection parameters from the config.json format
    if 'database' in config_data and 'connection' in config_data['database']:
        connection_config = config_data['database']['connection']
        if 'supabase' in connection_config:
            config = connection_config['supabase']
            # Add pooling settings if available
            if 'pooling' in config_data['database']:
                pooling = config_data['database']['pooling']
                if 'min_connections' in pooling:
                    config['min_connections'] = pooling['min_connections']
                if 'max_connections' in pooling:
                    config['max_connections'] = pooling['max_connections']
        else:
            print("Error: No supabase configuration found in config.json")
            return False
    else:
        print("Error: No database connection configuration found in config.json")
        return False
    
    # Display configuration (without password)
    display_config = config.copy()
    if 'password' in display_config:
        display_config['password'] = '********'
    print(f"Configuration: {json.dumps(display_config, indent=2)}")
    
    # Initialize db module
    print("\nInitializing database modules with the same configuration...")
    
    # Stop all current pools
    print("Closing any existing connection pools...")
    db.close_all_connections()
    db_service_role.close_all_connections()
    db_public.close_all_connections()
    
    # Initialize each module with the same configuration
    success_db = db.init_connection_pool(config=config)
    success_service = db_service_role.init_connection_pool(config=config)
    success_public = db_public.init_connection_pool(config=config)
    
    print(f"✓ Main db module initialized: {success_db}")
    print(f"✓ Service role module initialized: {success_service}")
    print(f"✓ Public module initialized: {success_public}")
    
    # Test the connections
    print("\nTesting database connections...")
    
    # Test main db module
    try:
        conn = db.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db")
            result = cursor.fetchone()
            db_name = result['db']
            print(f"✓ Main db connection successful (database: {db_name})")
        db.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with main db connection: {e}")
        traceback.print_exc()
    
    # Test service role module
    try:
        conn = db_service_role.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db")
            result = cursor.fetchone()
            db_name = result['db']
            print(f"✓ Service role connection successful (database: {db_name})")
        db_service_role.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with service role connection: {e}")
        traceback.print_exc()
    
    # Test public module
    try:
        conn = db_public.get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database() as db")
            result = cursor.fetchone()
            db_name = result['db']
            print(f"✓ Public connection successful (database: {db_name})")
        db_public.release_connection(conn)
    except Exception as e:
        print(f"✗ Error with public connection: {e}")
        traceback.print_exc()
    
    print("\nDatabase modules synchronized!")
    return True

def fix_config_usage():
    """Fix how db_service_role and db_public use the configuration."""
    # Path to modules
    service_role_path = Path('database/db_service_role.py')
    public_path = Path('database/db_public.py')
    
    # Backup files
    print("Backing up database module files...")
    if service_role_path.exists():
        with open(service_role_path, 'r') as f:
            service_role_content = f.read()
        with open(str(service_role_path) + '.bak', 'w') as f:
            f.write(service_role_content)
        print(f"✓ Backed up {service_role_path}")
    
    if public_path.exists():
        with open(public_path, 'r') as f:
            public_content = f.read()
        with open(str(public_path) + '.bak', 'w') as f:
            f.write(public_content)
        print(f"✓ Backed up {public_path}")
    
    # Fix db_service_role.py if it exists
    if service_role_path.exists():
        print(f"\nFixing configuration usage in {service_role_path}...")
        
        if "missing_keys = [key for key in required_keys if not config.get(key)]" in service_role_content:
            # Update the validation code
            new_content = service_role_content.replace(
                "missing_keys = [key for key in required_keys if not config.get(key)]",
                "missing_keys = [key for key in required_keys if key not in config or config[key] is None or config[key] == '']"
            )
            
            # Write the updated content
            with open(service_role_path, 'w') as f:
                f.write(new_content)
            print(f"✓ Fixed config validation in {service_role_path}")
        else:
            print(f"Config validation already fixed in {service_role_path}")
    
    # Fix db_public.py if it exists
    if public_path.exists():
        print(f"\nFixing configuration usage in {public_path}...")
        
        if "missing_keys = [key for key in required_keys if not config.get(key)]" in public_content:
            # Update the validation code
            new_content = public_content.replace(
                "missing_keys = [key for key in required_keys if not config.get(key)]",
                "missing_keys = [key for key in required_keys if key not in config or config[key] is None or config[key] == '']"
            )
            
            # Write the updated content
            with open(public_path, 'w') as f:
                f.write(new_content)
            print(f"✓ Fixed config validation in {public_path}")
        else:
            print(f"Config validation already fixed in {public_path}")
    
    return True

if __name__ == "__main__":
    print("=" * 60)
    print("Database Modules Sync")
    print("=" * 60)
    
    # Fix how config is used in each module
    fix_config_usage()
    
    # Sync the modules to use the same configuration
    sync_modules()
    
    print("\nAll fixes applied successfully!")
    print("\nRun the debug_connection_flow.py script to verify everything works:")
    print("python debug_connection_flow.py --config-path=config/config.json")
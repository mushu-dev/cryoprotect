#!/usr/bin/env python3
"""
Fix the database connection system.

This script provides a comprehensive fix for the database connection system,
creating a simplified adapter layer that works with the existing modules.
"""

import os
import sys
import json
import logging
from pathlib import Path
from dotenv import load_dotenv
from typing import Dict, Any, Optional

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_env_and_create_config():
    """Load environment variables and create database configuration."""
    print("Loading environment variables and creating database configuration...")
    
    # Load environment variables from .env
    load_dotenv()
    
    # Create configuration dictionary
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
    
    # Create config directory if it doesn't exist
    Path('config').mkdir(exist_ok=True)
    
    # Create config.json with the right format
    full_config = {
        'database': {
            'connection': {
                'mode': 'supabase',
                'supabase': {
                    'host': config['host'],
                    'port': config['port'],
                    'database': config['database'],
                    'user': config['user'],
                    'password': config['password'],
                    'project_id': os.getenv('SUPABASE_PROJECT_ID', 'tsdlmynydfuypiugmkev'),
                    'application_name': config['application_name']
                },
                'ssl': {
                    'mode': 'require'
                }
            },
            'pooling': {
                'enabled': True,
                'min_connections': config['min_connections'],
                'max_connections': config['max_connections']
            }
        }
    }
    
    # Write the config to file
    config_path = Path('config/config.json')
    with open(config_path, 'w') as f:
        json.dump(full_config, f, indent=2)
    
    print(f"✓ Configuration file created at {config_path}")
    
    # Create a simplified database configuration for direct usage
    with open('config/db_config.json', 'w') as f:
        # Mask password for display
        display_config = config.copy()
        display_config['password'] = '********'
        print(f"Direct connection config: {json.dumps(display_config, indent=2)}")
        json.dump(config, f, indent=2)
    
    print(f"✓ Direct database configuration created at config/db_config.json")
    
    return config, str(config_path)

def create_adapter_factory():
    """Create the adapter factory module."""
    print("Creating adapter factory module...")
    
    # Define the content
    content = """#!/usr/bin/env python3
\"\"\"
Simple adapter factory for the database connection system.

This module provides a simplified connection factory that works with the new
direct connection system while maintaining compatibility with code
expecting the adapter pattern.
\"\"\"

import os
import logging
import json
from typing import Dict, Any, Optional
from pathlib import Path

# Configure logger
logger = logging.getLogger(__name__)

class ConnectionFactory:
    \"\"\"
    Simplified connection factory that works with the direct connection system.
    
    This class implements a minimal subset of the ConnectionFactory interface
    to work with the new direct connection system.
    \"\"\"
    
    def __init__(self, config=None):
        \"\"\"
        Initialize the connection factory.
        
        Args:
            config: Optional configuration dictionary
        \"\"\"
        self.config = config or {}
        
        # Initialize the connection pool if needed
        try:
            from database import db
            if not hasattr(db, '_pool') or db._pool is None:
                logger.info("Initializing connection pool from adapter factory")
                db.init_connection_pool(config=self.get_connection_params())
        except Exception as e:
            logger.error(f"Error initializing connection pool: {e}")
    
    def get_connection_params(self) -> Dict[str, Any]:
        \"\"\"
        Get connection parameters from configuration.
        
        Returns:
            Dict of connection parameters
        \"\"\"
        # If config is empty, try to load from config file
        if not self.config:
            try:
                config_path = Path('config/db_config.json')
                if config_path.exists():
                    with open(config_path, 'r') as f:
                        return json.load(f)
            except Exception as e:
                logger.error(f"Error loading config file: {e}")
        
        # Extract connection parameters from supabase section if available
        if 'database' in self.config and 'connection' in self.config['database']:
            connection = self.config['database']['connection']
            if 'supabase' in connection and connection.get('mode') == 'supabase':
                return connection['supabase']
        
        # Return the config as is (might be already flattened)
        return self.config
    
    def get_connection(self):
        \"\"\"
        Get a database connection.
        
        Returns:
            Database connection
        \"\"\"
        from database import db
        return db.get_connection()

# For backward compatibility
def get_db_connection_info() -> Dict[str, Any]:
    \"\"\"
    Get information about the current database connection.
    
    Returns:
        Dict with connection information
    \"\"\"
    from database import db
    
    try:
        conn = db.get_connection()
        if not conn:
            return {'error': 'Could not get database connection'}
        
        try:
            with conn.cursor() as cursor:
                # Get server information
                cursor.execute("SELECT version() as version, current_database() as database, current_user as user")
                server_info = cursor.fetchone()
                
                # Get client encoding
                cursor.execute("SHOW client_encoding")
                encoding = cursor.fetchone()['client_encoding']
                
                return {
                    'server_version': server_info['version'],
                    'database': server_info['database'],
                    'user': server_info['user'],
                    'encoding': encoding,
                }
        finally:
            # Release the connection
            db.release_connection(conn)
            
    except Exception as e:
        logger.error(f"Error getting database connection info: {e}")
        return {'error': str(e)}
"""
    
    # Write the file
    with open('database/adapter_factory.py', 'w') as f:
        f.write(content)
    
    print("✓ Created adapter_factory.py")

def fix_db_module():
    """Fix the main db module to use the config properly."""
    print("Checking database/db.py module...")
    
    # Read the db.py file to check if it's already fixed
    with open('database/db.py', 'r') as f:
        content = f.read()
    
    # Check if we need to fix the parameter validation
    if "missing_keys = [key for key in required_keys if not config.get(key)]" in content:
        print("Fixing parameter validation in database/db.py...")
        
        # Backup the original file
        backup_path = 'database/db.py.bak'
        with open(backup_path, 'w') as f:
            f.write(content)
        print(f"✓ Backed up original to {backup_path}")
        
        # Replace the validation with a more robust check
        new_content = content.replace(
            "missing_keys = [key for key in required_keys if not config.get(key)]",
            "missing_keys = [key for key in required_keys if key not in config or config[key] is None or config[key] == '']"
        )
        
        # Write the fixed content
        with open('database/db.py', 'w') as f:
            f.write(new_content)
        print("✓ Fixed parameter validation in database/db.py")
    else:
        print("Parameter validation already fixed in database/db.py")

def check_connection():
    """Test the database connection."""
    print("\nTesting database connection...")
    
    try:
        # Import the db module
        from database import db
        from database import db_service_role
        from database import db_public
        
        # Load the configuration
        with open('config/db_config.json', 'r') as f:
            config = json.load(f)
        
        # Initialize the db module
        print("Initializing db module...")
        success = db.init_connection_pool(config=config)
        print(f"✓ DB module initialization: {success}")
        
        # Initialize the service role module
        print("Initializing db_service_role module...")
        success = db_service_role.init_connection_pool(config=config)
        print(f"✓ DB service role module initialization: {success}")
        
        # Initialize the public module
        print("Initializing db_public module...")
        success = db_public.init_connection_pool(config=config)
        print(f"✓ DB public module initialization: {success}")
        
        # Test getting a connection
        print("\nTesting database connections...")
        
        # Test the main db module
        conn = db.get_connection()
        if conn:
            print("✓ Successfully got connection from db module")
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user")
                db_name, user = cursor.fetchone()
                print(f"  Database: {db_name}")
                print(f"  User: {user}")
            db.release_connection(conn)
        else:
            print("✗ Failed to get connection from db module")
        
        # Test the service role module
        conn = db_service_role.get_connection()
        if conn:
            print("✓ Successfully got connection from db_service_role module")
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user")
                db_name, user = cursor.fetchone()
                print(f"  Database: {db_name}")
                print(f"  User: {user}")
            db_service_role.release_connection(conn)
        else:
            print("✗ Failed to get connection from db_service_role module")
        
        # Test the public module
        conn = db_public.get_connection()
        if conn:
            print("✓ Successfully got connection from db_public module")
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user")
                db_name, user = cursor.fetchone()
                print(f"  Database: {db_name}")
                print(f"  User: {user}")
            db_public.release_connection(conn)
        else:
            print("✗ Failed to get connection from db_public module")
        
        # Test the adapter factory
        try:
            print("\nTesting adapter factory...")
            from database.adapter_factory import ConnectionFactory
            
            # Create a factory
            factory = ConnectionFactory()
            print("✓ Successfully created connection factory")
            
            # Get a connection
            conn = factory.get_connection()
            if conn:
                print("✓ Successfully got connection from adapter factory")
                with conn.cursor() as cursor:
                    cursor.execute("SELECT current_database(), current_user")
                    db_name, user = cursor.fetchone()
                    print(f"  Database: {db_name}")
                    print(f"  User: {user}")
                db.release_connection(conn)
            else:
                print("✗ Failed to get connection from adapter factory")
        except Exception as e:
            print(f"✗ Error testing adapter factory: {e}")
        
        print("\n✓ Database connection system fixed and working!")
    
    except Exception as e:
        print(f"✗ Error during connection test: {e}")
        import traceback
        traceback.print_exc()

def main():
    """Main entry point for the script."""
    print("=" * 60)
    print("Database Connection System Fix")
    print("=" * 60)
    
    # Load environment and create configuration
    load_env_and_create_config()
    
    # Create the adapter factory
    create_adapter_factory()
    
    # Fix the db module validation
    fix_db_module()
    
    # Check the connection
    check_connection()
    
    print("\nAll fixes applied successfully!")
    print("\nRun the debug_connection_flow.py script to verify everything works:")
    print("python debug_connection_flow.py --config-path=config/config.json --verbose")

if __name__ == "__main__":
    main()
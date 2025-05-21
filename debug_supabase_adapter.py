#\!/usr/bin/env python3
"""
Debug Supabase Adapter

This script tests and debugs specifically the Supabase adapter in the CryoProtect 
database abstraction layer. It includes detailed logging of each step in the 
connection process and helps identify issues with the Supabase connection.
"""

import os
import sys
import json
import logging
import argparse
import importlib
import traceback
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('debug_supabase')

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Debug Supabase adapter')
    parser.add_argument('--config', '-c', dest='config_path', 
                        help='Path to the configuration file (overrides environment variable)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose output with detailed tracing')
    parser.add_argument('--direct', '-d', action='store_true',
                        help='Test direct connection using psycopg2 instead of the adapter')
    return parser.parse_args()

def find_config_path(args):
    """Find the configuration path."""
    # Check command-line arguments
    if args.config_path:
        if os.path.exists(args.config_path):
            return args.config_path
        else:
            logger.error(f"Configuration file not found at {args.config_path}")
            sys.exit(1)
    
    # Check environment variable
    env_config_path = os.environ.get("CRYOPROTECT_CONFIG_PATH")
    if env_config_path and os.path.exists(env_config_path):
        return env_config_path
    
    # Check default locations
    default_paths = [
        os.path.join(os.path.dirname(__file__), "config", "test_config.json"),
        os.path.join(os.path.dirname(__file__), "config.json"),
        os.path.join(os.path.dirname(__file__), ".config", "config.json")
    ]
    
    for path in default_paths:
        if os.path.exists(path):
            return path
    
    # If we get here, no configuration file was found
    logger.error("No configuration file found. Please specify a path with --config or set CRYOPROTECT_CONFIG_PATH.")
    sys.exit(1)

def get_supabase_config(config_path):
    """Get Supabase configuration from config file."""
    try:
        with open(config_path, 'r') as f:
            config_data = json.load(f)
        
        supabase_config = config_data.get('database', {}).get('connection', {}).get('supabase', {})
        if not supabase_config:
            logger.error("No Supabase configuration found in config file")
            return None
        
        # Validate essential fields
        required_fields = ['host', 'port', 'database', 'user', 'password']
        missing_fields = [field for field in required_fields if field not in supabase_config]
        
        if missing_fields:
            logger.error(f"Missing required fields in Supabase configuration: {', '.join(missing_fields)}")
            return None
        
        return supabase_config
    except Exception as e:
        logger.error(f"Error loading Supabase configuration: {str(e)}")
        return None

def test_direct_connection(config):
    """Test direct connection using psycopg2."""
    logger.info("Testing direct connection using psycopg2")
    
    try:
        import psycopg2
        from psycopg2.extras import RealDictCursor
        
        logger.info(f"Connecting to {config['host']}:{config['port']} as {config['user']}...")
        
        connection = psycopg2.connect(
            host=config['host'],
            port=config['port'],
            database=config['database'],
            user=config['user'],
            password=config['password'],
            application_name="CryoProtect-SupabaseTest"
        )
        
        logger.info("✅ Connection successful")
        
        # Test a simple query
        with connection.cursor(cursor_factory=RealDictCursor) as cursor:
            logger.info("Executing test query: SELECT COUNT(*) FROM molecules")
            cursor.execute("SELECT COUNT(*) FROM molecules;")
            result = cursor.fetchone()
            logger.info(f"Query result: {result['count']} molecules found")
        
        # Close connection
        connection.close()
        logger.info("Connection closed")
        
        return True
    except Exception as e:
        logger.error(f"Error in direct connection: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def test_supabase_adapter(config):
    """Test the Supabase adapter."""
    logger.info("Testing Supabase adapter")
    
    # Add project root to path
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    
    try:
        # Import adapter class
        from database.adapters.supabase import SupabaseAdapter
        
        # Create adapter
        logger.info("Creating Supabase adapter")
        adapter = SupabaseAdapter(config)
        
        # Test connection
        logger.info("Testing connection")
        connection_result = adapter.connect()
        
        if connection_result:
            logger.info("✅ Connection successful")
            
            # Test query
            logger.info("Executing test query: SELECT COUNT(*) FROM molecules")
            result = adapter.execute_query("SELECT COUNT(*) FROM molecules;")
            
            if result:
                logger.info(f"Query result: {result[0]['count']} molecules found")
            else:
                logger.warning("Query returned no results")
            
            # Test transaction
            logger.info("Testing transaction")
            with adapter.transaction():
                transaction_result = adapter.execute_query("SELECT COUNT(*) FROM molecules;")
                logger.info(f"Transaction query result: {transaction_result[0]['count']} molecules found")
            
            logger.info("Transaction completed successfully")
            
            # Close connection
            adapter.close()
            logger.info("Connection closed")
            
            return True
        else:
            logger.error("❌ Connection failed")
            return False
        
    except ImportError as e:
        logger.error(f"Error importing Supabase adapter: {str(e)}")
        logger.error("The SupabaseAdapter class may not exist in database/adapters/supabase.py")
        logger.error("Checking if supabase.py exists...")
        
        adapter_path = os.path.join(os.path.dirname(__file__), "database", "adapters", "supabase.py")
        if os.path.exists(adapter_path):
            logger.info(f"Supabase adapter file exists at {adapter_path}")
            
            # Try to find the adapter directory structure
            with open(adapter_path, 'r') as f:
                contents = f.read()
                logger.info(f"File size: {len(contents)} bytes")
                logger.info(f"First 100 characters: {contents[:100]}...")
            
            # Check if there are any classes defined
            import re
            classes = re.findall(r'class\s+(\w+)', contents)
            logger.info(f"Classes defined in file: {classes}")
            
            # Import a different way
            spec = importlib.util.spec_from_file_location("supabase_adapter", adapter_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Check what's in the module
            module_items = dir(module)
            logger.info(f"Items in module: {[item for item in module_items if not item.startswith('__')]}")
            
        else:
            logger.error(f"Supabase adapter file does not exist at {adapter_path}")
            
            # Try to find adapter files
            adapter_dir = os.path.join(os.path.dirname(__file__), "database", "adapters")
            if os.path.exists(adapter_dir):
                logger.info(f"Adapter directory exists at {adapter_dir}")
                adapter_files = os.listdir(adapter_dir)
                logger.info(f"Files in adapter directory: {adapter_files}")
            else:
                logger.error(f"Adapter directory does not exist at {adapter_dir}")
        
        return False
    except Exception as e:
        logger.error(f"Error in Supabase adapter test: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def main():
    """Main function."""
    # Parse arguments
    args = parse_arguments()
    
    # Set log level based on verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    print("=" * 60)
    print("Supabase Adapter Debugger")
    print("=" * 60)
    
    # Find configuration path
    config_path = find_config_path(args)
    logger.info(f"Using configuration file: {config_path}")
    
    # Load Supabase configuration
    supabase_config = get_supabase_config(config_path)
    
    if not supabase_config:
        print("❌ Failed to load Supabase configuration")
        return
    
    # Log configuration without password
    safe_config = {k: v for k, v in supabase_config.items() if k \!= 'password'}
    logger.info(f"Supabase configuration: {json.dumps(safe_config, indent=2)}")
    
    # Test connection
    if args.direct:
        success = test_direct_connection(supabase_config)
    else:
        success = test_supabase_adapter(supabase_config)
    
    if success:
        print("\n✅ Supabase connection test completed successfully")
    else:
        print("\n❌ Supabase connection test failed")
    
    print("=" * 60)

if __name__ == "__main__":
    main()

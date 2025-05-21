#\!/usr/bin/env python3
"""
Debug Connection Adapter

This script debugs the database connection adapter by adding detailed logging
and tracing the execution flow through the adapter classes. It helps identify
where connection issues might be occurring in the abstraction layer.
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
logger = logging.getLogger('debug_connection')

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Debug database connection adapter')
    parser.add_argument('--config', '-c', dest='config_path', 
                        help='Path to the configuration file (overrides environment variable)')
    parser.add_argument('--adapter', '-a', dest='adapter_type', default='supabase',
                        choices=['local', 'supabase', 'pooler'],
                        help='Specific adapter type to test (default: supabase)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose output with detailed tracing')
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

def setup_debug_environment(args):
    """Set up debug environment with monkey patching."""
    # Add project root to path
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    
    # First, import the modules we need to monkey patch
    try:
        from database.connection import ConnectionFactory
        from database.connection_config import validate_config, get_connection_config
        modules_imported = True
    except ImportError as e:
        logger.error(f"Error importing modules: {str(e)}")
        modules_imported = False
        
    if not modules_imported:
        return False
    
    # Only proceed if modules were successfully imported
    if modules_imported and args.verbose:
        # Save original methods
        original_get_connection = ConnectionFactory.get_connection
        
        # Create a debug wrapper for get_connection method
        def debug_get_connection(self, adapter_type=None, config=None):
            logger.debug(f"=== ConnectionFactory.get_connection called ===")
            logger.debug(f"adapter_type: {adapter_type}")
            logger.debug(f"config: {config if config else 'Using default config'}")
            
            try:
                logger.debug("Calling original get_connection method")
                result = original_get_connection(self, adapter_type, config)
                logger.debug(f"get_connection result: {result}")
                return result
            except Exception as e:
                logger.error(f"Exception in get_connection: {str(e)}")
                logger.error(traceback.format_exc())
                raise
        
        # Apply the monkey patches
        ConnectionFactory.get_connection = debug_get_connection
        logger.debug("Applied debug wrappers to ConnectionFactory methods")
        
    return True

def debug_connection_creation(args, config_path):
    """Debug the connection creation process."""
    logger.info(f"Using configuration file: {config_path}")
    
    # Set the environment variable
    os.environ["CRYOPROTECT_CONFIG_PATH"] = config_path
    
    try:
        # Load configuration
        with open(config_path, 'r') as f:
            config_data = json.load(f)
        logger.debug(f"Loaded configuration from {config_path}")
        
        # Display connection mode
        connection_mode = config_data.get('database', {}).get('connection', {}).get('mode', 'unknown')
        logger.info(f"Connection mode from config: {connection_mode}")
        
        # Import required modules
        from database.connection import ConnectionFactory
        from database.connection_config import validate_config, get_connection_config
        
        # Validate configuration
        logger.info("Validating configuration...")
        validate_config()
        logger.info("Configuration validation successful")
        
        # Get configuration for adapter
        adapter_type = args.adapter_type
        logger.info(f"Getting configuration for adapter type: {adapter_type}")
        adapter_config = get_connection_config(adapter_type)
        
        if not adapter_config:
            logger.error(f"No configuration found for adapter type: {adapter_type}")
            return False
        
        # Log adapter configuration (excluding password)
        safe_config = {k: v for k, v in adapter_config.items() if k \!= 'password'}
        logger.info(f"Adapter configuration: {json.dumps(safe_config, indent=2)}")
        
        # Create connection factory
        logger.info("Creating connection factory")
        factory = ConnectionFactory()
        
        # Try to get a connection
        logger.info(f"Trying to get connection with adapter type: {adapter_type}")
        connection = factory.get_connection(adapter_type=adapter_type)
        
        if connection:
            logger.info("✅ Connection established successfully")
            logger.info(f"Connection type: {type(connection).__name__}")
            
            # Try a simple query
            logger.info("Executing test query...")
            try:
                results = connection.execute_query("SELECT COUNT(*) FROM molecules;")
                if results:
                    count = results[0].get('count', 0)
                    logger.info(f"✅ Query successful: {count} molecules found")
                else:
                    logger.warning("Query returned no results")
            except Exception as e:
                logger.error(f"Error executing query: {str(e)}")
            
            return True
        else:
            logger.error("❌ Failed to establish connection")
            return False
        
    except Exception as e:
        logger.error(f"Exception during connection debugging: {str(e)}")
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
    print("Database Connection Adapter Debugger")
    print("=" * 60)
    
    # Find configuration path
    config_path = find_config_path(args)
    
    # Set up debug environment
    if not setup_debug_environment(args):
        print("Failed to set up debug environment")
        return
    
    # Debug connection creation
    success = debug_connection_creation(args, config_path)
    
    if success:
        print("\n✅ Connection debugging completed successfully")
    else:
        print("\n❌ Connection debugging failed")
    
    print("=" * 60)

if __name__ == "__main__":
    main()

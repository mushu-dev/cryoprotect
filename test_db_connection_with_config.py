#\!/usr/bin/env python3
"""
Test Database Connection With Configuration System

This script tests the database connection using the new configuration system.
It can use the CRYOPROTECT_CONFIG_PATH environment variable or a path 
specified as a command-line argument.
"""

import os
import sys
import json
import argparse
from pathlib import Path

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Test database connection with configuration system')
    parser.add_argument('--config', '-c', dest='config_path', 
                        help='Path to the configuration file (overrides environment variable)')
    return parser.parse_args()

def find_config_path():
    """Find the configuration path from environment variable or default location."""
    # Check command-line arguments
    args = parse_arguments()
    if args.config_path:
        if os.path.exists(args.config_path):
            return args.config_path
        else:
            print(f"Error: Configuration file not found at {args.config_path}")
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
    print("Error: No configuration file found. Please specify a path with --config or set CRYOPROTECT_CONFIG_PATH.")
    print("\nYou can create a test configuration file with:")
    print("  python config/create_test_config.py")
    sys.exit(1)

# Find configuration path
config_path = find_config_path()

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Set the environment variable for the configuration path
os.environ["CRYOPROTECT_CONFIG_PATH"] = config_path

try:
    # Import the configuration system
    from database.connection_config import (
        validate_config,
        get_connection_config,
        is_adapter_enabled,
        get_adapter_order,
        test_adapter_configuration
    )
    
    # Import the connection factory
    from database.connection import ConnectionFactory
    
    # Import the database modules
    from database import db, db_service_role, db_public
    
except ImportError as e:
    print(f"Error importing database modules: {str(e)}")
    print("Make sure you run this script from the project root.")
    sys.exit(1)

def test_configuration():
    """Test the configuration system."""
    print("\n=== Testing Configuration System ===")
    
    try:
        # Validate configuration
        print("Validating configuration...")
        validate_config()
        print("✅ Configuration validated successfully")
        
        # Check adapter order
        print("\nChecking adapter order...")
        adapter_order = get_adapter_order()
        print(f"Adapter order: {adapter_order}")
        
        # Check enabled adapters
        print("\nChecking enabled adapters...")
        for adapter in ['local', 'supabase', 'pooler']:
            enabled = is_adapter_enabled(adapter)
            status = "✅ Enabled" if enabled else "❌ Disabled"
            print(f"  - {adapter}: {status}")
        
        # Test adapter configurations
        print("\nTesting adapter configurations...")
        for adapter in adapter_order:
            success, message = test_adapter_configuration(adapter)
            status = "✅ Valid" if success else "❌ Invalid"
            print(f"  - {adapter}: {status} - {message}")
        
        return True
    except Exception as e:
        print(f"❌ Error testing configuration: {str(e)}")
        return False

def test_connection_factory():
    """Test the connection factory."""
    print("\n=== Testing Connection Factory ===")
    
    try:
        # Create connection factory
        print("Creating connection factory...")
        factory = ConnectionFactory()
        print("✅ Connection factory created successfully")
        
        # Get connection
        print("\nGetting connection...")
        connection = factory.get_connection()
        
        if connection:
            print("✅ Connection established successfully")
            print(f"Connection type: {type(connection).__name__}")
            
            # Test a simple query
            print("\nExecuting test query...")
            results = connection.execute_query("SELECT COUNT(*) FROM molecules;")
            if results:
                count = results[0].get('count', 0)
                print(f"✅ Query successful: {count} molecules found")
            else:
                print("❌ Query returned no results")
            
            return True
        else:
            print("❌ Failed to establish connection")
            return False
    except Exception as e:
        print(f"❌ Error testing connection factory: {str(e)}")
        return False

def test_database_modules():
    """Test the database modules."""
    print("\n=== Testing Database Modules ===")
    
    try:
        # Test db module
        print("Testing db module...")
        result = db.execute_query("SELECT COUNT(*) FROM molecules;")
        if result:
            count = result[0].get('count', 0)
            print(f"✅ db module query successful: {count} molecules found")
        else:
            print("❌ db module query returned no results")
        
        # Test db_service_role module
        print("\nTesting db_service_role module...")
        result = db_service_role.execute_query("SELECT COUNT(*) FROM molecules;")
        if result:
            count = result[0].get('count', 0)
            print(f"✅ db_service_role module query successful: {count} molecules found")
        else:
            print("❌ db_service_role module query returned no results")
        
        # Test db_public module
        print("\nTesting db_public module...")
        result = db_public.execute_query("SELECT COUNT(*) FROM molecules;")
        if result:
            count = result[0].get('count', 0)
            print(f"✅ db_public module query successful: {count} molecules found")
        else:
            print("❌ db_public module query returned no results")
        
        return True
    except Exception as e:
        print(f"❌ Error testing database modules: {str(e)}")
        return False

def main():
    """Main function."""
    print("=" * 60)
    print("Database Connection Test With Configuration System")
    print("=" * 60)
    print(f"Using configuration file: {config_path}")
    
    # Load and display configuration
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        connection_mode = config.get('database', {}).get('connection', {}).get('mode', 'unknown')
        print(f"Connection mode: {connection_mode}")
    except Exception as e:
        print(f"Error loading configuration: {str(e)}")
    
    # Run tests
    config_success = test_configuration()
    if config_success:
        factory_success = test_connection_factory()
        if factory_success:
            db_success = test_database_modules()
            if db_success:
                print("\n✅ All tests passed successfully\!")
            else:
                print("\n❌ Database module tests failed")
        else:
            print("\n❌ Connection factory tests failed")
    else:
        print("\n❌ Configuration tests failed")
    
    print("=" * 60)

if __name__ == "__main__":
    main()

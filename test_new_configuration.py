#!/usr/bin/env python3
"""
CryoProtect v2 - Test New Configuration System

This script tests the new configuration system by attempting to establish connections
to the database using the different connection methods.

Tests:
1. Configuration validation
2. Local database connection
3. Supabase direct connection
4. Connection pooler
5. Service role connection
6. Public (authenticated) connection
"""

import os
import sys
import logging
import time
import json
from typing import Dict, Any, Tuple, List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("configuration_test")

# Try to import from configuration system
try:
    from database.connection_config import (
        validate_config,
        get_connection_config,
        test_adapter_configuration,
        get_adapter_order,
        is_adapter_enabled,
        get_adapter_configs
    )
    logger.info("Successfully imported configuration system modules")
except ImportError as e:
    logger.error(f"Failed to import configuration system modules: {str(e)}")
    sys.exit(1)

# Try to import connection factory
try:
    from database.connection import get_db_connection, get_db_connection_info
    logger.info("Successfully imported connection factory")
except ImportError as e:
    logger.error(f"Failed to import connection factory: {str(e)}")
    sys.exit(1)

# Try to import specific database modules
try:
    from database.db import get_connection as get_db_connection
    from database.db import test_connection as test_db_connection
    from database.db import close_all_connections as close_db_connections
    from database.db_service_role import get_connection as get_service_role_connection
    from database.db_service_role import test_connection as test_service_role_connection
    from database.db_service_role import close_all_connections as close_service_role_connections
    from database.db_public import get_connection as get_public_connection
    from database.db_public import test_connection as test_public_connection
    from database.db_public import close_all_connections as close_public_connections
    logger.info("Successfully imported database modules")
except ImportError as e:
    logger.error(f"Failed to import database modules: {str(e)}")
    sys.exit(1)

# Try to import Supabase client
try:
    from api.supabase_client import create_client
    logger.info("Successfully imported Supabase client")
    HAS_SUPABASE_CLIENT = True
except ImportError:
    logger.warning("Supabase client not available, skipping related tests")
    HAS_SUPABASE_CLIENT = False

def test_configuration_validation():
    """Test configuration validation."""
    logger.info("Testing configuration validation...")
    try:
        validate_config()
        logger.info("✅ Configuration validation successful")
        return True
    except Exception as e:
        logger.error(f"❌ Configuration validation failed: {str(e)}")
        return False

def test_get_connection_config():
    """Test getting connection configurations."""
    logger.info("Testing get_connection_config...")
    
    results = {}
    
    # Test getting all adapters
    try:
        adapter_order = get_adapter_order()
        logger.info(f"Adapter order: {adapter_order}")
        results["adapter_order"] = adapter_order
        
        # Test each adapter
        for adapter in adapter_order:
            try:
                is_enabled = is_adapter_enabled(adapter)
                logger.info(f"Adapter {adapter} enabled: {is_enabled}")
                
                if is_enabled:
                    config = get_connection_config(adapter)
                    # Filter out sensitive information for logging
                    safe_config = {k: v for k, v in config.items() 
                                 if k not in ['password', 'key', 'service_key']}
                    logger.info(f"Got configuration for {adapter}: {safe_config}")
                    results[f"{adapter}_config"] = True
                else:
                    results[f"{adapter}_config"] = False
            except Exception as e:
                logger.error(f"Failed to get configuration for {adapter}: {str(e)}")
                results[f"{adapter}_config"] = False
        
        return True, results
    except Exception as e:
        logger.error(f"Failed to get connection configurations: {str(e)}")
        return False, {}

def test_connection_factory():
    """Test the connection factory."""
    logger.info("Testing connection factory...")
    
    try:
        # Get the connection
        db = get_db_connection()
        if not db:
            logger.error("❌ Failed to get connection from factory")
            return False
        
        # Test the connection
        success, message = db.test_connection()
        if success:
            logger.info(f"✅ Connection test successful: {message}")
            
            # Get connection info
            try:
                info = get_db_connection_info()
                adapter_type = info.get('adapter_type', 'unknown')
                logger.info(f"Connection is using adapter type: {adapter_type}")
            except Exception as e:
                logger.warning(f"Failed to get connection info: {str(e)}")
            
            return True
        else:
            logger.error(f"❌ Connection test failed: {message}")
            return False
    except Exception as e:
        logger.error(f"❌ Connection factory test failed: {str(e)}")
        return False
        
def test_db_module():
    """Test the database module."""
    logger.info("Testing db.py module...")
    
    try:
        # Test connection
        success, message = test_db_connection()
        if success:
            logger.info(f"✅ db.py connection test successful: {message}")
            
            # Get a connection and test executing a query
            conn = get_db_connection()
            if conn:
                logger.info("Successfully got connection from db module")
                # Make sure to close connections
                close_db_connections()
                return True
            else:
                logger.error("Failed to get connection from db module")
                return False
        else:
            logger.error(f"❌ db.py connection test failed: {message}")
            return False
    except Exception as e:
        logger.error(f"❌ db.py module test failed: {str(e)}")
        return False

def test_service_role_module():
    """Test the service role module."""
    logger.info("Testing db_service_role.py module...")
    
    try:
        # Test connection
        success, message = test_service_role_connection()
        if success:
            logger.info(f"✅ db_service_role.py connection test successful: {message}")
            
            # Get a connection and test executing a query
            conn = get_service_role_connection()
            if conn:
                logger.info("Successfully got connection from service role module")
                # Make sure to close connections
                close_service_role_connections()
                return True
            else:
                logger.error("Failed to get connection from service role module")
                return False
        else:
            logger.error(f"❌ db_service_role.py connection test failed: {message}")
            return False
    except Exception as e:
        logger.error(f"❌ db_service_role.py module test failed: {str(e)}")
        return False

def test_public_module():
    """Test the public module."""
    logger.info("Testing db_public.py module...")
    
    try:
        # Test connection
        success, message = test_public_connection()
        if success:
            logger.info(f"✅ db_public.py connection test successful: {message}")
            
            # Get a connection and test executing a query
            conn = get_public_connection()
            if conn:
                logger.info("Successfully got connection from public module")
                # Make sure to close connections
                close_public_connections()
                return True
            else:
                logger.error("Failed to get connection from public module")
                return False
        else:
            logger.error(f"❌ db_public.py connection test failed: {message}")
            return False
    except Exception as e:
        logger.error(f"❌ db_public.py module test failed: {str(e)}")
        return False

def test_supabase_client():
    """Test the Supabase client."""
    if not HAS_SUPABASE_CLIENT:
        logger.warning("Supabase client not available, skipping test")
        return None
        
    logger.info("Testing Supabase client...")
    
    try:
        client = create_client()
        if client:
            # Test a simple query
            try:
                response = client.table('property_types').select('count', count='exact').limit(1).execute()
                if hasattr(response, 'count'):
                    logger.info(f"✅ Supabase client successfully executed query, found {response.count} rows")
                    return True
                else:
                    logger.error("❌ Supabase client query did not return expected response")
                    return False
            except Exception as e:
                logger.error(f"❌ Supabase client failed to execute query: {str(e)}")
                return False
        else:
            logger.error("❌ Failed to create Supabase client")
            return False
    except Exception as e:
        logger.error(f"❌ Supabase client test failed: {str(e)}")
        return False

def generate_report(results):
    """
    Generate a report of the test results.
    
    Args:
        results: Dictionary of test results
        
    Returns:
        Report as a string
    """
    report = "# Database Configuration Test Report\n\n"
    report += f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    # Overall status
    overall_success = all(result for result in results.values() if result is not None)
    report += f"## Overall Status: {'✅ Success' if overall_success else '❌ Failed'}\n\n"
    
    # Individual test results
    report += "## Test Results\n\n"
    for test_name, result in results.items():
        status = "✅ Success" if result else "❌ Failed" if result is not None else "⏩ Skipped"
        report += f"- {test_name}: {status}\n"
    
    return report

def main():
    """Main entry point."""
    results = {}
    
    # Test configuration validation
    results["Configuration Validation"] = test_configuration_validation()
    
    # Test getting connection configurations
    config_success, config_results = test_get_connection_config()
    results["Get Connection Configurations"] = config_success
    
    # Test connection factory
    results["Connection Factory"] = test_connection_factory()
    
    # Test database modules
    results["DB Module"] = test_db_module()
    results["Service Role Module"] = test_service_role_module()
    results["Public Module"] = test_public_module()
    
    # Test Supabase client
    results["Supabase Client"] = test_supabase_client()
    
    # Generate and display report
    report = generate_report(results)
    print("\n" + report)
    
    # Save report to file
    report_path = "database_configuration_test_report.md"
    with open(report_path, "w") as f:
        f.write(report)
    logger.info(f"Report saved to {report_path}")
    
    # Return success or failure
    return 0 if all(result for result in results.values() if result is not None) else 1

if __name__ == "__main__":
    sys.exit(main())
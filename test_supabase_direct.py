#!/usr/bin/env python3
"""
Test script for supabase_direct.py

This script tests the direct PostgreSQL connection to Supabase
using the SupabaseDirectConnection class.
"""

import os
import sys
import logging
from supabase_direct import SupabaseDirectConnection

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_connection():
    """Test the connection to Supabase."""
    try:
        # Get connection instance
        db = SupabaseDirectConnection.get_instance()
        
        # Test connection health
        if db.check_connection_health():
            logger.info("Connection health check passed")
        else:
            logger.error("Connection health check failed")
            return False
        
        # Execute simple query
        result = db.execute_query("SELECT 1 as test")
        logger.info(f"Query result: {result}")
        
        if result and result[0]['test'] == 1:
            logger.info("Basic query test passed")
        else:
            logger.error("Basic query test failed")
            return False
        
        # Test parameterized query
        result = db.execute_query("SELECT %s as param_test", {"param_test": "success"})
        logger.info(f"Parameterized query result: {result}")
        
        if result and result[0]['param_test'] == "success":
            logger.info("Parameterized query test passed")
        else:
            logger.error("Parameterized query test failed")
            return False
        
        # Test batch execution
        queries = [
            "CREATE TEMPORARY TABLE test_batch (id SERIAL PRIMARY KEY, name TEXT)",
            "INSERT INTO test_batch (name) VALUES ('test1')",
            "INSERT INTO test_batch (name) VALUES ('test2')"
        ]
        
        if db.execute_batch(queries):
            logger.info("Batch execution test passed")
        else:
            logger.error("Batch execution test failed")
            return False
        
        # Verify batch execution
        result = db.execute_query("SELECT * FROM test_batch ORDER BY id")
        logger.info(f"Batch verification result: {result}")
        
        if result and len(result) == 2:
            logger.info("Batch verification passed")
        else:
            logger.error("Batch verification failed")
            return False
        
        # Get stats
        stats = db.get_stats()
        logger.info(f"Connection stats: {stats}")
        
        return True
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        return False
    finally:
        # Close all connections
        if 'db' in locals():
            db.close_all()
            logger.info("Closed all connections")

def setup_test_environment():
    """Set up the test environment with required environment variables."""
    # Check if environment variables are already set
    if os.environ.get('SUPABASE_DB_HOST') and os.environ.get('SUPABASE_DB_PASSWORD'):
        return True
    
    # Prompt for credentials if not set
    print("Supabase database credentials not found in environment variables.")
    print("Please enter the following information for testing:")
    
    host = input("Supabase DB Host: ")
    port = input("Supabase DB Port [5432]: ") or "5432"
    name = input("Supabase DB Name [postgres]: ") or "postgres"
    user = input("Supabase DB User [postgres]: ") or "postgres"
    password = input("Supabase DB Password: ")
    
    if not host or not password:
        print("Error: Host and password are required.")
        return False
    
    # Set environment variables
    os.environ['SUPABASE_DB_HOST'] = host
    os.environ['SUPABASE_DB_PORT'] = port
    os.environ['SUPABASE_DB_NAME'] = name
    os.environ['SUPABASE_DB_USER'] = user
    os.environ['SUPABASE_DB_PASSWORD'] = password
    
    return True

if __name__ == "__main__":
    print("Testing SupabaseDirectConnection...")
    
    if not setup_test_environment():
        print("Failed to set up test environment.")
        sys.exit(1)
    
    if test_connection():
        print("All tests passed successfully!")
        sys.exit(0)
    else:
        print("Tests failed. See log for details.")
        sys.exit(1)
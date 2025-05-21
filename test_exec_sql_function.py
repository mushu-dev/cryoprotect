#!/usr/bin/env python3
"""
CryoProtect v2 - Test exec_sql Function

This script tests the exec_sql function in the Supabase database to ensure it works correctly.
It executes a simple SQL query using the exec_sql function and verifies the result.
"""

import os
import json
import logging
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def get_supabase_client():
    """Get a Supabase client with service role key."""
    try:
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY")
        
        if not SUPABASE_URL or not SUPABASE_KEY:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        raise

def test_exec_sql_function():
    """Test the exec_sql function with various SQL queries."""
    try:
        # Get Supabase client
        supabase = get_supabase_client()
        
        # Test 1: Simple SELECT query
        logger.info("Test 1: Simple SELECT query")
        sql = "SELECT 1 as test"
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False
        
        logger.info(f"Result: {response.data}")
        
        # Test 2: CREATE TABLE query
        logger.info("\nTest 2: CREATE TABLE query")
        sql = """
        CREATE TABLE IF NOT EXISTS public.exec_sql_test (
            id SERIAL PRIMARY KEY,
            name TEXT,
            created_at TIMESTAMPTZ DEFAULT NOW()
        )
        """
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False
        
        logger.info(f"Result: {response.data}")
        
        # Test 3: INSERT query
        logger.info("\nTest 3: INSERT query")
        sql = """
        INSERT INTO public.exec_sql_test (name)
        VALUES ('Test 1'), ('Test 2')
        """
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False
        
        logger.info(f"Result: {response.data}")
        
        # Test 4: SELECT query with results
        logger.info("\nTest 4: SELECT query with results")
        sql = "SELECT * FROM public.exec_sql_test"
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False
        
        logger.info(f"Result: {response.data}")
        
        # Test 5: Error handling (intentional error)
        logger.info("\nTest 5: Error handling (intentional error)")
        sql = "SELECT * FROM non_existent_table"
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Expected error: {response.error}")
        else:
            logger.info(f"Result: {response.data}")
        
        # Test 6: Clean up test table
        logger.info("\nTest 6: Clean up test table")
        sql = "DROP TABLE IF EXISTS public.exec_sql_test"
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False
        
        logger.info(f"Result: {response.data}")
        
        logger.info("\nAll tests completed successfully!")
        return True
    
    except Exception as e:
        logger.error(f"Error testing exec_sql function: {str(e)}")
        return False

if __name__ == "__main__":
    logger.info("Testing exec_sql function...")
    success = test_exec_sql_function()
    
    if success:
        logger.info("exec_sql function is working correctly!")
    else:
        logger.error("exec_sql function test failed!")
#!/usr/bin/env python3
"""
Test Migration Tracking Table Creation (Direct Implementation)

This script tests that the migration tracking table can be created directly
without relying on the has_table RPC function.
"""

import os
import sys
import logging
import json
from dotenv import load_dotenv
from supabase import create_client, Client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("test_migration_tracking.log")
    ]
)
logger = logging.getLogger(__name__)

def get_supabase_client():
    """
    Get a Supabase client using environment variables.
    
    Returns:
        Client: Supabase client
    """
    supabase_url = os.getenv("SUPABASE_URL")
    supabase_key = os.getenv("SUPABASE_KEY")
    
    if not supabase_url or not supabase_key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
    
    return create_client(supabase_url, supabase_key)

def create_migrations_table_direct(supabase):
    """
    Create the migrations table directly using a PostgreSQL query.
    
    Args:
        supabase: Supabase client
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Use PostgreSQL query to create the table if it doesn't exist
        query = """
        CREATE TABLE IF NOT EXISTS migrations (
            id SERIAL PRIMARY KEY,
            version TEXT NOT NULL,
            name TEXT NOT NULL,
            applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        -- Create unique index if it doesn't exist
        DO $$
        BEGIN
            IF NOT EXISTS (
                SELECT 1 FROM pg_indexes 
                WHERE indexname = 'migrations_version_idx'
            ) THEN
                CREATE UNIQUE INDEX migrations_version_idx ON migrations(version);
            END IF;
        END $$;
        
        -- Return success status
        SELECT 'success' as status;
        """
        
        # Execute the query
        response = supabase.rpc('exec_sql', {'query': query}).execute()
        
        # Log the response for debugging
        logger.info(f"Create table response: {response}")
        
        # Check if the query was successful
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error creating migrations table: {response.error}")
            return False
        
        logger.info("Created migrations table")
        return True
    except Exception as e:
        logger.error(f"Error creating migrations table: {str(e)}")
        return False

def check_migrations_table_exists(supabase):
    """
    Check if the migrations table exists.
    
    Args:
        supabase: Supabase client
        
    Returns:
        bool: True if the table exists, False otherwise
    """
    try:
        # Try to query the migrations table
        query = "SELECT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_schema = 'public' AND table_name = 'migrations') AS table_exists;"
        response = supabase.rpc('exec_sql', {'query': query}).execute()
        
        # Log the response for debugging
        logger.info(f"Check table exists response: {response}")
        
        # Check if the query was successful
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error checking if migrations table exists: {response.error}")
            return False
        
        # Check if the table exists
        if hasattr(response, 'data') and response.data:
            # The response data might be in different formats depending on the implementation
            # Try different approaches to extract the boolean value
            
            # If data is a list of dictionaries
            if isinstance(response.data, list) and len(response.data) > 0:
                if 'table_exists' in response.data[0]:
                    return response.data[0]['table_exists']
                elif 'exists' in response.data[0]:
                    return response.data[0]['exists']
            
            # If data is a string (JSON)
            if isinstance(response.data, str):
                try:
                    data = json.loads(response.data)
                    if isinstance(data, list) and len(data) > 0:
                        if 'table_exists' in data[0]:
                            return data[0]['table_exists']
                        elif 'exists' in data[0]:
                            return data[0]['exists']
                except:
                    pass
            
            # If we can't parse the response, try a direct query to the table
            try:
                direct_query = "SELECT COUNT(*) FROM migrations;"
                direct_response = supabase.rpc('exec_sql', {'query': direct_query}).execute()
                logger.info(f"Direct query response: {direct_response}")
                
                # If we get here without an error, the table exists
                return True
            except Exception as e:
                logger.warning(f"Error in direct query: {str(e)}")
                return False
        
        return False
    except Exception as e:
        logger.error(f"Error checking if migrations table exists: {str(e)}")
        return False

def run_test():
    """
    Run the migration tracking table creation test.
    
    Returns:
        dict: Test results with status, summary, and details
    """
    results = {
        "status": "ERROR",
        "summary": "Test not completed",
        "details": []
    }
    
    try:
        logger.info("Starting migration tracking table creation test")
        
        # Get Supabase client
        logger.info("Getting Supabase client")
        supabase = get_supabase_client()
        
        # Create the migrations table (using IF NOT EXISTS so it's safe to run multiple times)
        logger.info("Creating migrations table")
        success = create_migrations_table_direct(supabase)
        
        if success:
            logger.info("Successfully created migrations table")
            
            # Try a direct query to verify the table exists
            try:
                query = "SELECT COUNT(*) FROM migrations;"
                response = supabase.rpc('exec_sql', {'query': query}).execute()
                logger.info(f"Verification query response: {response}")
                
                # If we get here without an error, the table exists
                logger.info("Confirmed migrations table exists")
                results["status"] = "SUCCESS"
                results["summary"] = "Migration tracking table created successfully"
                results["details"].append("Successfully created migrations table")
                results["details"].append("Confirmed table exists by querying it")
            except Exception as e:
                logger.error(f"Error verifying migrations table: {str(e)}")
                results["status"] = "ERROR"
                results["summary"] = "Failed to verify migrations table exists"
                results["details"].append(f"Error verifying migrations table: {str(e)}")
        else:
            logger.error("Failed to create migrations table")
            results["status"] = "ERROR"
            results["summary"] = "Failed to create migrations table"
            results["details"].append("Failed to create migrations table")
    
    except Exception as e:
        logger.error(f"Error during test: {str(e)}")
        results["status"] = "ERROR"
        results["summary"] = f"Error during test: {str(e)}"
        results["details"].append(f"Exception: {str(e)}")
        
    return results

def main():
    """Main entry point for the script."""
    # Load environment variables
    load_dotenv()
    
    # Run the test
    results = run_test()
    
    # Print results
    print("\n=== Test Results ===")
    print(f"Status: {results['status']}")
    print(f"Summary: {results['summary']}")
    print("\nDetails:")
    for detail in results["details"]:
        print(f"- {detail}")
    
    # Return exit code based on status
    return 0 if results["status"] == "SUCCESS" else 1

if __name__ == "__main__":
    sys.exit(main())
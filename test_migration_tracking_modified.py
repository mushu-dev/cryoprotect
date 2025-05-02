#!/usr/bin/env python3
"""
Test Migration Tracking Table Creation (Modified)

This script tests that the migration tracking table can be created by the migration management module.
It uses the table() method instead of sql() to interact with the database.
"""

import os
import sys
import logging
from dotenv import load_dotenv

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

def check_table_exists(supabase, table_name):
    """
    Check if a table exists using information_schema query.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table to check
        
    Returns:
        bool: True if the table exists, False otherwise
    """
    try:
        # Query information_schema.tables to check if the table exists
        response = supabase.table("information_schema.tables").select("table_name").eq("table_schema", "public").eq("table_name", table_name).execute()
        return len(response.data) > 0 if response.data else False
    except Exception as e:
        logger.error(f"Error checking if table exists: {str(e)}")
        return False

def create_migrations_table(supabase):
    """
    Create the migrations table directly.
    
    Args:
        supabase: Supabase client
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # First check if the table already exists
        if check_table_exists(supabase, 'migrations'):
            logger.info("Migrations table already exists")
            return True
        
        # If we're here, we need to create the table
        # Since we can't use raw SQL directly, we'll use the database.migrations.initialize_migration_tracking function
        from database.migrations import initialize_migration_tracking
        success = initialize_migration_tracking(supabase)
        
        if success:
            logger.info("Created migrations table")
            return True
        else:
            logger.error("Failed to create migrations table")
            return False
    except Exception as e:
        logger.error(f"Error creating migrations table: {str(e)}")
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
        
        # Import required modules
        from service_role_helper import get_supabase_client
        from database.migrations import initialize_migration_tracking
        
        # Get Supabase client
        logger.info("Getting Supabase client")
        supabase = get_supabase_client()
        
        # Call initialize_migration_tracking
        logger.info("Calling initialize_migration_tracking")
        success = initialize_migration_tracking(supabase)
        
        # Check result
        if success:
            logger.info("Migration tracking table created or already exists")
            results["status"] = "SUCCESS"
            results["summary"] = "Migration tracking table created or already exists"
            results["details"].append("initialize_migration_tracking returned True")
            
            # Verify table exists by querying it
            try:
                logger.info("Verifying migrations table exists")
                table_exists = check_table_exists(supabase, 'migrations')
                
                if table_exists:
                    logger.info("Confirmed migrations table exists")
                    results["details"].append("Confirmed migrations table exists via information_schema query")
                    
                    # Try to query the migrations table
                    try:
                        response = supabase.table("migrations").select("*").limit(5).execute()
                        logger.info(f"Successfully queried migrations table, found {len(response.data) if response.data else 0} records")
                        results["details"].append(f"Successfully queried migrations table, found {len(response.data) if response.data else 0} records")
                    except Exception as e:
                        logger.warning(f"Error querying migrations table: {str(e)}")
                        results["details"].append(f"Warning: Error querying migrations table: {str(e)}")
                else:
                    logger.warning("Could not confirm migrations table exists")
                    results["details"].append("Warning: Could not confirm migrations table exists via information_schema query")
            except Exception as e:
                logger.warning(f"Error verifying migrations table: {str(e)}")
                results["details"].append(f"Warning: Error verifying migrations table: {str(e)}")
        else:
            logger.error("Failed to create migration tracking table")
            results["status"] = "ERROR"
            results["summary"] = "Failed to create migration tracking table"
            results["details"].append("initialize_migration_tracking returned False")
    
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
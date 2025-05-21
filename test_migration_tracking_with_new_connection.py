#!/usr/bin/env python3
"""
Test Migration Tracking Table Creation with New Connection Utility

This script tests that the migration tracking table can be created by the migration management module
using the newly implemented create_connection function from database.utils.connection.
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
        logging.FileHandler("test_migration_tracking_new_connection.log")
    ]
)
logger = logging.getLogger(__name__)

def run_test():
    """
    Run the migration tracking table creation test using the new connection utility.
    
    Returns:
        dict: Test results with status, summary, and details
    """
    results = {
        "status": "ERROR",
        "summary": "Test not completed",
        "details": []
    }
    
    try:
        logger.info("Starting migration tracking table creation test with new connection utility")
        
        # Import required modules
        from database.utils.connection import create_connection
        from database.migrations.runner import initialize_migration_tracking
        
        # Create database connection using our new utility
        logger.info("Creating database connection using new connection utility")
        supabase = create_connection()
        results["details"].append("Successfully created connection using new connection utility")
        
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
                response = supabase.rpc('has_table', {'table_name': 'migrations'}).execute()
                table_exists = response.data[0] if response.data else False
                
                if table_exists:
                    logger.info("Confirmed migrations table exists")
                    results["details"].append("Confirmed migrations table exists via has_table RPC")
                else:
                    logger.warning("Could not confirm migrations table exists")
                    results["details"].append("Warning: Could not confirm migrations table exists via has_table RPC")
                    
                # Try a direct query to verify the table exists
                try:
                    logger.info("Verifying migrations table with direct query")
                    response = supabase.table('migrations').select('count').limit(1).execute()
                    logger.info(f"Direct query response: {response}")
                    results["details"].append("Successfully queried migrations table directly")
                except Exception as e:
                    logger.warning(f"Error in direct query: {str(e)}")
                    results["details"].append(f"Warning: Error in direct query: {str(e)}")
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
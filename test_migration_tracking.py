#!/usr/bin/env python3
"""
Test Migration Tracking Table Creation

This script tests that the migration tracking table can be created by the migration management module.
It calls initialize_migration_tracking with a valid database connection and confirms that the
migrations table is created (or already exists) and that the function returns True.
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
        from database.utils.connection import create_connection
        from database.migrations.runner import initialize_migration_tracking
        
        # Get Supabase client using our new connection utility
        logger.info("Creating database connection using new connection utility")
        supabase = create_connection()
        
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
#!/usr/bin/env python3
"""
CryoProtect v2 - Test Fix Relationships Script

This script tests the fix_relationships.py script to ensure it correctly
identifies and fixes relationship design issues in the CryoProtect Supabase project.

Usage:
    python test_fix_relationships.py [--full-test]
"""

import os
import sys
import json
import argparse
import logging
import subprocess
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("test_fix_relationships.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase project ID
SUPABASE_PROJECT_ID = "tsdlmynydfuypiugmkev"

def run_command(command):
    """Run a command and return the output."""
    logger.info(f"Running command: {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Command output: {result.stdout}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with error: {e.stderr}")
        return None

def verify_supabase_connection():
    """Verify that we can connect to the Supabase project."""
    logger.info("Verifying Supabase connection...")
    
    # Check if the required environment variables are set
    supabase_url = os.getenv("SUPABASE_URL")
    supabase_key = os.getenv("SUPABASE_KEY")
    
    if not supabase_url or not supabase_key:
        logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        return False
    
    # Run a simple Python script to test the connection
    test_script = """
import os
from supabase import create_client, Client

supabase_url = os.getenv("SUPABASE_URL")
supabase_key = os.getenv("SUPABASE_KEY")

supabase = create_client(supabase_url, supabase_key)
response = supabase.table("molecules").select("*").limit(1).execute()
print(f"Connection successful. Found {len(response.data)} records in molecules table.")
"""
    
    with open("test_connection.py", "w") as f:
        f.write(test_script)
    
    result = run_command("python test_connection.py")
    os.remove("test_connection.py")
    
    if result and "Connection successful" in result:
        logger.info("Supabase connection verified")
        return True
    else:
        logger.error("Failed to connect to Supabase")
        return False

def test_dry_run():
    """Test the fix_relationships.py script with --dry-run flag."""
    logger.info("Testing fix_relationships.py with --dry-run flag...")
    
    result = run_command("python fix_relationships.py --dry-run")
    
    if result:
        logger.info("Dry run test completed successfully")
        return True
    else:
        logger.error("Dry run test failed")
        return False

def test_verify_only():
    """Test the fix_relationships.py script with --verify-only flag."""
    logger.info("Testing fix_relationships.py with --verify-only flag...")
    
    result = run_command("python fix_relationships.py --verify-only")
    
    if result:
        logger.info("Verify-only test completed successfully")
        return True
    else:
        logger.error("Verify-only test failed")
        return False

def run_full_test():
    """Run a full test of the fix_relationships.py script."""
    logger.info("Running full test of fix_relationships.py...")
    
    # First, run with --dry-run to see what would be done
    dry_run_result = test_dry_run()
    if not dry_run_result:
        logger.error("Dry run test failed. Aborting full test.")
        return False
    
    # Ask for confirmation before proceeding
    print("\n" + "=" * 80)
    print("WARNING: This will make actual changes to the database.")
    print("Make sure you have a backup before proceeding.")
    print("=" * 80)
    
    confirmation = input("\nDo you want to proceed with the full test? (yes/no): ")
    if confirmation.lower() != "yes":
        logger.info("Full test aborted by user")
        print("Full test aborted.")
        return False
    
    # Run the script without any flags to make actual changes
    logger.info("Running fix_relationships.py to make actual changes...")
    result = run_command("python fix_relationships.py")
    
    if not result:
        logger.error("Failed to run fix_relationships.py")
        return False
    
    # Verify the changes
    verify_result = test_verify_only()
    if not verify_result:
        logger.error("Verification after changes failed")
        return False
    
    logger.info("Full test completed successfully")
    return True

def main():
    """Main function to test the fix_relationships.py script."""
    parser = argparse.ArgumentParser(description='Test the fix_relationships.py script')
    parser.add_argument('--full-test', action='store_true', help='Run a full test including making actual changes')
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print(f"CryoProtect v2 - Test Fix Relationships Script (Project ID: {SUPABASE_PROJECT_ID})")
    print("=" * 80)
    
    try:
        # Verify Supabase connection
        print("\nVerifying Supabase connection...")
        if not verify_supabase_connection():
            print("Failed to connect to Supabase. Check the logs for details.")
            return 1
        
        # Test with --dry-run flag
        print("\nTesting with --dry-run flag...")
        if not test_dry_run():
            print("Dry run test failed. Check the logs for details.")
            return 1
        
        # Test with --verify-only flag
        print("\nTesting with --verify-only flag...")
        if not test_verify_only():
            print("Verify-only test failed. Check the logs for details.")
            return 1
        
        # Run full test if requested
        if args.full_test:
            print("\nRunning full test...")
            if not run_full_test():
                print("Full test failed. Check the logs for details.")
                return 1
            print("\nFull test completed successfully.")
        
        print("\n" + "=" * 60)
        print("Test Completed Successfully")
        print("=" * 60)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error testing fix_relationships.py: {str(e)}")
        print(f"\nError testing fix_relationships.py: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
CryoProtect v2 - API Fixes Runner

This script runs the necessary fixes to make the API fully functional:
1. Fixes the database table name mismatch
2. Verifies that the API endpoints are working correctly

Usage:
    python run_api_fixes.py
"""

import os
import sys
import json
import logging
import subprocess
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("api_fixes.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_command(command, description):
    """Run a command and log the output."""
    logger.info(f"Running {description}...")
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"{description} completed successfully")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"{description} failed: {e}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return False, e.stderr

def fix_database_tables():
    """Run the database table fix script."""
    success, output = run_command("python fix_database_tables.py", "Database table fix script")
    return success

def run_api_verification():
    """Run the API verification script."""
    success, output = run_command("python tests/verify_api_endpoints.py", "API verification script")
    return success, output

def main():
    print("\n" + "=" * 80)
    print("CryoProtect v2 - API Fixes Runner")
    print("=" * 80)
    
    # Step 1: Fix database tables
    print("\nStep 1: Fixing database tables...")
    if not fix_database_tables():
        print("Failed to fix database tables. Check the logs for details.")
        return 1
    
    # Step 2: Run API verification
    print("\nStep 2: Verifying API endpoints...")
    success, output = run_api_verification()
    if not success:
        print("API verification failed. Check the logs for details.")
        return 1
    
    # Step 3: Generate summary report
    print("\nStep 3: Generating summary report...")
    
    # Try to parse the verification report
    report_path = None
    for line in output.splitlines():
        if "Report saved to:" in line:
            report_path = line.split("Report saved to:")[1].strip()
            break
    
    if report_path and os.path.exists(report_path):
        try:
            with open(report_path, 'r') as f:
                report = json.load(f)
            
            print("\n" + "=" * 60)
            print("API Fixes Summary Report")
            print("=" * 60)
            print(f"Status: {report['status']}")
            print(f"Total Endpoints: {report['summary']['total']}")
            print(f"Implemented: {report['summary']['implemented']} ({report['summary']['implemented']/report['summary']['total']*100:.1f}%)")
            print(f"Functional: {report['summary']['functional']} ({report['summary']['functional']/report['summary']['total']*100:.1f}%)")
            print(f"Discrepancies: {report['summary']['discrepancies']}")
            
            if report['status'] == "SUCCESS":
                print("\nAll API endpoints are now functional!")
            else:
                print("\nSome API endpoints still have issues. Check the verification report for details.")
        except Exception as e:
            logger.error(f"Error parsing verification report: {e}")
            print("\nCould not parse verification report. Check the logs for details.")
    else:
        print("\nCould not find verification report. Check the logs for details.")
    
    print("\n" + "=" * 60)
    print("API Fixes Complete")
    print("=" * 60)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
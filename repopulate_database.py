#!/usr/bin/env python3
"""
CryoProtect v2 - Database Repopulation Script

This script runs all population scripts in sequence to repopulate the database
with the correct plural table names after the update_table_names_comprehensive.py
script has been run.

Usage:
    python repopulate_database.py [--dry-run]

Author: Claude
Date: April 18, 2025
"""

import os
import sys
import time
import argparse
import logging
import subprocess
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("repopulate_database.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Population scripts in the order they should be run
POPULATION_SCRIPTS = [
    "populate_molecules.py",
    "populate_mixtures.py",
    "populate_predictions.py",
    "populate_experiments.py"
]

def run_script(script_path, dry_run=False):
    """Run a script and return True if successful, False otherwise."""
    try:
        if not os.path.exists(script_path):
            logger.error(f"Script not found: {script_path}")
            return False
        
        # Build the command
        if dry_run:
            cmd = [sys.executable, script_path, "--dry-run"]
        else:
            cmd = [sys.executable, script_path]
        
        logger.info(f"Running: {' '.join(cmd)}")
        print(f"\nRunning: {' '.join(cmd)}")
        print("-" * 60)
        
        # Run the script
        start_time = time.time()
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        elapsed_time = time.time() - start_time
        
        # Log the output
        if result.stdout:
            for line in result.stdout.split('\n'):
                if line.strip():
                    logger.info(f"[{script_path}] {line}")
        
        if result.stderr:
            for line in result.stderr.split('\n'):
                if line.strip():
                    logger.warning(f"[{script_path}] {line}")
        
        logger.info(f"Successfully ran {script_path} in {elapsed_time:.2f} seconds")
        print(f"Successfully ran {script_path} in {elapsed_time:.2f} seconds")
        return True
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running {script_path}: {e}")
        if e.stdout:
            logger.info(f"STDOUT: {e.stdout}")
        if e.stderr:
            logger.error(f"STDERR: {e.stderr}")
        print(f"Error running {script_path}: {e}")
        return False
    
    except Exception as e:
        logger.error(f"Unexpected error running {script_path}: {str(e)}")
        print(f"Unexpected error running {script_path}: {str(e)}")
        return False

def verify_api_endpoints():
    """Verify API endpoints are working after population."""
    try:
        # Check if the verification script exists
        verify_script = "test_all_api_endpoints.py"
        if not os.path.exists(verify_script):
            logger.warning(f"API verification script not found: {verify_script}")
            print(f"Warning: API verification script not found: {verify_script}")
            return False
        
        # Run the verification script
        print("\nVerifying API endpoints...")
        logger.info("Running API endpoint verification")
        
        cmd = [sys.executable, verify_script]
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        
        # Log the output
        if result.stdout:
            for line in result.stdout.split('\n'):
                if line.strip():
                    logger.info(f"[{verify_script}] {line}")
        
        if result.stderr:
            for line in result.stderr.split('\n'):
                if line.strip():
                    logger.warning(f"[{verify_script}] {line}")
        
        logger.info("API endpoint verification complete")
        print("API endpoint verification complete")
        return True
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error verifying API endpoints: {e}")
        if e.stdout:
            logger.info(f"STDOUT: {e.stdout}")
        if e.stderr:
            logger.error(f"STDERR: {e.stderr}")
        print(f"Error verifying API endpoints: {e}")
        return False
    
    except Exception as e:
        logger.error(f"Unexpected error verifying API endpoints: {str(e)}")
        print(f"Unexpected error verifying API endpoints: {str(e)}")
        return False

def main():
    """Run all population scripts in sequence."""
    parser = argparse.ArgumentParser(description="Repopulate the database with the correct plural table names.")
    parser.add_argument("--dry-run", action="store_true", help="Run population scripts in dry run mode")
    parser.add_argument("--skip-verify", action="store_true", help="Skip API endpoint verification")
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Database Repopulation")
    print("=" * 80)
    print("\nStarting repopulation at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Starting database repopulation")
    
    if args.dry_run:
        print("\nRunning in DRY RUN mode (no data will be changed)")
        logger.info("Running in DRY RUN mode")
    
    # Run each script in sequence
    successful_scripts = []
    failed_scripts = []
    
    for script in POPULATION_SCRIPTS:
        if run_script(script, args.dry_run):
            successful_scripts.append(script)
        else:
            failed_scripts.append(script)
            logger.error(f"Failed to run {script}, stopping sequence")
            break
    
    # Verify API endpoints if not in dry run mode and no failures
    api_verification = None
    if not args.dry_run and not failed_scripts and not args.skip_verify:
        api_verification = verify_api_endpoints()
    
    # Summary
    print("\n" + "=" * 60)
    print("Repopulation Summary")
    print("=" * 60)
    
    if successful_scripts:
        print(f"\nSuccessfully ran {len(successful_scripts)} scripts:")
        for script in successful_scripts:
            print(f"- {script}")
    
    if failed_scripts:
        print(f"\nFailed to run {len(failed_scripts)} scripts:")
        for script in failed_scripts:
            print(f"- {script}")
    
    if api_verification is not None:
        if api_verification:
            print("\nAPI endpoint verification: SUCCESS")
        else:
            print("\nAPI endpoint verification: FAILED")
    
    print("\nFinished at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Database repopulation complete")
    
    if failed_scripts:
        return 1
    elif api_verification is False:
        return 2
    else:
        return 0

if __name__ == "__main__":
    sys.exit(main())

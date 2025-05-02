#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Table Names and Repopulate Database

This script runs the comprehensive table name update followed by database repopulation
to fix the issue with singular vs. plural table names.

Usage:
    python fix_table_names_and_repopulate.py [--dry-run]

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
        logging.FileHandler("fix_table_names_and_repopulate.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_script(script_path, args=None):
    """Run a script and return True if successful, False otherwise."""
    try:
        if not os.path.exists(script_path):
            logger.error(f"Script not found: {script_path}")
            return False
        
        # Build the command
        cmd = [sys.executable, script_path]
        if args:
            cmd.extend(args)
        
        logger.info(f"Running: {' '.join(cmd)}")
        print(f"\nRunning: {' '.join(cmd)}")
        print("-" * 60)
        
        # Run the script
        start_time = time.time()
        result = subprocess.run(cmd, check=True, text=True)
        elapsed_time = time.time() - start_time
        
        logger.info(f"Successfully ran {script_path} in {elapsed_time:.2f} seconds")
        print(f"Successfully ran {script_path} in {elapsed_time:.2f} seconds")
        return True
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running {script_path}: {e}")
        print(f"Error running {script_path}: {e}")
        return False
    
    except Exception as e:
        logger.error(f"Unexpected error running {script_path}: {str(e)}")
        print(f"Unexpected error running {script_path}: {str(e)}")
        return False

def main():
    """Run the table name update and database repopulation scripts."""
    parser = argparse.ArgumentParser(description="Fix table names and repopulate the database.")
    parser.add_argument("--dry-run", action="store_true", help="Simulate changes without making them")
    parser.add_argument("--update-only", action="store_true", help="Only update table names, don't repopulate")
    parser.add_argument("--repopulate-only", action="store_true", help="Only repopulate, don't update table names")
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Fix Table Names and Repopulate Database")
    print("=" * 80)
    print("\nStarting process at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Starting fix and repopulate process")
    
    if args.dry_run:
        print("\nRunning in DRY RUN mode (no data will be changed)")
        logger.info("Running in DRY RUN mode")
    
    # Step 1: Update table names in scripts
    if not args.repopulate_only:
        print("\n" + "=" * 60)
        print("STEP 1: Update Table Names in Scripts")
        print("=" * 60)
        
        update_script = "update_table_names_comprehensive.py"
        update_args = ["--dry-run"] if args.dry_run else []
        
        if not run_script(update_script, update_args):
            logger.error("Failed to update table names")
            print("\nError: Failed to update table names. Check the logs for details.")
            return 1
    
    # Step 2: Repopulate database
    if not args.update_only and not args.dry_run:
        print("\n" + "=" * 60)
        print("STEP 2: Repopulate Database with Correct Table Names")
        print("=" * 60)
        
        repopulate_script = "repopulate_database.py"
        repopulate_args = ["--dry-run"] if args.dry_run else []
        
        if not run_script(repopulate_script, repopulate_args):
            logger.error("Failed to repopulate database")
            print("\nError: Failed to repopulate database. Check the logs for details.")
            return 2
    
    # Summary
    print("\n" + "=" * 60)
    print("Process Summary")
    print("=" * 60)
    
    if args.repopulate_only:
        print("\nSkipped table name updates (--repopulate-only)")
    elif args.dry_run:
        print("\nTable name updates were simulated (--dry-run)")
    else:
        print("\nTable names were successfully updated in scripts")
    
    if args.update_only:
        print("\nSkipped database repopulation (--update-only)")
    elif args.dry_run:
        print("\nDatabase repopulation was simulated (--dry-run)")
    else:
        print("\nDatabase was successfully repopulated with correct table names")
    
    print("\nFinished at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Fix and repopulate process complete")
    
    # Next steps
    print("\n" + "=" * 60)
    print("Next Steps")
    print("=" * 60)
    
    if args.dry_run:
        print("\n1. Review the logs to ensure everything looks correct")
        print("2. Run the script again without the --dry-run flag to apply the changes")
    else:
        print("\n1. Test the API endpoints to ensure they're working correctly")
        print("2. Check the database to confirm data is in the correct tables")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

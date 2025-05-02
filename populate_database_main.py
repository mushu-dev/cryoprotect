#!/usr/bin/env python3
"""
CryoProtect v2 - Database Population Main Script

This script runs all the database population scripts in the correct order to ensure
proper relationships between tables. It first audits the database to identify empty
tables, then populates them with scientifically accurate data.

Usage:
    python populate_database_main.py [--dry-run] [--skip-audit]

Environment variables required (from .env):
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
"""

import os
import sys
import argparse
import logging
import subprocess
import time
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_population_main.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Population scripts in dependency order
POPULATION_SCRIPTS = [
    "supabase_database_audit.py",
    "populate_molecules.py",
    "populate_mixtures.py",
    "populate_experiments.py",
    "populate_predictions.py"
]

def run_script(script_path, dry_run=False):
    """Run a Python script and return its exit code."""
    cmd = [sys.executable, script_path]
    if dry_run:
        cmd.append("--dry-run")
    
    logger.info(f"Running {' '.join(cmd)}")
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - start_time
        logger.info(f"Script {script_path} completed successfully in {elapsed:.2f} seconds")
        return 0
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        logger.error(f"Script {script_path} failed after {elapsed:.2f} seconds with exit code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return e.returncode

def main():
    parser = argparse.ArgumentParser(description="Run all database population scripts in the correct order.")
    parser.add_argument("--dry-run", action="store_true", help="Run all scripts in dry-run mode")
    parser.add_argument("--skip-audit", action="store_true", help="Skip the database audit step")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Database Population")
    start_time = time.time()
    
    # Run scripts in order
    for i, script in enumerate(POPULATION_SCRIPTS):
        # Skip audit if requested
        if args.skip_audit and script == "supabase_database_audit.py":
            logger.info("Skipping database audit as requested")
            continue
        
        logger.info(f"Step {i+1}/{len(POPULATION_SCRIPTS)}: Running {script}")
        exit_code = run_script(script, args.dry_run)
        
        if exit_code != 0:
            logger.error(f"Script {script} failed with exit code {exit_code}. Stopping.")
            return exit_code
        
        logger.info(f"Step {i+1}/{len(POPULATION_SCRIPTS)} completed successfully")
    
    total_elapsed = time.time() - start_time
    logger.info(f"All database population steps completed successfully in {total_elapsed:.2f} seconds")
    
    # Print summary
    print("\n" + "="*60)
    print("CryoProtect v2 Database Population Summary")
    print("="*60)
    print(f"Total time: {total_elapsed:.2f} seconds")
    print(f"Mode: {'Dry run (no changes made)' if args.dry_run else 'Live run (database updated)'}")
    print(f"Audit: {'Skipped' if args.skip_audit else 'Performed'}")
    print("\nThe following scripts were executed:")
    for script in POPULATION_SCRIPTS:
        if args.skip_audit and script == "supabase_database_audit.py":
            print(f"  - {script} (skipped)")
        else:
            print(f"  - {script}")
    print("\nCheck individual log files for detailed information:")
    print("  - database_audit.log")
    print("  - molecule_population.log")
    print("  - mixture_population.log")
    print("  - experiment_population.log")
    print("  - prediction_population.log")
    print("="*60)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
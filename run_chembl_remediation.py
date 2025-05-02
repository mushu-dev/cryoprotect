#!/usr/bin/env python3
"""
ChEMBL Remediation Orchestration Script

This script runs the entire ChEMBL remediation process:
1. Add chembl_id column to molecules table
2. Import ChEMBL compounds
3. Reconcile property values
4. Verify the remediation
"""

import os
import sys
import logging
import argparse
import subprocess
import time
from datetime import datetime

# Set up logging
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"logs/chembl_remediation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_script(script_path, args=None, desc=None):
    """
    Run a Python script and return success status.
    
    Args:
        script_path: Path to the script
        args: List of command-line arguments
        desc: Description of the step
        
    Returns:
        True if script executed successfully, False otherwise
    """
    if desc:
        logger.info(f"Running: {desc}")
    
    cmd = [sys.executable, script_path]
    if args:
        cmd.extend(args)
    
    logger.info(f"Executing: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Script completed successfully: {script_path}")
        if result.stdout:
            for line in result.stdout.splitlines():
                logger.info(f"STDOUT: {line}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Script failed: {script_path}")
        logger.error(f"Return code: {e.returncode}")
        if e.stdout:
            for line in e.stdout.splitlines():
                logger.info(f"STDOUT: {line}")
        if e.stderr:
            for line in e.stderr.splitlines():
                logger.error(f"STDERR: {line}")
        return False
    except Exception as e:
        logger.error(f"Error running script {script_path}: {str(e)}")
        return False

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run ChEMBL remediation process")
    parser.add_argument("--skip-migration", action="store_true", help="Skip the migration step")
    parser.add_argument("--skip-import", action="store_true", help="Skip the import step")
    parser.add_argument("--skip-reconciliation", action="store_true", help="Skip the reconciliation step")
    parser.add_argument("--skip-verification", action="store_true", help="Skip the verification step")
    parser.add_argument("--dry-run", action="store_true", help="Run reconciliation in dry-run mode")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info("Starting ChEMBL remediation process")
    
    # Record start time
    start_time = time.time()
    
    # Step 1: Run the fixed orchestration script for migration
    if not args.skip_migration:
        logger.info("STEP 1: Running database migration to add chembl_id column")
        if not run_script(
            "chembl_remediation_main_fixed.py",
            desc="Database migration"
        ):
            logger.error("Migration step failed - stopping")
            return 1
    else:
        logger.info("Skipping migration step")
    
    # Step 2: Import ChEMBL compounds
    if not args.skip_import:
        logger.info("STEP 2: Importing ChEMBL compounds")
        if not run_script(
            "ChEMBL_Integrated_Import.py",
            args=["--limit", "2500"],
            desc="ChEMBL data import"
        ):
            logger.error("Import step failed - stopping")
            return 1
    else:
        logger.info("Skipping import step")
    
    # Step 3: Reconcile property values
    if not args.skip_reconciliation:
        logger.info("STEP 3: Reconciling property values")
        reconcile_args = []
        if args.dry_run:
            reconcile_args.append("--dry-run")
        if not run_script(
            "reconcile_chembl_properties.py",
            args=reconcile_args,
            desc="Property reconciliation"
        ):
            logger.error("Reconciliation step failed - stopping")
            return 1
    else:
        logger.info("Skipping reconciliation step")
    
    # Step 4: Verify the remediation
    if not args.skip_verification:
        logger.info("STEP 4: Verifying the remediation")
        if not run_script(
            "verify_chembl_reconciliation.py",
            args=["--generate-html"],
            desc="Verification"
        ):
            logger.error("Verification step failed")
            # We continue even if verification fails, so we can see the report
    else:
        logger.info("Skipping verification step")
    
    # Record end time and calculate duration
    end_time = time.time()
    duration_seconds = int(end_time - start_time)
    minutes, seconds = divmod(duration_seconds, 60)
    
    logger.info("="*50)
    logger.info("CHEMBL REMEDIATION PROCESS COMPLETE")
    logger.info("="*50)
    logger.info(f"Total runtime: {minutes} minutes, {seconds} seconds")
    logger.info("Check the 'reports' directory for verification results")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
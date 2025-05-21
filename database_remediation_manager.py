#!/usr/bin/env python3
"""
CryoProtect v2 - Database Remediation Manager

This script orchestrates the complete database remediation process for CryoProtect v2,
following a 5-phase approach to address all critical issues:

Phase 1: Security Lockdown (HIGHEST PRIORITY)
- Enable RLS on all tables
- Restrict anonymous access
- Create basic RLS policies
- Add indexes for policy columns

Phase 2: Schema Standardization
- Standardize table names to plural form
- Update all foreign key references
- Document all changes in migration scripts

Phase 3: Relationship Remediation
- Resolve fan traps by creating junction tables
- Update application code to use new structure
- Add appropriate indexes to junction tables

Phase 4: Data Migration
- Consolidate duplicate tables
- Migrate data to normalized structure
- Validate data integrity after migration

Phase 5: Performance Optimization
- Add missing indexes
- Create optimized views for common queries
- Implement connection pooling optimizations

Usage:
    python database_remediation_manager.py [--phase PHASE] [--dry-run] [--verify-only]
"""

import os
import sys
import json
import time
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
        logging.FileHandler(f"database_remediation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Define the phases and their corresponding scripts
PHASES = {
    1: {
        "name": "Security Lockdown",
        "description": "Enable RLS on all tables, restrict anonymous access, create basic RLS policies",
        "scripts": [
            {"script": "implement_security.py", "args": []},
            {"script": "apply_service_role_rls.py", "args": []}
        ],
        "verification": "verify_rls_policies"
    },
    2: {
        "name": "Schema Standardization",
        "description": "Standardize table names to plural form, update all foreign key references",
        "scripts": [
            {"script": "standardize_schema.py", "args": []}
        ],
        "verification": "verify_schema_standardization"
    },
    3: {
        "name": "Relationship Remediation",
        "description": "Resolve fan traps by creating junction tables, add appropriate indexes",
        "scripts": [
            {"script": "fix_relationships.py", "args": []}
        ],
        "verification": "verify_relationships"
    },
    4: {
        "name": "Data Migration",
        "description": "Consolidate duplicate tables, migrate data to normalized structure",
        "scripts": [
            {"script": "complete_database_remediation.py", "args": []}
        ],
        "verification": "verify_database_integrity"
    },
    5: {
        "name": "Performance Optimization",
        "description": "Add missing indexes, create optimized views for common queries",
        "scripts": [
            {"script": "apply_performance_indexes.bat" if os.name == "nt" else "apply_performance_indexes.sh", "args": []}
        ],
        "verification": "verify_performance_indexes"
    }
}

def run_script(script, args=None, dry_run=False):
    """Run a script with the given arguments."""
    if args is None:
        args = []
    
    # Add dry-run argument if applicable
    if dry_run and script.endswith(".py"):
        args.append("--dry-run")
    
    # Construct the command
    if script.endswith(".py"):
        cmd = [sys.executable, script] + args
    elif script.endswith(".bat") and os.name == "nt":
        cmd = [script] + args
    elif script.endswith(".sh") and os.name != "nt":
        cmd = ["bash", script] + args
    else:
        logger.error(f"Unsupported script type: {script}")
        return False
    
    # Log the command
    logger.info(f"Running: {' '.join(cmd)}")
    
    if dry_run:
        logger.info(f"DRY RUN: Would execute {' '.join(cmd)}")
        return True
    
    # Execute the command
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Script output: {result.stdout}")
        if result.stderr:
            logger.warning(f"Script errors: {result.stderr}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Script failed with exit code {e.returncode}")
        logger.error(f"Output: {e.stdout}")
        logger.error(f"Error: {e.stderr}")
        return False

def verify_rls_policies(dry_run=False):
    """Verify that RLS policies are correctly applied."""
    if dry_run:
        logger.info("DRY RUN: Would verify RLS policies")
        return True
    
    return run_script("implement_security.py", ["--verify-only"])

def verify_schema_standardization(dry_run=False):
    """Verify that schema standardization was successful."""
    if dry_run:
        logger.info("DRY RUN: Would verify schema standardization")
        return True
    
    # Check if tables have been renamed to plural form
    try:
        from supabase import create_client
        
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY")
        
        if not supabase_url or not supabase_key:
            logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
            return False
        
        supabase = create_client(supabase_url, supabase_key)
        
        # Check for plural table names
        plural_tables = ["molecules", "mixtures", "predictions", "experiments", 
                         "experiment_properties", "mixture_components", 
                         "calculation_methods", "property_types", "projects", "teams"]
        
        all_exist = True
        for table in plural_tables:
            try:
                # Try to select a single row from the table
                response = supabase.table(table).select("*").limit(1).execute()
                logger.info(f"Table '{table}' exists")
            except Exception as e:
                logger.error(f"Table '{table}' does not exist: {str(e)}")
                all_exist = False
        
        return all_exist
    except Exception as e:
        logger.error(f"Error verifying schema standardization: {str(e)}")
        return False

def verify_relationships(dry_run=False):
    """Verify that relationship fixes were successful."""
    if dry_run:
        logger.info("DRY RUN: Would verify relationships")
        return True
    
    return run_script("fix_relationships.py", ["--verify-only"])

def verify_database_integrity(dry_run=False):
    """Verify database integrity after data migration."""
    if dry_run:
        logger.info("DRY RUN: Would verify database integrity")
        return True
    
    return run_script("verify_database_integrity.py")

def verify_performance_indexes(dry_run=False):
    """Verify that performance indexes were created."""
    if dry_run:
        logger.info("DRY RUN: Would verify performance indexes")
        return True
    
    script = "verify_performance_indexes.bat" if os.name == "nt" else "verify_performance_indexes.sh"
    return run_script(script)

def run_phase(phase_number, dry_run=False, verify_only=False):
    """Run a specific phase of the remediation process."""
    if phase_number not in PHASES:
        logger.error(f"Invalid phase number: {phase_number}")
        return False
    
    phase = PHASES[phase_number]
    
    logger.info(f"=== Phase {phase_number}: {phase['name']} ===")
    logger.info(phase['description'])
    
    if verify_only:
        logger.info(f"Verifying Phase {phase_number}...")
        verification_func = globals()[phase['verification']]
        return verification_func(dry_run)
    
    # Run each script in the phase
    success = True
    for script_info in phase['scripts']:
        script = script_info['script']
        args = script_info['args']
        
        logger.info(f"Running {script}...")
        if not run_script(script, args, dry_run):
            logger.error(f"Failed to run {script}")
            success = False
            break
    
    # Verify the phase
    if success and not dry_run:
        logger.info(f"Verifying Phase {phase_number}...")
        verification_func = globals()[phase['verification']]
        if not verification_func(dry_run):
            logger.error(f"Verification failed for Phase {phase_number}")
            success = False
    
    return success

def run_all_phases(start_phase=1, end_phase=5, dry_run=False):
    """Run all phases of the remediation process."""
    for phase_number in range(start_phase, end_phase + 1):
        logger.info(f"\n{'=' * 80}")
        logger.info(f"Starting Phase {phase_number}: {PHASES[phase_number]['name']}")
        logger.info(f"{'=' * 80}\n")
        
        if not run_phase(phase_number, dry_run):
            logger.error(f"Phase {phase_number} failed. Stopping remediation process.")
            return False
        
        logger.info(f"Phase {phase_number} completed successfully.")
    
    return True

def main():
    """Main function to run the database remediation process."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 Database Remediation Manager")
    parser.add_argument("--phase", type=int, choices=range(1, 6), help="Run a specific phase (1-5)")
    parser.add_argument("--start-phase", type=int, choices=range(1, 6), default=1, help="Start from a specific phase (1-5)")
    parser.add_argument("--end-phase", type=int, choices=range(1, 6), default=5, help="End at a specific phase (1-5)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    parser.add_argument("--verify-only", action="store_true", help="Only verify the specified phase(s)")
    args = parser.parse_args()
    
    # Print banner
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Database Remediation Manager")
    print("=" * 80 + "\n")
    
    if args.dry_run:
        logger.info("DRY RUN MODE: No changes will be made to the database")
    
    if args.verify_only:
        logger.info("VERIFY ONLY MODE: Only verification will be performed")
    
    # Run a specific phase or all phases
    if args.phase:
        success = run_phase(args.phase, args.dry_run, args.verify_only)
    else:
        if args.verify_only:
            # Verify all phases from start_phase to end_phase
            success = True
            for phase_number in range(args.start_phase, args.end_phase + 1):
                phase_success = run_phase(phase_number, args.dry_run, True)
                if not phase_success:
                    logger.error(f"Verification failed for Phase {phase_number}")
                    success = False
        else:
            # Run all phases from start_phase to end_phase
            success = run_all_phases(args.start_phase, args.end_phase, args.dry_run)
    
    # Print summary
    print("\n" + "=" * 80)
    print("CryoProtect v2 Database Remediation Summary")
    print("=" * 80)
    
    if success:
        print("\nStatus: SUCCESS")
        print("All specified phases completed successfully.")
    else:
        print("\nStatus: FAILED")
        print("One or more phases failed. Check the log for details.")
    
    print("\nPhases:")
    for phase_number, phase in PHASES.items():
        status = "Not Run"
        if args.phase:
            if phase_number == args.phase:
                status = "Success" if success else "Failed"
        else:
            if args.start_phase <= phase_number <= args.end_phase:
                status = "Success" if success else "Failed"
        
        print(f"  {phase_number}. {phase['name']}: {status}")
    
    print("\nFor detailed information, check the log file.")
    print("=" * 80 + "\n")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
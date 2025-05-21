#!/usr/bin/env python3
"""
This script runs the database cleanup operations in the correct order.
It uses environment variables for database connection and provides
a simple CLI interface with safety checks.
"""

import os
import sys
import argparse
import subprocess
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def run_script(script_name, dry_run=False):
    """Run a Python script with proper environment variables."""
    # Get environment variables
    db_host = os.getenv('SUPABASE_DB_HOST')
    db_port = os.getenv('SUPABASE_DB_PORT', '5432')
    db_name = os.getenv('SUPABASE_DB_NAME', 'postgres')
    db_user = os.getenv('SUPABASE_DB_USER')
    db_password = os.getenv('SUPABASE_DB_PASSWORD')
    
    if not all([db_host, db_port, db_name, db_user, db_password]):
        print("ERROR: Missing database connection details in environment variables.")
        print("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        return False
    
    # Build command
    cmd = [sys.executable, script_name]
    
    # Add arguments
    if dry_run:
        cmd.append('--dry-run')
    
    # Set environment variables for the subprocess
    env = os.environ.copy()
    
    # Run the script
    print(f"Running {script_name}...")
    try:
        result = subprocess.run(cmd, env=env, check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Script {script_name} failed with exit code {e.returncode}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Run database cleanup operations")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show what would be done without making changes")
    parser.add_argument("--step", choices=["chembl", "mixtures", "classification", "all"],
                        default="all", help="Specific step to run")
    args = parser.parse_args()
    
    print("CryoProtect Database Remediation Tool")
    print("======================================")
    
    success = True
    
    # First step: Clean ChEMBL identifiers
    if args.step in ["chembl", "all"]:
        print("\n### Step 1: Cleaning ChEMBL identifiers ###")
        if not run_script("clean_chembl_identifiers.py", args.dry_run):
            print("WARNING: ChEMBL identifier cleanup failed or was skipped.")
            success = False
    
    # Second step: Complete mixtures
    if args.step in ["mixtures", "all"] and success:
        print("\n### Step 2: Completing or removing placeholder mixtures ###")
        if not run_script("complete_mixtures.py", args.dry_run):
            print("WARNING: Mixture completion failed or was skipped.")
            success = False
    
    # Third step: Add cryoprotectant classifications
    if args.step in ["classification", "all"] and success:
        print("\n### Step 3: Adding cryoprotectant classifications ###")
        if not run_script("add_cryoprotectant_classification.py", args.dry_run):
            print("WARNING: Classification addition failed or was skipped.")
            success = False
    
    # Summary
    print("\n### Summary ###")
    if success:
        print("All database cleanup steps completed successfully!")
    else:
        print("Database cleanup completed with warnings or errors.")
        print("Please check the output above for details.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
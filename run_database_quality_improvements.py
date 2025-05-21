#!/usr/bin/env python3
"""
This script runs all database quality improvement tools in sequence:
1. Consolidate duplicate molecules
2. Standardize SMILES notation
3. Assign users to teams

It provides a unified interface for improving database quality.
"""

import os
import sys
import argparse
import subprocess
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"database_quality_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_command(cmd, description):
    """Run a command and log its output."""
    logger.info(f"Running: {description}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Success: {description}")
        logger.debug(f"Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed: {description}")
        logger.error(f"Error: {e}")
        logger.error(f"Output: {e.stdout}")
        logger.error(f"Error Output: {e.stderr}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Run all database quality improvement tools.")
    parser.add_argument("--dry-run", action="store_true", help="Run all tools in dry-run mode")
    parser.add_argument("--skip-consolidation", action="store_true", help="Skip molecule consolidation")
    parser.add_argument("--skip-smiles", action="store_true", help="Skip SMILES standardization")
    parser.add_argument("--skip-teams", action="store_true", help="Skip team assignments")
    parser.add_argument("--team-role", default="member", choices=["member", "admin", "owner"], 
                        help="Role to assign to users (default: member)")
    parser.add_argument("--batch-size", type=int, default=500, 
                        help="Batch size for SMILES standardization (default: 500)")
    parser.add_argument("--max-molecules", type=int, default=None, 
                        help="Maximum number of molecules to process (default: all)")
    args = parser.parse_args()
    
    print("\n" + "=" * 70)
    print("DATABASE QUALITY IMPROVEMENT PROCESS")
    print("=" * 70)
    print(f"Mode: {'DRY RUN (no changes will be made)' if args.dry_run else 'EXECUTION (changes will be made)'}")
    print("=" * 70 + "\n")
    
    # List of scripts to run
    scripts = []
    
    # 1. Molecule Consolidation
    if not args.skip_consolidation:
        cmd = ["python", "run_duplicate_molecule_consolidation.py"]
        if args.dry_run:
            cmd.append("--dry-run")
        scripts.append((cmd, "Molecule Consolidation"))
    
    # 2. SMILES Standardization
    if not args.skip_smiles:
        cmd = ["python", "standardize_smiles_notation.py"]
        if args.dry_run:
            cmd.append("--dry-run")
        if args.batch_size:
            cmd.extend(["--batch-size", str(args.batch_size)])
        if args.max_molecules:
            cmd.extend(["--max-molecules", str(args.max_molecules)])
        scripts.append((cmd, "SMILES Standardization"))
    
    # 3. Team Assignments
    if not args.skip_teams:
        cmd = ["python", "assign_users_to_teams.py"]
        if args.dry_run:
            cmd.append("--dry-run")
        cmd.extend(["--role", args.team_role])
        scripts.append((cmd, "Team Assignment"))
    
    # Run all scripts
    results = {}
    for cmd, description in scripts:
        success = run_command(cmd, description)
        results[description] = "Success" if success else "Failed"
    
    # Print summary
    print("\n" + "=" * 70)
    print("DATABASE QUALITY IMPROVEMENT SUMMARY")
    print("=" * 70)
    for description, status in results.items():
        print(f"{description}: {status}")
    print("=" * 70)
    
    # Return success only if all scripts succeeded
    return all(status == "Success" for status in results.values())

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
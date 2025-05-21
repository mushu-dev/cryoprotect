#!/usr/bin/env python3
"""
This script runs the duplicate molecule consolidation process.
It can be run in dry-run mode to see what would be done without making changes.
"""

import os
import sys
import logging
import argparse
from datetime import datetime
import consolidate_duplicate_molecules

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the duplicate molecule consolidation process.")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry-run mode without making changes")
    args = parser.parse_args()
    
    print(f"Running duplicate molecule consolidation {'in dry-run mode' if args.dry_run else 'in execution mode'}")
    
    # Initialize a todo list
    if args.dry_run:
        print("\nDRY RUN - No changes will be made to the database")
        
    # Run the main function from the consolidation script
    consolidate_duplicate_molecules.main(dry_run=args.dry_run)
    
    print("\nDuplicate molecule consolidation process complete!")
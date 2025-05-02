#!/usr/bin/env python3
"""
CryoProtect v2 - Run Database Remediation

This script runs the database integrity remediation process and displays the results.

Usage:
    python run_database_remediation.py [--dry-run]
"""

import os
import sys
import json
import argparse
import subprocess
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Run database integrity remediation and display results.")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry-run mode (no changes made)")
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Database Integrity Remediation")
    print("=" * 80)
    
    # Check if the remediate_database_integrity.py script exists
    if not os.path.exists("remediate_database_integrity.py"):
        print("ERROR: remediate_database_integrity.py script not found.")
        return 1
    
    # Run the remediation script
    print("\nRunning database remediation...")
    cmd = [sys.executable, "remediate_database_integrity.py"]
    if args.dry_run:
        cmd.append("--dry-run")
        print("(DRY RUN MODE - No changes will be made)")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Print the output from the remediation script
        if result.stdout:
            print("\nRemediation Output:")
            print(result.stdout)
        
        if result.stderr:
            print("\nRemediation Errors:")
            print(result.stderr)
        
        # Check if the remediation summary file exists
        if os.path.exists("database_remediation_summary.json"):
            try:
                with open("database_remediation_summary.json", "r") as f:
                    summary = json.load(f)
                
                print("\n" + "=" * 80)
                print("Remediation Summary")
                print("=" * 80)
                print(f"Status: {summary.get('status', 'UNKNOWN')}")
                print(f"Timestamp: {summary.get('timestamp', datetime.now().isoformat())}")
                
                if "actions" in summary:
                    print("\nActions Taken:")
                    for action in summary["actions"]:
                        print(f"  {action['table']}:")
                        print(f"    - Orphaned records found: {action['orphaned_records_found']}")
                        print(f"    - Orphaned records removed: {action['orphaned_records_removed']}")
                
                if "verification_result" in summary:
                    print(f"\nVerification Result: {summary['verification_result']}")
                
                # Check if there's a database integrity report after remediation
                if os.path.exists("database_integrity_report.json") and not args.dry_run:
                    try:
                        with open("database_integrity_report.json", "r") as f:
                            integrity_report = json.load(f)
                        
                        print("\nDatabase Integrity Report:")
                        print(f"  Status: {integrity_report.get('status', 'UNKNOWN')}")
                        print(f"  Empty Tables: {len(integrity_report.get('empty_tables', []))}")
                        print(f"  Foreign Key Issues: {len(integrity_report.get('foreign_key_issues', []))}")
                        print(f"  Scientific Issues: {len(integrity_report.get('scientific_issues', []))}")
                        
                        if integrity_report.get('empty_tables'):
                            print("\n  Empty Tables:")
                            for table in integrity_report['empty_tables']:
                                print(f"    - {table}")
                        
                        if integrity_report.get('foreign_key_issues'):
                            print("\n  Foreign Key Issues:")
                            for issue in integrity_report['foreign_key_issues']:
                                print(f"    - {issue}")
                        
                        if integrity_report.get('scientific_issues'):
                            print("\n  Scientific Issues:")
                            for issue in integrity_report['scientific_issues']:
                                print(f"    - {issue}")
                    except Exception as e:
                        print(f"\nError parsing database integrity report: {e}")
            except Exception as e:
                print(f"\nError parsing remediation summary: {e}")
        else:
            print("\nNo remediation summary file found.")
        
        print("\n" + "=" * 80)
        if args.dry_run:
            print("Remediation completed in DRY RUN mode. No changes were made.")
            print("Run without --dry-run to apply the changes.")
        elif result.returncode == 0:
            print("Remediation completed successfully!")
        else:
            print("Remediation completed with errors. Check the logs for details.")
        print("=" * 80 + "\n")
        
        return result.returncode
    except Exception as e:
        print(f"\nError running remediation script: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
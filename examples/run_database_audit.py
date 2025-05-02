#!/usr/bin/env python3
"""
CryoProtect v2 - Database Audit Example

This script demonstrates how to use the database audit functionality to verify
database resources before data population. It shows how to:

1. Run a complete database audit
2. Handle verification failures
3. Generate and interpret audit reports

Usage:
    python examples/run_database_audit.py

Requirements:
    - Supabase credentials in .env file or environment variables
    - Database schema already applied (migrations/001_initial_schema.sql)
"""

import os
import sys
import json
from datetime import datetime
from pathlib import Path

# Add parent directory to path to import supabase_database_audit
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the database audit module
import supabase_database_audit

def main():
    """Run the database audit and handle the results."""
    print("\n" + "="*60)
    print("CryoProtect v2 - Database Audit Example")
    print("="*60)
    
    # Ensure required directories exist
    Path("logs").mkdir(exist_ok=True)
    Path("reports").mkdir(exist_ok=True)
    
    try:
        # Run the database audit
        print("\nRunning database audit...")
        supabase_database_audit.main()
        
        # If we get here, the audit was successful
        print("\nDatabase audit completed successfully!")
        
        # Load and display the latest audit report
        try:
            with open("reports/database_audit_report_latest.json", "r") as f:
                report = json.load(f)
                
            print("\nAudit Report Summary:")
            print(f"- Timestamp: {report['verification']['timestamp']}")
            print(f"- Tables Verified: {report['verification']['tables_verified']}")
            print(f"- Integrity Verified: {report['verification']['integrity_verified']}")
            print(f"- Total Tables: {report['summary']['total_tables']}")
            print(f"- Populated Tables: {report['summary']['populated_tables']}")
            
            if report['verification']['quota_warnings']:
                print("\nQuota Warnings:")
                for warning in report['verification']['quota_warnings']:
                    print(f"- {warning}")
            
            if report['population_recommendations']:
                print("\nPopulation Recommendations:")
                for rec in sorted(report['population_recommendations'], key=lambda x: x['priority'])[:5]:
                    print(f"- {rec['recommendation']}")
        except Exception as e:
            print(f"\nError reading audit report: {str(e)}")
    
    except SystemExit as e:
        # The audit script called sys.exit()
        if e.code != 0:
            print("\nDatabase audit failed with verification errors.")
            print("Please check the logs for details and fix the issues before proceeding with data population.")
            
            # Try to read the latest report if it exists
            try:
                with open("reports/database_audit_report_latest.json", "r") as f:
                    report = json.load(f)
                
                if 'verification' in report:
                    if not report['verification']['tables_verified']:
                        print("\nMissing Tables:")
                        for table in report['verification']['missing_tables']:
                            print(f"- {table}")
                        print("\nPlease run the initial schema migration to create these tables.")
            except Exception:
                pass
    
    except Exception as e:
        print(f"\nUnexpected error during database audit: {str(e)}")
    
    print("\n" + "="*60)

if __name__ == "__main__":
    main()
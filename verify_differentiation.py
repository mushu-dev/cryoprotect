#!/usr/bin/env python3
"""
Verify differentiation results.

This script:
1. Checks that molecules in differentiate groups have proper differentiation properties
2. Verifies that the differentiation information correctly identifies each molecule
3. Generates a verification report
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Constants
PLAN_FILE = "duplicate_consolidation_plan.json"
DIFFERENTIATION_LOG_FILE = "differentiation_log.json"
VERIFICATION_REPORT_FILE = "differentiation_verification.json"

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    return psycopg2.connect(**db_params)

def load_plan_and_log():
    """Load the consolidation plan and differentiation log."""
    try:
        with open(PLAN_FILE, 'r') as f:
            plan = json.load(f)
        
        if os.path.exists(DIFFERENTIATION_LOG_FILE):
            with open(DIFFERENTIATION_LOG_FILE, 'r') as f:
                log = json.load(f)
        else:
            log = []
        
        return plan, log
    except Exception as e:
        print(f"Error loading plan and log: {e}")
        sys.exit(1)

def verify_differentiation(conn, plan, log):
    """Verify that differentiation was completed successfully."""
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "groups_verified": 0,
        "groups_total": len(plan.get('differentiate', [])),
        "molecules_differentiated": 0,
        "issues": [],
        "groups": {}
    }
    
    differentiate_groups = plan.get('differentiate', [])
    differentiated_group_ids = [entry.get('group_id') for entry in log]
    
    print(f"Verifying {len(differentiate_groups)} differentiate groups")
    print(f"Found {len(differentiated_group_ids)} differentiated groups in log")
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        for group in differentiate_groups:
            group_id = group['group_id']
            print(f"Verifying group {group_id} with {group['molecule_count']} molecules")
            
            # Check if group was differentiated
            if group_id not in differentiated_group_ids:
                verification_report["issues"].append({
                    "group_id": group_id,
                    "issue": "Group not found in differentiation log"
                })
                continue
            
            verification_report["groups"][group_id] = {
                "molecule_count": group['molecule_count'],
                "molecules": [],
                "successfully_differentiated": True
            }
            
            # Check each molecule in the group
            for mol in group['molecules_detail']:
                mol_id = mol['id']
                mol_name = mol['name']
                
                # Check if molecule exists and has differentiation property
                cursor.execute("""
                    SELECT id, name, pubchem_cid, properties
                    FROM molecules
                    WHERE id = %s
                """, (mol_id,))
                
                molecule = cursor.fetchone()
                if not molecule:
                    verification_report["issues"].append({
                        "group_id": group_id,
                        "issue": f"Molecule {mol_id} not found"
                    })
                    verification_report["groups"][group_id]["successfully_differentiated"] = False
                    continue
                
                # Check for differentiation property
                has_diff_property = False
                diff_info = None
                
                if molecule['properties'] and 'differentiation' in molecule['properties']:
                    has_diff_property = True
                    diff_info = molecule['properties']['differentiation']
                
                mol_report = {
                    "id": mol_id,
                    "name": mol_name,
                    "has_differentiation_property": has_diff_property,
                    "differentiation_info": diff_info
                }
                
                verification_report["groups"][group_id]["molecules"].append(mol_report)
                
                if not has_diff_property:
                    verification_report["issues"].append({
                        "group_id": group_id,
                        "issue": f"Molecule {mol_id} doesn't have differentiation property"
                    })
                    verification_report["groups"][group_id]["successfully_differentiated"] = False
                else:
                    verification_report["molecules_differentiated"] += 1
            
            if verification_report["groups"][group_id]["successfully_differentiated"]:
                verification_report["groups_verified"] += 1
    
    return verification_report

def main():
    """Verify differentiation and generate a report."""
    # Load plan and log
    plan, log = load_plan_and_log()
    
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        # Verify differentiation
        verification_report = verify_differentiation(conn, plan, log)
        
        # Save verification report
        with open(VERIFICATION_REPORT_FILE, 'w') as f:
            json.dump(verification_report, f, indent=2, default=str)
        
        # Print summary
        print("\nVerification Summary:")
        print(f"Groups verified: {verification_report['groups_verified']} of {verification_report['groups_total']}")
        print(f"Molecules differentiated: {verification_report['molecules_differentiated']}")
        print(f"Issues found: {len(verification_report['issues'])}")
        
        if verification_report['issues']:
            print("\nIssues:")
            for issue in verification_report['issues']:
                print(f"  Group {issue['group_id']}: {issue['issue']}")
        
        print(f"\nVerification report saved to {VERIFICATION_REPORT_FILE}")
        
    except Exception as e:
        print(f"Error during verification: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
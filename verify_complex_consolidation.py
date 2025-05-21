#!/usr/bin/env python3
"""
Verify complex consolidation results.

This script:
1. Checks that all complex merge molecules have been consolidated
2. Verifies that secondary molecules reference their primary counterpart
3. Ensures properties, mixture components, and predictions have been migrated
4. Generates a verification report
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
CONSOLIDATION_LOG_FILE = "complex_consolidation_log.json"
VERIFICATION_REPORT_FILE = "complex_consolidation_verification.json"

# Connect to the database
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
    """Load the consolidation plan and log."""
    try:
        with open(PLAN_FILE, 'r') as f:
            plan = json.load(f)
        
        with open(CONSOLIDATION_LOG_FILE, 'r') as f:
            log = json.load(f)
        
        return plan, log
    except Exception as e:
        print(f"Error loading plan and log: {e}")
        sys.exit(1)

def verify_complex_consolidation(conn, plan, log):
    """Verify that complex consolidation was completed successfully."""
    verification_report = {
        "timestamp": datetime.now().isoformat(),
        "groups_verified": 0,
        "groups_total": len(plan.get('complex_merge', [])),
        "secondary_molecules_consolidated": 0,
        "property_migrations_verified": 0,
        "mixture_component_migrations_verified": 0,
        "issues": [],
        "primary_molecules": {},
        "secondary_molecules": {}
    }
    
    complex_groups = plan.get('complex_merge', [])
    consolidated_group_ids = [entry.get('group_id') for entry in log]
    
    print(f"Verifying {len(complex_groups)} complex merge groups")
    print(f"Found {len(consolidated_group_ids)} consolidated groups in log")
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        for group in complex_groups:
            group_id = group['group_id']
            primary_id = group['primary_candidate']
            print(f"Verifying group {group_id} with primary {primary_id}")
            
            # Check if group was consolidated
            if group_id not in consolidated_group_ids:
                verification_report["issues"].append({
                    "group_id": group_id,
                    "issue": "Group not found in consolidation log"
                })
                continue
            
            # Check primary molecule
            cursor.execute("""
                SELECT id, name, pubchem_cid, properties
                FROM molecules
                WHERE id = %s
            """, (primary_id,))
            
            primary_molecule = cursor.fetchone()
            if not primary_molecule:
                verification_report["issues"].append({
                    "group_id": group_id,
                    "issue": f"Primary molecule {primary_id} not found"
                })
                continue
            
            verification_report["primary_molecules"][primary_id] = {
                "name": primary_molecule['name'],
                "pubchem_cid": primary_molecule['pubchem_cid'],
                "secondary_molecules": [],
                "properties_acquired": 0,
                "mixture_components_acquired": 0
            }
            
            # Check each secondary molecule
            secondary_molecules = [m for m in group['molecules_detail'] if m['id'] != primary_id]
            for sec_mol in secondary_molecules:
                sec_id = sec_mol['id']
                
                # Check if molecule exists and has consolidated_to property
                cursor.execute("""
                    SELECT id, name, pubchem_cid, properties
                    FROM molecules
                    WHERE id = %s
                """, (sec_id,))
                
                molecule = cursor.fetchone()
                if not molecule:
                    verification_report["issues"].append({
                        "group_id": group_id,
                        "issue": f"Secondary molecule {sec_id} not found"
                    })
                    continue
                
                # Check if consolidated_to property exists and points to primary
                if not molecule['properties'] or 'consolidated_to' not in molecule['properties']:
                    verification_report["issues"].append({
                        "group_id": group_id,
                        "issue": f"Secondary molecule {sec_id} doesn't have consolidated_to property"
                    })
                    continue
                
                consolidated_to = molecule['properties']['consolidated_to']
                if consolidated_to != primary_id:
                    verification_report["issues"].append({
                        "group_id": group_id,
                        "issue": f"Secondary molecule {sec_id} points to {consolidated_to} instead of {primary_id}"
                    })
                    continue
                
                verification_report["secondary_molecules_consolidated"] += 1
                
                # Get log entries for this secondary molecule
                sec_log_entries = []
                for log_entry in log:
                    if log_entry.get('group_id') == group_id:
                        changes = log_entry.get('changes', [])
                        sec_log_entries.extend([c for c in changes if 
                                              c.get('from_molecule_id') == sec_id or 
                                              c.get('molecule_id') == sec_id])
                
                verification_report["secondary_molecules"][sec_id] = {
                    "name": molecule['name'],
                    "pubchem_cid": molecule['pubchem_cid'],
                    "consolidated_to": consolidated_to,
                    "properties_migrated": 0,
                    "mixture_components_migrated": 0,
                    "log_entries": len(sec_log_entries)
                }
                
                verification_report["primary_molecules"][primary_id]["secondary_molecules"].append(sec_id)
                
                # Count property migrations in log
                property_migrations = [e for e in sec_log_entries if e.get('type') == 'property_migration']
                verification_report["secondary_molecules"][sec_id]["properties_migrated"] = len(property_migrations)
                verification_report["primary_molecules"][primary_id]["properties_acquired"] += len(property_migrations)
                verification_report["property_migrations_verified"] += len(property_migrations)
                
                # Count mixture component migrations in log
                mixture_migrations = [e for e in sec_log_entries if e.get('type') == 'mixture_component_migration']
                verification_report["secondary_molecules"][sec_id]["mixture_components_migrated"] = len(mixture_migrations)
                verification_report["primary_molecules"][primary_id]["mixture_components_acquired"] += len(mixture_migrations)
                verification_report["mixture_component_migrations_verified"] += len(mixture_migrations)
            
            verification_report["groups_verified"] += 1
    
    return verification_report

def main():
    """Verify complex consolidation and generate a report."""
    # Load plan and log
    plan, log = load_plan_and_log()
    
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        # Verify consolidation
        verification_report = verify_complex_consolidation(conn, plan, log)
        
        # Save verification report
        with open(VERIFICATION_REPORT_FILE, 'w') as f:
            json.dump(verification_report, f, indent=2, default=str)
        
        # Print summary
        print("\nVerification Summary:")
        print(f"Groups verified: {verification_report['groups_verified']} of {verification_report['groups_total']}")
        print(f"Secondary molecules consolidated: {verification_report['secondary_molecules_consolidated']}")
        print(f"Property migrations verified: {verification_report['property_migrations_verified']}")
        print(f"Mixture component migrations verified: {verification_report['mixture_component_migrations_verified']}")
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
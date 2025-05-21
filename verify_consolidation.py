#!/usr/bin/env python3
"""
Verify the results of duplicate molecule consolidation.

This script:
1. Verifies that all consolidated molecules have proper references
2. Checks that properties were correctly migrated
3. Ensures relationships were preserved
4. Reports on the status of the consolidation
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Constants
CONSOLIDATION_LOG_FILE = "duplicate_consolidation_log.json"

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

def load_consolidation_log():
    """Load the consolidation log."""
    if not os.path.exists(CONSOLIDATION_LOG_FILE):
        print(f"Consolidation log file {CONSOLIDATION_LOG_FILE} not found.")
        sys.exit(1)
        
    with open(CONSOLIDATION_LOG_FILE, 'r') as f:
        return json.load(f)

def verify_consolidation(conn, log_entries):
    """Verify the consolidation results."""
    verification_results = []
    
    for entry in log_entries:
        # Collect the group information
        group_id = entry['group_id']
        consolidation_type = entry['consolidation_type']
        primary_id = entry['primary_molecule']['id']
        primary_name = entry['primary_molecule']['name']
        secondary_ids = [m['id'] for m in entry['secondary_molecules']]
        
        print(f"\nVerifying group {group_id} ({consolidation_type}):")
        print(f"  Primary: {primary_name} ({primary_id})")
        print(f"  Secondary: {len(secondary_ids)} molecules")
        
        # Check each secondary molecule
        all_valid = True
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check primary molecule exists
            cursor.execute("SELECT id, name FROM molecules WHERE id = %s", (primary_id,))
            if not cursor.fetchone():
                print(f"  ERROR: Primary molecule {primary_id} not found!")
                all_valid = False
            
            # Check each secondary molecule
            for sec_id in secondary_ids:
                # Verify molecule exists
                cursor.execute("SELECT id, name, properties FROM molecules WHERE id = %s", (sec_id,))
                sec_mol = cursor.fetchone()
                
                if not sec_mol:
                    print(f"  ERROR: Secondary molecule {sec_id} not found!")
                    all_valid = False
                    continue
                
                # Verify consolidated_to property
                if not sec_mol['properties'] or 'consolidated_to' not in sec_mol['properties']:
                    print(f"  ERROR: {sec_mol['name']} ({sec_id}) is missing consolidated_to property!")
                    all_valid = False
                elif sec_mol['properties']['consolidated_to'] != primary_id:
                    print(f"  ERROR: {sec_mol['name']} ({sec_id}) points to wrong primary: {sec_mol['properties']['consolidated_to']}")
                    all_valid = False
                else:
                    print(f"  OK: {sec_mol['name']} ({sec_id}) correctly points to primary")
                
                # Check migrated properties
                for change in entry['changes']:
                    if change['type'] == 'property_migration' and change['from_molecule_id'] == sec_id:
                        prop_type = change['property_type']
                        new_prop_id = change['new_property_id']
                        
                        # Verify property was migrated
                        cursor.execute("""
                            SELECT id FROM molecular_properties 
                            WHERE id = %s AND molecule_id = %s
                        """, (new_prop_id, primary_id))
                        
                        if not cursor.fetchone():
                            print(f"  ERROR: Property {prop_type} was not migrated from {sec_id} to {primary_id}!")
                            all_valid = False
                        else:
                            print(f"  OK: Property {prop_type} successfully migrated from {sec_id} to {primary_id}")
            
            # Verification summary
            result = {
                'group_id': group_id,
                'consolidation_type': consolidation_type,
                'primary_id': primary_id,
                'primary_name': primary_name,
                'secondary_count': len(secondary_ids),
                'verified': all_valid
            }
            verification_results.append(result)
        
    # Print summary
    print("\nVerification Summary:")
    verified_count = sum(1 for r in verification_results if r['verified'])
    print(f"Total groups: {len(verification_results)}")
    print(f"Successfully verified: {verified_count}")
    print(f"Failed verification: {len(verification_results) - verified_count}")
    
    # Return verification results
    return verification_results

def main():
    """Verify the consolidation results and report status."""
    # Load consolidation log
    print("Loading consolidation log...")
    log_entries = load_consolidation_log()
    
    # Group entries by group_id to avoid duplicate checks
    unique_entries = {}
    for entry in log_entries:
        group_id = entry['group_id']
        # Keep the most recent entry for each group
        unique_entries[group_id] = entry
    
    print(f"Found {len(unique_entries)} consolidated groups to verify")
    
    # Connect to database
    print("Connecting to database...")
    conn = connect_to_db()
    
    try:
        # Verify consolidation
        verification_results = verify_consolidation(conn, unique_entries.values())
        
        # Save verification results
        verification_file = "consolidation_verification.json"
        with open(verification_file, 'w') as f:
            json.dump(verification_results, f, indent=2)
        
        print(f"\nVerification results saved to {verification_file}")
        
    finally:
        conn.close()

if __name__ == "__main__":
    main()
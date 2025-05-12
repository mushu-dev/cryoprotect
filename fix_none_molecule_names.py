#!/usr/bin/env python3
"""
Fix molecules with 'None' names by retrieving proper names from PubChem.

This script:
1. Finds all molecules with 'None' as their name and a valid PubChem CID
2. Retrieves the preferred name from PubChem using the PubChem API
3. Updates the molecule name in the database
4. Records the original name (None) in the properties field
5. Updates modification_history
"""

import os
import sys
import json
import time
import requests
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Constants
BATCH_SIZE = 20  # Process this many molecules at a time to avoid rate limiting
RETRY_DELAY = 0.5  # Delay between retries in seconds
MAX_RETRIES = 3
PUBCHEM_CACHE = {}

# Database connection
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

def get_none_molecules_with_cid(conn):
    """Get all molecules with 'None' as name and a valid PubChem CID."""
    molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                   data_source, properties
            FROM molecules 
            WHERE name = 'None' 
            AND pubchem_cid IS NOT NULL
            ORDER BY pubchem_cid
        """)
        
        molecules = cursor.fetchall()
    
    return molecules

def get_name_from_pubchem(pubchem_cid):
    """Get the preferred name for a molecule from PubChem."""
    if not pubchem_cid:
        return None
        
    # Check cache first
    if pubchem_cid in PUBCHEM_CACHE:
        return PUBCHEM_CACHE[pubchem_cid]
    
    # Call PubChem API
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/synonyms/JSON"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                if 'InformationList' in data and 'Information' in data['InformationList']:
                    if len(data['InformationList']['Information']) > 0:
                        if 'Synonym' in data['InformationList']['Information'][0]:
                            # Get the first (preferred) synonym
                            name = data['InformationList']['Information'][0]['Synonym'][0]
                            PUBCHEM_CACHE[pubchem_cid] = name
                            return name
            elif response.status_code == 429:  # Rate limit
                # Wait before trying again
                time.sleep(RETRY_DELAY * (attempt + 1))
                continue
            else:
                break
                
        except Exception as e:
            print(f"  Error fetching PubChem data: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
            
    # If we couldn't get a name, cache None to avoid further attempts
    PUBCHEM_CACHE[pubchem_cid] = None
    return None

def update_molecule_name(conn, molecule_id, old_name, new_name, dry_run=False):
    """Update a molecule's name in the database."""
    try:
        if not dry_run:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules
                    SET name = %s,
                        properties = COALESCE(properties, '{}'::jsonb) || %s::jsonb,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                            'timestamp', NOW(),
                            'action', 'fix_none_molecule_names',
                            'details', 'Updated None name with PubChem preferred name',
                            'old_name', %s,
                            'new_name', %s
                        )::jsonb
                    WHERE id = %s
                """, (
                    new_name,
                    json.dumps({"original_name": old_name}),
                    old_name,
                    new_name,
                    molecule_id
                ))
                
                return cursor.rowcount
        else:
            # Return what would be updated in dry run mode
            return 1
    except Exception as e:
        print(f"  Error updating molecule name: {e}")
        return 0

def main():
    """Fix molecules with 'None' names by retrieving proper names from PubChem."""
    import argparse
    parser = argparse.ArgumentParser(description="Fix None molecule names")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    parser.add_argument("--limit", type=int, default=0, help="Limit the number of molecules to process (0 for all)")
    args = parser.parse_args()
    
    print("Fixing 'None' molecule names using PubChem data...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get molecules with 'None' names and PubChem CIDs
        molecules = get_none_molecules_with_cid(conn)
        print(f"Found {len(molecules)} molecules with 'None' names and valid PubChem CIDs")
        
        # Process molecules in batches to respect API rate limits
        to_update = []
        molecules_processed = 0
        current_batch = 0
        
        total_batches = (len(molecules) + BATCH_SIZE - 1) // BATCH_SIZE
        
        for i in range(0, len(molecules), BATCH_SIZE):
            current_batch += 1
            batch = molecules[i:i+BATCH_SIZE]
            print(f"\nProcessing batch {current_batch}/{total_batches} ({len(batch)} molecules)")
            
            for molecule in batch:
                molecule_id = molecule['id']
                pubchem_cid = molecule['pubchem_cid']
                
                print(f"  Processing molecule ID: {molecule_id}, PubChem CID: {pubchem_cid}")
                
                # Get name from PubChem
                pubchem_name = get_name_from_pubchem(pubchem_cid)
                
                if pubchem_name:
                    to_update.append({
                        'id': molecule_id,
                        'old_name': 'None',
                        'new_name': pubchem_name,
                        'pubchem_cid': pubchem_cid
                    })
                    print(f"    Found name: {pubchem_name}")
                else:
                    print(f"    No name found from PubChem")
                
                molecules_processed += 1
                if args.limit > 0 and molecules_processed >= args.limit:
                    break
            
            # Take a short break between batches to avoid hitting rate limits
            if not current_batch == total_batches and not args.dry_run:
                time.sleep(1)
                
            if args.limit > 0 and molecules_processed >= args.limit:
                break
                
        # Update molecules in the database
        if to_update:
            print(f"\nFound {len(to_update)} molecules with names from PubChem")
            updated_count = 0
            
            for molecule in to_update:
                molecule_id = molecule['id']
                old_name = molecule['old_name']
                new_name = molecule['new_name']
                
                if not args.dry_run:
                    print(f"  Updating {old_name} to {new_name}")
                else:
                    print(f"  Would update {old_name} to {new_name}")
                
                updated = update_molecule_name(conn, molecule_id, old_name, new_name, args.dry_run)
                if updated:
                    updated_count += 1
            
            # Save results to file
            results = {
                "total_none_molecules": len(molecules),
                "molecules_processed": molecules_processed,
                "molecules_with_pubchem_names": len(to_update),
                "molecules_updated": updated_count,
                "updates": to_update,
                "dry_run": args.dry_run
            }
            
            with open("none_molecule_name_fixes_report.json", "w") as f:
                json.dump(results, f, indent=2, default=str)
            
            print("\nSummary:")
            print(f"Total molecules with 'None' names and PubChem CIDs: {len(molecules)}")
            print(f"Molecules processed: {molecules_processed}")
            print(f"Molecules with names from PubChem: {len(to_update)}")
            print(f"Molecules successfully updated: {updated_count}")
            print(f"Report saved to 'none_molecule_name_fixes_report.json'")
            
            # Commit changes if not dry run
            if not args.dry_run:
                if args.auto_commit:
                    conn.commit()
                    print("Changes automatically committed to the database")
                else:
                    confirm = input("\nCommit these changes to the database? (y/n): ")
                    if confirm.lower() == 'y':
                        conn.commit()
                        print("Changes committed to the database")
                    else:
                        conn.rollback()
                        print("Changes rolled back, no updates were made")
        else:
            print("No molecules could be updated with names from PubChem")
    
    except Exception as e:
        conn.rollback()
        print(f"Error fixing None molecule names: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
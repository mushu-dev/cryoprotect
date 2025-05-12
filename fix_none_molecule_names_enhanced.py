#!/usr/bin/env python3
"""
Fix molecules with 'None' names by retrieving proper names from PubChem using
a more comprehensive approach that queries the full compound details.

This script:
1. Finds all molecules with 'None' as their name and a valid PubChem CID
2. Retrieves compound information from PubChem using multiple endpoints
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
BATCH_SIZE = 10  # Process this many molecules at a time to avoid rate limiting
RETRY_DELAY = 1  # Delay between retries in seconds
MAX_RETRIES = 3
PUBCHEM_CACHE = {}
USER_AGENT = "CryoProtect/1.0 (CryoProtect Database Maintenance; mailto:your-email@example.com)"

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

def get_name_from_pubchem_comprehensive(pubchem_cid):
    """
    Get the name for a molecule from PubChem using multiple endpoints
    to maximize the chance of finding a name.
    """
    if not pubchem_cid:
        return None
        
    # Check cache first
    if pubchem_cid in PUBCHEM_CACHE:
        return PUBCHEM_CACHE[pubchem_cid]
    
    # Try multiple PubChem API endpoints
    
    # 1. Try Compound Summary first (most reliable for names)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{pubchem_cid}/JSON"
    headers = {'User-Agent': USER_AGENT}
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = response.json()
                
                # Navigate the complex PubChem data structure to find the name
                if 'Record' in data and 'Section' in data['Record']:
                    for section in data['Record']['Section']:
                        if section.get('TOCHeading') == 'Names and Identifiers':
                            if 'Section' in section:
                                for subsection in section['Section']:
                                    if subsection.get('TOCHeading') == 'Record Title':
                                        if 'Information' in subsection:
                                            for info in subsection['Information']:
                                                if 'StringValue' in info and info['StringValue']:
                                                    name = info['StringValue']
                                                    PUBCHEM_CACHE[pubchem_cid] = name
                                                    return name
            
            # If we hit a rate limit or other error, wait and retry
            if response.status_code == 429:
                time.sleep(RETRY_DELAY * (attempt + 1))
                continue
                
        except Exception as e:
            print(f"  Error fetching PubChem data from summary: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
    
    # 2. Try Property endpoint if summary failed
    time.sleep(0.5)  # Brief pause before trying another endpoint
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/property/Title,IUPACName,MolecularFormula/JSON"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = response.json()
                
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    props = data['PropertyTable']['Properties']
                    if props and len(props) > 0:
                        # Try Title first
                        if 'Title' in props[0] and props[0]['Title']:
                            name = props[0]['Title']
                            PUBCHEM_CACHE[pubchem_cid] = name
                            return name
                        
                        # Try IUPAC name as fallback
                        if 'IUPACName' in props[0] and props[0]['IUPACName']:
                            name = props[0]['IUPACName']
                            PUBCHEM_CACHE[pubchem_cid] = name
                            return name
            
            # If we hit a rate limit or other error, wait and retry
            if response.status_code == 429:
                time.sleep(RETRY_DELAY * (attempt + 1))
                continue
                
        except Exception as e:
            print(f"  Error fetching PubChem property data: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
    
    # 3. Try Synonyms endpoint as last resort
    time.sleep(0.5)  # Brief pause before trying another endpoint
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/synonyms/JSON"
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = response.json()
                
                if 'InformationList' in data and 'Information' in data['InformationList']:
                    if len(data['InformationList']['Information']) > 0:
                        if 'Synonym' in data['InformationList']['Information'][0]:
                            synonyms = data['InformationList']['Information'][0]['Synonym']
                            if synonyms and len(synonyms) > 0:
                                # Get the first (preferred) synonym
                                name = synonyms[0]
                                PUBCHEM_CACHE[pubchem_cid] = name
                                return name
            
            # If we hit a rate limit or other error, wait and retry
            if response.status_code == 429:
                time.sleep(RETRY_DELAY * (attempt + 1))
                continue
                
        except Exception as e:
            print(f"  Error fetching PubChem synonym data: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
    
    # 4. Generate a fallback name using Compound- prefix + CID if no other name is found
    fallback_name = f"Compound-{pubchem_cid}"
    PUBCHEM_CACHE[pubchem_cid] = fallback_name
    return fallback_name

def update_molecule_name(conn, molecule_id, old_name, new_name, name_source, dry_run=False):
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
                            'action', 'fix_none_molecule_names_enhanced',
                            'details', %s,
                            'old_name', %s,
                            'new_name', %s
                        )::jsonb
                    WHERE id = %s
                """, (
                    new_name,
                    json.dumps({
                        "original_name": old_name,
                        "name_source": name_source
                    }),
                    f"Updated None name with {name_source}",
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
    
    print("Fixing 'None' molecule names using enhanced PubChem lookup...")
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
                
                # Get name from PubChem using comprehensive approach
                pubchem_name = get_name_from_pubchem_comprehensive(pubchem_cid)
                
                if pubchem_name:
                    name_source = "PubChem API"
                    # For fallback-generated names, note the source
                    if pubchem_name.startswith("Compound-"):
                        name_source = "fallback using CID"
                        
                    to_update.append({
                        'id': molecule_id,
                        'old_name': 'None',
                        'new_name': pubchem_name,
                        'pubchem_cid': pubchem_cid,
                        'name_source': name_source
                    })
                    print(f"    Found name: {pubchem_name} (Source: {name_source})")
                else:
                    print(f"    No name found from PubChem")
                
                molecules_processed += 1
                if args.limit > 0 and molecules_processed >= args.limit:
                    break
            
            # Take a short break between batches to avoid hitting rate limits
            if not current_batch == total_batches and not args.dry_run:
                time.sleep(2)
                
            if args.limit > 0 and molecules_processed >= args.limit:
                break
                
        # Update molecules in the database
        if to_update:
            print(f"\nFound {len(to_update)} molecules with names from PubChem")
            updated_count = 0
            fallback_count = 0
            
            for molecule in to_update:
                molecule_id = molecule['id']
                old_name = molecule['old_name']
                new_name = molecule['new_name']
                name_source = molecule['name_source']
                
                if name_source == "fallback using CID":
                    fallback_count += 1
                
                if not args.dry_run:
                    print(f"  Updating {old_name} to {new_name} (Source: {name_source})")
                else:
                    print(f"  Would update {old_name} to {new_name} (Source: {name_source})")
                
                updated = update_molecule_name(conn, molecule_id, old_name, new_name, name_source, args.dry_run)
                if updated:
                    updated_count += 1
            
            # Save results to file
            results = {
                "total_none_molecules": len(molecules),
                "molecules_processed": molecules_processed,
                "molecules_with_names": len(to_update),
                "molecules_with_fallback_names": fallback_count,
                "molecules_updated": updated_count,
                "updates": to_update,
                "dry_run": args.dry_run
            }
            
            with open("none_molecule_name_fixes_enhanced_report.json", "w") as f:
                json.dump(results, f, indent=2, default=str)
            
            print("\nSummary:")
            print(f"Total molecules with 'None' names and PubChem CIDs: {len(molecules)}")
            print(f"Molecules processed: {molecules_processed}")
            print(f"Molecules with retrieved names: {len(to_update)}")
            print(f"Molecules using fallback naming: {fallback_count}")
            print(f"Molecules successfully updated: {updated_count}")
            print(f"Report saved to 'none_molecule_name_fixes_enhanced_report.json'")
            
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
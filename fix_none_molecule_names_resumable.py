#!/usr/bin/env python3
"""
Fix molecules with 'None' names by retrieving proper names from PubChem using
a resumable approach that can be interrupted and continued later.

This script:
1. Finds all molecules with 'None' as their name and a valid PubChem CID
2. Retrieves compound information from PubChem using multiple endpoints
3. Updates the molecule name in the database
4. Records the original name (None) in the properties field
5. Updates modification_history
6. Saves progress after each batch of molecules
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
PROGRESS_FILE = "none_name_fix_progress.json"
RESULTS_FILE = "none_molecule_name_fixes_report.json"

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

def load_progress():
    """Load progress from file if it exists."""
    if os.path.exists(PROGRESS_FILE):
        try:
            with open(PROGRESS_FILE, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading progress file: {e}")
    
    # Default progress data
    return {
        "processed_ids": [],
        "updated_ids": [],
        "updates": [],
        "total_processed": 0,
        "total_updated": 0,
        "total_fallback": 0
    }

def save_progress(progress_data):
    """Save progress to file."""
    try:
        with open(PROGRESS_FILE, 'w') as f:
            json.dump(progress_data, f, indent=2, default=str)
    except Exception as e:
        print(f"Error saving progress file: {e}")

def save_results(results):
    """Save results to file."""
    try:
        with open(RESULTS_FILE, 'w') as f:
            json.dump(results, f, indent=2, default=str)
    except Exception as e:
        print(f"Error saving results file: {e}")

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
                            'action', 'fix_none_molecule_names',
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
    parser = argparse.ArgumentParser(description="Fix None molecule names (resumable)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    parser.add_argument("--limit", type=int, default=0, help="Limit the number of molecules to process (0 for all)")
    parser.add_argument("--reset-progress", action="store_true", help="Reset progress and start from beginning")
    args = parser.parse_args()
    
    print("Fixing 'None' molecule names using resumable PubChem lookup...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Load progress from previous run
    progress = load_progress()
    
    if args.reset_progress:
        print("Resetting progress and starting from beginning")
        progress = {
            "processed_ids": [],
            "updated_ids": [],
            "updates": [],
            "total_processed": 0,
            "total_updated": 0,
            "total_fallback": 0
        }
    else:
        print(f"Resuming from previous run: {progress['total_processed']} molecules processed, {progress['total_updated']} updated")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get molecules with 'None' names and PubChem CIDs
        all_molecules = get_none_molecules_with_cid(conn)
        print(f"Found {len(all_molecules)} molecules with 'None' names and valid PubChem CIDs")
        
        # Filter out already processed molecules
        molecules = [m for m in all_molecules if str(m['id']) not in progress['processed_ids']]
        print(f"{len(molecules)} molecules remaining to process")
        
        # Process molecules in batches to respect API rate limits
        current_batch = 0
        
        total_batches = (len(molecules) + BATCH_SIZE - 1) // BATCH_SIZE
        
        for i in range(0, len(molecules), BATCH_SIZE):
            current_batch += 1
            batch = molecules[i:i+BATCH_SIZE]
            print(f"\nProcessing batch {current_batch}/{total_batches} ({len(batch)} molecules)")
            
            batch_updates = []
            
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
                        progress['total_fallback'] += 1
                        
                    update_info = {
                        'id': molecule_id,
                        'old_name': 'None',
                        'new_name': pubchem_name,
                        'pubchem_cid': pubchem_cid,
                        'name_source': name_source
                    }
                    
                    batch_updates.append(update_info)
                    print(f"    Found name: {pubchem_name} (Source: {name_source})")
                else:
                    print(f"    No name found from PubChem")
                
                # Add to processed IDs
                progress['processed_ids'].append(str(molecule_id))
                progress['total_processed'] += 1
                
                if args.limit > 0 and progress['total_processed'] >= args.limit:
                    break
            
            # Update molecules in the database
            if batch_updates:
                print(f"\nUpdating {len(batch_updates)} molecules with names from PubChem")
                
                for update_info in batch_updates:
                    molecule_id = update_info['id']
                    old_name = update_info['old_name']
                    new_name = update_info['new_name']
                    name_source = update_info['name_source']
                    
                    if not args.dry_run:
                        print(f"  Updating {old_name} to {new_name} (Source: {name_source})")
                    else:
                        print(f"  Would update {old_name} to {new_name} (Source: {name_source})")
                    
                    updated = update_molecule_name(conn, molecule_id, old_name, new_name, name_source, args.dry_run)
                    if updated:
                        progress['updated_ids'].append(str(molecule_id))
                        progress['total_updated'] += 1
                        progress['updates'].append(update_info)
                
                # Commit changes after each batch if not dry run
                if not args.dry_run:
                    if args.auto_commit:
                        conn.commit()
                        print("  Changes automatically committed to the database")
                    else:
                        confirm = input("\n  Commit these changes to the database? (y/n): ")
                        if confirm.lower() == 'y':
                            conn.commit()
                            print("  Changes committed to the database")
                        else:
                            conn.rollback()
                            print("  Changes rolled back, no updates were made")
                            
                            # Remove rolled back changes from progress
                            for update_info in batch_updates:
                                if str(update_info['id']) in progress['updated_ids']:
                                    progress['updated_ids'].remove(str(update_info['id']))
                                    progress['total_updated'] -= 1
                                    
                                    # Remove from updates list
                                    progress['updates'] = [u for u in progress['updates'] if u['id'] != update_info['id']]
            
            # Save progress after each batch
            save_progress(progress)
            print(f"  Progress saved: {progress['total_processed']}/{len(all_molecules)} processed, {progress['total_updated']} updated")
            
            # Take a short break between batches to avoid hitting rate limits
            if not current_batch == total_batches and not args.dry_run:
                time.sleep(2)
                
            if args.limit > 0 and progress['total_processed'] >= args.limit:
                break
                
        # Save final results
        results = {
            "total_none_molecules": len(all_molecules),
            "molecules_processed": progress['total_processed'],
            "molecules_with_names": progress['total_updated'],
            "molecules_with_fallback_names": progress['total_fallback'],
            "updates": progress['updates'],
            "dry_run": args.dry_run
        }
        
        save_results(results)
        
        print("\nFinal Summary:")
        print(f"Total molecules with 'None' names and PubChem CIDs: {len(all_molecules)}")
        print(f"Molecules processed: {progress['total_processed']}")
        print(f"Molecules with names from PubChem: {progress['total_updated']}")
        print(f"Molecules using fallback naming: {progress['total_fallback']}")
        print(f"Results saved to '{RESULTS_FILE}'")
        
    except KeyboardInterrupt:
        print("\nProcess interrupted by user")
        print("Progress has been saved and can be resumed later")
        conn.rollback()
    
    except Exception as e:
        conn.rollback()
        print(f"Error fixing None molecule names: {e}")
        print("Progress has been saved and can be resumed later")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
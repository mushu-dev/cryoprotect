#!/usr/bin/env python3
"""
Standardize molecule names according to established naming conventions.

This script:
1. Applies name standardization rules from MOLECULE_NAME_STANDARDIZATION_CONVENTIONS.md
2. Updates molecule names in the database
3. Records the original name in the properties
4. Updates modification_history
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
COMMON_CRYOPROTECTANTS = {
    # Original name patterns -> Standardized name
    "glycerol": "Glycerol",
    "propane-1,2,3-triol": "Glycerol",
    "ethylene glycol": "Ethylene glycol",
    "ethane-1,2-diol": "Ethylene glycol", 
    "propylene glycol": "Propylene glycol",
    "propane-1,2-diol": "Propylene glycol",
    "dimethyl sulfoxide": "Dimethyl sulfoxide",
    "dmso": "Dimethyl sulfoxide",
    "methylsulfinylmethane": "Dimethyl sulfoxide",
    "trehalose": "Trehalose",
    "alpha,alpha-trehalose": "Trehalose",
    "α,α-trehalose": "Trehalose",
    "methanol": "Methanol",
    "formamide": "Formamide",
    "acetamide": "Acetamide"
}

PUBCHEM_CACHE = {}
RETRY_DELAY = 0.5  # Delay between retries in seconds
MAX_RETRIES = 3

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

def get_molecules_for_standardization(conn):
    """Get all molecules that need name standardization."""
    molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get all molecules
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, smiles, 
                   data_source, properties
            FROM molecules
            WHERE name IS NOT NULL AND name NOT LIKE 'TEST\\_%'
            ORDER BY name
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

def standardize_name(molecule):
    """Apply name standardization rules to a molecule."""
    original_name = molecule['name']
    
    # Skip already standardized names
    if molecule.get('properties') and molecule['properties'].get('original_name'):
        return original_name, False
    
    # Handle None cases
    if original_name.lower() == 'none':
        # Try to get a name from PubChem
        if molecule['pubchem_cid']:
            pubchem_name = get_name_from_pubchem(molecule['pubchem_cid'])
            if pubchem_name:
                return pubchem_name, True
        
        # Fallback if no PubChem name
        if molecule['molecular_formula']:
            return f"Compound-{molecule['molecular_formula']}", True
        
        # No good options, keep as is
        return original_name, False
    
    # Check for common cryoprotectants
    lower_name = original_name.lower()
    for pattern, standardized in COMMON_CRYOPROTECTANTS.items():
        if lower_name == pattern:
            if original_name != standardized:
                return standardized, True
            else:
                return original_name, False
    
    # Apply first letter capitalization rule
    if len(original_name) > 1:
        standardized = original_name[0].upper() + original_name[1:].lower()
        if standardized != original_name:
            return standardized, True
    
    # No changes needed
    return original_name, False

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
                            'action', 'standardize_molecule_names',
                            'details', 'Standardized molecule name according to conventions',
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
    """Standardize molecule names according to established conventions."""
    import argparse
    parser = argparse.ArgumentParser(description="Standardize molecule names")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    parser.add_argument("--limit", type=int, default=0, help="Limit the number of molecules to process (0 for all)")
    args = parser.parse_args()
    
    print("Standardizing molecule names...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get molecules that need standardization
        molecules = get_molecules_for_standardization(conn)
        print(f"Found {len(molecules)} molecules to analyze")
        
        # Apply name standardization rules
        to_standardize = []
        molecules_processed = 0
        
        for molecule in molecules:
            molecule_id = molecule['id']
            old_name = molecule['name']
            
            # Apply standardization
            new_name, changed = standardize_name(molecule)
            
            if changed:
                to_standardize.append({
                    'id': molecule_id,
                    'old_name': old_name,
                    'new_name': new_name
                })
            
            molecules_processed += 1
            if args.limit > 0 and molecules_processed >= args.limit:
                break
                
        # Update molecules in the database
        if to_standardize:
            print(f"Found {len(to_standardize)} molecules that need name standardization")
            updated_count = 0
            
            for molecule in to_standardize:
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
                "total_molecules": molecules_processed,
                "standardized_molecules": updated_count,
                "changes": to_standardize,
                "dry_run": args.dry_run
            }
            
            with open("molecule_name_standardization_report.json", "w") as f:
                json.dump(results, f, indent=2, default=str)
            
            print("\nSummary:")
            print(f"Molecules analyzed: {molecules_processed}")
            print(f"Molecules requiring standardization: {len(to_standardize)}")
            print(f"Molecules successfully updated: {updated_count}")
            print(f"Report saved to 'molecule_name_standardization_report.json'")
            
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
            print("No molecules require name standardization")
    
    except Exception as e:
        conn.rollback()
        print(f"Error standardizing molecule names: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
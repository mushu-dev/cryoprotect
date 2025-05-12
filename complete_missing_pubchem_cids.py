#!/usr/bin/env python3
"""
Complete missing PubChem CID values in molecules table.

This script:
1. Identifies molecules with missing pubchem_cid values
2. Attempts to find PubChem CIDs using name/formula via PubChem API
3. Updates the database with found CID values
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

def get_missing_cid_molecules(conn):
    """Get molecules with missing pubchem_cid values."""
    molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, smiles, inchi, inchikey
            FROM molecules
            WHERE pubchem_cid IS NULL
            ORDER BY name
        """)
        molecules = cursor.fetchall()
    
    return molecules

def search_pubchem_by_name(name):
    """Search PubChem for a compound by name."""
    if not name:
        return None
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                return data['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"Error searching PubChem by name: {e}")
    
    return None

def search_pubchem_by_formula(formula):
    """Search PubChem for a compound by formula."""
    if not formula:
        return None
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula}/cids/JSON"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                # Return the first CID, might not be the exact compound
                return data['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"Error searching PubChem by formula: {e}")
    
    return None

def search_pubchem_by_inchi(inchi):
    """Search PubChem for a compound by InChI."""
    if not inchi:
        return None
    
    # URL encode the InChI
    import urllib.parse
    encoded_inchi = urllib.parse.quote(inchi)
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/{encoded_inchi}/cids/JSON"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                return data['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"Error searching PubChem by InChI: {e}")
    
    return None

def search_pubchem_by_inchikey(inchikey):
    """Search PubChem for a compound by InChIKey."""
    if not inchikey:
        return None
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                return data['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"Error searching PubChem by InChIKey: {e}")
    
    return None

def search_for_cid(molecule):
    """Search for PubChem CID using multiple methods."""
    results = {
        "molecule_id": molecule["id"],
        "name": molecule["name"],
        "found_cid": None,
        "search_method": None
    }
    
    # Try different search methods
    methods = [
        {"func": search_pubchem_by_name, "key": "name", "delay": 0.5},
        {"func": search_pubchem_by_inchikey, "key": "inchikey", "delay": 0.5},
        {"func": search_pubchem_by_inchi, "key": "inchi", "delay": 0.5},
        {"func": search_pubchem_by_formula, "key": "molecular_formula", "delay": 0.5}
    ]
    
    for method in methods:
        if molecule.get(method["key"]):
            print(f"  Searching for {molecule['name']} by {method['key']}...")
            cid = method["func"](molecule[method["key"]])
            
            if cid:
                results["found_cid"] = cid
                results["search_method"] = method["key"]
                print(f"  Found CID {cid} for {molecule['name']} using {method['key']}")
                break
            
            # Respect PubChem API rate limits
            time.sleep(method["delay"])
    
    return results

def check_cid_exists(conn, cid):
    """Check if a PubChem CID already exists in the database."""
    with conn.cursor() as cursor:
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid = %s", (cid,))
        count = cursor.fetchone()[0]
        return count > 0

def update_molecule_cid(conn, molecule_id, cid):
    """Update molecule with found PubChem CID."""
    # First check if this CID already exists
    if check_cid_exists(conn, cid):
        print(f"  WARNING: CID {cid} already exists in the database, skipping update")
        return 0

    with conn.cursor() as cursor:
        cursor.execute("""
            UPDATE molecules
            SET pubchem_cid = %s,
                updated_at = NOW(),
                modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                    'timestamp', NOW(),
                    'action', 'complete_missing_pubchem_cids',
                    'details', 'Added missing PubChem CID'
                )::jsonb
            WHERE id = %s
        """, (cid, molecule_id))

        return cursor.rowcount

def main():
    """Complete missing PubChem CID values."""
    print("Completing missing PubChem CID values...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get molecules with missing pubchem_cid
        molecules = get_missing_cid_molecules(conn)
        print(f"Found {len(molecules)} molecules with missing PubChem CID values")
        
        if not molecules:
            print("No molecules to update.")
            return
        
        results = {
            "updated": 0,
            "not_found": 0,
            "details": []
        }
        
        # Process each molecule
        for i, molecule in enumerate(molecules):
            print(f"Processing {i+1}/{len(molecules)}: {molecule['name']}")
            
            search_result = search_for_cid(molecule)
            results["details"].append(search_result)
            
            if search_result["found_cid"]:
                # Update molecule in database
                updated = update_molecule_cid(conn, molecule["id"], search_result["found_cid"])
                if updated:
                    results["updated"] += 1
                    print(f"  Updated molecule {molecule['name']} with CID {search_result['found_cid']}")
            else:
                results["not_found"] += 1
                print(f"  No PubChem CID found for {molecule['name']}")
        
        # Commit changes
        conn.commit()
        
        # Summary
        print("\nResults Summary:")
        print(f"Total molecules processed: {len(molecules)}")
        print(f"Updated with found CIDs: {results['updated']}")
        print(f"CIDs not found: {results['not_found']}")
        
        # Save results to file
        with open("pubchem_cid_completion_report.json", "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        print("Report saved to 'pubchem_cid_completion_report.json'")
    
    except Exception as e:
        conn.rollback()
        print(f"Error completing PubChem CIDs: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Analyze remaining molecules with missing PubChem CIDs.

This script:
1. Identifies molecules still missing pubchem_cid values
2. Analyzes their properties to determine why they're missing
3. Categorizes them (test molecules, duplicates, etc.)
4. Provides recommendations for handling each category
"""

import os
import sys
import json
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
    """Get detailed information about molecules with missing pubchem_cid values."""
    molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, smiles, inchi, inchikey, 
                   data_source, created_at, is_public, properties
            FROM molecules
            WHERE pubchem_cid IS NULL
            ORDER BY name
        """)
        molecules = cursor.fetchall()
    
    return molecules

def check_potential_duplicates(conn, molecules):
    """Check if molecules missing CIDs are potential duplicates of existing ones."""
    results = []
    
    for molecule in molecules:
        duplicates = []
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check by name
            if molecule["name"]:
                cursor.execute("""
                    SELECT id, name, pubchem_cid, molecular_formula
                    FROM molecules
                    WHERE lower(name) = lower(%s) AND pubchem_cid IS NOT NULL
                    LIMIT 5
                """, (molecule["name"],))
                name_matches = cursor.fetchall()
                if name_matches:
                    for match in name_matches:
                        duplicates.append({
                            "match_type": "name",
                            "molecule_id": match["id"],
                            "name": match["name"],
                            "pubchem_cid": match["pubchem_cid"]
                        })
            
            # Check by formula
            if molecule["molecular_formula"]:
                cursor.execute("""
                    SELECT id, name, pubchem_cid, molecular_formula
                    FROM molecules
                    WHERE molecular_formula = %s AND pubchem_cid IS NOT NULL
                    AND id != %s
                    LIMIT 5
                """, (molecule["molecular_formula"], molecule["id"]))
                formula_matches = cursor.fetchall()
                if formula_matches:
                    for match in formula_matches:
                        duplicates.append({
                            "match_type": "formula",
                            "molecule_id": match["id"],
                            "name": match["name"],
                            "pubchem_cid": match["pubchem_cid"]
                        })
        
        results.append({
            "molecule_id": molecule["id"],
            "name": molecule["name"],
            "formula": molecule["molecular_formula"],
            "potential_duplicates": duplicates
        })
    
    return results

def categorize_molecules(molecules, duplicate_info):
    """Categorize molecules based on their properties and potential duplicates."""
    categories = {
        "test_molecules": [],
        "duplicates": [],
        "incomplete_data": [],
        "unknown": []
    }
    
    duplicate_map = {item["molecule_id"]: item["potential_duplicates"] for item in duplicate_info}
    
    for molecule in molecules:
        # Check if it's a test molecule
        if molecule["name"] and any(test_term in molecule["name"].lower() for test_term in ["test", "dummy", "example"]):
            categories["test_molecules"].append(molecule)
        
        # Check if it has potential duplicates
        elif molecule["id"] in duplicate_map and duplicate_map[molecule["id"]]:
            molecule["potential_duplicates"] = duplicate_map[molecule["id"]]
            categories["duplicates"].append(molecule)
        
        # Check if it has incomplete data
        elif not molecule["molecular_formula"] or not molecule["smiles"]:
            categories["incomplete_data"].append(molecule)
        
        # Otherwise, unknown reason
        else:
            categories["unknown"].append(molecule)
    
    return categories

def analyze_missing_cids():
    """Analyze molecules with missing PubChem CIDs."""
    print("Analyzing molecules with missing PubChem CIDs...")
    
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Get molecules with missing CIDs
        molecules = get_missing_cid_molecules(conn)
        print(f"Found {len(molecules)} molecules with missing PubChem CIDs")
        
        if not molecules:
            print("No molecules to analyze.")
            return
        
        # Check for potential duplicates
        print("Checking for potential duplicates...")
        duplicate_info = check_potential_duplicates(conn, molecules)
        
        # Categorize molecules
        print("Categorizing molecules...")
        categories = categorize_molecules(molecules, duplicate_info)
        
        # Print summary
        print("\nAnalysis Summary:")
        print(f"Total molecules with missing CIDs: {len(molecules)}")
        print(f"Test molecules: {len(categories['test_molecules'])}")
        print(f"Potential duplicates: {len(categories['duplicates'])}")
        print(f"Incomplete data: {len(categories['incomplete_data'])}")
        print(f"Unknown reason: {len(categories['unknown'])}")
        
        # Print details of each category
        if categories["test_molecules"]:
            print("\nTest Molecules:")
            for molecule in categories["test_molecules"]:
                print(f"  - {molecule['name']} (ID: {molecule['id']})")
        
        if categories["duplicates"]:
            print("\nPotential Duplicates:")
            for molecule in categories["duplicates"]:
                print(f"  - {molecule['name']} (ID: {molecule['id']})")
                for dup in molecule["potential_duplicates"]:
                    print(f"    * Match by {dup['match_type']}: {dup['name']} (CID: {dup['pubchem_cid']})")
        
        if categories["incomplete_data"]:
            print("\nIncomplete Data:")
            for molecule in categories["incomplete_data"]:
                print(f"  - {molecule['name']} (ID: {molecule['id']})")
                print(f"    * Formula: {molecule['molecular_formula'] or 'Missing'}")
                print(f"    * SMILES: {molecule['smiles'] or 'Missing'}")
        
        if categories["unknown"]:
            print("\nUnknown Reason:")
            for molecule in categories["unknown"]:
                print(f"  - {molecule['name']} (ID: {molecule['id']})")
        
        # Save analysis to file
        analysis_results = {
            "total_count": len(molecules),
            "categories": {
                "test_molecules": [{"id": m["id"], "name": m["name"]} for m in categories["test_molecules"]],
                "duplicates": [{
                    "id": m["id"], 
                    "name": m["name"],
                    "duplicates": m["potential_duplicates"]
                } for m in categories["duplicates"]],
                "incomplete_data": [{
                    "id": m["id"], 
                    "name": m["name"],
                    "formula": m["molecular_formula"],
                    "smiles": m["smiles"]
                } for m in categories["incomplete_data"]],
                "unknown": [{"id": m["id"], "name": m["name"]} for m in categories["unknown"]]
            },
            "recommendations": {
                "test_molecules": "Keep as is or add 'TEST_' prefix to name for clarity",
                "duplicates": "Consider updating with CID from matched molecules",
                "incomplete_data": "Complete molecular information to enable CID lookup",
                "unknown": "Perform manual lookup in PubChem or leave as is"
            }
        }
        
        with open("missing_cid_analysis.json", "w") as f:
            json.dump(analysis_results, f, indent=2, default=str)
        
        print("\nAnalysis saved to 'missing_cid_analysis.json'")
        
        # Generate recommendations
        print("\nRecommendations:")
        print("1. Test Molecules: Keep as is or add 'TEST_' prefix to names for clarity")
        print("2. Duplicates: Update with CID from matched molecules")
        print("3. Incomplete Data: Complete molecular information to enable CID lookup")
        print("4. Unknown: Perform manual lookup in PubChem or leave as is")
    
    except Exception as e:
        print(f"Error analyzing missing CIDs: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    analyze_missing_cids()
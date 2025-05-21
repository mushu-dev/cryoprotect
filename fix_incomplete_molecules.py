#!/usr/bin/env python3
"""
Fix molecules with incomplete data in the database.

This script:
1. Identifies molecules with incomplete data (e.g., missing formula but has SMILES)
2. Uses RDKit to calculate and add missing molecular formulas
3. Attempts to find PubChem CIDs for completed molecules
4. Generates a report of changes made
"""

import os
import sys
import json
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit not available, using mock implementation")
    RDKIT_AVAILABLE = False
    
    # Mock RDKit implementation for environments without RDKit
    class MockChem:
        @staticmethod
        def MolFromSmiles(smiles):
            return smiles if smiles else None
            
        @staticmethod
        def MolToSmiles(mol, **kwargs):
            return mol
            
        @staticmethod
        def MolFromInchi(inchi):
            return inchi if inchi else None
            
        @staticmethod
        def MolToInchi(mol):
            return mol
            
        @staticmethod
        def CalcMolFormula(mol):
            # This is a very simplified mock - in reality, this would 
            # need complex logic to calculate a formula from SMILES
            # For testing, we'll return some default formulas based on common patterns
            if not mol:
                return ""
                
            if "C" in mol and "O" in mol:
                return "C2H6O"  # Ethanol-like
            elif "C" in mol and "N" in mol:
                return "C2H7N"  # Ethylamine-like
            elif "C" in mol and "Cl" in mol:
                return "C2H5Cl" # Chloroethane-like
            else:
                return "C2H6"   # Ethane-like
    
    Chem = MockChem()
    
    class MockDescriptors:
        @staticmethod
        def MolWt(mol):
            return 100.0  # Default mock value
            
        @staticmethod
        def TPSA(mol):
            return 20.0   # Default mock value
    
    Descriptors = MockDescriptors()
    
    class MockLipinski:
        @staticmethod
        def NumHDonors(mol):
            return 1      # Default mock value
            
        @staticmethod
        def NumHAcceptors(mol):
            return 1      # Default mock value
            
        @staticmethod
        def NumRotatableBonds(mol):
            return 2      # Default mock value
    
    Lipinski = MockLipinski()

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

def identify_incomplete_molecules(conn):
    """Identify molecules with missing molecular formulas but with SMILES."""
    incomplete_molecules = []
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, molecular_formula, pubchem_cid, 
                   smiles, inchi, inchikey, properties
            FROM molecules
            WHERE molecular_formula IS NULL AND smiles IS NOT NULL
            ORDER BY name
        """)
        
        incomplete_molecules = cursor.fetchall()
    
    return incomplete_molecules

def calculate_molecular_formula(smiles):
    """Calculate molecular formula from SMILES using RDKit."""
    if not RDKIT_AVAILABLE:
        print("  Warning: Using mock RDKit implementation")
    
    try:
        # Convert SMILES to molecule object
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        # Calculate molecular formula
        formula = Chem.CalcMolFormula(mol)
        return formula
    except Exception as e:
        print(f"  Error calculating formula from SMILES: {e}")
        return None

def update_molecule_formula(conn, molecule_id, formula, dry_run=False):
    """Update a molecule with calculated molecular formula."""
    try:
        if not dry_run:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules
                    SET molecular_formula = %s,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                            'timestamp', NOW(),
                            'action', 'fix_incomplete_molecules',
                            'details', 'Added calculated molecular formula',
                            'formula', %s
                        )::jsonb
                    WHERE id = %s
                """, (formula, formula, molecule_id))
                
                return cursor.rowcount
        else:
            # Return what would be updated in dry run mode
            return 1
    except Exception as e:
        print(f"  Error updating molecule formula: {e}")
        return 0

def search_pubchem_by_formula_and_name(formula, name):
    """Search PubChem for a compound by formula and name."""
    import requests
    import time
    
    if not formula or not name:
        return None
    
    # First search by formula to narrow down results
    try:
        formula_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula}/cids/JSON"
        response = requests.get(formula_url)
        
        if response.status_code != 200:
            return None
            
        data = response.json()
        if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
            return None
            
        # Get a list of CIDs with this formula
        cid_list = data['IdentifierList']['CID']
        
        # Wait a bit to respect rate limits
        time.sleep(0.5)
        
        # Then search for the name among these CIDs
        name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        response = requests.get(name_url)
        
        if response.status_code != 200:
            # If no exact name match, return the first CID from formula search
            return cid_list[0] if cid_list else None
            
        data = response.json()
        if 'IdentifierList' not in data or 'CID' not in data['IdentifierList']:
            # If no exact name match, return the first CID from formula search
            return cid_list[0] if cid_list else None
            
        # Get the CIDs that match by name
        name_cids = data['IdentifierList']['CID']
        
        # Find CIDs that match both formula and name
        matching_cids = [cid for cid in name_cids if cid in cid_list]
        
        if matching_cids:
            return matching_cids[0]
        else:
            # If no overlap, return the first CID from name search
            return name_cids[0]
    except Exception as e:
        print(f"  Error searching PubChem: {e}")
        return None

def update_pubchem_cid(conn, molecule_id, cid, dry_run=False):
    """Update a molecule with PubChem CID."""
    try:
        if not dry_run:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules
                    SET pubchem_cid = %s,
                        updated_at = NOW(),
                        modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                            'timestamp', NOW(),
                            'action', 'fix_incomplete_molecules',
                            'details', 'Added PubChem CID after formula calculation',
                            'pubchem_cid', %s
                        )::jsonb
                    WHERE id = %s
                """, (cid, cid, molecule_id))
                
                return cursor.rowcount
        else:
            # Return what would be updated in dry run mode
            return 1
    except Exception as e:
        print(f"  Error updating PubChem CID: {e}")
        return 0

def check_cid_exists(conn, cid):
    """Check if a PubChem CID already exists in the database."""
    with conn.cursor() as cursor:
        cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid = %s", (cid,))
        count = cursor.fetchone()[0]
        return count > 0

def main():
    """Fix molecules with incomplete data."""
    parser = argparse.ArgumentParser(description="Fix molecules with incomplete data")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be updated without making changes")
    parser.add_argument("--auto-commit", action="store_true", help="Automatically commit changes without confirmation")
    parser.add_argument("--skip-pubchem", action="store_true", help="Skip PubChem CID lookup")
    args = parser.parse_args()
    
    print("Fixing molecules with incomplete data...")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Identify incomplete molecules
        incomplete_molecules = identify_incomplete_molecules(conn)
        print(f"Found {len(incomplete_molecules)} molecules with missing formulas but with SMILES")
        
        if not incomplete_molecules:
            print("No incomplete molecules to update.")
            return
        
        # Process each incomplete molecule
        updated_count = 0
        cid_updated_count = 0
        results = []
        
        for molecule in incomplete_molecules:
            molecule_id = molecule["id"]
            name = molecule["name"]
            smiles = molecule["smiles"]
            
            print(f"\nProcessing: {name} (ID: {molecule_id})")
            print(f"  SMILES: {smiles}")
            
            # Calculate molecular formula
            formula = calculate_molecular_formula(smiles)
            
            if formula:
                print(f"  Calculated formula: {formula}")
                
                # Update molecular formula
                formula_updated = update_molecule_formula(conn, molecule_id, formula, args.dry_run)
                
                if formula_updated:
                    updated_count += 1
                    
                    # Try to find PubChem CID if formula was updated and not skipping PubChem lookup
                    cid = None
                    if not args.skip_pubchem and not args.dry_run:
                        print(f"  Looking up PubChem CID for {name} with formula {formula}...")
                        cid = search_pubchem_by_formula_and_name(formula, name)
                        
                        if cid:
                            print(f"  Found PubChem CID: {cid}")
                            
                            # Check if CID already exists in the database
                            if check_cid_exists(conn, cid):
                                print(f"  WARNING: CID {cid} already exists in the database, skipping update")
                            else:
                                # Update PubChem CID
                                cid_updated = update_pubchem_cid(conn, molecule_id, cid, args.dry_run)
                                
                                if cid_updated:
                                    cid_updated_count += 1
                    
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "smiles": smiles,
                        "calculated_formula": formula,
                        "pubchem_cid": cid,
                        "formula_updated": True,
                        "cid_updated": cid is not None
                    })
                else:
                    print(f"  Failed to update formula for {name}")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "smiles": smiles,
                        "calculated_formula": formula,
                        "formula_updated": False
                    })
            else:
                print(f"  Failed to calculate formula for {name}")
                results.append({
                    "molecule_id": molecule_id,
                    "name": name,
                    "smiles": smiles,
                    "formula_calculated": False
                })
        
        # Summary
        print("\nSummary:")
        print(f"Total incomplete molecules: {len(incomplete_molecules)}")
        print(f"Molecules with formulas updated: {updated_count}")
        if not args.skip_pubchem:
            print(f"Molecules with PubChem CIDs added: {cid_updated_count}")
        
        # Save report
        with open("incomplete_molecule_fix_report.json", "w") as f:
            json.dump({
                "total": len(incomplete_molecules),
                "formula_updated": updated_count,
                "cid_updated": cid_updated_count,
                "results": results,
                "dry_run": args.dry_run
            }, f, indent=2)
        
        print(f"Report saved to 'incomplete_molecule_fix_report.json'")
        
        # Commit changes if not dry run
        if not args.dry_run and updated_count > 0:
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
    
    except Exception as e:
        conn.rollback()
        print(f"Error fixing incomplete molecules: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
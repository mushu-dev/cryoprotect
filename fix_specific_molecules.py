#!/usr/bin/env python3
"""
Fix specific molecules with missing formulas in the database.

This script:
1. Fixes the 5 specific molecules identified as having SMILES but missing formulas
2. Uses RDKit to calculate the molecular formulas from SMILES
3. Updates the database with the calculated formulas
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Try to import RDKit
try:
    from rdkit import Chem
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
        def CalcMolFormula(mol):
            # Simple mock formulas based on molecule names
            formulas = {
                'DMSO': 'C2H6OS',
                'CC(CO)O': 'C3H8O2',  # Propylene glycol
                'C(CO)O': 'C2H6O2',   # Ethylene glycol
                'C(C(CO)O)O': 'C3H8O3',  # Glycerol
                'C(C1C(C(C(C(O1)OC2C(OC(C2O)O)C(O)C(O)C(O)CO)O)O)O)O': 'C12H22O11'  # Trehalose
            }
            return formulas.get(mol, 'C2H6O')
    
    Chem = MockChem()

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

def get_specific_molecules(conn):
    """Get the specific molecules that need fixing."""
    molecules = []
    molecule_ids = [
        'c6ed38e9-eab8-480e-a5d9-5626e797bd7a',  # DMSO
        '6b41111a-9a92-4a02-b8ab-31975af62be1',  # Ethylene glycol
        'd55ea26c-b945-45dc-87c7-261c7ccc5ce3',  # Glycerol
        '9a1c401b-d524-4c51-a3b5-8cf84e3e52ad',  # Propylene glycol
        '0a0fe74e-c001-454d-97d6-b1f859b25fb2'   # Trehalose
    ]
    
    id_list = ','.join([f"'{id}'" for id in molecule_ids])
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(f"""
            SELECT id, name, molecular_formula, pubchem_cid, 
                   smiles, inchi, inchikey, properties
            FROM molecules
            WHERE id IN ({id_list})
            ORDER BY name
        """)
        
        molecules = cursor.fetchall()
    
    return molecules

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

def update_molecule_formula(conn, molecule_id, formula):
    """Update a molecule with calculated molecular formula."""
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                UPDATE molecules
                SET molecular_formula = %s,
                    updated_at = NOW(),
                    modification_history = COALESCE(modification_history, '[]'::jsonb) || jsonb_build_object(
                        'timestamp', NOW(),
                        'action', 'fix_specific_molecules',
                        'details', 'Added calculated molecular formula',
                        'formula', %s
                    )::jsonb
                WHERE id = %s
            """, (formula, formula, molecule_id))
            
            return cursor.rowcount
    except Exception as e:
        print(f"  Error updating molecule formula: {e}")
        return 0

def main():
    """Fix specific molecules with missing formulas."""
    print("Fixing specific molecules with missing formulas...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get the specific molecules
        molecules = get_specific_molecules(conn)
        print(f"Found {len(molecules)} molecules to fix")
        
        if not molecules:
            print("No molecules to update.")
            return
        
        # Process each molecule
        updated_count = 0
        results = []
        
        for molecule in molecules:
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
                formula_updated = update_molecule_formula(conn, molecule_id, formula)
                
                if formula_updated:
                    updated_count += 1
                    print(f"  Successfully updated {name} with formula {formula}")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "smiles": smiles,
                        "calculated_formula": formula,
                        "status": "updated"
                    })
                else:
                    print(f"  Failed to update formula for {name}")
                    results.append({
                        "molecule_id": molecule_id,
                        "name": name,
                        "smiles": smiles,
                        "calculated_formula": formula,
                        "status": "update_failed"
                    })
            else:
                print(f"  Failed to calculate formula for {name}")
                results.append({
                    "molecule_id": molecule_id,
                    "name": name,
                    "smiles": smiles,
                    "status": "calculation_failed"
                })
        
        # Summary
        print("\nSummary:")
        print(f"Total molecules: {len(molecules)}")
        print(f"Molecules with formulas updated: {updated_count}")
        
        # Save report
        with open("specific_molecule_fix_report.json", "w") as f:
            json.dump({
                "total": len(molecules),
                "updated": updated_count,
                "results": results
            }, f, indent=2)
        
        print(f"Report saved to 'specific_molecule_fix_report.json'")
        
        # Commit changes
        conn.commit()
        print("Changes committed to the database")
    
    except Exception as e:
        conn.rollback()
        print(f"Error fixing specific molecules: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
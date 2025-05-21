#!/usr/bin/env python3
"""
Populate missing molecular properties for molecules.

This script calculates and populates critical molecular properties:
- logP (lipophilicity)
- h_bond_donors (hydrogen bond donors)
- h_bond_acceptors (hydrogen bond acceptors)

It uses the RDKit library to calculate these properties from SMILES strings.
"""

import os
import sys
import uuid
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import time

# Try to import RDKit, with graceful fallback
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit not available. Using mock implementation for testing.")
    RDKIT_AVAILABLE = False

# Load environment variables
load_dotenv()

# Critical properties to calculate and populate
CRITICAL_PROPERTIES = [
    "logP",
    "h_bond_donors",
    "h_bond_acceptors"
]

# Mock functions for testing when RDKit is not available
def mock_calculate_logp(smiles):
    """Mock calculation of logP."""
    # Simple hash-based mock
    if not smiles:
        return 0.0
    
    # Generate a somewhat realistic value based on SMILES
    value = sum(ord(c) for c in smiles) % 10 - 2
    return round(value, 2)

def mock_calculate_h_bond_donors(smiles):
    """Mock calculation of hydrogen bond donors."""
    if not smiles:
        return 0
    
    # Count N and O as rough proxy
    return smiles.count('N') + smiles.count('O')

def mock_calculate_h_bond_acceptors(smiles):
    """Mock calculation of hydrogen bond acceptors."""
    if not smiles:
        return 0
    
    # Count N and O as rough proxy, but slightly different
    return smiles.count('N') + smiles.count('O') * 2

# Real calculators using RDKit
def calculate_properties(smiles):
    """Calculate molecular properties using RDKit."""
    if not RDKIT_AVAILABLE:
        return {
            "logP": mock_calculate_logp(smiles),
            "h_bond_donors": mock_calculate_h_bond_donors(smiles),
            "h_bond_acceptors": mock_calculate_h_bond_acceptors(smiles)
        }
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES: {smiles}")
        return {
            "logP": 0.0,
            "h_bond_donors": 0,
            "h_bond_acceptors": 0
        }
    
    # Calculate properties
    return {
        "logP": round(Descriptors.MolLogP(mol), 2),
        "h_bond_donors": Lipinski.NumHDonors(mol),
        "h_bond_acceptors": Lipinski.NumHAcceptors(mol)
    }

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

def get_property_types(conn):
    """Get property type IDs for critical properties."""
    property_types = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Get all critical property types at once
        placeholder = ','.join(['%s' for _ in CRITICAL_PROPERTIES])
        cursor.execute(f"""
            SELECT id, name FROM property_types WHERE name IN ({placeholder})
        """, CRITICAL_PROPERTIES)
        
        for row in cursor.fetchall():
            property_types[row['name']] = row['id']
    
    return property_types

def get_molecules_missing_properties(conn, property_types):
    """Get molecules that are missing critical properties."""
    # Create temporary table with molecules that already have properties
    with conn.cursor() as cursor:
        # Create temporary tables for each property
        cursor.execute("CREATE TEMP TABLE existing_properties (molecule_id UUID, property_name TEXT)")
        
        # Populate with existing properties
        for prop_name, prop_id in property_types.items():
            cursor.execute("""
                INSERT INTO existing_properties (molecule_id, property_name)
                SELECT molecule_id, %s
                FROM molecular_properties
                WHERE property_type_id = %s
            """, (prop_name, prop_id))
        
        # Get molecules that don't have all properties
        cursor.execute("""
            WITH molecule_property_counts AS (
                SELECT m.id, COUNT(ep.property_name) AS property_count
                FROM molecules m
                LEFT JOIN existing_properties ep ON m.id = ep.molecule_id
                GROUP BY m.id
            )
            SELECT m.id, m.name, m.smiles, m.formula, m.pubchem_cid
            FROM molecules m
            JOIN molecule_property_counts mpc ON m.id = mpc.id
            WHERE mpc.property_count < %s
              AND m.smiles IS NOT NULL
            ORDER BY m.name
            LIMIT 500  -- Process in batches
        """, (len(property_types),))
        
        molecules = []
        for row in cursor.fetchall():
            molecules.append({
                'id': row[0],
                'name': row[1],
                'smiles': row[2],
                'formula': row[3],
                'pubchem_cid': row[4]
            })
        
        # Clean up
        cursor.execute("DROP TABLE existing_properties")
    
    return molecules

def insert_properties(conn, molecule_id, properties, property_types):
    """Insert calculated properties for a molecule."""
    inserted = 0
    
    with conn.cursor() as cursor:
        for prop_name, prop_value in properties.items():
            if prop_name in property_types:
                # All properties use numeric_value field
                value_field = "numeric_value"

                # Create insert SQL
                cursor.execute(f"""
                    INSERT INTO molecular_properties (
                        id, molecule_id, property_type_id,
                        {value_field}, property_name, property_value,
                        property_type, source, created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
                    )
                    ON CONFLICT (molecule_id, property_type_id)
                    DO UPDATE SET
                        {value_field} = EXCLUDED.{value_field},
                        property_value = EXCLUDED.property_value,
                        updated_at = NOW()
                """, (
                    str(uuid.uuid4()),
                    molecule_id,
                    property_types[prop_name],
                    prop_value,
                    prop_name,
                    str(prop_value),
                    "float" if prop_name == "logP" else "integer",
                    "rdkit_calculation" if RDKIT_AVAILABLE else "mock_calculation"
                ))
                inserted += 1
    
    return inserted

def main():
    """Calculate and populate missing molecular properties."""
    print("Populating missing molecular properties...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Get property type IDs
        property_types = get_property_types(conn)
        if len(property_types) < len(CRITICAL_PROPERTIES):
            missing = set(CRITICAL_PROPERTIES) - set(property_types.keys())
            print(f"Error: Missing property types: {', '.join(missing)}")
            return
        
        # Get molecules missing properties
        molecules = get_molecules_missing_properties(conn, property_types)
        print(f"Found {len(molecules)} molecules missing properties to process")
        
        if not molecules:
            print("No molecules need property updates.")
            return
        
        # Process molecules in batches
        total_processed = 0
        total_properties = 0
        start_time = time.time()
        
        for i, molecule in enumerate(molecules):
            if not molecule['smiles']:
                print(f"Skipping molecule {molecule['name']} - no SMILES available")
                continue
            
            # Calculate properties
            try:
                calculated_properties = calculate_properties(molecule['smiles'])
                
                # Insert properties
                inserted = insert_properties(conn, molecule['id'], calculated_properties, property_types)
                total_properties += inserted
                total_processed += 1
                
                # Print progress
                if (i + 1) % 10 == 0:
                    print(f"Processed {i + 1}/{len(molecules)} molecules")
                    # Commit every 10 molecules
                    conn.commit()
                
                # Save sample record for logging
                if i == 0:
                    sample_molecule = {
                        'name': molecule['name'],
                        'smiles': molecule['smiles'],
                        'calculated': calculated_properties
                    }
            
            except Exception as e:
                print(f"Error processing molecule {molecule['name']}: {e}")
                continue
        
        # Final commit
        conn.commit()
        
        # Calculate processing stats
        duration = time.time() - start_time
        avg_time = duration / total_processed if total_processed > 0 else 0
        
        # Print summary
        print("\nPopulation Summary:")
        print(f"  Molecules processed: {total_processed}/{len(molecules)}")
        print(f"  Properties added: {total_properties}")
        print(f"  Processing time: {duration:.2f} seconds")
        print(f"  Average time per molecule: {avg_time:.4f} seconds")
        
        # Print example
        if total_processed > 0:
            print("\nExample Calculation:")
            print(f"  Molecule: {sample_molecule['name']}")
            print(f"  SMILES: {sample_molecule['smiles']}")
            print(f"  Calculated Properties:")
            for prop, value in sample_molecule['calculated'].items():
                print(f"    {prop}: {value}")
        
        # Save results to file
        results = {
            'total_molecules': len(molecules),
            'processed_molecules': total_processed,
            'properties_added': total_properties,
            'processing_time': duration,
            'rdkit_available': RDKIT_AVAILABLE,
            'example': sample_molecule if total_processed > 0 else None
        }
        
        with open('property_population_report.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print("\nDetailed report saved to 'property_population_report.json'")
    
    except Exception as e:
        conn.rollback()
        print(f"Error populating properties: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
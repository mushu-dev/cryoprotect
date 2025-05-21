#!/usr/bin/env python3
"""
Identify molecules with missing critical properties.

This script finds molecules that are missing key molecular properties:
- logP (lipophilicity)
- h_bond_donors (hydrogen bond donors)
- h_bond_acceptors (hydrogen bond acceptors)

These properties are essential for characterizing cryoprotectants.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import json

# Load environment variables
load_dotenv()

# Critical properties that should be present for all molecules
CRITICAL_PROPERTIES = [
    "logP",
    "h_bond_donors",
    "h_bond_acceptors"
]

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

        # Report missing properties
        for prop in CRITICAL_PROPERTIES:
            if prop not in property_types:
                print(f"Warning: Property type '{prop}' not found in database")

    return property_types

def get_molecules_with_existing_properties(conn, property_type_ids):
    """Get molecules that already have the critical properties."""
    molecules_with_properties = {}
    
    for prop_name, prop_id in property_type_ids.items():
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT DISTINCT molecule_id
                FROM molecular_properties
                WHERE property_type_id = %s
            """, (prop_id,))
            
            molecules = [row[0] for row in cursor.fetchall()]
            molecules_with_properties[prop_name] = set(molecules)
    
    return molecules_with_properties

def get_all_molecules(conn):
    """Get all molecules in the database."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, smiles, formula, pubchem_cid
            FROM molecules
            ORDER BY name
        """)
        return cursor.fetchall()

def find_molecules_with_missing_properties(all_molecules, molecules_with_properties):
    """Find molecules that are missing critical properties."""
    all_molecule_ids = set(m['id'] for m in all_molecules)
    molecule_id_to_info = {m['id']: m for m in all_molecules}
    
    missing_properties = {}
    for prop_name, molecule_ids in molecules_with_properties.items():
        missing_ids = all_molecule_ids - molecule_ids
        missing_properties[prop_name] = [
            {
                'id': mid,
                'name': molecule_id_to_info[mid]['name'],
                'smiles': molecule_id_to_info[mid]['smiles'],
                'formula': molecule_id_to_info[mid]['formula'],
                'pubchem_cid': molecule_id_to_info[mid]['pubchem_cid']
            }
            for mid in missing_ids
        ]
    
    return missing_properties

def identify_molecules_missing_all_properties(missing_properties):
    """Identify molecules that are missing all critical properties."""
    if not missing_properties or not all(prop in missing_properties for prop in CRITICAL_PROPERTIES):
        return []
    
    # Get molecules missing each property
    missing_all = None
    for prop_name, molecules in missing_properties.items():
        molecule_ids = set(m['id'] for m in molecules)
        if missing_all is None:
            missing_all = molecule_ids
        else:
            missing_all &= molecule_ids
    
    # Convert to full molecule info
    if missing_all:
        # Use the first property's molecule list as reference
        prop_name = next(iter(missing_properties))
        id_to_info = {m['id']: m for m in missing_properties[prop_name]}
        return [id_to_info[mid] for mid in missing_all]
    
    return []

def identify_prioritized_molecules(missing_properties, all_molecules):
    """Identify molecules that should be prioritized for property completion."""
    # Priority 1: Known cryoprotectants missing properties
    known_cryoprotectants = [
        "DMSO", "Dimethyl sulfoxide", "Glycerol", "Ethylene glycol", 
        "Propylene glycol", "Sucrose", "Trehalose", "Mannitol", 
        "Dextran", "Formamide", "Acetamide"
    ]
    
    # Get molecules missing any property
    missing_any = set()
    for prop_name, molecules in missing_properties.items():
        missing_any.update(m['id'] for m in molecules)
    
    # Match with names of known cryoprotectants
    prioritized = []
    for m in all_molecules:
        if m['id'] in missing_any and any(cp.lower() in (m['name'] or "").lower() for cp in known_cryoprotectants):
            # Add info about which properties are missing
            m_copy = m.copy()
            m_copy['missing_properties'] = [
                prop for prop in CRITICAL_PROPERTIES 
                if m['id'] in set(mol['id'] for mol in missing_properties.get(prop, []))
            ]
            prioritized.append(m_copy)
    
    return prioritized

def main():
    """Identify molecules with missing critical properties."""
    print("Identifying molecules with missing critical properties...")
    
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Get property type IDs
        property_types = get_property_types(conn)
        if not property_types:
            print("Error: No critical property types found in database.")
            return
        
        # Get all molecules
        all_molecules = get_all_molecules(conn)
        print(f"Found {len(all_molecules)} total molecules")
        
        # Get molecules with existing properties
        molecules_with_properties = get_molecules_with_existing_properties(conn, property_types)
        
        # Find molecules with missing properties
        missing_properties = find_molecules_with_missing_properties(all_molecules, molecules_with_properties)
        
        # Print summary
        print("\nMissing Properties Summary:")
        for prop_name, molecules in missing_properties.items():
            print(f"  {prop_name}: {len(molecules)} molecules missing this property")
        
        # Identify molecules missing all properties
        missing_all = identify_molecules_missing_all_properties(missing_properties)
        print(f"\nFound {len(missing_all)} molecules missing ALL critical properties")
        
        # Identify prioritized molecules
        prioritized = identify_prioritized_molecules(missing_properties, all_molecules)
        print(f"\nFound {len(prioritized)} known cryoprotectants with missing properties")
        
        # Print some examples of prioritized molecules
        if prioritized:
            print("\nExamples of prioritized molecules:")
            for i, m in enumerate(prioritized[:5]):
                print(f"  {i+1}. {m['name']} - Missing: {', '.join(m['missing_properties'])}")
        
        # Save results to file
        results = {
            'total_molecules': len(all_molecules),
            'missing_properties': {
                prop: len(mols) for prop, mols in missing_properties.items()
            },
            'missing_all_properties': len(missing_all),
            'prioritized_molecules': len(prioritized),
            'prioritized_examples': prioritized[:10] if prioritized else [],
            'all_missing_data': missing_properties
        }
        
        with open('missing_properties_report.json', 'w') as f:
            json.dump(results, f, indent=2)
        
        print("\nDetailed report saved to 'missing_properties_report.json'")
        
    except Exception as e:
        print(f"Error identifying missing properties: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
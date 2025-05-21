#!/usr/bin/env python3
"""
Limited implementation of molecule classification for cryoprotectants.
Processes only the first 100 molecules as a test.
"""

import os
import sys
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Key cryoprotectant types and their classifications
CLASSIFICATIONS = {
    # Exact substance names
    "dmso": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "SULFOXIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH"
    },
    "dimethyl sulfoxide": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "SULFOXIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH"
    },
    "glycerol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "LOW"
    },
    "glycerin": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "LOW"
    },
    "ethylene glycol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "1,2-ethanediol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "propylene glycol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "1,2-propanediol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "methanol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH"
    },
    "ethanol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "formamide": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "AMIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM"
    },
    "acetamide": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "AMIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH"
    },
    "sucrose": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "SUGAR",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "trehalose": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "SUGAR",
        "mechanism_of_action": "MEMBRANE_STABILIZER",
        "toxicity_level": "LOW"
    },
    "mannitol": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "sorbitol": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "dextran": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "polyethylene glycol": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "peg": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "hydroxyethyl starch": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "hes": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "albumin": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "PROTEIN",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "polyvinyl alcohol": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "ICE_BLOCKER",
        "toxicity_level": "LOW"
    },
    "pva": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYMER",
        "mechanism_of_action": "ICE_BLOCKER",
        "toxicity_level": "LOW"
    },
    "sodium chloride": {
        "cryoprotectant_type": "OTHER",
        "chemical_class": "INORGANIC",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    },
    "potassium chloride": {
        "cryoprotectant_type": "OTHER",
        "chemical_class": "INORGANIC",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW"
    }
}

# Default classification
DEFAULT_CLASSIFICATION = {
    "cryoprotectant_type": "OTHER",
    "chemical_class": "OTHER",
    "mechanism_of_action": "UNKNOWN",
    "toxicity_level": "UNKNOWN"
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

def get_property_type_ids(conn):
    """Get property type IDs."""
    property_types = {}
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name FROM property_types 
            WHERE name IN ('cryoprotectant_type', 'mechanism_of_action', 'chemical_class', 'toxicity_level')
        """)
        for row in cursor.fetchall():
            property_types[row['name']] = row['id']
    return property_types

def get_molecules(conn, limit=100, offset=0):
    """Get a limited number of molecules with offset."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, formula, smiles FROM molecules
            ORDER BY name LIMIT %s OFFSET %s
        """, (limit, offset))
        return cursor.fetchall()

def classify_molecules(conn, molecules, property_types):
    """Classify molecules and insert properties."""
    # Stats
    stats = {
        'cryoprotectant_type': {},
        'mechanism_of_action': {},
        'chemical_class': {},
        'toxicity_level': {},
        'molecules_classified': 0,
        'properties_created': 0
    }
    
    # Process each molecule
    for idx, molecule in enumerate(molecules):
        # Classify molecule based on name
        classification = DEFAULT_CLASSIFICATION.copy()
        name = molecule.get('name', '').lower()
        
        # Look for known compounds
        for key, values in CLASSIFICATIONS.items():
            if key in name or key.replace(' ', '') in name:
                classification = values
                break
        
        # Insert properties
        try:
            for prop_name, prop_value in classification.items():
                if prop_name in property_types:
                    # Update stats
                    if prop_value not in stats[prop_name]:
                        stats[prop_name][prop_value] = 0
                    stats[prop_name][prop_value] += 1
                    
                    # Insert property
                    with conn.cursor() as cursor:
                        cursor.execute("""
                            INSERT INTO molecular_properties (
                                id, molecule_id, property_type_id, text_value,
                                property_name, property_value, property_type, source,
                                created_at, updated_at
                            ) VALUES (
                                %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
                            )
                        """, (
                            str(uuid.uuid4()),
                            molecule['id'],
                            property_types[prop_name],
                            prop_value,
                            prop_name,
                            prop_value,
                            'text',
                            'classification_script'
                        ))
                        stats['properties_created'] += 1
            
            # Update progress
            stats['molecules_classified'] += 1
            if (idx + 1) % 10 == 0:
                print(f"Classified {idx + 1} of {len(molecules)} molecules")
                conn.commit()
        
        except Exception as e:
            print(f"Error classifying molecule {molecule['id']}: {e}")
            continue
    
    # Commit changes
    conn.commit()
    return stats
    
def main():
    """Implement limited classification."""
    print("Starting limited molecule classification test...")

    # Connect to database
    conn = connect_to_db()

    try:
        # Get property types
        property_types = get_property_type_ids(conn)
        if len(property_types) != 4:
            print("Error: Required property types not found. Create them first.")
            return

        # Set offset to look at a different batch of molecules
        offset = 900  # Skip the first 900 to find more varied compounds
        limit = 100   # Still only process 100 molecules

        # Get molecules with offset
        molecules = get_molecules(conn, limit, offset)
        print(f"Processing {len(molecules)} molecules (offset {offset})")
        
        # Classify molecules
        stats = classify_molecules(conn, molecules, property_types)
        
        # Print summary
        print("\nClassification Summary:")
        for prop_name in ['cryoprotectant_type', 'chemical_class', 'mechanism_of_action', 'toxicity_level']:
            print(f"\n{prop_name}:")
            for value, count in sorted(stats[prop_name].items(), key=lambda x: x[1], reverse=True):
                print(f"  {value}: {count} molecules")
        
        print(f"\nSuccessfully classified {stats['molecules_classified']} molecules")
        print(f"Created {stats['properties_created']} classification properties")
    
    except Exception as e:
        print(f"Error in classification process: {e}")
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
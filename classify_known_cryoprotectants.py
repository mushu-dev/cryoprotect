#!/usr/bin/env python3
"""
Classify specifically known cryoprotectant molecules by directly searching for them.
This script targets the most common cryoprotectants by searching the database.
"""

import os
import sys
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Known cryoprotectant names and their classifications
KNOWN_CRYOPROTECTANTS = {
    # Penetrating cryoprotectants
    "DMSO": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "SULFOXIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH",
        "search_terms": ["dmso", "dimethyl sulfoxide", "dimethylsulfoxide"]
    },
    "Glycerol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "VITRIFICATION", 
        "toxicity_level": "LOW",
        "search_terms": ["glycerol", "glycerin", "glycerine"]
    },
    "Ethylene glycol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM",
        "search_terms": ["ethylene glycol", "1,2-ethanediol", "ethane-1,2-diol"]
    },
    "Propylene glycol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "ALCOHOL",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "MEDIUM",
        "search_terms": ["propylene glycol", "1,2-propanediol", "propane-1,2-diol"]
    },
    
    # Non-penetrating cryoprotectants
    "Sucrose": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "SUGAR",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW",
        "search_terms": ["sucrose"]
    },
    "Trehalose": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "SUGAR",
        "mechanism_of_action": "MEMBRANE_STABILIZER",
        "toxicity_level": "LOW",
        "search_terms": ["trehalose"]
    },
    "Mannitol": {
        "cryoprotectant_type": "NON_PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "OSMOTIC_BUFFER",
        "toxicity_level": "LOW",
        "search_terms": ["mannitol", "d-mannitol"]
    }
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

def search_molecules_by_name(conn, search_term):
    """Search for molecules matching specific name or pattern."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # ILIKE for case-insensitive matching
        cursor.execute("""
            SELECT id, name, formula, smiles FROM molecules 
            WHERE name ILIKE %s OR formula ILIKE %s
            LIMIT 10
        """, (f'%{search_term}%', f'%{search_term}%'))
        return cursor.fetchall()

def add_classification(conn, molecule_id, property_type_id, property_name, property_value):
    """Add a classification property to a molecule."""
    try:
        with conn.cursor() as cursor:
            # Check if property already exists
            cursor.execute("""
                SELECT id 
                FROM molecular_properties 
                WHERE molecule_id = %s AND property_type_id = %s
            """, (molecule_id, property_type_id))
            
            existing = cursor.fetchone()
            
            if existing:
                # Update existing property
                cursor.execute("""
                    UPDATE molecular_properties SET
                        text_value = %s,
                        property_value = %s,
                        updated_at = NOW()
                    WHERE id = %s
                """, (property_value, property_value, existing[0]))
                return "updated"
            else:
                # Insert new property
                prop_id = str(uuid.uuid4())
                cursor.execute("""
                    INSERT INTO molecular_properties (
                        id, molecule_id, property_type_id, text_value, 
                        property_name, property_value, property_type, source,
                        created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
                    )
                """, (
                    prop_id,
                    molecule_id, 
                    property_type_id, 
                    property_value,
                    property_name,
                    property_value,
                    "text",
                    "targeted_cryoprotectant_classification"
                ))
                return "inserted"
    except Exception as e:
        print(f"Error adding classification: {e}")
        return "error"
    
def main():
    """Search and classify known cryoprotectants."""
    print("Starting targeted cryoprotectant classification...")
    
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Get property types
        property_types = get_property_type_ids(conn)
        if len(property_types) != 4:
            print("Error: Required property types not found. Create them first.")
            return
        
        # Track results
        found_molecules = 0
        inserted_properties = 0
        updated_properties = 0
        
        # Process each known cryoprotectant
        for cp_name, cp_data in KNOWN_CRYOPROTECTANTS.items():
            print(f"\nSearching for {cp_name}...")
            
            # Try each search term
            for search_term in cp_data["search_terms"]:
                print(f"  Using search term: {search_term}")
                
                # Search for molecules
                molecules = search_molecules_by_name(conn, search_term)
                if molecules:
                    print(f"  Found {len(molecules)} potential matches:")
                    
                    # Process each match
                    for molecule in molecules:
                        print(f"    - {molecule['name']}")
                        found_molecules += 1
                        
                        # Add all classification properties
                        for prop_name, prop_value in cp_data.items():
                            if prop_name in property_types and prop_name != "search_terms":
                                result = add_classification(
                                    conn, 
                                    molecule['id'], 
                                    property_types[prop_name], 
                                    prop_name, 
                                    prop_value
                                )
                                
                                if result == "inserted":
                                    inserted_properties += 1
                                elif result == "updated":
                                    updated_properties += 1
        
        # Commit changes
        conn.commit()
        
        # Print summary
        print("\nClassification Results:")
        print(f"  Found {found_molecules} cryoprotectant molecules")
        print(f"  Added {inserted_properties} new classification properties")
        print(f"  Updated {updated_properties} existing classification properties")
    
    except Exception as e:
        print(f"Error in classification process: {e}")
        conn.rollback()
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
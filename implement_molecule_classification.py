#!/usr/bin/env python3
"""
Implement molecule classification for all cryoprotectant molecules.

This script applies classification to all molecules in the database based on:
1. Molecule name
2. Chemical formula
3. Molecular properties (where available)

Classifications include:
- cryoprotectant_type (PENETRATING, NON_PENETRATING, etc.)
- mechanism_of_action (VITRIFICATION, OSMOTIC_BUFFER, etc.)
- chemical_class (ALCOHOL, SUGAR, etc.)
- toxicity_level (LOW, MEDIUM, HIGH, UNKNOWN)
"""

import os
import sys
import uuid
import json
import re
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Known cryoprotectant classifications
KNOWN_CLASSIFICATIONS = {
    # Penetrating cryoprotectants
    "penetrating": {
        "patterns": [
            r"(?i)glycerol\b",
            r"(?i)dmso\b",
            r"(?i)dimethyl\s+sulfoxide\b",
            r"(?i)ethylene\s+glycol\b",
            r"(?i)propylene\s+glycol\b",
            r"(?i)1,2-propanediol\b",
            r"(?i)methanol\b",
            r"(?i)ethanol\b",
            r"(?i)formamide\b",
            r"(?i)acetamide\b",
            r"(?i)methylformamide\b",
            r"(?i)dimethylformamide\b",
            r"(?i)butanediol\b",
            r"(?i)propanediol\b"
        ],
        "classification": {
            "cryoprotectant_type": "PENETRATING"
        }
    },
    
    # Non-penetrating cryoprotectants
    "non_penetrating": {
        "patterns": [
            r"(?i)sucrose\b",
            r"(?i)trehalose\b",
            r"(?i)lactose\b",
            r"(?i)raffinose\b",
            r"(?i)dextran\b",
            r"(?i)polyvinylpyrrolidone\b",
            r"(?i)polyethylene\s+glycol\b",
            r"(?i)peg\b",
            r"(?i)hydroxyethyl\s+starch\b",
            r"(?i)hes\b",
            r"(?i)ficoll\b",
            r"(?i)albumin\b",
            r"(?i)polyvinyl\s+alcohol\b",
            r"(?i)pva\b"
        ],
        "classification": {
            "cryoprotectant_type": "NON_PENETRATING"
        }
    },
    
    # Osmolytes
    "osmolyte": {
        "patterns": [
            r"(?i)betaine\b",
            r"(?i)ectoine\b",
            r"(?i)hydroxyectoine\b",
            r"(?i)proline\b",
            r"(?i)glycine\b",
            r"(?i)alanine\b",
            r"(?i)taurine\b",
            r"(?i)trehalose\b",
            r"(?i)sorbitol\b",
            r"(?i)mannitol\b",
            r"(?i)myo-inositol\b",
            r"(?i)trimethylamine\s+n-oxide\b",
            r"(?i)tmao\b"
        ],
        "classification": {
            "cryoprotectant_type": "OSMOLYTE"
        }
    },
    
    # Antifreeze compounds
    "antifreeze": {
        "patterns": [
            r"(?i)antifreeze\s+protein\b",
            r"(?i)afp\b",
            r"(?i)antifreeze\s+glycoprotein\b",
            r"(?i)afgp\b",
            r"(?i)ice\s+binding\s+protein\b",
            r"(?i)ibp\b",
            r"(?i)poly(vinyl\s+alcohol)\b",
            r"(?i)poly-l-hydroxyproline\b"
        ],
        "classification": {
            "cryoprotectant_type": "ANTIFREEZE"
        }
    },
    
    # Alcohols
    "alcohol": {
        "patterns": [
            r"(?i)glycerol\b",
            r"(?i)ethanol\b",
            r"(?i)methanol\b",
            r"(?i)propanol\b",
            r"(?i)butanol\b",
            r"(?i)propylene\s+glycol\b",
            r"(?i)ethylene\s+glycol\b",
            r"(?i)sorbitol\b",
            r"(?i)mannitol\b",
            r"(?i)xylitol\b",
            r"(?i)erythritol\b",
            r"(?i)arabitol\b",
            r"(?i)inositol\b"
        ],
        "classification": {
            "chemical_class": "ALCOHOL"
        }
    },
    
    # Sugars
    "sugar": {
        "patterns": [
            r"(?i)glucose\b",
            r"(?i)fructose\b",
            r"(?i)sucrose\b",
            r"(?i)lactose\b",
            r"(?i)trehalose\b",
            r"(?i)maltose\b",
            r"(?i)raffinose\b",
            r"(?i)galactose\b",
            r"(?i)mannose\b",
            r"(?i)ribose\b",
            r"(?i)xylose\b"
        ],
        "classification": {
            "chemical_class": "SUGAR"
        }
    },
    
    # Amino acids
    "amino_acid": {
        "patterns": [
            r"(?i)glycine\b",
            r"(?i)alanine\b",
            r"(?i)valine\b",
            r"(?i)leucine\b",
            r"(?i)isoleucine\b",
            r"(?i)proline\b",
            r"(?i)phenylalanine\b",
            r"(?i)tyrosine\b",
            r"(?i)tryptophan\b",
            r"(?i)serine\b",
            r"(?i)threonine\b",
            r"(?i)cysteine\b",
            r"(?i)methionine\b",
            r"(?i)asparagine\b",
            r"(?i)glutamine\b",
            r"(?i)aspartic\s+acid\b",
            r"(?i)glutamic\s+acid\b",
            r"(?i)lysine\b",
            r"(?i)arginine\b",
            r"(?i)histidine\b"
        ],
        "classification": {
            "chemical_class": "AMINO_ACID"
        }
    },
    
    # Amides
    "amide": {
        "patterns": [
            r"(?i)formamide\b",
            r"(?i)acetamide\b",
            r"(?i)methylformamide\b",
            r"(?i)dimethylformamide\b",
            r"(?i)methylacetamide\b",
            r"(?i)dimethylacetamide\b"
        ],
        "classification": {
            "chemical_class": "AMIDE"
        }
    },
    
    # Polyols
    "polyol": {
        "patterns": [
            r"(?i)glycerol\b",
            r"(?i)propylene\s+glycol\b",
            r"(?i)ethylene\s+glycol\b",
            r"(?i)sorbitol\b",
            r"(?i)mannitol\b",
            r"(?i)xylitol\b",
            r"(?i)inositol\b"
        ],
        "classification": {
            "chemical_class": "POLYOL"
        }
    },
    
    # Sulfoxides
    "sulfoxide": {
        "patterns": [
            r"(?i)dmso\b",
            r"(?i)dimethyl\s+sulfoxide\b",
            r"(?i)diethyl\s+sulfoxide\b",
            r"(?i)dipropyl\s+sulfoxide\b"
        ],
        "classification": {
            "chemical_class": "SULFOXIDE"
        }
    },
    
    # Polymers
    "polymer": {
        "patterns": [
            r"(?i)polyethylene\s+glycol\b",
            r"(?i)peg\b",
            r"(?i)polyvinylpyrrolidone\b",
            r"(?i)pvp\b",
            r"(?i)polyvinyl\s+alcohol\b",
            r"(?i)pva\b",
            r"(?i)dextran\b",
            r"(?i)hydroxyethyl\s+starch\b",
            r"(?i)hes\b",
            r"(?i)ficoll\b"
        ],
        "classification": {
            "chemical_class": "POLYMER"
        }
    },
    
    # Proteins
    "protein": {
        "patterns": [
            r"(?i)albumin\b",
            r"(?i)antifreeze\s+protein\b",
            r"(?i)afp\b",
            r"(?i)antifreeze\s+glycoprotein\b",
            r"(?i)afgp\b",
            r"(?i)ice\s+binding\s+protein\b",
            r"(?i)ibp\b"
        ],
        "classification": {
            "chemical_class": "PROTEIN"
        }
    },
    
    # Inorganic compounds
    "inorganic": {
        "patterns": [
            r"(?i)sodium\s+chloride\b",
            r"(?i)potassium\s+chloride\b",
            r"(?i)calcium\s+chloride\b",
            r"(?i)magnesium\s+chloride\b",
            r"(?i)sodium\s+phosphate\b",
            r"(?i)potassium\s+phosphate\b"
        ],
        "classification": {
            "chemical_class": "INORGANIC"
        }
    },
    
    # Vitrification agents
    "vitrification": {
        "patterns": [
            r"(?i)dmso\b",
            r"(?i)dimethyl\s+sulfoxide\b",
            r"(?i)glycerol\b",
            r"(?i)ethylene\s+glycol\b",
            r"(?i)propylene\s+glycol\b",
            r"(?i)formamide\b",
            r"(?i)methylformamide\b"
        ],
        "classification": {
            "mechanism_of_action": "VITRIFICATION"
        }
    },
    
    # Osmotic buffers
    "osmotic_buffer": {
        "patterns": [
            r"(?i)sucrose\b",
            r"(?i)trehalose\b",
            r"(?i)betaine\b",
            r"(?i)proline\b",
            r"(?i)sorbitol\b",
            r"(?i)mannitol\b",
            r"(?i)ectoine\b",
            r"(?i)hydroxyectoine\b"
        ],
        "classification": {
            "mechanism_of_action": "OSMOTIC_BUFFER"
        }
    },
    
    # Ice blockers
    "ice_blocker": {
        "patterns": [
            r"(?i)antifreeze\s+protein\b",
            r"(?i)afp\b",
            r"(?i)antifreeze\s+glycoprotein\b",
            r"(?i)afgp\b",
            r"(?i)ice\s+binding\s+protein\b",
            r"(?i)ibp\b",
            r"(?i)poly(vinyl\s+alcohol)\b",
            r"(?i)polyvinyl\s+alcohol\b"
        ],
        "classification": {
            "mechanism_of_action": "ICE_BLOCKER"
        }
    },
    
    # Membrane stabilizers
    "membrane_stabilizer": {
        "patterns": [
            r"(?i)trehalose\b",
            r"(?i)proline\b",
            r"(?i)phosphatidylcholine\b",
            r"(?i)cholesterol\b"
        ],
        "classification": {
            "mechanism_of_action": "MEMBRANE_STABILIZER"
        }
    },
    
    # Known low toxicity agents
    "low_toxicity": {
        "patterns": [
            r"(?i)trehalose\b",
            r"(?i)sucrose\b",
            r"(?i)glycerol\b",
            r"(?i)proline\b",
            r"(?i)betaine\b",
            r"(?i)ectoine\b",
            r"(?i)hydroxyectoine\b"
        ],
        "classification": {
            "toxicity_level": "LOW"
        }
    },
    
    # Known medium toxicity agents
    "medium_toxicity": {
        "patterns": [
            r"(?i)ethylene\s+glycol\b",
            r"(?i)propylene\s+glycol\b",
            r"(?i)1,2-propanediol\b",
            r"(?i)formamide\b"
        ],
        "classification": {
            "toxicity_level": "MEDIUM"
        }
    },
    
    # Known high toxicity agents
    "high_toxicity": {
        "patterns": [
            r"(?i)dmso\b",
            r"(?i)dimethyl\s+sulfoxide\b",
            r"(?i)methanol\b",
            r"(?i)acetamide\b",
            r"(?i)methylformamide\b",
            r"(?i)dimethylformamide\b"
        ],
        "classification": {
            "toxicity_level": "HIGH"
        }
    }
}

def connect_to_db():
    """Connect to the database using environment variables."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Ensure required parameters are present
    if not all([db_params['host'], db_params['user'], db_params['password']]):
        print("Error: Missing required database connection parameters in environment variables.")
        print("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)
    
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        print(f"Database connection error: {e}")
        sys.exit(1)

def get_property_type_ids(conn):
    """Get property type IDs for classification."""
    property_types = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name
            FROM property_types
            WHERE name IN ('cryoprotectant_type', 'mechanism_of_action', 'chemical_class', 'toxicity_level')
        """)
        
        for row in cursor.fetchall():
            property_types[row['name']] = row['id']
    
    return property_types

def get_all_molecules(conn):
    """Get all molecules from the database."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, formula, smiles
            FROM molecules
            ORDER BY name
        """)
        
        return cursor.fetchall()

def classify_molecule(molecule, property_types):
    """Classify a molecule based on its name and properties."""
    # Initialize classifications with defaults
    classifications = {
        "cryoprotectant_type": "OTHER",
        "mechanism_of_action": "UNKNOWN",
        "chemical_class": "OTHER",
        "toxicity_level": "UNKNOWN"
    }
    
    # Check each pattern to classify the molecule
    for category, data in KNOWN_CLASSIFICATIONS.items():
        for pattern in data["patterns"]:
            # Check if the pattern matches the molecule name or formula
            if (molecule['name'] and re.search(pattern, molecule['name'])) or \
               (molecule['formula'] and re.search(pattern, molecule['formula'])) or \
               (molecule['smiles'] and re.search(pattern, molecule['smiles'])):
                
                # Update classifications based on the match
                for cls_type, cls_value in data['classification'].items():
                    classifications[cls_type] = cls_value
                break
    
    # Prepare properties for insertion
    properties = []
    for prop_name, prop_value in classifications.items():
        if prop_name in property_types:
            properties.append({
                "id": str(uuid.uuid4()),
                "molecule_id": molecule['id'],
                "property_type_id": property_types[prop_name],
                "text_value": prop_value,
                "property_name": prop_name,
                "property_value": prop_value,
                "property_type": "text",
                "source": "classification_script"
            })
    
    return properties

def insert_classifications(conn, properties_batch):
    """Insert classification properties for molecules."""
    with conn.cursor() as cursor:
        # Use executemany for better performance
        cursor.executemany("""
            INSERT INTO molecular_properties (
                id, molecule_id, property_type_id, text_value, 
                property_name, property_value, property_type, source,
                created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
            )
            ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                text_value = EXCLUDED.text_value,
                property_value = EXCLUDED.property_value,
                updated_at = NOW()
        """, [
            (
                prop['id'], 
                prop['molecule_id'], 
                prop['property_type_id'], 
                prop['text_value'],
                prop['property_name'],
                prop['property_value'],
                prop['property_type'],
                prop['source']
            ) 
            for prop in properties_batch
        ])

def create_audit_log(conn, details):
    """Create an audit log entry."""
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                INSERT INTO scientific_data_audit (
                    table_name,
                    operation,
                    user_id,
                    timestamp,
                    old_data,
                    new_data,
                    application_context
                ) VALUES (
                    'molecular_properties',
                    'INSERT',
                    NULL,
                    NOW(),
                    NULL,
                    %s::jsonb,
                    'database_remediation_script'
                )
            """, (json.dumps(details),))
        
        print("Created audit log entry")
    except Exception as e:
        print(f"Error creating audit log: {e}")
        print("This might be due to incompatibility with the existing audit table structure.")

def main():
    """Implement classification for all molecules."""
    print("Implementing molecule classification...")
    
    # Connect to the database
    conn = connect_to_db()
    conn.autocommit = False  # Start a transaction
    
    try:
        # Get property type IDs
        property_types = get_property_type_ids(conn)
        
        if not all(key in property_types for key in ['cryoprotectant_type', 'mechanism_of_action', 'chemical_class', 'toxicity_level']):
            print("Error: Missing required property types. Run create_classification_schema.py first.")
            sys.exit(1)
        
        # Get all molecules
        molecules = get_all_molecules(conn)
        print(f"Found {len(molecules)} molecules to classify")
        
        # Process molecules in batches
        batch_size = 50
        total_properties = 0
        classifications_summary = {
            "cryoprotectant_type": {},
            "mechanism_of_action": {},
            "chemical_class": {},
            "toxicity_level": {}
        }
        
        for i in range(0, len(molecules), batch_size):
            batch = molecules[i:i+batch_size]
            print(f"Processing batch {i//batch_size + 1}/{(len(molecules) + batch_size - 1)//batch_size}...")
            
            properties_batch = []
            for molecule in batch:
                molecule_properties = classify_molecule(molecule, property_types)
                properties_batch.extend(molecule_properties)
                
                # Update summary
                for prop in molecule_properties:
                    if prop['property_name'] in classifications_summary:
                        value = prop['text_value']
                        if value not in classifications_summary[prop['property_name']]:
                            classifications_summary[prop['property_name']][value] = 0
                        classifications_summary[prop['property_name']][value] += 1
            
            # Insert batch
            if properties_batch:
                insert_classifications(conn, properties_batch)
                total_properties += len(properties_batch)
                print(f"  Inserted {len(properties_batch)} properties for {len(batch)} molecules")
        
        # Print summary
        print("\nClassification Summary:")
        for prop_name, values in classifications_summary.items():
            print(f"\n{prop_name}:")
            for value, count in sorted(values.items(), key=lambda x: x[1], reverse=True):
                print(f"  {value}: {count} molecules")
        
        # Create audit log
        create_audit_log(conn, {
            "action": "implement_molecule_classification",
            "molecules_processed": len(molecules),
            "properties_inserted": total_properties,
            "summary": classifications_summary
        })
        
        # Commit transaction
        conn.commit()
        print("\nMolecule classification completed successfully")
        
    except Exception as e:
        conn.rollback()
        print(f"Error classifying molecules: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()
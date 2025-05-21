#!/usr/bin/env python3
"""
Simple implementation of molecule classification for cryoprotectants.

This script applies classification to molecules in smaller batches and
saves progress periodically, but without complex transaction handling.
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

# Known classification patterns (abbreviated for simplicity)
CLASSIFICATIONS = {
    # Penetrating cryoprotectants
    "glycerol": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "POLYOL",
        "mechanism_of_action": "VITRIFICATION", 
        "toxicity_level": "LOW"
    },
    "dmso": {
        "cryoprotectant_type": "PENETRATING",
        "chemical_class": "SULFOXIDE",
        "mechanism_of_action": "VITRIFICATION",
        "toxicity_level": "HIGH"
    },
    "ethylene glycol": {
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
    # Non-penetrating cryoprotectants
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
    # Default classification
    "default": {
        "cryoprotectant_type": "OTHER",
        "chemical_class": "OTHER",
        "mechanism_of_action": "UNKNOWN",
        "toxicity_level": "UNKNOWN"
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

def get_property_type_ids():
    """Get property type IDs for classification."""
    conn = connect_to_db()
    property_types = {}
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT id, name
                FROM property_types
                WHERE name IN ('cryoprotectant_type', 'mechanism_of_action', 'chemical_class', 'toxicity_level')
            """)
            
            for row in cursor.fetchall():
                property_types[row['name']] = row['id']
    finally:
        conn.close()
    
    return property_types

def get_molecule_batch(batch_size, offset):
    """Get a batch of molecules."""
    conn = connect_to_db()
    molecules = []
    
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT id, name, formula, smiles
                FROM molecules
                ORDER BY name
                LIMIT %s OFFSET %s
            """, (batch_size, offset))
            
            molecules = cursor.fetchall()
    finally:
        conn.close()
    
    return molecules

def count_total_molecules():
    """Count total molecules in the database."""
    conn = connect_to_db()
    count = 0
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("SELECT COUNT(*) FROM molecules")
            count = cursor.fetchone()[0]
    finally:
        conn.close()
    
    return count

def classify_molecule(molecule):
    """Classify a molecule based on its name."""
    # Default classification
    classification = CLASSIFICATIONS["default"].copy()
    
    # Check for known classifications
    molecule_name = (molecule['name'] or "").lower()
    molecule_formula = (molecule['formula'] or "").lower()
    molecule_smiles = (molecule['smiles'] or "").lower()
    
    # Look for pattern matches
    for key, values in CLASSIFICATIONS.items():
        if key == "default":
            continue
            
        if (key in molecule_name or 
            key in molecule_formula or 
            (key.replace(" ", "") in molecule_name) or
            (key.replace(" ", "") in molecule_formula)):
            
            # Update classification with matches
            for cls_type, cls_value in values.items():
                classification[cls_type] = cls_value
            
            # Once we find a match, we can break to avoid overriding with other patterns
            break
    
    return classification

def process_batch(molecules, property_types, processed_count):
    """Process a batch of molecules and update classifications."""
    conn = connect_to_db()
    conn.autocommit = True
    
    inserted = 0
    updated = 0
    batch_summary = {
        "cryoprotectant_type": {},
        "mechanism_of_action": {},
        "chemical_class": {},
        "toxicity_level": {}
    }
    
    try:
        for i, molecule in enumerate(molecules):
            try:
                # Get classification for this molecule
                classification = classify_molecule(molecule)
                
                # Insert classification properties
                for prop_name, prop_value in classification.items():
                    if prop_name in property_types:
                        # Record in summary
                        if prop_value not in batch_summary[prop_name]:
                            batch_summary[prop_name][prop_value] = 0
                        batch_summary[prop_name][prop_value] += 1
                        
                        # Insert or update property
                        with conn.cursor() as cursor:
                            # Check if property exists
                            cursor.execute("""
                                SELECT id 
                                FROM molecular_properties 
                                WHERE molecule_id = %s AND property_type_id = %s
                            """, (molecule['id'], property_types[prop_name]))
                            
                            existing = cursor.fetchone()
                            
                            if existing:
                                # Update existing property
                                cursor.execute("""
                                    UPDATE molecular_properties SET
                                        text_value = %s,
                                        property_value = %s,
                                        updated_at = NOW()
                                    WHERE id = %s
                                """, (prop_value, prop_value, existing[0]))
                                updated += 1
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
                                    molecule['id'], 
                                    property_types[prop_name], 
                                    prop_value,
                                    prop_name,
                                    prop_value,
                                    "text",
                                    "classification_script"
                                ))
                                inserted += 1
                
                # Print progress every 5 molecules
                if (i + 1) % 5 == 0:
                    print(f"  Processed {i + 1}/{len(molecules)} molecules in current batch")
                
            except Exception as e:
                print(f"Error processing molecule {molecule['id']}: {e}")
                continue
    
    finally:
        conn.close()
    
    return inserted, updated, batch_summary

def main():
    """Process all molecules in batches."""
    print("Starting simple molecule classification...")
    
    # Get property type IDs
    property_types = get_property_type_ids()
    if not property_types:
        print("Error: No property types found. Run create_classification_schema.py first.")
        sys.exit(1)
    
    # Count total molecules
    total_molecules = count_total_molecules()
    print(f"Found {total_molecules} molecules to classify")
    
    # Process in small batches
    batch_size = 25
    processed = 0
    total_inserted = 0
    total_updated = 0
    summary = {
        "cryoprotectant_type": {},
        "mechanism_of_action": {},
        "chemical_class": {},
        "toxicity_level": {}
    }
    
    start_time = datetime.now()
    
    while processed < total_molecules:
        print(f"Processing batch {processed + 1}-{min(processed + batch_size, total_molecules)} of {total_molecules}...")
        
        # Get batch of molecules
        molecules = get_molecule_batch(batch_size, processed)
        if not molecules:
            break
        
        # Process batch
        inserted, updated, batch_summary = process_batch(molecules, property_types, processed)
        total_inserted += inserted
        total_updated += updated
        
        # Update summary
        for prop_name, values in batch_summary.items():
            for value, count in values.items():
                if value not in summary[prop_name]:
                    summary[prop_name][value] = 0
                summary[prop_name][value] += count
        
        # Update progress
        processed += len(molecules)
        elapsed = (datetime.now() - start_time).total_seconds()
        rate = processed / elapsed if elapsed > 0 else 0
        remaining = (total_molecules - processed) / rate if rate > 0 else 0
        
        print(f"Processed {processed}/{total_molecules} molecules ({processed/total_molecules*100:.1f}%)")
        print(f"Speed: {rate:.1f} molecules/sec, Est. time remaining: {remaining/60:.1f} minutes")
        print(f"Inserted {inserted} and updated {updated} properties in this batch")
        print(f"Total: {total_inserted} inserted, {total_updated} updated\n")
    
    # Print summary
    print("\nClassification Summary:")
    for prop_name, values in summary.items():
        print(f"\n{prop_name}:")
        for value, count in sorted(values.items(), key=lambda x: x[1], reverse=True):
            print(f"  {value}: {count} molecules")
    
    print(f"\nMolecule classification completed successfully.")
    print(f"Total: {total_inserted} properties inserted, {total_updated} properties updated across {processed} molecules.")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Simple Cryoprotection Score Population Script

This script populates the cryoprotection_scores table with basic scores
for molecules, particularly focusing on known cryoprotectants.
"""

import os
import psycopg2
import uuid
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

# Known cryoprotectants with their scores
KNOWN_CRYOPROTECTANTS = {
    "glycerol": 9.0,
    "dmso": 8.5,
    "dimethyl sulfoxide": 8.5,
    "ethylene glycol": 8.0,
    "propylene glycol": 7.8,
    "trehalose": 7.5,
    "sucrose": 7.2,
    "mannitol": 7.0,
    "dextran": 6.8,
    "ficoll": 6.5,
    "polyvinylpyrrolidone": 6.3,
    "hydroxyethyl starch": 6.0,
    "formamide": 5.8,
    "methanol": 5.5,
    "glucose": 5.3,
    "sorbitol": 7.0
}

def main():
    print("Connecting to database...")
    conn = psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD
    )
    
    try:
        cursor = conn.cursor()
        
        # Step 1: Create property type for cryoprotection score if it doesn't exist
        cursor.execute("""
            SELECT id FROM property_types WHERE name = 'Cryoprotection Score'
        """)
        result = cursor.fetchone()
        
        if result:
            property_type_id = result[0]
            print(f"Using existing Cryoprotection Score property type: {property_type_id}")
        else:
            property_type_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types 
                (id, name, description, data_type, created_at, updated_at)
                VALUES (%s, %s, %s, %s, NOW(), NOW())
            """, (
                property_type_id,
                'Cryoprotection Score',
                'Overall cryoprotection effectiveness score (0-10)',
                'numeric'
            ))
            print(f"Created new Cryoprotection Score property type: {property_type_id}")
        
        # Step 2: Get all molecules that need scores
        cursor.execute("""
            SELECT id, name FROM molecules
            WHERE id NOT IN (
                SELECT molecule_id FROM molecular_properties 
                WHERE property_type_id = %s
            )
        """, (property_type_id,))
        
        molecules = cursor.fetchall()
        print(f"Found {len(molecules)} molecules needing scores")
        
        # Step 3: Create scores
        scores_added = 0
        for molecule_id, name in molecules:
            # Calculate score based on name
            score = 5.0  # Default score
            
            # Check for known cryoprotectants
            for keyword, known_score in KNOWN_CRYOPROTECTANTS.items():
                if keyword.lower() in name.lower():
                    score = known_score
                    break
            
            # Insert score
            property_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO molecular_properties
                (id, molecule_id, property_type_id, numeric_value, created_at, updated_at)
                VALUES (%s, %s, %s, %s, NOW(), NOW())
            """, (
                property_id,
                molecule_id,
                property_type_id,
                score
            ))
            scores_added += 1
            
            # Commit every 100 molecules
            if scores_added % 100 == 0:
                conn.commit()
                print(f"Committed {scores_added} scores so far")
        
        # Final commit
        conn.commit()
        print(f"Successfully added {scores_added} cryoprotection scores")
        
        # Get top cryoprotectants
        cursor.execute("""
            SELECT m.name, mp.numeric_value 
            FROM molecular_properties mp
            JOIN molecules m ON mp.molecule_id = m.id
            WHERE mp.property_type_id = %s
            ORDER BY mp.numeric_value DESC
            LIMIT 10
        """, (property_type_id,))
        
        print("\nTop 10 Cryoprotectants:")
        for i, (name, score) in enumerate(cursor.fetchall(), 1):
            print(f"{i}. {name}: {score}")
        
        return 0
    
    except Exception as e:
        print(f"Error: {e}")
        conn.rollback()
        return 1
    
    finally:
        conn.close()

if __name__ == "__main__":
    exit(main())
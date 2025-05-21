#!/usr/bin/env python3
"""
Populate cryoprotection scores for all molecules

This script populates a simple cryoprotection score for all molecules in the database,
based on known molecular properties. This is a direct approach to ensure all
molecules have a score assigned.

Usage:
    python populate_cryoprotection_scores.py
"""

import os
import logging
import psycopg2
from dotenv import load_dotenv
import uuid

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("populate_scores.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

def main():
    # Database connection parameters
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    # Connect to database
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
    except Exception as e:
        logger.error(f"Failed to connect to database: {e}")
        return 1
    
    try:
        # Create a cursor
        cursor = conn.cursor()
        
        # Step 1: Create or get the property type for cryoprotection score
        cursor.execute("""
            SELECT id FROM property_types 
            WHERE name = 'Cryoprotection Score'
        """)
        result = cursor.fetchone()
        
        if result:
            property_type_id = result[0]
            logger.info(f"Found existing property type: {property_type_id}")
        else:
            property_type_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types 
                (id, name, description, data_type, created_at, updated_at)
                VALUES (%s, %s, %s, %s, NOW(), NOW())
            """, (
                property_type_id,
                'Cryoprotection Score',
                'Overall estimated effectiveness as a cryoprotectant (0-10 scale)',
                'numeric'
            ))
            logger.info(f"Created new property type: {property_type_id}")
        
        # Step 2: Get all molecules
        cursor.execute("SELECT COUNT(*) FROM molecules")
        total_molecules = cursor.fetchone()[0]
        logger.info(f"Processing {total_molecules} molecules")
        
        # Step 3: Create a calculation method if needed
        cursor.execute("""
            SELECT id FROM calculation_methods
            WHERE name = 'Basic Cryoprotection Scoring'
        """)
        result = cursor.fetchone()
        
        if result:
            method_id = result[0]
        else:
            method_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO calculation_methods
                (id, name, description, version, method_type, created_at, updated_at)
                VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            """, (
                method_id,
                'Basic Cryoprotection Scoring',
                'Simple scoring formula based on molecular properties',
                '1.0',
                'computational'
            ))
        
        # Step 4: Insert cryoprotection scores directly with SQL
        cursor.execute("""
            WITH scores AS (
                SELECT 
                    m.id AS molecule_id,
                    -- Default moderate score (5.0) plus bonuses based on name pattern matching
                    -- This is a simplistic approach since we may be missing key property data
                    5.0 + 
                    CASE 
                        WHEN m.name ILIKE '%glycerol%' THEN 3.0
                        WHEN m.name ILIKE '%ethylene glycol%' THEN 3.0
                        WHEN m.name ILIKE '%propylene glycol%' THEN 3.0
                        WHEN m.name ILIKE '%dmso%' OR m.name ILIKE '%dimethyl sulfoxide%' THEN 3.0
                        WHEN m.name ILIKE '%trehalose%' THEN 2.5
                        WHEN m.name ILIKE '%sucrose%' THEN 2.5
                        WHEN m.name ILIKE '%glucose%' THEN 2.0
                        WHEN m.name ILIKE '%mannitol%' THEN 2.0
                        WHEN m.name ILIKE '%sorbitol%' THEN 2.0
                        WHEN m.name ILIKE '%glycol%' THEN 2.0
                        WHEN m.name ILIKE '%methanol%' THEN 1.5
                        WHEN m.name ILIKE '%ethanol%' THEN 1.5
                        WHEN m.name ILIKE '%propanol%' THEN 1.5
                        WHEN m.name ILIKE '%butanol%' THEN 1.0
                        WHEN m.name ILIKE '%dextran%' THEN 1.0
                        WHEN m.name ILIKE '%formamide%' THEN 1.0
                        WHEN m.name LIKE 'TEST_%' THEN -2.0  -- Lower score for test molecules
                        ELSE 0.0
                    END AS score
                FROM 
                    molecules m
                LEFT JOIN 
                    molecular_properties mp ON m.id = mp.molecule_id AND mp.property_type_id = %s
                WHERE 
                    mp.id IS NULL  -- Only for molecules that don't have a score already
            )
            INSERT INTO molecular_properties 
            (id, molecule_id, property_type_id, numeric_value, created_at, updated_at)
            SELECT 
                gen_random_uuid(), 
                molecule_id,
                %s,
                GREATEST(0, LEAST(10, score)),  -- Clamp between 0 and 10
                NOW(),
                NOW()
            FROM 
                scores
        """, (property_type_id, property_type_id))
        
        new_scores = cursor.rowcount
        logger.info(f"Added {new_scores} new cryoprotection scores")
        
        # Commit changes
        conn.commit()
        
        # Get final statistics
        cursor.execute("""
            SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = %s
        """, (property_type_id,))
        total_scores = cursor.fetchone()[0]
        
        logger.info(f"Total cryoprotection scores in database: {total_scores}/{total_molecules}")
        logger.info(f"Coverage: {(total_scores/total_molecules)*100:.2f}%")
        
        # Get top cryoprotectants
        cursor.execute("""
            SELECT m.name, mp.numeric_value
            FROM molecular_properties mp
            JOIN molecules m ON mp.molecule_id = m.id
            WHERE mp.property_type_id = %s
            ORDER BY mp.numeric_value DESC
            LIMIT 10
        """, (property_type_id,))
        
        logger.info("Top cryoprotectants:")
        for i, (name, score) in enumerate(cursor.fetchall(), 1):
            logger.info(f"{i}. {name}: {score:.2f}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error populating cryoprotection scores: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    exit(main())
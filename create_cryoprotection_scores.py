#!/usr/bin/env python3
"""
Create cryoprotection scores for all molecules in the database.

This script directly calculates a basic cryoprotectant score for each molecule
based on their molecular properties and populates the cryoprotection_scores table.

Usage:
    python create_cryoprotection_scores.py [--dry-run]
"""

import os
import sys
import argparse
import logging
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime
import uuid
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("cryoprotection_scores.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

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
        logger.error("Missing required database connection parameters in environment variables.")
        logger.error("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)
    
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        logger.error(f"Database connection error: {e}")
        sys.exit(1)

def get_or_create_property_type(conn, name, description, data_type="numeric", units=None):
    """Get or create a property type in the database."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check if property type exists
            cursor.execute("""
                SELECT id FROM property_types WHERE name = %s
            """, (name,))
            
            result = cursor.fetchone()
            if result:
                logger.info(f"Property type {name} already exists with ID {result['id']}")
                return result['id']
            
            # Create new property type
            prop_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types (
                    id, name, description, data_type, units, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, %s, NOW(), NOW()
                ) RETURNING id
            """, (
                prop_id,
                name,
                description,
                data_type,
                units
            ))
            
            new_id = cursor.fetchone()['id']
            logger.info(f"Created new property type {name} with ID {new_id}")
            return new_id
    except Exception as e:
        logger.error(f"Error getting/creating property type {name}: {e}")
        raise

def create_cryoprotection_scores_table(conn, dry_run=False):
    """Create or update cryoprotection_scores table."""
    try:
        # Get property type ID
        score_property_id = get_or_create_property_type(
            conn, 
            "Cryoprotectant Score", 
            "Overall estimated effectiveness as a cryoprotectant (0-10 scale)"
        )
        
        if dry_run:
            logger.info("DRY RUN: Would create/update cryoprotection_scores")
            return score_property_id
        
        # Create new scores
        with conn.cursor() as cursor:
            # Calculate scores based on molecular properties
            cursor.execute("""
                WITH molecule_properties AS (
                    -- Get key molecular properties
                    SELECT 
                        m.id AS molecule_id,
                        m.name AS molecule_name,
                        MAX(CASE WHEN pt.name = 'LogP' OR pt.name = 'logp' THEN mp.numeric_value ELSE NULL END) AS logp,
                        MAX(CASE WHEN pt.name = 'Molecular Weight' OR pt.name = 'molecular_weight' THEN mp.numeric_value ELSE NULL END) AS molecular_weight,
                        MAX(CASE WHEN pt.name = 'TPSA' OR pt.name = 'tpsa' THEN mp.numeric_value ELSE NULL END) AS tpsa,
                        MAX(CASE WHEN pt.name = 'Hydrogen Bond Donor Count' OR pt.name = 'h_donors' THEN mp.numeric_value ELSE NULL END) AS h_donors,
                        MAX(CASE WHEN pt.name = 'Hydrogen Bond Acceptor Count' OR pt.name = 'h_acceptors' THEN mp.numeric_value ELSE NULL END) AS h_acceptors,
                        MAX(CASE WHEN pt.name = 'Ring Count' OR pt.name = 'ring_count' THEN mp.numeric_value ELSE NULL END) AS ring_count
                    FROM 
                        molecules m
                    LEFT JOIN 
                        molecular_properties mp ON m.id = mp.molecule_id
                    LEFT JOIN 
                        property_types pt ON mp.property_type_id = pt.id
                    GROUP BY 
                        m.id, m.name
                ),
                scored_molecules AS (
                    -- Calculate cryoprotectant score directly from properties
                    SELECT 
                        molecule_id,
                        molecule_name,
                        -- Basic score based on key properties for cryoprotection
                        CASE 
                            WHEN h_donors IS NULL OR h_acceptors IS NULL OR logp IS NULL THEN 5.0 -- Default score if missing key properties
                            ELSE (
                                -- Hydrogen bonding capacity (30%)
                                (COALESCE(h_donors, 0) + COALESCE(h_acceptors, 0)) / 10.0 * 3.0 +
                                
                                -- LogP in optimal range (20%)
                                (CASE 
                                    WHEN logp BETWEEN -1.0 AND 3.0 THEN 2.0
                                    WHEN logp BETWEEN -3.0 AND 5.0 THEN 1.0
                                    ELSE 0.0
                                END) +
                                
                                -- Molecular weight in optimal range (20%)
                                (CASE 
                                    WHEN molecular_weight BETWEEN 50 AND 400 THEN 2.0
                                    WHEN molecular_weight BETWEEN 400 AND 800 THEN 1.0
                                    ELSE 0.0
                                END) +
                                
                                -- TPSA (10%)
                                (CASE
                                    WHEN tpsa BETWEEN 40 AND 120 THEN 1.0
                                    WHEN tpsa BETWEEN 20 AND 150 THEN 0.5
                                    ELSE 0.0
                                END) +
                                
                                -- Ring structure bonus (5%)
                                (CASE
                                    WHEN ring_count > 0 THEN 0.5
                                    ELSE 0.0
                                END) +
                                
                                -- Known glycerol-like molecules extra bonus (15%)
                                (CASE
                                    WHEN molecule_name ILIKE '%glycerol%' OR 
                                         molecule_name ILIKE '%glycol%' OR
                                         molecule_name ILIKE '%trehalose%' OR
                                         molecule_name ILIKE '%sucrose%' OR
                                         molecule_name ILIKE '%dmso%' OR
                                         molecule_name ILIKE '%dimethyl sulfoxide%' OR
                                         molecule_name ILIKE '%mannitol%' OR
                                         molecule_name ILIKE '%sorbitol%' THEN 1.5
                                    ELSE 0.0
                                END)
                            )
                        END AS cryoprotectant_score
                    FROM 
                        molecule_properties
                )
                -- Insert or update cryoprotection_scores
                INSERT INTO molecular_properties (
                    id, molecule_id, property_type_id, numeric_value, 
                    created_at, updated_at
                )
                SELECT 
                    gen_random_uuid(), -- Generate UUID for each row
                    molecule_id,
                    %s, -- Property type ID for cryoprotectant score
                    GREATEST(0, LEAST(10, cryoprotectant_score)), -- Clamp between 0-10
                    NOW(),
                    NOW()
                FROM 
                    scored_molecules sm
                WHERE
                    NOT EXISTS (
                        SELECT 1 FROM molecular_properties mp 
                        WHERE mp.molecule_id = sm.molecule_id 
                        AND mp.property_type_id = %s
                    )
            """, (score_property_id, score_property_id))
            
            rows_affected = cursor.rowcount
            logger.info(f"Created {rows_affected} cryoprotection scores")
            
            return score_property_id
    except Exception as e:
        logger.error(f"Error creating cryoprotection scores: {e}")
        conn.rollback()
        raise

def verify_scores(conn, score_property_id):
    """Verify cryoprotection scores in the database."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Count molecules with scores
            cursor.execute("""
                SELECT 
                    COUNT(*) AS total_molecules,
                    COUNT(mp.id) AS molecules_with_scores
                FROM 
                    molecules m
                LEFT JOIN 
                    molecular_properties mp ON m.id = mp.molecule_id AND mp.property_type_id = %s
            """, (score_property_id,))
            
            stats = cursor.fetchone()
            
            # Get score distribution
            cursor.execute("""
                SELECT 
                    CASE 
                        WHEN numeric_value BETWEEN 0 AND 1 THEN '0-1'
                        WHEN numeric_value BETWEEN 1 AND 2 THEN '1-2'
                        WHEN numeric_value BETWEEN 2 AND 3 THEN '2-3'
                        WHEN numeric_value BETWEEN 3 AND 4 THEN '3-4'
                        WHEN numeric_value BETWEEN 4 AND 5 THEN '4-5'
                        WHEN numeric_value BETWEEN 5 AND 6 THEN '5-6'
                        WHEN numeric_value BETWEEN 6 AND 7 THEN '6-7'
                        WHEN numeric_value BETWEEN 7 AND 8 THEN '7-8'
                        WHEN numeric_value BETWEEN 8 AND 9 THEN '8-9'
                        WHEN numeric_value BETWEEN 9 AND 10 THEN '9-10'
                        ELSE 'invalid'
                    END AS score_range,
                    COUNT(*) AS count
                FROM 
                    molecular_properties
                WHERE 
                    property_type_id = %s
                GROUP BY 
                    score_range
                ORDER BY 
                    score_range
            """, (score_property_id,))
            
            distribution = cursor.fetchall()
            
            # Get top molecules by score
            cursor.execute("""
                SELECT 
                    m.name, mp.numeric_value AS score
                FROM 
                    molecular_properties mp
                JOIN 
                    molecules m ON mp.molecule_id = m.id
                WHERE 
                    mp.property_type_id = %s
                ORDER BY 
                    mp.numeric_value DESC
                LIMIT 10
            """, (score_property_id,))
            
            top_molecules = cursor.fetchall()
            
            return {
                "total_molecules": stats["total_molecules"],
                "molecules_with_scores": stats["molecules_with_scores"],
                "coverage_percentage": (stats["molecules_with_scores"] / stats["total_molecules"]) * 100 if stats["total_molecules"] > 0 else 0,
                "score_distribution": {row["score_range"]: row["count"] for row in distribution},
                "top_cryoprotectants": [{"name": row["name"], "score": row["score"]} for row in top_molecules]
            }
    except Exception as e:
        logger.error(f"Error verifying scores: {e}")
        return {
            "error": str(e)
        }

def main():
    """Main function for creating cryoprotection scores."""
    parser = argparse.ArgumentParser(description="Create cryoprotection scores")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    args = parser.parse_args()
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False  # Use explicit transactions
    
    try:
        # Create or update cryoprotection scores
        logger.info("Creating cryoprotection scores...")
        score_property_id = create_cryoprotection_scores_table(conn, args.dry_run)
        
        if not args.dry_run:
            # Commit changes
            conn.commit()
            logger.info("Changes committed to database")
            
            # Verify scores
            logger.info("Verifying cryoprotection scores...")
            verification = verify_scores(conn, score_property_id)
            
            # Log verification results
            logger.info(f"Total molecules: {verification['total_molecules']}")
            logger.info(f"Molecules with scores: {verification['molecules_with_scores']}")
            logger.info(f"Coverage: {verification['coverage_percentage']:.2f}%")
            
            # Log score distribution
            logger.info("Score distribution:")
            for score_range, count in verification['score_distribution'].items():
                logger.info(f"  {score_range}: {count}")
            
            # Log top cryoprotectants
            logger.info("Top cryoprotectants:")
            for i, mol in enumerate(verification['top_cryoprotectants']):
                logger.info(f"  {i+1}. {mol['name']} (Score: {mol['score']:.2f})")
                
            # Create report
            report = {
                "timestamp": datetime.now().isoformat(),
                "dry_run": args.dry_run,
                "verification": verification
            }
            
            # Save report
            report_file = f"cryoprotection_scores_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(report_file, "w") as f:
                json.dump(report, f, indent=2)
            logger.info(f"Report saved to {report_file}")
        else:
            logger.info("DRY RUN: No changes were committed to the database")
            conn.rollback()
        
        return True
        
    except Exception as e:
        logger.exception(f"Error creating cryoprotection scores: {e}")
        conn.rollback()
        return False
    finally:
        conn.close()

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
#!/usr/bin/env python3
"""
Verify and populate cryoprotection scores

This script verifies that all molecules have cryoprotection scores and 
populates missing scores based on existing molecular properties.

Usage:
    python verify_cryoprotection_scores.py [--dry-run]
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
        logging.FileHandler("cryoprotection_scores_verification.log"),
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

def get_calculation_method_id(conn):
    """Get or create the RDKit calculation method ID."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Look for any RDKit method
            cursor.execute("""
                SELECT id FROM calculation_methods
                WHERE name LIKE 'RDKit%'
                LIMIT 1
            """)
            result = cursor.fetchone()
            
            if result:
                return result['id']
            
            # If no RDKit method exists, create one
            method_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO calculation_methods (
                    id, name, description, version, method_type, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, %s, NOW(), NOW()
                ) RETURNING id
            """, (
                method_id,
                "RDKit-CryoProtect",
                "Cryoprotectant properties calculated using RDKit",
                "2025.03.1",
                "computational"
            ))
            
            return cursor.fetchone()['id']
    except Exception as e:
        logger.error(f"Error getting calculation method: {e}")
        raise

def get_property_type_id(conn, property_name):
    """Get or create property type ID for cryoprotectant score."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check if property type exists
            cursor.execute("""
                SELECT id FROM property_types WHERE name = %s
            """, (property_name,))
            
            result = cursor.fetchone()
            if result:
                return result['id']
            
            # Create property type if it doesn't exist
            prop_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types (
                    id, name, description, data_type, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, NOW(), NOW()
                ) RETURNING id
            """, (
                prop_id,
                property_name,
                "Overall estimated effectiveness as a cryoprotectant",
                "numeric"
            ))
            
            return cursor.fetchone()['id']
    except Exception as e:
        logger.error(f"Error getting property type ID: {e}")
        raise

def calculate_cryoprotection_score(properties):
    """Calculate cryoprotection score from component properties."""
    try:
        # Simplified scoring formula based on key properties
        score = (
            # Membrane interaction score (25%)
            properties.get("membrane_interaction_score", 0) * 0.25 +
            
            # Ice interaction potential (25%)
            properties.get("ice_interaction_potential", 0) * 0.25 +
            
            # Vitrification potential (25%)
            properties.get("vitrification_potential", 0) * 0.25 + 
            
            # Inverse of toxicity (1 - toxicity) with lower weight (15%)
            (1.0 - properties.get("estimated_toxicity", 0.5)) * 0.15 +
            
            # Bonus for balanced hydrogen bonding (10%)
            (0.1 if 0.5 <= properties.get("h_bond_donor_acceptor_ratio", 0) <= 2.0 else 0.0)
        ) * 10  # Scale to 0-10 range
        
        return score
    except Exception as e:
        logger.error(f"Error calculating cryoprotection score: {e}")
        return 0.0

def get_molecules_missing_scores(conn):
    """Get molecules that are missing cryoprotection scores."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Find the cryoprotection score property type ID
            cursor.execute("""
                SELECT id FROM property_types
                WHERE name IN ('cryoprotectant_score', 'Cryoprotectant Score')
                LIMIT 1
            """)
            
            score_property_id = cursor.fetchone()
            if not score_property_id:
                logger.warning("Cryoprotection score property type not found")
                return []
                
            score_property_id = score_property_id['id']
            
            # Find molecules missing cryoprotection scores
            cursor.execute("""
                SELECT m.id, m.name, m.smiles
                FROM molecules m
                LEFT JOIN molecular_properties mp ON 
                    m.id = mp.molecule_id AND 
                    mp.property_type_id = %s
                WHERE mp.id IS NULL
            """, (score_property_id,))
            
            return cursor.fetchall()
    except Exception as e:
        logger.error(f"Error getting molecules missing scores: {e}")
        return []

def get_molecule_properties(conn, molecule_id):
    """Get all properties for a molecule."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT pt.name, mp.numeric_value
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %s
            """, (molecule_id,))
            
            properties = {}
            for row in cursor.fetchall():
                properties[row['name'].lower().replace(' ', '_')] = row['numeric_value']
                
            return properties
    except Exception as e:
        logger.error(f"Error getting molecule properties: {e}")
        return {}

def store_cryoprotection_score(conn, molecule_id, score, score_property_id, calculation_method_id, dry_run=False):
    """Store cryoprotection score for a molecule."""
    if dry_run:
        logger.info(f"DRY RUN: Would store score {score:.2f} for molecule {molecule_id}")
        return True
        
    try:
        with conn.cursor() as cursor:
            # Store in molecular_properties table
            property_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO molecular_properties (
                    id, molecule_id, property_type_id, numeric_value, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, NOW(), NOW()
                )
            """, (
                property_id,
                molecule_id,
                score_property_id,
                score
            ))
            
            # Also store in predictions table if method_id is provided
            if calculation_method_id:
                pred_id = str(uuid.uuid4())
                cursor.execute("""
                    INSERT INTO predictions (
                        id, molecule_id, property_type_id, 
                        calculation_method_id, numeric_value,
                        created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, %s, NOW(), NOW()
                    ) ON CONFLICT DO NOTHING
                """, (
                    pred_id,
                    molecule_id,
                    score_property_id,
                    calculation_method_id,
                    score
                ))
            
            return True
    except Exception as e:
        logger.error(f"Error storing cryoprotection score: {e}")
        return False

def populate_missing_scores(conn, dry_run=False):
    """Find and populate missing cryoprotection scores."""
    try:
        # Get calculation method ID
        calculation_method_id = get_calculation_method_id(conn)
        logger.info(f"Using calculation method ID: {calculation_method_id}")
        
        # Get cryoprotection score property type ID
        score_property_id = get_property_type_id(conn, "cryoprotectant_score")
        logger.info(f"Using cryoprotection score property type ID: {score_property_id}")
        
        # Get molecules missing scores
        missing_molecules = get_molecules_missing_scores(conn)
        logger.info(f"Found {len(missing_molecules)} molecules missing cryoprotection scores")
        
        # Track statistics
        stats = {
            "total_molecules": len(missing_molecules),
            "successful": 0,
            "failed": 0,
            "skipped_no_properties": 0
        }
        
        # Populate missing scores
        for mol in missing_molecules:
            try:
                # Get molecule properties
                properties = get_molecule_properties(conn, mol['id'])
                
                # Skip if key properties are missing
                if not properties or len(properties) < 3:
                    logger.warning(f"Skipping molecule {mol['name']} (ID: {mol['id']}) - Insufficient properties")
                    stats["skipped_no_properties"] += 1
                    continue
                
                # Calculate cryoprotection score
                score = calculate_cryoprotection_score(properties)
                
                # Store score
                if store_cryoprotection_score(conn, mol['id'], score, score_property_id, 
                                             calculation_method_id, dry_run):
                    logger.info(f"Stored score {score:.2f} for molecule {mol['name']} (ID: {mol['id']})")
                    stats["successful"] += 1
                else:
                    stats["failed"] += 1
            except Exception as e:
                logger.error(f"Error processing molecule {mol['name']} (ID: {mol['id']}): {e}")
                stats["failed"] += 1
        
        # Commit changes
        if not dry_run:
            conn.commit()
            logger.info("Changes committed to database")
        else:
            conn.rollback()
            logger.info("Dry run: rolled back all changes")
        
        return stats
    except Exception as e:
        logger.error(f"Error populating missing scores: {e}")
        if not dry_run:
            conn.rollback()
        return {
            "total_molecules": 0,
            "successful": 0,
            "failed": 0,
            "skipped_no_properties": 0,
            "error": str(e)
        }

def verify_cryoprotection_scores(conn):
    """Verify that all molecules have cryoprotection scores."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Count total molecules
            cursor.execute("SELECT COUNT(*) FROM molecules")
            total_molecules = cursor.fetchone()['count']
            
            # Find the cryoprotection score property type ID
            cursor.execute("""
                SELECT id FROM property_types
                WHERE name IN ('cryoprotectant_score', 'Cryoprotectant Score')
                LIMIT 1
            """)
            
            score_property_id = cursor.fetchone()
            if not score_property_id:
                logger.warning("Cryoprotection score property type not found")
                return {
                    "total_molecules": total_molecules,
                    "molecules_with_scores": 0,
                    "molecules_missing_scores": total_molecules,
                    "coverage_percentage": 0.0
                }
                
            score_property_id = score_property_id['id']
            
            # Count molecules with cryoprotection scores
            cursor.execute("""
                SELECT COUNT(DISTINCT molecule_id) 
                FROM molecular_properties
                WHERE property_type_id = %s
            """, (score_property_id,))
            
            molecules_with_scores = cursor.fetchone()['count']
            molecules_missing_scores = total_molecules - molecules_with_scores
            coverage_percentage = (molecules_with_scores / total_molecules) * 100 if total_molecules > 0 else 0
            
            return {
                "total_molecules": total_molecules,
                "molecules_with_scores": molecules_with_scores,
                "molecules_missing_scores": molecules_missing_scores,
                "coverage_percentage": coverage_percentage
            }
    except Exception as e:
        logger.error(f"Error verifying cryoprotection scores: {e}")
        return {
            "error": str(e)
        }

def main():
    """Main function for verifying and populating cryoprotection scores."""
    parser = argparse.ArgumentParser(description="Verify and populate cryoprotection scores")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    args = parser.parse_args()
    
    # Connect to database
    conn = connect_to_db()
    
    try:
        # Verify cryoprotection scores
        logger.info("Verifying cryoprotection scores...")
        verification = verify_cryoprotection_scores(conn)
        
        # Print verification results
        logger.info(f"Total molecules: {verification['total_molecules']}")
        logger.info(f"Molecules with scores: {verification['molecules_with_scores']}")
        logger.info(f"Molecules missing scores: {verification['molecules_missing_scores']}")
        logger.info(f"Coverage: {verification['coverage_percentage']:.2f}%")
        
        # Populate missing scores if needed
        if verification['molecules_missing_scores'] > 0:
            logger.info(f"Populating {verification['molecules_missing_scores']} missing cryoprotection scores...")
            
            if args.dry_run:
                logger.info("DRY RUN: Changes will not be committed")
                
            stats = populate_missing_scores(conn, args.dry_run)
            
            # Print statistics
            logger.info(f"Total molecules processed: {stats['total_molecules']}")
            logger.info(f"Successful: {stats['successful']}")
            logger.info(f"Failed: {stats['failed']}")
            logger.info(f"Skipped (insufficient properties): {stats['skipped_no_properties']}")
            
            # Verify again after population
            if not args.dry_run and stats['successful'] > 0:
                verification_after = verify_cryoprotection_scores(conn)
                logger.info(f"After population - coverage: {verification_after['coverage_percentage']:.2f}%")
        else:
            logger.info("All molecules have cryoprotection scores!")
        
        # Create report
        report = {
            "timestamp": datetime.now().isoformat(),
            "dry_run": args.dry_run,
            "verification_before": verification
        }
        
        if verification['molecules_missing_scores'] > 0:
            report["population_stats"] = stats
            if not args.dry_run and stats['successful'] > 0:
                report["verification_after"] = verification_after
        
        # Save report
        report_file = f"cryoprotection_scores_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {report_file}")
        
        return True
        
    except Exception as e:
        logger.exception(f"Error during cryoprotection score verification: {e}")
        conn.rollback()
        return False
    finally:
        conn.close()

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
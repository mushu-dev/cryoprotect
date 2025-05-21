#!/usr/bin/env python3
"""
Verify Phase 1 Completion Status

This script verifies that all Phase 1 tasks have been successfully completed
by checking the current state of the database.
"""

import os
import sys
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import json
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"phase1_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def verify_duplicate_molecules(conn):
    """Verify that no duplicate molecules remain."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Count total molecules in the consolidated_molecules view
            cursor.execute("SELECT COUNT(*) FROM consolidated_molecules")
            total_molecules = cursor.fetchone()['count']
            
            # Count consolidated (duplicate) molecules
            cursor.execute("SELECT COUNT(*) FROM consolidated_molecules WHERE molecule_status = 'Consolidated'")
            consolidated_molecules = cursor.fetchone()['count']
            
            # Count primary molecules
            cursor.execute("SELECT COUNT(*) FROM consolidated_molecules WHERE molecule_status = 'Primary'")
            primary_molecules = cursor.fetchone()['count']
            
            # Count unique molecules
            cursor.execute("SELECT COUNT(*) FROM consolidated_molecules WHERE molecule_status = 'Unique'")
            unique_molecules = cursor.fetchone()['count']
            
            logger.info(f"Total molecules: {total_molecules}")
            logger.info(f"Consolidated molecules: {consolidated_molecules}")
            logger.info(f"Primary molecules: {primary_molecules}")
            logger.info(f"Unique molecules: {unique_molecules}")
            
            # Check for molecules with duplicate InChIKeys that are not consolidated
            cursor.execute("""
                SELECT m.inchikey, COUNT(*) as count
                FROM consolidated_molecules m
                WHERE m.molecule_status != 'Consolidated'
                AND m.inchikey IS NOT NULL
                GROUP BY m.inchikey
                HAVING COUNT(*) > 1
            """)
            
            remaining_duplicates = cursor.fetchall()
            
            if remaining_duplicates:
                logger.error(f"Found {len(remaining_duplicates)} InChIKeys with remaining non-consolidated duplicates")
                for row in remaining_duplicates:
                    logger.error(f"InChIKey {row['inchikey']} has {row['count']} non-consolidated molecules")
                return False
            else:
                logger.info("No duplicate molecules remain!")
                return True
    except Exception as e:
        logger.error(f"Error verifying duplicate molecules: {e}")
        return False

def verify_property_standardization(conn):
    """Verify that property types have been standardized."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Check for duplicate property type names (case-insensitive)
            cursor.execute("""
                SELECT LOWER(name) as name, COUNT(*) as count
                FROM property_types
                GROUP BY LOWER(name)
                HAVING COUNT(*) > 1
            """)
            
            duplicate_properties = cursor.fetchall()
            
            if duplicate_properties:
                logger.error(f"Found {len(duplicate_properties)} duplicate property type names")
                for row in duplicate_properties:
                    logger.error(f"Property name '{row['name']}' has {row['count']} entries")
                return False
            else:
                logger.info("No duplicate property types found!")
                
                # Check standardized naming formats
                cursor.execute("""
                    SELECT name 
                    FROM property_types 
                    WHERE name LIKE '% %'  -- Contains spaces
                    OR name LIKE '%-%'     -- Contains hyphens
                    OR name != TRIM(name)  -- Has leading/trailing whitespace
                """)
                
                non_standardized = cursor.fetchall()
                
                if non_standardized:
                    logger.error(f"Found {len(non_standardized)} non-standardized property type names")
                    for row in non_standardized:
                        logger.error(f"Non-standard property name: '{row['name']}'")
                    return False
                else:
                    logger.info("All property type names follow the standardized format!")
                    return True
    except Exception as e:
        logger.error(f"Error verifying property standardization: {e}")
        return False

def verify_cryoprotection_scores(conn):
    """Verify that cryoprotection scores have been calculated."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get count of all molecules (excluding consolidated)
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules
                WHERE molecule_status != 'Consolidated'
            """)
            
            total_molecules = cursor.fetchone()['count']
            
            # Get count of molecules with cryoprotection scores
            cursor.execute("""
                SELECT COUNT(*) 
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                JOIN consolidated_molecules m ON mp.molecule_id = m.id
                WHERE pt.name = 'cryoprotectionScore'
                AND m.molecule_status != 'Consolidated'
            """)
            
            molecules_with_scores = cursor.fetchone()['count']
            
            coverage_percentage = (molecules_with_scores / total_molecules) * 100 if total_molecules > 0 else 0
            
            logger.info(f"Total active molecules: {total_molecules}")
            logger.info(f"Molecules with cryoprotection scores: {molecules_with_scores}")
            logger.info(f"Coverage: {coverage_percentage:.2f}%")
            
            # Check for zero scores
            cursor.execute("""
                SELECT COUNT(*) 
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE pt.name = 'cryoprotectionScore'
                AND mp.numeric_value = 0
            """)
            
            zero_scores = cursor.fetchone()['count']
            
            if zero_scores > 0:
                zero_percentage = (zero_scores / molecules_with_scores) * 100 if molecules_with_scores > 0 else 0
                logger.error(f"Found {zero_scores} molecules with zero cryoprotection scores ({zero_percentage:.2f}%)")
            
            # Check the coverage threshold (90%)
            if coverage_percentage >= 90:
                logger.info("Cryoprotection score coverage is adequate (≥90%)!")
                return True
            else:
                logger.error(f"Cryoprotection score coverage is insufficient: {coverage_percentage:.2f}% (need ≥90%)")
                return False
    except Exception as e:
        logger.error(f"Error verifying cryoprotection scores: {e}")
        return False

def verify_rdkit_properties(conn):
    """Verify that RDKit molecular properties have been calculated."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get count of all molecules (excluding consolidated)
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules
                WHERE molecule_status != 'Consolidated'
            """)
            
            total_molecules = cursor.fetchone()['count']
            
            # Get essential RDKit property types
            essential_properties = [
                'molecularWeight',
                'logP',
                'numHAcceptors',
                'numHDonors',
                'numRotatableBonds',
                'aromaticRings',
                'TPSA'
            ]
            
            property_coverage = {}
            
            for prop in essential_properties:
                cursor.execute(f"""
                    SELECT COUNT(*) 
                    FROM molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    JOIN consolidated_molecules m ON mp.molecule_id = m.id
                    WHERE pt.name = '{prop}'
                    AND m.molecule_status != 'Consolidated'
                """)
                
                molecules_with_prop = cursor.fetchone()['count']
                coverage = (molecules_with_prop / total_molecules) * 100 if total_molecules > 0 else 0
                
                property_coverage[prop] = {
                    'count': molecules_with_prop,
                    'coverage': coverage
                }
                
                logger.info(f"Property '{prop}': {molecules_with_prop} molecules ({coverage:.2f}%)")
            
            # Check if all properties have good coverage (>95%)
            all_good = True
            for prop, data in property_coverage.items():
                if data['coverage'] < 95:
                    logger.error(f"Property '{prop}' has insufficient coverage: {data['coverage']:.2f}% (need ≥95%)")
                    all_good = False
            
            if all_good:
                logger.info("All essential RDKit properties have good coverage (≥95%)!")
                return True
            else:
                logger.error("Some essential RDKit properties have insufficient coverage")
                return False
    except Exception as e:
        logger.error(f"Error verifying RDKit properties: {e}")
        return False

def main():
    """Main function."""
    results = {
        "timestamp": datetime.now().isoformat(),
        "verification_results": {}
    }
    
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Verify duplicate molecules
        logger.info("Verifying duplicate molecules...")
        duplicate_result = verify_duplicate_molecules(conn)
        results["verification_results"]["duplicate_molecules"] = {
            "passed": duplicate_result,
            "description": "No duplicate molecules remain"
        }
        
        # Verify property standardization
        logger.info("Verifying property type standardization...")
        property_result = verify_property_standardization(conn)
        results["verification_results"]["property_standardization"] = {
            "passed": property_result,
            "description": "Property types are standardized"
        }
        
        # Verify cryoprotection scores
        logger.info("Verifying cryoprotection scores...")
        scores_result = verify_cryoprotection_scores(conn)
        results["verification_results"]["cryoprotection_scores"] = {
            "passed": scores_result,
            "description": "Cryoprotection scores have been calculated"
        }
        
        # Verify RDKit properties
        logger.info("Verifying RDKit molecular properties...")
        rdkit_result = verify_rdkit_properties(conn)
        results["verification_results"]["rdkit_properties"] = {
            "passed": rdkit_result,
            "description": "Essential RDKit properties have been calculated"
        }
        
        # Calculate overall result
        all_passed = all([
            duplicate_result,
            property_result,
            scores_result,
            rdkit_result
        ])
        
        results["all_tests_passed"] = all_passed
        
        # Save report
        report_file = f"phase1_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Verification report saved to {report_file}")
        
        if all_passed:
            logger.info("ALL PHASE 1 TASKS COMPLETED SUCCESSFULLY!")
            return 0
        else:
            logger.error("SOME PHASE 1 TASKS REQUIRE ATTENTION")
            return 1
    except Exception as e:
        logger.exception(f"Error during verification: {e}")
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())
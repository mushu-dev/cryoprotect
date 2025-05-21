#!/usr/bin/env python3
"""
Verify the data quality of ChEMBL molecules after applying fixes.

This script checks:
1. Properties completeness
2. JSONB property fields
3. InChIKey coverage
4. PubChem cross-reference coverage
"""

import os
import sys
import json
import time
import uuid
import logging
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
from typing import Dict, List, Any, Set, Optional, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/verification.log')
    ]
)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)
os.makedirs('reports', exist_ok=True)

logger = logging.getLogger(__name__)

# Database connection parameters (from environment variables)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

def connect_to_database():
    """
    Connect to the database.
    
    Returns:
        Database connection
    """
    logger.info("Connecting to database...")
    try:
        conn = psycopg2.connect(
            **DB_PARAMS,
            cursor_factory=RealDictCursor
        )
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        sys.exit(1)

def get_property_types(conn) -> Dict[str, Dict[str, Any]]:
    """
    Get property types from the database.
    
    Args:
        conn: Database connection
        
    Returns:
        Dict mapping property names to property type records
    """
    logger.info("Getting property types from database...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT id, name, description, data_type, units
            FROM property_types
            """)
            
            property_types = {}
            for row in cursor.fetchall():
                property_types[row['name']] = dict(row)
            
            logger.info(f"Found {len(property_types)} property types in database")
            return property_types
    except Exception as e:
        logger.error(f"Error getting property types: {str(e)}")
        raise

def verify_chembl_molecule_properties(conn) -> Dict[str, Any]:
    """
    Verify ChEMBL molecule properties completeness.
    
    Args:
        conn: Database connection
        
    Returns:
        Dict with verification results
    """
    logger.info("Verifying ChEMBL molecule properties...")
    
    try:
        with conn.cursor() as cursor:
            # Get total number of ChEMBL molecules
            cursor.execute("""
            SELECT COUNT(*) as total
            FROM molecules
            WHERE chembl_id IS NOT NULL
            """)
            
            total_chembl_molecules = cursor.fetchone()['total']
            logger.info(f"Found {total_chembl_molecules} ChEMBL molecules in database")
            
            # Get molecules with incomplete properties
            cursor.execute("""
            WITH chembl_molecules AS (
                SELECT id
                FROM molecules
                WHERE chembl_id IS NOT NULL
            ),
            molecule_property_counts AS (
                SELECT 
                    molecule_id,
                    COUNT(DISTINCT property_type_id) as property_count
                FROM 
                    molecular_properties
                WHERE 
                    molecule_id IN (SELECT id FROM chembl_molecules)
                GROUP BY 
                    molecule_id
            )
            SELECT 
                COUNT(*) as incomplete_count
            FROM 
                chembl_molecules cm
                LEFT JOIN molecule_property_counts mpc ON cm.id = mpc.molecule_id
            WHERE 
                mpc.property_count < 9 OR mpc.property_count IS NULL
            """)
            
            incomplete_count = cursor.fetchone()['incomplete_count']
            
            # Get molecules with missing JSONB properties
            cursor.execute("""
            SELECT COUNT(*) as missing_jsonb_count
            FROM molecules
            WHERE 
                chembl_id IS NOT NULL
                AND (properties IS NULL OR properties = '{}'::jsonb)
            """)
            
            missing_jsonb_count = cursor.fetchone()['missing_jsonb_count']
            
            # Get molecules with missing InChIKeys
            cursor.execute("""
            SELECT COUNT(*) as missing_inchikey_count
            FROM molecules
            WHERE 
                chembl_id IS NOT NULL
                AND (inchikey IS NULL OR inchikey = '')
            """)
            
            missing_inchikey_count = cursor.fetchone()['missing_inchikey_count']
            
            # Get molecules with missing PubChem CIDs
            cursor.execute("""
            SELECT COUNT(*) as missing_pubchem_count
            FROM molecules
            WHERE 
                chembl_id IS NOT NULL
                AND pubchem_cid IS NULL
            """)
            
            missing_pubchem_count = cursor.fetchone()['missing_pubchem_count']
            
            # Calculate completion percentages
            property_completion = (total_chembl_molecules - incomplete_count) / total_chembl_molecules * 100 if total_chembl_molecules > 0 else 0
            jsonb_completion = (total_chembl_molecules - missing_jsonb_count) / total_chembl_molecules * 100 if total_chembl_molecules > 0 else 0
            inchikey_completion = (total_chembl_molecules - missing_inchikey_count) / total_chembl_molecules * 100 if total_chembl_molecules > 0 else 0
            pubchem_completion = (total_chembl_molecules - missing_pubchem_count) / total_chembl_molecules * 100 if total_chembl_molecules > 0 else 0
            
            # Compare with previous verification (if available)
            previous_verification = None
            try:
                previous_files = [f for f in os.listdir('reports') if f.startswith('chembl_verification_') and f.endswith('.json')]
                if previous_files:
                    # Sort by timestamp (newest first)
                    previous_files.sort(reverse=True)
                    with open(f"reports/{previous_files[0]}", 'r') as f:
                        previous_data = json.load(f)
                        logger.info(f"Loaded previous verification from {previous_files[0]}")
                        
                        # Ensure proper structure
                        if all(key in previous_data for key in ["properties_completion", "jsonb_properties_completion", "inchikey_completion", "pubchem_completion"]):
                            previous_verification = previous_data
                        else:
                            logger.warning(f"Previous verification file has incompatible format")
            except Exception as e:
                logger.warning(f"Could not load previous verification: {str(e)}")
            
            # Prepare results
            results = {
                "timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
                "total_chembl_molecules": total_chembl_molecules,
                "properties_completion": {
                    "molecules_with_incomplete_properties": incomplete_count,
                    "completion_percentage": property_completion,
                    "change_since_last_verification": None if not previous_verification else property_completion - previous_verification["properties_completion"]["completion_percentage"]
                },
                "jsonb_properties_completion": {
                    "molecules_missing_jsonb_properties": missing_jsonb_count,
                    "completion_percentage": jsonb_completion,
                    "change_since_last_verification": None if not previous_verification else jsonb_completion - previous_verification["jsonb_properties_completion"]["completion_percentage"]
                },
                "inchikey_completion": {
                    "molecules_missing_inchikey": missing_inchikey_count,
                    "completion_percentage": inchikey_completion,
                    "change_since_last_verification": None if not previous_verification else inchikey_completion - previous_verification["inchikey_completion"]["completion_percentage"]
                },
                "pubchem_completion": {
                    "molecules_missing_pubchem_cid": missing_pubchem_count,
                    "completion_percentage": pubchem_completion,
                    "change_since_last_verification": None if not previous_verification else pubchem_completion - previous_verification["pubchem_completion"]["completion_percentage"]
                }
            }
            
            return results
    except Exception as e:
        logger.error(f"Error verifying ChEMBL molecule properties: {str(e)}")
        raise

def main():
    """Main function."""
    try:
        # Connect to database
        conn = connect_to_database()
        
        try:
            # Verify ChEMBL molecule properties
            results = verify_chembl_molecule_properties(conn)
            
            # Save results
            timestamp = time.strftime('%Y%m%d_%H%M%S')
            results_file = f"reports/chembl_verification_{timestamp}.json"
            
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            logger.info(f"Results saved to {results_file}")
            
            # Print summary
            print("\n=== ChEMBL Data Quality Verification ===")
            print(f"Total ChEMBL molecules: {results['total_chembl_molecules']}")
            print(f"Properties completion: {results['properties_completion']['completion_percentage']:.2f}% ({results['total_chembl_molecules'] - results['properties_completion']['molecules_with_incomplete_properties']}/{results['total_chembl_molecules']} complete)")
            print(f"JSONB properties completion: {results['jsonb_properties_completion']['completion_percentage']:.2f}% ({results['total_chembl_molecules'] - results['jsonb_properties_completion']['molecules_missing_jsonb_properties']}/{results['total_chembl_molecules']} complete)")
            print(f"InChIKey completion: {results['inchikey_completion']['completion_percentage']:.2f}% ({results['total_chembl_molecules'] - results['inchikey_completion']['molecules_missing_inchikey']}/{results['total_chembl_molecules']} complete)")
            print(f"PubChem CID completion: {results['pubchem_completion']['completion_percentage']:.2f}% ({results['total_chembl_molecules'] - results['pubchem_completion']['molecules_missing_pubchem_cid']}/{results['total_chembl_molecules']} complete)")
            
            if results['properties_completion']['change_since_last_verification'] is not None:
                print("\nChanges since last verification:")
                print(f"Properties completion: {results['properties_completion']['change_since_last_verification']:.2f}%")
                print(f"JSONB properties completion: {results['jsonb_properties_completion']['change_since_last_verification']:.2f}%")
                print(f"InChIKey completion: {results['inchikey_completion']['change_since_last_verification']:.2f}%")
                print(f"PubChem CID completion: {results['pubchem_completion']['change_since_last_verification']:.2f}%")
            
            print("==========================================\n")
            
            return 0
        
        except Exception as e:
            logger.error(f"Error in main function: {str(e)}")
            return 1
        
        finally:
            conn.close()
            logger.info("Database connection closed")
    
    except Exception as e:
        logger.error(f"Unhandled exception: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
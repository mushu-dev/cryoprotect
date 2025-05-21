#!/usr/bin/env python3
"""
Identify molecules with data quality issues in the database.

This script identifies:
- Molecules with incomplete properties
- Molecules with missing JSONB properties
- Molecules without PubChem cross-references
- Molecules without InChIKeys
"""

import os
import sys
import json
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from typing import Dict, List, Any, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/identify_issues.log')
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

# Required property types
REQUIRED_PROPERTIES = [
    'LogP',
    'TPSA',
    'Molecular Weight',
    'Heavy Atom Count',
    'Hydrogen Bond Donor Count',
    'Hydrogen Bond Acceptor Count',
    'Rotatable Bond Count',
    'Ring Count',
    'Aromatic Ring Count'
]

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

def get_property_type_ids(conn) -> Dict[str, str]:
    """
    Get mapping of property type names to IDs.
    
    Args:
        conn: Database connection
        
    Returns:
        Dict mapping property type names to IDs
    """
    logger.info("Getting property type IDs...")
    property_ids = {}
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT id, name FROM property_types
            WHERE name = ANY(%s)
            """, (REQUIRED_PROPERTIES,))
            
            for row in cursor:
                property_ids[row['name']] = row['id']
            
            logger.info(f"Found {len(property_ids)} property types")
            return property_ids
    except Exception as e:
        logger.error(f"Error getting property type IDs: {str(e)}")
        conn.rollback()
        raise

def find_molecules_with_incomplete_properties(conn, property_ids: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Find molecules with incomplete properties.
    
    Args:
        conn: Database connection
        property_ids: Dict mapping property types to IDs
        
    Returns:
        List of molecules with incomplete properties
    """
    logger.info("Finding molecules with incomplete properties...")
    
    try:
        with conn.cursor() as cursor:
            # Get all ChEMBL molecules
            cursor.execute("""
            SELECT 
                id, 
                name, 
                chembl_id, 
                pubchem_cid, 
                inchikey,
                smiles
            FROM 
                molecules
            WHERE 
                chembl_id IS NOT NULL
            """)
            
            all_molecules = cursor.fetchall()
            logger.info(f"Found {len(all_molecules)} total ChEMBL molecules")
            
            # For each molecule, get its properties
            molecules_with_incomplete_properties = []
            
            for molecule in all_molecules:
                cursor.execute("""
                SELECT 
                    pt.name AS property_name
                FROM 
                    molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    mp.molecule_id = %s
                    AND pt.name = ANY(%s)
                """, (molecule['id'], REQUIRED_PROPERTIES))
                
                properties = cursor.fetchall()
                property_names = {prop['property_name'] for prop in properties}
                
                # Check if all required properties exist
                missing_properties = [prop for prop in REQUIRED_PROPERTIES if prop not in property_names]
                
                if missing_properties:
                    molecule_with_issues = dict(molecule)
                    molecule_with_issues['missing_properties'] = missing_properties
                    molecule_with_issues['existing_properties'] = list(property_names)
                    molecules_with_incomplete_properties.append(molecule_with_issues)
            
            logger.info(f"Found {len(molecules_with_incomplete_properties)} molecules with incomplete properties")
            return molecules_with_incomplete_properties
    except Exception as e:
        logger.error(f"Error finding molecules with incomplete properties: {str(e)}")
        conn.rollback()
        raise

def find_molecules_missing_jsonb_properties(conn) -> List[Dict[str, Any]]:
    """
    Find molecules with missing JSONB properties.
    
    Args:
        conn: Database connection
        
    Returns:
        List of molecules missing JSONB properties
    """
    logger.info("Finding molecules missing JSONB properties...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT 
                id, 
                name, 
                chembl_id, 
                pubchem_cid, 
                inchikey,
                smiles,
                properties
            FROM 
                molecules
            WHERE 
                chembl_id IS NOT NULL
                AND (properties IS NULL OR properties = '{}' OR properties = '[]')
            ORDER BY name
            """)
            
            results = cursor.fetchall()
            logger.info(f"Found {len(results)} molecules missing JSONB properties")
            return results
    except Exception as e:
        logger.error(f"Error finding molecules missing JSONB properties: {str(e)}")
        conn.rollback()
        raise

def find_molecules_missing_pubchem_cid(conn) -> List[Dict[str, Any]]:
    """
    Find molecules missing PubChem CIDs.
    
    Args:
        conn: Database connection
        
    Returns:
        List of molecules missing PubChem CIDs
    """
    logger.info("Finding molecules missing PubChem CIDs...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT 
                id, 
                name, 
                chembl_id, 
                inchikey,
                smiles
            FROM 
                molecules
            WHERE 
                chembl_id IS NOT NULL
                AND pubchem_cid IS NULL
            ORDER BY name
            """)
            
            results = cursor.fetchall()
            logger.info(f"Found {len(results)} molecules missing PubChem CIDs")
            return results
    except Exception as e:
        logger.error(f"Error finding molecules missing PubChem CIDs: {str(e)}")
        conn.rollback()
        raise

def find_molecules_missing_inchikey(conn) -> List[Dict[str, Any]]:
    """
    Find molecules missing InChIKeys.
    
    Args:
        conn: Database connection
        
    Returns:
        List of molecules missing InChIKeys
    """
    logger.info("Finding molecules missing InChIKeys...")
    
    try:
        with conn.cursor() as cursor:
            cursor.execute("""
            SELECT 
                id, 
                name, 
                chembl_id, 
                pubchem_cid,
                smiles
            FROM 
                molecules
            WHERE 
                chembl_id IS NOT NULL
                AND (inchikey IS NULL OR inchikey = '')
            ORDER BY name
            """)
            
            results = cursor.fetchall()
            logger.info(f"Found {len(results)} molecules missing InChIKeys")
            return results
    except Exception as e:
        logger.error(f"Error finding molecules missing InChIKeys: {str(e)}")
        conn.rollback()
        raise

def save_results(results: Dict[str, Any], output_file: str) -> None:
    """
    Save results to a JSON file.
    
    Args:
        results: Results dict
        output_file: Output file path
    """
    logger.info(f"Saving results to {output_file}...")
    
    try:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        logger.info(f"Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Error saving results: {str(e)}")
        raise

def main():
    """Main function."""
    # Connect to database
    conn = connect_to_database()
    
    try:
        # Get property type IDs
        property_ids = get_property_type_ids(conn)
        
        # Find molecules with issues
        incomplete_properties = find_molecules_with_incomplete_properties(conn, property_ids)
        missing_jsonb = find_molecules_missing_jsonb_properties(conn)
        missing_pubchem = find_molecules_missing_pubchem_cid(conn)
        missing_inchikey = find_molecules_missing_inchikey(conn)
        
        # Prepare results
        results = {
            "molecules_with_incomplete_properties": incomplete_properties,
            "molecules_missing_jsonb_properties": missing_jsonb,
            "molecules_missing_pubchem_cid": missing_pubchem,
            "molecules_missing_inchikey": missing_inchikey,
            "summary": {
                "total_with_incomplete_properties": len(incomplete_properties),
                "total_missing_jsonb_properties": len(missing_jsonb),
                "total_missing_pubchem_cid": len(missing_pubchem),
                "total_missing_inchikey": len(missing_inchikey)
            }
        }
        
        # Save results
        save_results(results, "reports/molecule_issues.json")
        
        # Print summary
        print("\n=== Molecule Data Issues Summary ===")
        print(f"Molecules with incomplete properties: {len(incomplete_properties)}")
        print(f"Molecules missing JSONB properties: {len(missing_jsonb)}")
        print(f"Molecules missing PubChem CIDs: {len(missing_pubchem)}")
        print(f"Molecules missing InChIKeys: {len(missing_inchikey)}")
        print(f"Detailed report saved to: reports/molecule_issues.json")
        print("====================================\n")
        
        # Return success
        return 0
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return 1
    finally:
        # Close database connection
        conn.close()
        logger.info("Database connection closed")

if __name__ == "__main__":
    sys.exit(main())
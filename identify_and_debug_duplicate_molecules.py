#!/usr/bin/env python3
"""
Script to identify and debug duplicate molecule consolidation issues.
"""

import os
import sys
import logging
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv
from datetime import datetime
import json

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

def get_db_connection():
    """Get a direct database connection using psycopg2."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        logger.info("Connected to database using direct PostgreSQL connection")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def identify_duplicate_molecules():
    """Identify duplicate molecules by InChIKey."""
    logger.info("Identifying duplicate molecules by InChIKey...")
    
    conn = get_db_connection()
    if not conn:
        return []
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    duplicate_groups = []
    
    try:
        # Find all InChIKeys with multiple molecule entries
        cursor.execute("""
            SELECT inchikey, COUNT(*) as count, array_agg(id) as molecule_ids, 
                   array_agg(name) as names, array_agg(pubchem_cid) as pubchem_cids
            FROM molecules
            WHERE inchikey IS NOT NULL
            GROUP BY inchikey
            HAVING COUNT(*) > 1
        """)
        
        duplicate_groups = cursor.fetchall()
        logger.info(f"Found {len(duplicate_groups)} InChIKeys with duplicate molecules")
        
    except Exception as e:
        logger.error(f"Error identifying duplicate molecules: {e}")
        conn.rollback()
    
    cursor.close()
    conn.close()
    return duplicate_groups

def debug_select_primary_molecule(group):
    """Debug version of function to select the primary molecule from a group of duplicates."""
    # Parse PostgreSQL array string into a Python list of UUIDs
    molecule_ids_str = group['molecule_ids']
    # Remove curly braces and split by comma
    if molecule_ids_str.startswith('{') and molecule_ids_str.endswith('}'): 
        molecule_ids_str = molecule_ids_str[1:-1]
    molecule_ids = [id.strip() for id in molecule_ids_str.split(',')]
    
    logger.info(f"Processing molecules with IDs: {molecule_ids}")
    
    conn = get_db_connection()
    if not conn:
        return None
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    primary_molecule = None
    
    try:
        # Check if any of these molecules actually exist in the database
        for molecule_id in molecule_ids:
            cursor.execute("""
                SELECT id, name, pubchem_cid, created_at, data_source, 
                       molecular_weight, smiles, inchi, inchikey
                FROM molecules
                WHERE id = %s
            """, (molecule_id,))
            
            result = cursor.fetchone()
            if result:
                logger.info(f"Found molecule: {result['name']} (ID: {result['id']}, PubChem CID: {result['pubchem_cid']})")
            else:
                logger.warning(f"Molecule with ID {molecule_id} not found in database")
                
        # Now try to get all molecules with the group's InChIKey directly
        inchikey = group['inchikey']
        logger.info(f"Trying to get molecules with InChIKey: {inchikey}")
        logger.info(f"InChIKey type: {type(inchikey)}, value: {inchikey}")
        
        # Use direct string substitution for safety (since we're debugging)
        query = f"""
            SELECT id, name, pubchem_cid, created_at, data_source, 
                   molecular_weight, smiles, inchi, inchikey
            FROM molecules
            WHERE inchikey = '{inchikey}'
            ORDER BY 
                CASE WHEN pubchem_cid IS NOT NULL THEN 0 ELSE 1 END,
                CASE WHEN name LIKE 'TEST_%' THEN 1 ELSE 0 END,
                created_at
        """
        logger.info(f"Query: {query}")
        cursor.execute(query)
        
        molecules = cursor.fetchall()
        logger.info(f"Found {len(molecules)} molecules with InChIKey {inchikey}")
        
        if molecules and len(molecules) > 0:
            # Debug output
            for i, molecule in enumerate(molecules):
                logger.info(f"  Molecule {i+1}: {molecule['name']} (ID: {molecule['id']}, PubChem CID: {molecule['pubchem_cid']})")
            
            # First molecule in sorted order is our primary
            primary_molecule = molecules[0]
            logger.info(f"Would select as primary: {primary_molecule['name']} (ID: {primary_molecule['id']}, PubChem CID: {primary_molecule['pubchem_cid']})")
            
            return primary_molecule
        else:
            logger.error(f"No molecules found with InChIKey {inchikey}")
            return None
        
    except Exception as e:
        logger.error(f"Error in debug_select_primary_molecule: {e}")
        conn.rollback()
        import traceback
        logger.error(traceback.format_exc())
        return None
    finally:
        cursor.close()
        conn.close()

def debug_first_duplicate_group():
    """Debug the first duplicate group."""
    duplicate_groups = identify_duplicate_molecules()
    
    if not duplicate_groups:
        logger.info("No duplicate molecules found")
        return
        
    # Debug just the first group
    first_group = duplicate_groups[0]
    logger.info(f"Debugging first group with InChIKey: {first_group['inchikey']} ({first_group['count']} molecules)")
    
    # Try to select the primary molecule
    primary_molecule = debug_select_primary_molecule(first_group)
    
    if primary_molecule:
        logger.info(f"Successfully selected primary molecule: {primary_molecule['name']}")
    else:
        logger.error("Failed to select primary molecule")

if __name__ == "__main__":
    debug_first_duplicate_group()
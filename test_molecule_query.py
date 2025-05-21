#!/usr/bin/env python3
"""
Test script to check if we can query molecules by UUID.
"""

import os
import sys
import logging
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv
import json
from datetime import datetime
from decimal import Decimal

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

def test_molecule_query():
    """Test querying molecules by UUID."""
    conn = get_db_connection()
    if not conn:
        logger.error("Failed to connect to database")
        return
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    
    try:
        # Get a few molecule UUIDs for testing
        cursor.execute("""
            SELECT id, name, inchikey, pubchem_cid
            FROM molecules
            LIMIT 3
        """)
        
        test_molecules = cursor.fetchall()
        logger.info(f"Found {len(test_molecules)} test molecules")
        
        if not test_molecules:
            logger.error("No molecules found in database")
            return
            
        # Test querying a single molecule by UUID
        test_uuid = test_molecules[0]['id']
        logger.info(f"Testing query for molecule with UUID: {test_uuid}")
        
        cursor.execute("""
            SELECT id, name, pubchem_cid, created_at, data_source, 
                   molecular_weight, smiles, inchi, inchikey
            FROM molecules
            WHERE id = %s
        """, (test_uuid,))
        
        result = cursor.fetchall()
        logger.info(f"Query returned {len(result)} rows")
        
        if result:
            # Instead of trying to JSON serialize, just print key info directly
            first_result = result[0]
            logger.info(f"First result: ID={first_result['id']}, Name={first_result['name']}, PubChem CID={first_result['pubchem_cid']}")
            logger.info(f"InChIKey={first_result['inchikey']}, Molecular Weight={first_result['molecular_weight']}")
            logger.info(f"Created at={first_result['created_at']}")
        else:
            logger.error("Query returned no results")
            
        # Test querying multiple molecules by UUID
        # Force UUIDs to be strings for psycopg2
        test_uuids = [str(m['id']) for m in test_molecules]
        placeholders = ','.join(['%s'] * len(test_uuids))
        logger.info(f"Testing query for multiple molecules with UUIDs: {test_uuids}")
        
        query = f"""
            SELECT id, name, pubchem_cid, created_at, data_source, 
                   molecular_weight, smiles, inchi, inchikey
            FROM molecules
            WHERE id IN ({placeholders})
        """
        
        cursor.execute(query, test_uuids)
        multi_result = cursor.fetchall()
        logger.info(f"Multi-UUID query returned {len(multi_result)} rows")
        
        if multi_result:
            for i, row in enumerate(multi_result):
                logger.info(f"Result {i+1}: {row['name']} (ID: {row['id']}, PubChem CID: {row['pubchem_cid']})")
                
            # Get duplicate molecules with the InChIKey we were looking at earlier
            cursor.execute("""
                SELECT id, name, pubchem_cid, created_at, data_source, 
                       molecular_weight, smiles, inchi, inchikey
                FROM molecules
                WHERE inchikey = 'CZMRCDWAGMRECN-UGDNZRGBSA-N'
            """)
            dup_result = cursor.fetchall()
            logger.info(f"Found {len(dup_result)} molecules with test InChIKey")
            
            if dup_result:
                for i, row in enumerate(dup_result):
                    logger.info(f"Duplicate {i+1}: {row['name']} (ID: {row['id']}, PubChem CID: {row['pubchem_cid']})")
        else:
            logger.error("Multi-UUID query returned no results")
            
    except Exception as e:
        logger.error(f"Error testing molecule query: {e}")
        conn.rollback()
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    test_molecule_query()
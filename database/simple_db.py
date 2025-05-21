"""
Simplified database module for direct database operations.
This module provides simple functions for database connection and operations,
designed to work in container environments without external dependencies.
"""

import os
import psycopg2
from psycopg2.extras import RealDictCursor, execute_batch
from typing import List, Dict, Any, Optional, Tuple, Union
import logging
import time
import random

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/simple_db.log')
    ]
)

logger = logging.getLogger(__name__)

# Database connection parameters
DB_HOST = os.environ.get('SUPABASE_DB_HOST')
DB_PORT = os.environ.get('SUPABASE_DB_PORT', '5432')
DB_NAME = os.environ.get('SUPABASE_DB_NAME', 'postgres')
DB_USER = os.environ.get('SUPABASE_DB_USER')
DB_PASSWORD = os.environ.get('SUPABASE_DB_PASSWORD')

# Connection pools are not used in this simplified version
# Each function will create and close its own connection

def get_db_connection():
    """Create a new database connection"""
    try:
        connection = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        connection.autocommit = False
        return connection
    except Exception as e:
        logger.error(f"Error connecting to database: {e}")
        raise

def execute_query(query: str, params: Optional[tuple] = None, fetch: bool = True):
    """Execute a database query with optional parameters"""
    connection = None
    try:
        connection = get_db_connection()
        cursor = connection.cursor(cursor_factory=RealDictCursor)
        cursor.execute(query, params)
        
        if fetch:
            result = cursor.fetchall()
        else:
            result = None
            connection.commit()
            
        cursor.close()
        return result
    except Exception as e:
        if connection:
            connection.rollback()
        logger.error(f"Error executing query: {e}")
        raise
    finally:
        if connection:
            connection.close()

def execute_with_retry(query: str, params: Optional[tuple] = None, 
                      fetch: bool = True, max_retries: int = 3):
    """Execute a query with retry logic for transient errors"""
    retries = 0
    while retries < max_retries:
        try:
            return execute_query(query, params, fetch)
        except (psycopg2.OperationalError, psycopg2.InterfaceError) as e:
            retries += 1
            if retries >= max_retries:
                logger.error(f"Max retries reached. Query failed: {e}")
                raise
            
            # Exponential backoff with jitter
            sleep_time = (2 ** retries) + random.uniform(0, 1)
            logger.warning(f"Database error, retrying in {sleep_time:.2f} seconds: {e}")
            time.sleep(sleep_time)
        except Exception as e:
            # Don't retry other errors
            logger.error(f"Query failed with non-retriable error: {e}")
            raise

def insert_molecule(molecule: Dict[str, Any]) -> int:
    """Insert a molecule into the database"""
    connection = None
    try:
        connection = get_db_connection()
        cursor = connection.cursor(cursor_factory=RealDictCursor)
        
        # Check if molecule already exists by CID
        cursor.execute(
            "SELECT id FROM molecules WHERE pubchem_cid = %s",
            (molecule.get('pubchem_cid'),)
        )
        existing = cursor.fetchone()
        
        if existing:
            # Update existing molecule
            update_query = """
            UPDATE molecules SET
                name = %s,
                formula = %s,
                smiles = %s,
                inchi = %s,
                inchi_key = %s,
                cryoprotectant = %s
            WHERE id = %s
            RETURNING id
            """
            cursor.execute(update_query, (
                molecule.get('name'),
                molecule.get('formula'),
                molecule.get('smiles'),
                molecule.get('inchi'),
                molecule.get('inchi_key'),
                molecule.get('cryoprotectant', True),
                existing['id']
            ))
            result = cursor.fetchone()
            molecule_id = result['id']
        else:
            # Insert new molecule
            insert_query = """
            INSERT INTO molecules (
                pubchem_cid, chembl_id, name, formula, smiles, inchi, inchi_key, cryoprotectant
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id
            """
            cursor.execute(insert_query, (
                molecule.get('pubchem_cid'),
                molecule.get('chembl_id'),
                molecule.get('name'),
                molecule.get('formula'),
                molecule.get('smiles'),
                molecule.get('inchi'),
                molecule.get('inchi_key'),
                molecule.get('cryoprotectant', True)
            ))
            result = cursor.fetchone()
            molecule_id = result['id']
        
        connection.commit()
        return molecule_id
    except Exception as e:
        if connection:
            connection.rollback()
        logger.error(f"Error inserting molecule {molecule.get('pubchem_cid')}: {e}")
        raise
    finally:
        if connection:
            connection.close()

def batch_insert_molecules(molecules: List[Dict[str, Any]]) -> Dict[str, int]:
    """Insert multiple molecules in a batch transaction"""
    connection = None
    try:
        connection = get_db_connection()
        cursor = connection.cursor(cursor_factory=RealDictCursor)
        
        results = {}
        for molecule in molecules:
            # Check if molecule already exists by CID
            cursor.execute(
                "SELECT id FROM molecules WHERE pubchem_cid = %s",
                (molecule.get('pubchem_cid'),)
            )
            existing = cursor.fetchone()
            
            if existing:
                # Update existing molecule
                update_query = """
                UPDATE molecules SET
                    name = %s,
                    formula = %s,
                    smiles = %s,
                    inchi = %s,
                    inchi_key = %s,
                    cryoprotectant = %s
                WHERE id = %s
                RETURNING id
                """
                cursor.execute(update_query, (
                    molecule.get('name'),
                    molecule.get('formula'),
                    molecule.get('smiles'),
                    molecule.get('inchi'),
                    molecule.get('inchi_key'),
                    molecule.get('cryoprotectant', True),
                    existing['id']
                ))
                result = cursor.fetchone()
                results[molecule.get('pubchem_cid')] = result['id']
            else:
                # Insert new molecule
                insert_query = """
                INSERT INTO molecules (
                    pubchem_cid, chembl_id, name, formula, smiles, inchi, inchi_key, cryoprotectant
                ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                RETURNING id
                """
                cursor.execute(insert_query, (
                    molecule.get('pubchem_cid'),
                    molecule.get('chembl_id'),
                    molecule.get('name'),
                    molecule.get('formula'),
                    molecule.get('smiles'),
                    molecule.get('inchi'),
                    molecule.get('inchi_key'),
                    molecule.get('cryoprotectant', True)
                ))
                result = cursor.fetchone()
                results[molecule.get('pubchem_cid')] = result['id']
        
        connection.commit()
        return results
    except Exception as e:
        if connection:
            connection.rollback()
        logger.error(f"Error batch inserting molecules: {e}")
        raise
    finally:
        if connection:
            connection.close()

def insert_property(molecule_id: int, property_data: Dict[str, Any]) -> int:
    """Insert molecular property data"""
    connection = None
    try:
        connection = get_db_connection()
        cursor = connection.cursor(cursor_factory=RealDictCursor)
        
        # Check if property already exists
        cursor.execute(
            "SELECT id FROM molecular_properties WHERE molecule_id = %s AND property_name = %s",
            (molecule_id, property_data.get('property_name'))
        )
        existing = cursor.fetchone()
        
        if existing:
            # Update existing property
            update_query = """
            UPDATE molecular_properties SET
                property_value = %s,
                units = %s,
                source = %s
            WHERE id = %s
            RETURNING id
            """
            cursor.execute(update_query, (
                property_data.get('property_value'),
                property_data.get('units'),
                property_data.get('source', 'PubChem'),
                existing['id']
            ))
            result = cursor.fetchone()
            property_id = result['id']
        else:
            # Insert new property
            insert_query = """
            INSERT INTO molecular_properties (
                molecule_id, property_name, property_value, units, source
            ) VALUES (%s, %s, %s, %s, %s)
            RETURNING id
            """
            cursor.execute(insert_query, (
                molecule_id,
                property_data.get('property_name'),
                property_data.get('property_value'),
                property_data.get('units'),
                property_data.get('source', 'PubChem')
            ))
            result = cursor.fetchone()
            property_id = result['id']
        
        connection.commit()
        return property_id
    except Exception as e:
        if connection:
            connection.rollback()
        logger.error(f"Error inserting property for molecule {molecule_id}: {e}")
        raise
    finally:
        if connection:
            connection.close()

def batch_insert_properties(properties_data: List[Dict[str, Any]]) -> List[int]:
    """Insert multiple properties in a batch transaction"""
    connection = None
    try:
        connection = get_db_connection()
        cursor = connection.cursor(cursor_factory=RealDictCursor)
        
        property_ids = []
        for prop_data in properties_data:
            molecule_id = prop_data.get('molecule_id')
            property_name = prop_data.get('property_name')
            
            # Check if property already exists
            cursor.execute(
                "SELECT id FROM molecular_properties WHERE molecule_id = %s AND property_name = %s",
                (molecule_id, property_name)
            )
            existing = cursor.fetchone()
            
            if existing:
                # Update existing property
                update_query = """
                UPDATE molecular_properties SET
                    property_value = %s,
                    units = %s,
                    source = %s
                WHERE id = %s
                RETURNING id
                """
                cursor.execute(update_query, (
                    prop_data.get('property_value'),
                    prop_data.get('units'),
                    prop_data.get('source', 'PubChem'),
                    existing['id']
                ))
                result = cursor.fetchone()
                property_ids.append(result['id'])
            else:
                # Insert new property
                insert_query = """
                INSERT INTO molecular_properties (
                    molecule_id, property_name, property_value, units, source
                ) VALUES (%s, %s, %s, %s, %s)
                RETURNING id
                """
                cursor.execute(insert_query, (
                    molecule_id,
                    property_name,
                    prop_data.get('property_value'),
                    prop_data.get('units'),
                    prop_data.get('source', 'PubChem')
                ))
                result = cursor.fetchone()
                property_ids.append(result['id'])
        
        connection.commit()
        return property_ids
    except Exception as e:
        if connection:
            connection.rollback()
        logger.error(f"Error batch inserting properties: {e}")
        raise
    finally:
        if connection:
            connection.close()

def get_molecule_count() -> int:
    """Get the total count of molecules in the database"""
    try:
        result = execute_query("SELECT COUNT(*) as count FROM molecules")
        return result[0]['count'] if result else 0
    except Exception as e:
        logger.error(f"Error getting molecule count: {e}")
        return 0

def get_property_count() -> int:
    """Get the total count of molecular properties in the database"""
    try:
        result = execute_query("SELECT COUNT(*) as count FROM molecular_properties")
        return result[0]['count'] if result else 0
    except Exception as e:
        logger.error(f"Error getting property count: {e}")
        return 0
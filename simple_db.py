#!/usr/bin/env python3
"""
Simplified database module for the container environment.
This module provides a simplified interface for connecting to the database
and performing basic operations without external dependencies.
"""

import os
import psycopg2
from psycopg2.extras import RealDictCursor
import json
import logging
import time
from typing import Dict, List, Any, Optional, Tuple

logger = logging.getLogger(__name__)

def get_connection():
    """Get a database connection using environment variables."""
    # Get connection parameters from environment
    db_host = os.getenv('SUPABASE_DB_HOST', 'localhost')
    db_port = os.getenv('SUPABASE_DB_PORT', '5432')
    db_name = os.getenv('SUPABASE_DB_NAME', 'postgres')
    db_user = os.getenv('SUPABASE_DB_USER', 'postgres')
    db_password = os.getenv('SUPABASE_DB_PASSWORD', '')
    
    # Log connection details (without password)
    logger.info(f"Connecting to database: {db_host}:{db_port}/{db_name} as {db_user}")
    
    try:
        # Connect to database
        conn = psycopg2.connect(
            host=db_host,
            port=db_port,
            dbname=db_name,
            user=db_user,
            password=db_password
        )
        return conn
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        raise

def execute_query(query: str, params: Optional[tuple] = None, fetchall: bool = True) -> List[Dict[str, Any]]:
    """Execute a query and return the results."""
    conn = None
    try:
        conn = get_connection()
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute(query, params)
            if fetchall:
                results = cursor.fetchall()
                return [dict(row) for row in results]
            else:
                return []
    except Exception as e:
        logger.error(f"Error executing query: {str(e)}")
        raise
    finally:
        if conn:
            conn.close()

def execute_transaction(statements: List[Tuple[str, Any]]) -> bool:
    """Execute multiple statements in a transaction."""
    conn = None
    try:
        conn = get_connection()
        with conn.cursor() as cursor:
            for query, params in statements:
                cursor.execute(query, params)
        conn.commit()
        return True
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error executing transaction: {str(e)}")
        raise
    finally:
        if conn:
            conn.close()

def insert_molecule(data: Dict[str, Any]) -> str:
    """Insert a molecule into the database and return its ID."""
    query = """
    INSERT INTO molecules 
    (name, smiles, inchi, inchikey, formula, molecular_weight, properties, data_source, created_by)
    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, '00000000-0000-0000-0000-000000000001')
    ON CONFLICT (inchikey) 
    DO UPDATE SET
        properties = COALESCE(molecules.properties, '{}') || %s
    RETURNING id;
    """
    
    # Convert properties to JSON string
    properties_json = json.dumps(data.get('properties', {}))
    
    params = (
        data.get('name', ''),
        data.get('smiles', ''),
        data.get('inchi', ''),
        data.get('inchikey', ''),
        data.get('formula', ''),
        data.get('molecular_weight', 0),
        properties_json,
        data.get('data_source', 'PubChem Import'),
        properties_json
    )
    
    conn = None
    try:
        conn = get_connection()
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            result = cursor.fetchone()
            conn.commit()
            return result[0]
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error inserting molecule: {str(e)}")
        raise
    finally:
        if conn:
            conn.close()

def batch_insert_molecules(molecules: List[Dict[str, Any]]) -> Dict[str, int]:
    """Insert multiple molecules in a batch transaction."""
    if not molecules:
        return {"inserted": 0, "errors": 0}
    
    conn = None
    inserted = 0
    errors = 0
    
    try:
        conn = get_connection()
        with conn.cursor() as cursor:
            for molecule in molecules:
                try:
                    # Convert properties to JSON string
                    properties_json = json.dumps(molecule.get('properties', {}))
                    
                    query = """
                    INSERT INTO molecules 
                    (name, smiles, inchi, inchikey, formula, molecular_weight, properties, data_source, created_by)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, '00000000-0000-0000-0000-000000000001')
                    ON CONFLICT (inchikey) 
                    DO UPDATE SET
                        properties = COALESCE(molecules.properties, '{}') || %s
                    RETURNING id;
                    """
                    
                    params = (
                        molecule.get('name', ''),
                        molecule.get('smiles', ''),
                        molecule.get('inchi', ''),
                        molecule.get('inchikey', ''),
                        molecule.get('formula', ''),
                        molecule.get('molecular_weight', 0),
                        properties_json,
                        molecule.get('data_source', 'PubChem Import'),
                        properties_json
                    )
                    
                    cursor.execute(query, params)
                    inserted += 1
                except Exception as e:
                    logger.error(f"Error inserting molecule {molecule.get('inchikey', 'unknown')}: {str(e)}")
                    errors += 1
            
            conn.commit()
            return {"inserted": inserted, "errors": errors}
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error in batch insert: {str(e)}")
        return {"inserted": inserted, "errors": len(molecules) - inserted}
    finally:
        if conn:
            conn.close()

def check_molecule_exists(inchikey: str) -> Optional[str]:
    """Check if a molecule exists in the database by InChIKey."""
    query = "SELECT id FROM molecules WHERE inchikey = %s LIMIT 1;"
    
    conn = None
    try:
        conn = get_connection()
        with conn.cursor() as cursor:
            cursor.execute(query, (inchikey,))
            result = cursor.fetchone()
            if result:
                return result[0]
            return None
    except Exception as e:
        logger.error(f"Error checking if molecule exists: {str(e)}")
        return None
    finally:
        if conn:
            conn.close()

def get_database_stats() -> Dict[str, Any]:
    """Get database statistics."""
    stats = {}
    
    queries = [
        ("SELECT COUNT(*) FROM molecules;", "total_molecules"),
        ("SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%PubChem%';", "pubchem_molecules"),
        ("SELECT COUNT(*) FROM molecules WHERE inchikey IS NOT NULL;", "molecules_with_inchikey"),
        ("SELECT COUNT(*) FROM molecules WHERE smiles IS NOT NULL;", "molecules_with_smiles"),
        ("SELECT COUNT(*) FROM molecules WHERE properties IS NOT NULL AND properties::text != '{}';", "molecules_with_properties")
    ]
    
    conn = None
    try:
        conn = get_connection()
        with conn.cursor() as cursor:
            for query, key in queries:
                try:
                    cursor.execute(query)
                    result = cursor.fetchone()
                    stats[key] = result[0] if result else 0
                except Exception as e:
                    logger.error(f"Error executing query '{query}': {str(e)}")
                    stats[key] = -1
        return stats
    except Exception as e:
        logger.error(f"Error getting database stats: {str(e)}")
        return {"error": str(e)}
    finally:
        if conn:
            conn.close()
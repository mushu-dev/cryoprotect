#!/usr/bin/env python3
"""
Simplified database module that ensures all records are public.

This module provides direct connection to PostgreSQL and ensures that all
inserted records have is_public=true to prevent RLS issues.
"""

import os
import sys
import logging
import uuid
import psycopg2
import psycopg2.extras
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from contextlib import contextmanager
from typing import Dict, List, Any, Optional, Tuple, Union
from datetime import datetime

# Import from the connection configuration module
from .connection_config import (
    validate_config,
    get_connection_config,
    test_adapter_configuration
)

# Register UUID type adapter
psycopg2.extras.register_uuid()

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Global connection pool
_pool = None

def init_connection_pool(config: Dict[str, Any] = None, min_conn: int = None, max_conn: int = None) -> bool:
    """
    Initialize the PostgreSQL connection pool.
    
    Args:
        config: Connection parameters dictionary with these keys:
               host, port, database, user, password
        min_conn: Minimum number of connections in the pool
        max_conn: Maximum number of connections in the pool
        
    Returns:
        True if the connection pool was successfully initialized, False otherwise
    """
    global _pool
    
    if _pool is not None:
        logger.info("Connection pool already initialized")
        return True
    
    # Validate configuration
    validate_config()
    
    # If config wasn't provided, try to get it from the configuration system
    if config is None:
        config = get_connection_config('supabase')
        
    # Set min/max connections from config if not provided
    if min_conn is None:
        min_conn = config.get('min_connections', 1)
    if max_conn is None:
        max_conn = config.get('max_connections', 10)
        
    # Check for required config values
    required_keys = ['host', 'user', 'password']
    missing_keys = [key for key in required_keys if key not in config or config[key] is None or config[key] == '']
    
    if missing_keys:
        logger.error(f"Missing required database configuration: {', '.join(missing_keys)}")
        return False
    
    try:
        # Prepare connection parameters
        conn_params = {
            'host': config['host'],
            'port': config.get('port', 5432),
            'dbname': config.get('database', 'postgres'),
            'user': config['user'],
            'password': config['password'],
            'cursor_factory': RealDictCursor
        }
        
        # Add application name if provided
        if 'application_name' in config:
            conn_params['application_name'] = config['application_name']
            
        # Initialize the connection pool
        _pool = ThreadedConnectionPool(
            min_conn,
            max_conn,
            **conn_params
        )
        logger.info(f"Database connection pool initialized ({min_conn}-{max_conn} connections)")
        return True
    except Exception as e:
        logger.error(f"Failed to initialize connection pool: {str(e)}")
        return False

def get_connection():
    """
    Get a connection from the pool.
    
    Returns:
        A database connection or None if the pool is not initialized
    """
    global _pool
    
    if _pool is None:
        logger.warning("Attempting to get connection before pool initialization")
        # Try to initialize the pool with configuration from the configuration system
        if not init_connection_pool():
            logger.error("Failed to initialize connection pool on demand")
            return None
        
    try:
        return _pool.getconn()
    except Exception as e:
        logger.error(f"Error getting connection from pool: {str(e)}")
        return None

def release_connection(conn):
    """
    Release a connection back to the pool.
    
    Args:
        conn: The connection to release
    """
    global _pool
    
    if _pool is None:
        logger.error("Attempting to release connection before pool initialization")
        return
        
    try:
        _pool.putconn(conn)
    except Exception as e:
        logger.error(f"Error releasing connection to pool: {str(e)}")

def execute_query(query: str, params: tuple = None) -> List[Dict[str, Any]]:
    """
    Execute a SQL query and return the results.
    
    Args:
        query: SQL query string
        params: Query parameters tuple
        
    Returns:
        List of dictionaries where each dictionary is a row of results
    """
    conn = get_connection()
    if not conn:
        return []
    
    try:
        with conn.cursor() as cursor:
            # Try to make records viewable by settings jwt claims for a consistent user
            try:
                cursor.execute("""
                    SET LOCAL "request.jwt.claims" = '{"sub": "77777777-7777-7777-7777-777777777777", "role": "authenticated"}';
                """)
            except Exception as e:
                logger.warning(f"Could not set jwt claims: {str(e)}")
            
            # Execute the actual query
            cursor.execute(query, params)
            
            if cursor.description:  # If the query returns rows
                result = cursor.fetchall()
                return result
            else:
                conn.commit()
                return []
    except Exception as e:
        conn.rollback()
        logger.error(f"Error executing query: {str(e)}")
        return []
    finally:
        release_connection(conn)

def execute_batch(query: str, param_sets: List[tuple]) -> bool:
    """
    Execute a batch of SQL commands with different parameters.
    
    Args:
        query: SQL query template
        param_sets: List of parameter tuples
        
    Returns:
        True if the batch was executed successfully, False otherwise
    """
    from psycopg2.extras import execute_batch
    
    if not param_sets:
        return True
        
    conn = get_connection()
    if not conn:
        return False
        
    try:
        with conn.cursor() as cursor:
            # Try to make records viewable by settings jwt claims for a consistent user
            try:
                cursor.execute("""
                    SET LOCAL "request.jwt.claims" = '{"sub": "77777777-7777-7777-7777-777777777777", "role": "authenticated"}';
                """)
            except Exception as e:
                logger.warning(f"Could not set jwt claims: {str(e)}")
                
            execute_batch(cursor, query, param_sets)
            conn.commit()
            return True
    except Exception as e:
        conn.rollback()
        logger.error(f"Error executing batch: {str(e)}")
        return False
    finally:
        release_connection(conn)

def transaction():
    """
    Create a transaction context manager.

    Usage:
        with transaction() as cursor:
            cursor.execute("INSERT INTO ...")
            cursor.execute("UPDATE ...")

    Returns:
        Transaction context manager
    """
    class TransactionContextManager:
        def __init__(self):
            self.conn = None
            self.cursor = None

        def __enter__(self):
            self.conn = get_connection()
            if not self.conn:
                raise Exception("Could not get database connection")

            self.cursor = self.conn.cursor(cursor_factory=RealDictCursor)

            # Set the JWT claims for auth.uid()
            try:
                self.cursor.execute("""
                    SET LOCAL "request.jwt.claims" = '{"sub": "77777777-7777-7777-7777-777777777777", "role": "authenticated"}';
                """)
            except Exception as e:
                logger.warning(f"Could not set jwt claims: {str(e)}")

            return self.cursor

        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is None:
                # No exception, commit the transaction
                self.conn.commit()
            else:
                # Exception occurred, rollback
                self.conn.rollback()
                logger.error(f"Transaction rolled back due to error: {str(exc_val)}")

            self.cursor.close()
            release_connection(self.conn)
            return False  # Don't suppress exceptions

    return TransactionContextManager()

def insert_molecule(molecule_data: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """
    Insert a molecule and ensure it's public.
    
    Args:
        molecule_data: Dictionary of molecule data
        
    Returns:
        Tuple of (success, molecule_id)
    """
    # Ensure is_public is set to True
    molecule_data = {**molecule_data, "is_public": True}
    
    # Create a consistent created_by UUID for auth access
    # This will allow our connection to see records it creates
    if "created_by" not in molecule_data:
        # Create a predictable UUID based on a known value
        # This will be the same across all connections
        molecule_data["created_by"] = uuid.UUID("77777777-7777-7777-7777-777777777777")
    
    # Construct SQL query
    fields = []
    values = []
    placeholders = []
    
    for key, value in molecule_data.items():
        if value is not None:
            fields.append(key)
            values.append(value)
            placeholders.append("%s")
    
    # Add timestamps if not provided
    for field in ["created_at", "updated_at"]:
        if field not in fields:
            fields.append(field)
            values.append(datetime.now().isoformat())
            placeholders.append("%s")
    
    query = f"""
        INSERT INTO molecules 
        ({", ".join(fields)}) 
        VALUES ({", ".join(placeholders)})
        RETURNING id
    """
    
    try:
        result = execute_query(query, tuple(values))
        if result and len(result) > 0:
            return True, result[0]["id"]
        else:
            return False, None
    except Exception as e:
        logger.error(f"Error inserting molecule: {str(e)}")
        return False, None

def insert_property_type(name: str, data_type: str = "numeric") -> Optional[str]:
    """
    Insert or get a property type.

    Args:
        name: Property type name
        data_type: Data type (numeric, text, boolean)

    Returns:
        Property type ID or None if failed
    """
    query = """
        INSERT INTO property_types
        (name, data_type, created_at, updated_at)
        VALUES (%s, %s, NOW(), NOW())
        ON CONFLICT (name) DO UPDATE
        SET updated_at = NOW()
        RETURNING id
    """

    try:
        result = execute_query(query, (name, data_type))
        if result and len(result) > 0:
            return result[0]["id"]
        else:
            return None
    except Exception as e:
        logger.error(f"Error inserting property type: {str(e)}")
        return None

def insert_molecular_property(molecule_id: str, property_type_id: str,
                              value: Any, source: str = "Import") -> Optional[str]:
    """
    Insert a molecular property.

    Args:
        molecule_id: Molecule ID
        property_type_id: Property type ID
        value: Property value
        source: Data source

    Returns:
        Property ID or None if failed
    """
    # Determine value type and set the appropriate field
    if isinstance(value, (int, float)):
        value_field = "numeric_value"
    elif isinstance(value, bool):
        value_field = "boolean_value"
    else:
        value_field = "text_value"

    query = f"""
        INSERT INTO molecular_properties
        (molecule_id, property_type_id, {value_field}, source, created_at, updated_at)
        VALUES (%s, %s, %s, %s, NOW(), NOW())
        ON CONFLICT (molecule_id, property_type_id) DO UPDATE
        SET {value_field} = EXCLUDED.{value_field},
            updated_at = NOW()
        RETURNING id
    """

    try:
        result = execute_query(query, (molecule_id, property_type_id, value, source))
        if result and len(result) > 0:
            return result[0]["id"]
        else:
            return None
    except Exception as e:
        logger.error(f"Error inserting molecular property: {str(e)}")
        return None

def get_molecule_by_id(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule by ID.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        Molecule data or None if not found
    """
    query = """
        SELECT * FROM molecules 
        WHERE id = %s
    """
    
    try:
        result = execute_query(query, (molecule_id,))
        if result and len(result) > 0:
            return result[0]
        else:
            return None
    except Exception as e:
        logger.error(f"Error getting molecule: {str(e)}")
        return None

def get_molecular_properties(molecule_id: str) -> List[Dict[str, Any]]:
    """
    Get molecular properties for a molecule.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        List of properties
    """
    query = """
        SELECT mp.*, pt.name as property_name, pt.data_type 
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s
    """
    
    try:
        return execute_query(query, (molecule_id,))
    except Exception as e:
        logger.error(f"Error getting molecular properties: {str(e)}")
        return []

def test_connection() -> Tuple[bool, str]:
    """
    Test the database connection.
    
    Returns:
        Tuple of (success, message)
    """
    try:
        # Validate the configuration before testing
        validate_config()
        is_valid, message = test_adapter_configuration('supabase')
        if not is_valid:
            return False, f"Configuration error: {message}"
            
        # Test the actual connection
        result = execute_query("SELECT 1 AS test")
        if result and len(result) > 0 and result[0]["test"] == 1:
            return True, "Connection successful"
        else:
            return False, "Connection test failed"
    except Exception as e:
        logger.error(f"Connection test failed: {str(e)}")
        return False, f"Connection error: {str(e)}"

def close_all_connections():
    """Close all connections in the pool."""
    global _pool
    
    if _pool is not None:
        _pool.closeall()
        logger.info("All database connections closed")
        _pool = None
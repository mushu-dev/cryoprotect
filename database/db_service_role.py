#!/usr/bin/env python3
"""
Simplified database module with service role access.

This module provides direct connection to PostgreSQL using the service role,
which bypasses Row Level Security (RLS) policies.
"""

import os
import sys
import logging
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from contextlib import contextmanager
from typing import Dict, List, Any, Optional, Tuple, Union

# Import from the connection configuration module
from .connection_config import (
    validate_config,
    get_connection_config,
    test_adapter_configuration
)

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
        
    # Ensure options for service role
    if 'options' not in config:
        config['options'] = "-c role=service_role -c statement_timeout=60000"
    elif 'role=service_role' not in config['options']:
        config['options'] += " -c role=service_role"
    
    try:
        # Prepare connection parameters
        conn_params = {
            'host': config['host'],
            'port': config.get('port', 5432),
            'dbname': config.get('database', 'postgres'),
            'user': config['user'],
            'password': config['password'],
            'options': config.get('options'),
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
            # Try to disable RLS for this session
            try:
                cursor.execute("SET LOCAL role = 'service_role';")
                cursor.execute("SET LOCAL row_security = off;")
            except Exception as e:
                logger.warning(f"Could not set service role settings: {str(e)}")
            
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
            # Try to disable RLS for this session
            try:
                cursor.execute("SET LOCAL role = 'service_role';")
                cursor.execute("SET LOCAL row_security = off;")
            except Exception as e:
                logger.warning(f"Could not set service role settings: {str(e)}")
                
            execute_batch(cursor, query, param_sets)
            conn.commit()
            return True
    except Exception as e:
        conn.rollback()
        logger.error(f"Error executing batch: {str(e)}")
        return False
    finally:
        release_connection(conn)

@contextmanager
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
            
            # Try to disable RLS for this session
            try:
                self.cursor.execute("SET LOCAL role = 'service_role';")
                self.cursor.execute("SET LOCAL row_security = off;")
            except Exception as e:
                logger.warning(f"Could not set service role settings: {str(e)}")
            
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

def get_tables() -> List[str]:
    """
    Get a list of tables in the database.
    
    Returns:
        List of table names
    """
    result = execute_query("""
        SELECT table_name 
        FROM information_schema.tables 
        WHERE table_schema = 'public'
        ORDER BY table_name
    """)
    return [row["table_name"] for row in result]

def close_all_connections():
    """Close all connections in the pool."""
    global _pool
    
    if _pool is not None:
        _pool.closeall()
        logger.info("All database connections closed")
        _pool = None
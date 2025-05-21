"""
Simple database module for CryoProtect.

This module provides straightforward PostgreSQL connection management
for directly connecting to Supabase or local PostgreSQL databases.
"""

import os
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from psycopg2.pool import ThreadedConnectionPool
from dotenv import load_dotenv
from typing import Any, Dict, List, Optional, Union, Tuple

# Import from the connection configuration module
from .connection_config import (
    validate_config,
    get_connection_config,
    test_adapter_configuration
)

# Configure logger
logger = logging.getLogger(__name__)

# Global connection pool
_pool = None

def init_connection_pool(min_connections=None, max_connections=None, config=None):
    """
    Initialize the database connection pool.
    
    Args:
        min_connections: Minimum connections in pool (default: from config)
        max_connections: Maximum connections in pool (default: from config)
        config: Optional configuration dictionary
        
    Returns:
        bool: True if initialization successful, False otherwise
    """
    global _pool
    
    try:
        # Validate configuration
        validate_config()
        
        # Get configuration if not provided
        if not config:
            config = get_connection_config('local')
        
        # Set min/max connections from config if not provided
        if min_connections is None:
            min_connections = config.get('min_connections', 1)
        if max_connections is None:
            max_connections = config.get('max_connections', 10)
        
        # Check for required config values
        required_keys = ['host', 'user', 'password']
        missing_keys = [key for key in required_keys if key not in config or config[key] is None or config[key] == '']
        
        if missing_keys:
            logger.error(f"Missing required database configuration: {', '.join(missing_keys)}")
            return False
            
        # Prepare connection parameters
        conn_params = {
            'host': config['host'],
            'port': config.get('port', 5432),
            'dbname': config.get('database', 'postgres'),
            'user': config['user'],
            'password': config['password']
        }
        
        # Add application name if provided
        if 'application_name' in config:
            conn_params['application_name'] = config['application_name']
            
        # Initialize connection pool
        _pool = ThreadedConnectionPool(
            minconn=min_connections,
            maxconn=max_connections,
            **conn_params
        )
        
        logger.info(f"Database connection pool initialized ({min_connections}-{max_connections} connections)")
        return True
        
    except Exception as e:
        logger.error(f"Failed to initialize database connection pool: {str(e)}")
        return False

def get_connection():
    """
    Get a connection from the pool.
    
    Returns:
        Connection object or None if pool not initialized
    """
    global _pool
    
    if not _pool:
        logger.warning("Attempting to get connection before pool initialization")
        # Get configuration from the configuration system
        config = get_connection_config('local')
        if not init_connection_pool(config=config):
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
        conn: Connection to release
    """
    global _pool
    
    if not _pool:
        logger.warning("Attempting to release connection with no pool")
        return
        
    try:
        _pool.putconn(conn)
    except Exception as e:
        logger.error(f"Error releasing connection to pool: {str(e)}")

def close_all_connections():
    """
    Close all connections in the pool.
    
    Returns:
        bool: True if successful, False otherwise
    """
    global _pool
    
    if not _pool:
        logger.warning("Attempting to close connections with no pool")
        return True
        
    try:
        _pool.closeall()
        _pool = None
        logger.info("All database connections closed")
        return True
    except Exception as e:
        logger.error(f"Error closing all connections: {str(e)}")
        return False

def execute_query(query, params=None):
    """
    Execute a SQL query and return the results.
    
    Args:
        query: SQL query string
        params: Query parameters
        
    Returns:
        Query results or None on error
    """
    conn = None
    
    try:
        conn = get_connection()
        if not conn:
            return None
            
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute(query, params)
            
            if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                result = cursor.fetchall()
                return result
            else:
                conn.commit()
                return cursor.rowcount
                
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error executing query: {str(e)}")
        return None
    finally:
        if conn:
            release_connection(conn)

def execute_batch(queries):
    """
    Execute multiple SQL queries in a batch.
    
    Args:
        queries: List of SQL queries or tuple pairs of (query, params)
        
    Returns:
        List of results or None on error
    """
    conn = None
    results = []
    
    try:
        conn = get_connection()
        if not conn:
            return None
            
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            for item in queries:
                if isinstance(item, tuple) and len(item) == 2:
                    query, params = item
                else:
                    query, params = item, None
                    
                cursor.execute(query, params)
                
                if query.strip().upper().startswith('SELECT') or 'RETURNING' in query.upper():
                    results.append(cursor.fetchall())
                else:
                    conn.commit()
                    results.append(cursor.rowcount)
                    
        return results
        
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error executing batch: {str(e)}")
        return None
    finally:
        if conn:
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

def test_connection():
    """
    Test the database connection.
    
    Returns:
        (bool, str): Success status and message
    """
    try:
        # Validate the configuration before testing
        validate_config()
        is_valid, message = test_adapter_configuration('local')
        if not is_valid:
            return False, f"Configuration error: {message}"
            
        # Test the actual connection
        result = execute_query("SELECT 1 as test")
        if result and result[0]['test'] == 1:
            return True, "Connection successful"
        else:
            return False, "Connection test failed"
    except Exception as e:
        return False, f"Connection error: {str(e)}"

def get_tables():
    """
    Get list of tables in the database.
    
    Returns:
        List of table names or None on error
    """
    query = """
        SELECT table_name 
        FROM information_schema.tables 
        WHERE table_schema = 'public'
        ORDER BY table_name
    """
    
    try:
        result = execute_query(query)
        if result:
            return [row['table_name'] for row in result]
        return []
    except Exception as e:
        logger.error(f"Error getting table list: {str(e)}")
        return None

def get_table_info(table_name):
    """
    Get column information for a table.
    
    Args:
        table_name: Name of the table
        
    Returns:
        List of column details or None on error
    """
    query = """
        SELECT column_name, data_type, is_nullable
        FROM information_schema.columns
        WHERE table_name = %s AND table_schema = 'public'
        ORDER BY ordinal_position
    """
    
    try:
        return execute_query(query, (table_name,))
    except Exception as e:
        logger.error(f"Error getting table info for {table_name}: {str(e)}")
        return None
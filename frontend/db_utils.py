"""
Database utility functions for CryoProtect.
This module provides functions for working with the database.
"""

import sys
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from psycopg2 import pool
from contextlib import contextmanager
import db_config

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Global connection pool
connection_pool = None

def init_connection_pool():
    """Initialize the connection pool."""
    global connection_pool
    if connection_pool is not None:
        return
    
    try:
        # Get connection parameters
        params = db_config.get_db_connection_params()
        pool_config = db_config.get_connection_pool_config()
        
        # Create the connection pool
        connection_pool = pool.ThreadedConnectionPool(
            minconn=pool_config["minconn"],
            maxconn=pool_config["maxconn"],
            **params
        )
        logger.info("Database connection pool initialized")
    except Exception as e:
        logger.error(f"Failed to create connection pool: {e}")
        raise

@contextmanager
def get_db_connection():
    """
    Context manager for getting a database connection from the pool.
    Usage:
        with get_db_connection() as conn:
            # Use the connection
    """
    global connection_pool
    if connection_pool is None:
        init_connection_pool()
    
    conn = None
    try:
        conn = connection_pool.getconn()
        yield conn
    except Exception as e:
        logger.error(f"Error getting connection from pool: {e}")
        raise
    finally:
        if conn is not None:
            connection_pool.putconn(conn)

@contextmanager
def get_db_cursor(cursor_factory=None):
    """
    Context manager for getting a database cursor.
    Usage:
        with get_db_cursor() as cur:
            # Use the cursor
    """
    with get_db_connection() as conn:
        cursor = None
        try:
            if cursor_factory:
                cursor = conn.cursor(cursor_factory=cursor_factory)
            else:
                cursor = conn.cursor()
            yield cursor
            conn.commit()
        except Exception as e:
            conn.rollback()
            logger.error(f"Database error: {e}")
            raise
        finally:
            if cursor is not None:
                cursor.close()

def execute_sql_file(file_path, params=None):
    """
    Execute SQL from a file.
    Args:
        file_path: Path to the SQL file
        params: Optional parameters to pass to the SQL
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(file_path, 'r') as f:
            sql = f.read()
        
        with get_db_cursor() as cursor:
            if params:
                cursor.execute(sql, params)
            else:
                cursor.execute(sql)
        return True
    except Exception as e:
        logger.error(f"Error executing SQL file {file_path}: {e}")
        return False

def execute_query(query, params=None, fetch=True, cursor_factory=None):
    """
    Execute a SQL query.
    Args:
        query: SQL query to execute
        params: Optional parameters to pass to the query
        fetch: Whether to fetch and return results
        cursor_factory: Optional cursor factory (e.g., RealDictCursor)
    Returns:
        Query results if fetch=True, None otherwise
    """
    try:
        with get_db_cursor(cursor_factory=cursor_factory) as cursor:
            cursor.execute(query, params)
            if fetch:
                return cursor.fetchall()
            return None
    except Exception as e:
        logger.error(f"Error executing query: {e}")
        logger.error(f"Query: {query}")
        logger.error(f"Params: {params}")
        raise

def test_connection():
    """Test the database connection."""
    try:
        result = execute_query("SELECT 1 as test", cursor_factory=RealDictCursor)
        if result and result[0]['test'] == 1:
            logger.info("Database connection test successful")
            return True
        else:
            logger.error("Database connection test failed")
            return False
    except Exception as e:
        logger.error(f"Database connection test failed: {e}")
        return False

if __name__ == "__main__":
    # Test the database connection if this module is run directly
    if test_connection():
        print("Database connection test successful")
        sys.exit(0)
    else:
        print("Database connection test failed")
        sys.exit(1)
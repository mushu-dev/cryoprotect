#!/usr/bin/env python3
"""
CryoProtect v2 - Direct Supabase Connection

This module provides a direct PostgreSQL connection to Supabase using psycopg2
with connection pooling for improved performance over MCP.
"""

import os
import threading
import logging
import time
from typing import Dict, List, Any, Optional, Union
import psycopg2
import psycopg2.pool
import psycopg2.extras

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Default connection parameters
DEFAULT_DB_HOST = None  # Must be provided
DEFAULT_DB_PORT = 5432
DEFAULT_DB_NAME = "postgres"
DEFAULT_DB_USER = "postgres"
DEFAULT_DB_PASSWORD = None  # Must be provided
DEFAULT_DB_MIN_CONNECTIONS = 2
DEFAULT_DB_MAX_CONNECTIONS = 10
DEFAULT_DB_CONNECTION_TIMEOUT = 30


class SupabaseDirectConnection:
    """
    Singleton class for managing direct PostgreSQL connections to Supabase.
    Provides thread-safe connection pooling and SQL execution functions.
    """
    _instance = None
    _lock = threading.RLock()
    
    def __init__(self):
        """
        Initialize the connection pool. This should not be called directly.
        Use get_instance() instead.
        """
        if SupabaseDirectConnection._instance is not None:
            raise RuntimeError("Use SupabaseDirectConnection.get_instance() instead")
        
        # Read connection parameters from environment variables
        self.db_host = os.environ.get('SUPABASE_DB_HOST', DEFAULT_DB_HOST)
        self.db_port = int(os.environ.get('SUPABASE_DB_PORT', DEFAULT_DB_PORT))
        self.db_name = os.environ.get('SUPABASE_DB_NAME', DEFAULT_DB_NAME)
        self.db_user = os.environ.get('SUPABASE_DB_USER', DEFAULT_DB_USER)
        self.db_password = os.environ.get('SUPABASE_DB_PASSWORD', DEFAULT_DB_PASSWORD)
        self.min_connections = int(os.environ.get('SUPABASE_DB_MIN_CONNECTIONS', DEFAULT_DB_MIN_CONNECTIONS))
        self.max_connections = int(os.environ.get('SUPABASE_DB_MAX_CONNECTIONS', DEFAULT_DB_MAX_CONNECTIONS))
        
        # Validate required parameters
        if not self.db_host:
            raise ValueError("SUPABASE_DB_HOST environment variable is required")
        if not self.db_password:
            raise ValueError("SUPABASE_DB_PASSWORD environment variable is required")
        
        # Initialize connection pool
        self._initialize_pool()
        
        # Statistics
        self.query_count = 0
        self.error_count = 0
        self.last_error = None
        self.last_query_time = 0
        
        logger.info(f"Initialized SupabaseDirectConnection with pool size {self.min_connections}-{self.max_connections}")
    
    @classmethod
    def get_instance(cls) -> 'SupabaseDirectConnection':
        """
        Get the singleton instance of SupabaseDirectConnection.
        
        Returns:
            The singleton instance
        """
        with cls._lock:
            if cls._instance is None:
                cls._instance = SupabaseDirectConnection()
            return cls._instance
    
    def _initialize_pool(self):
        """Initialize the connection pool."""
        try:
            self.pool = psycopg2.pool.ThreadedConnectionPool(
                minconn=self.min_connections,
                maxconn=self.max_connections,
                host=self.db_host,
                port=self.db_port,
                dbname=self.db_name,
                user=self.db_user,
                password=self.db_password,
                connect_timeout=DEFAULT_DB_CONNECTION_TIMEOUT
            )
            logger.info("Connection pool initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing connection pool: {str(e)}")
            raise
    
    def get_connection(self):
        """
        Get a connection from the pool.
        
        Returns:
            A database connection
            
        Raises:
            psycopg2.pool.PoolError: If unable to get a connection from the pool
        """
        try:
            conn = self.pool.getconn()
            return conn
        except Exception as e:
            logger.error(f"Error getting connection from pool: {str(e)}")
            self.error_count += 1
            self.last_error = str(e)
            raise
    
    def release_connection(self, conn):
        """
        Release a connection back to the pool.
        
        Args:
            conn: The connection to release
        """
        try:
            self.pool.putconn(conn)
        except Exception as e:
            logger.error(f"Error releasing connection to pool: {str(e)}")
            # If we can't return it to the pool, try to close it
            try:
                conn.close()
            except:
                pass
    
    def close_all(self):
        """Close all connections in the pool."""
        try:
            self.pool.closeall()
            logger.info("All connections closed")
        except Exception as e:
            logger.error(f"Error closing connections: {str(e)}")
    
    def check_connection_health(self) -> bool:
        """
        Check if the connection to the database is healthy.
        
        Returns:
            True if the connection is healthy, False otherwise
        """
        conn = None
        try:
            conn = self.get_connection()
            with conn.cursor() as cursor:
                cursor.execute("SELECT 1")
                result = cursor.fetchone()
                return result[0] == 1
        except Exception as e:
            logger.error(f"Connection health check failed: {str(e)}")
            return False
        finally:
            if conn:
                self.release_connection(conn)
    
    def execute_query(self, query: str, params: Optional[Dict[str, Any]] = None) -> Optional[List[Dict[str, Any]]]:
        """
        Execute a SQL query and return the results as a list of dictionaries.
        
        Args:
            query: The SQL query to execute
            params: Optional parameters for the query
            
        Returns:
            A list of dictionaries representing the query results, or None for non-SELECT queries
            
        Raises:
            Exception: If an error occurs during query execution
        """
        conn = None
        start_time = time.time()
        
        try:
            conn = self.get_connection()
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                cursor.execute(query, params)
                
                # For SELECT queries, return the results
                if cursor.description:
                    result = cursor.fetchall()
                    # Convert from RealDictRow to regular dict
                    result = [dict(row) for row in result]
                else:
                    result = None
                
                conn.commit()
                self.query_count += 1
                self.last_query_time = time.time() - start_time
                return result
                
        except Exception as e:
            if conn:
                conn.rollback()
            self.error_count += 1
            self.last_error = str(e)
            logger.error(f"Error executing query: {str(e)}")
            raise
        finally:
            if conn:
                self.release_connection(conn)
    
    def execute_batch(self, queries: List[str], transaction: bool = True) -> bool:
        """
        Execute a batch of SQL queries, optionally in a transaction.
        
        Args:
            queries: A list of SQL queries to execute
            transaction: Whether to execute the queries in a transaction
            
        Returns:
            True if all queries were executed successfully, False otherwise
            
        Raises:
            Exception: If an error occurs during query execution and transaction=True
        """
        if not queries:
            return True
            
        conn = None
        start_time = time.time()
        
        try:
            conn = self.get_connection()
            with conn.cursor() as cursor:
                for query in queries:
                    cursor.execute(query)
                
                if transaction:
                    conn.commit()
                
                self.query_count += len(queries)
                self.last_query_time = time.time() - start_time
                return True
                
        except Exception as e:
            if conn and transaction:
                conn.rollback()
            self.error_count += 1
            self.last_error = str(e)
            logger.error(f"Error executing batch queries: {str(e)}")
            if transaction:
                raise
            return False
        finally:
            if conn:
                self.release_connection(conn)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the connection pool and query execution.
        
        Returns:
            A dictionary with statistics
        """
        return {
            "min_connections": self.min_connections,
            "max_connections": self.max_connections,
            "query_count": self.query_count,
            "error_count": self.error_count,
            "last_error": self.last_error,
            "last_query_time": self.last_query_time
        }


# Example usage
if __name__ == "__main__":
    # Set environment variables for testing
    # os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
    # os.environ['SUPABASE_DB_PASSWORD'] = 'your-password'
    
    # Get connection instance
    try:
        db = SupabaseDirectConnection.get_instance()
        
        # Execute query
        result = db.execute_query("SELECT 1 as test")
        print(f"Query Result: {result}")
        assert result[0]['test'] == 1
        print("Verification PASSED: SELECT 1 successful.")
        
        # Close connections when done
        db.close_all()
    except Exception as e:
        print(f"Verification FAILED: {e}")
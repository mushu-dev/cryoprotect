import os
import logging
import psycopg2
from psycopg2.pool import ThreadedConnectionPool
from psycopg2.extras import RealDictCursor
from typing import Any, Dict, List, Optional, Union, Tuple

from .adapter import DatabaseAdapter

logger = logging.getLogger(__name__)

class LocalPostgreSQLAdapter(DatabaseAdapter):
    """Local PostgreSQL database adapter implementation."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize adapter with configuration."""
        self.config = config
        self.connection_pool = None
        
    def connect(self) -> bool:
        """
        Establish connection pool to local PostgreSQL.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            # Extract connection parameters
            host = self.config.get('host', 'localhost')
            port = self.config.get('port', 5432)
            dbname = self.config.get('database', 'postgres')
            user = self.config.get('user', 'postgres')
            password = self.config.get('password', '')
            min_conn = int(self.config.get('min_connections', 1))
            max_conn = int(self.config.get('max_connections', 10))
            
            # Create connection pool
            self.connection_pool = ThreadedConnectionPool(
                minconn=min_conn,
                maxconn=max_conn,
                host=host,
                port=port,
                dbname=dbname,
                user=user,
                password=password
            )
            logger.info(f"Connected to local PostgreSQL at {host}:{port}/{dbname}")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to local PostgreSQL: {str(e)}")
            return False
            
    def disconnect(self) -> bool:
        """
        Close all connections in the pool.
        
        Returns:
            bool: True if disconnection successful, False otherwise
        """
        try:
            if self.connection_pool:
                self.connection_pool.closeall()
                logger.info("Disconnected from local PostgreSQL")
            return True
        except Exception as e:
            logger.error(f"Failed to disconnect from local PostgreSQL: {str(e)}")
            return False
            
    def execute_query(self, query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
        """
        Execute SQL query and return results.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            
        Returns:
            Query results
        """
        conn = None
        try:
            conn = self.connection_pool.getconn()
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
            raise
        finally:
            if conn:
                self.connection_pool.putconn(conn)
                
    def execute_batch(self, queries: List[str]) -> List[Any]:
        """
        Execute multiple SQL queries and return results.
        
        Args:
            queries: List of SQL queries to execute
            
        Returns:
            List of query results
        """
        conn = None
        results = []
        
        try:
            conn = self.connection_pool.getconn()
            with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                for query in queries:
                    cursor.execute(query)
                    
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
            raise
        finally:
            if conn:
                self.connection_pool.putconn(conn)
                
    def begin_transaction(self) -> Any:
        """
        Begin a database transaction.
        
        Returns:
            Connection object representing the transaction
        """
        try:
            conn = self.connection_pool.getconn()
            conn.autocommit = False
            return conn
        except Exception as e:
            logger.error(f"Error beginning transaction: {str(e)}")
            raise
            
    def commit_transaction(self, transaction: Any) -> bool:
        """
        Commit a database transaction.
        
        Args:
            transaction: Connection object representing the transaction
            
        Returns:
            bool: True if commit successful, False otherwise
        """
        try:
            transaction.commit()
            self.connection_pool.putconn(transaction)
            return True
        except Exception as e:
            logger.error(f"Error committing transaction: {str(e)}")
            return False
            
    def rollback_transaction(self, transaction: Any) -> bool:
        """
        Rollback a database transaction.
        
        Args:
            transaction: Connection object representing the transaction
            
        Returns:
            bool: True if rollback successful, False otherwise
        """
        try:
            transaction.rollback()
            self.connection_pool.putconn(transaction)
            return True
        except Exception as e:
            logger.error(f"Error rolling back transaction: {str(e)}")
            return False
            
    def get_connection_info(self) -> Dict[str, Any]:
        """
        Get connection information.
        
        Returns:
            Dict with connection information
        """
        return {
            'type': 'local_postgresql',
            'host': self.config.get('host', 'localhost'),
            'port': self.config.get('port', 5432),
            'database': self.config.get('database', 'postgres'),
            'user': self.config.get('user', 'postgres'),
            'pool_min_size': self.config.get('min_connections', 1),
            'pool_max_size': self.config.get('max_connections', 10),
            'connected': self.connection_pool is not None
        }
        
    def test_connection(self) -> Tuple[bool, str]:
        """
        Test database connection and return status with message.
        
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            result = self.execute_query("SELECT 1 as test")
            if result and result[0]['test'] == 1:
                return True, "Connection successful"
            else:
                return False, "Connection test failed"
        except Exception as e:
            return False, f"Connection error: {str(e)}"
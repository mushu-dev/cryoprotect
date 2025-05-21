# connection_pool_wrapper.py
import psycopg2
from psycopg2 import pool
import logging
import threading
import time
import socket
from typing import Dict, Any, Optional, List, Tuple

logger = logging.getLogger(__name__)

class ConnectionPoolWrapper:
    """
    Enhanced connection pool wrapper with lifecycle management, 
    health monitoring, and automatic reconnection.
    """
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the connection pool with configuration parameters.
        
        Args:
            config: Dictionary containing connection parameters
                - min_connections: Minimum connections in pool
                - max_connections: Maximum connections in pool
                - host: Database host
                - dbname: Database name
                - user: Database user
                - password: Database password
                - port: Database port
        """
        self.config = config
        self.min_conn = config.get('min_connections', 1)
        self.max_conn = config.get('max_connections', 10)
        self.pool = None
        self.health_check_interval = config.get('health_check_interval', 60)  # seconds
        self.connection_timeout = config.get('connection_timeout', 30)  # seconds
        self.pool_initialized = False
        self.active_connections = 0
        self.last_health_check = 0
        self._initialize_pool()
        
        # Start health check thread
        self._start_health_check()
    
    @classmethod
    def get_instance(cls, config: Optional[Dict[str, Any]] = None) -> 'ConnectionPoolWrapper':
        """
        Get singleton instance of connection pool wrapper.
        
        Args:
            config: Configuration dictionary (only used on first initialization)
            
        Returns:
            ConnectionPoolWrapper instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    if config is None:
                        raise ValueError("Config must be provided for initial pool creation")
                    cls._instance = cls(config)
        return cls._instance
    
    def _initialize_pool(self) -> None:
        """Initialize connection pool with hostname/IP fallback."""
        try:
            if self.pool is not None:
                # Close existing pool if it exists
                self.pool.closeall()
            
            # Create a copy of the config for connection
            conn_params = {
                'host': self.config.get('host'),
                'dbname': self.config.get('dbname'),
                'user': self.config.get('user'),
                'password': self.config.get('password'),
                'port': self.config.get('port'),
                'connect_timeout': self.connection_timeout
            }
            
            try:
                # First try connecting with hostname
                self.pool = psycopg2.pool.ThreadedConnectionPool(
                    self.min_conn,
                    self.max_conn,
                    **conn_params
                )
                self.pool_initialized = True
                self.active_connections = 0
                logger.info("Connection pool initialized with hostname %s", conn_params['host'])
            except (psycopg2.OperationalError, socket.gaierror) as e:
                # If hostname connection fails and we have an IP address, try with IP
                ip_address = self.config.get('ip_address')
                if ip_address:
                    logger.warning("Hostname connection failed (%s), trying with IP %s",
                                  str(e), ip_address)
                    conn_params['host'] = ip_address
                    self.pool = psycopg2.pool.ThreadedConnectionPool(
                        self.min_conn,
                        self.max_conn,
                        **conn_params
                    )
                    self.pool_initialized = True
                    self.active_connections = 0
                    logger.info("Connection pool initialized with IP %s", ip_address)
                else:
                    # Re-raise the exception if we don't have an IP fallback
                    logger.error("Hostname connection failed and no IP fallback available: %s", str(e))
                    raise
        except Exception as e:
            self.pool_initialized = False
            logger.error("Failed to initialize connection pool: %s", str(e))
            raise
    
    def _start_health_check(self) -> None:
        """Start a background thread for periodic health checks."""
        health_thread = threading.Thread(
            target=self._health_check_worker, 
            daemon=True
        )
        health_thread.start()
        logger.info("Health check thread started with interval %d seconds", 
                   self.health_check_interval)
    
    def _health_check_worker(self) -> None:
        """Worker function for periodic health checks."""
        while True:
            time.sleep(self.health_check_interval)
            try:
                self._check_pool_health()
            except Exception as e:
                logger.error("Health check failed: %s", str(e))
    
    def _check_pool_health(self) -> None:
        """Check health of connection pool and reinitialize if needed."""
        if not self.pool_initialized:
            logger.warning("Pool not initialized during health check, attempting to initialize")
            self._initialize_pool()
            return
            
        try:
            # Get a connection to test
            conn = self.pool.getconn()
            try:
                # Execute simple query to check connection
                with conn.cursor() as cursor:
                    cursor.execute("SELECT 1")
                    result = cursor.fetchone()
                    if result and result[0] == 1:
                        logger.debug("Connection pool health check passed")
                    else:
                        logger.warning("Connection returned unexpected result during health check")
                        self._initialize_pool()
            except Exception as e:
                logger.error("Connection test failed during health check: %s", str(e))
                self._initialize_pool()
            finally:
                # Return the connection to the pool
                self.pool.putconn(conn)
        except Exception as e:
            logger.error("Failed to get connection during health check: %s", str(e))
            self._initialize_pool()
        
        self.last_health_check = time.time()
    
    def get_connection(self) -> Tuple[Any, int]:
        """
        Get a connection from the pool with a unique identifier.
        
        Returns:
            Tuple of (connection, connection_id)
        """
        if not self.pool_initialized:
            self._initialize_pool()
            
        try:
            conn = self.pool.getconn()
            self.active_connections += 1
            conn_id = id(conn)
            logger.debug("Connection %d acquired from pool (active: %d)", 
                        conn_id, self.active_connections)
            return conn, conn_id
        except Exception as e:
            logger.error("Failed to get connection from pool: %s", str(e))
            # Try to reinitialize pool
            self._initialize_pool()
            # Retry once
            conn = self.pool.getconn()
            self.active_connections += 1
            conn_id = id(conn)
            logger.debug("Connection %d acquired from pool after retry (active: %d)", 
                        conn_id, self.active_connections)
            return conn, conn_id
    
    def return_connection(self, conn: Any, conn_id: int) -> None:
        """
        Return a connection to the pool.
        
        Args:
            conn: The connection to return
            conn_id: The connection identifier
        """
        if self.pool_initialized and conn is not None:
            try:
                self.pool.putconn(conn)
                self.active_connections -= 1
                logger.debug("Connection %d returned to pool (active: %d)", 
                            conn_id, self.active_connections)
            except Exception as e:
                logger.error("Failed to return connection %d to pool: %s", 
                            conn_id, str(e))
    
    def close_all(self) -> None:
        """Close all connections in the pool."""
        if self.pool_initialized:
            try:
                self.pool.closeall()
                logger.info("Closed all connections in the pool")
                self.pool_initialized = False
                self.active_connections = 0
            except Exception as e:
                logger.error("Failed to close all connections: %s", str(e))

# Helper functions for context manager usage
def get_db_connection():
    """Get a connection from the global connection pool."""
    from config import get_db_config
    pool = ConnectionPoolWrapper.get_instance(get_db_config())
    return pool.get_connection()

def return_db_connection(conn, conn_id):
    """Return a connection to the global connection pool."""
    from config import get_db_config
    pool = ConnectionPoolWrapper.get_instance(get_db_config())
    pool.return_connection(conn, conn_id)

class ConnectionManager:
    """Context manager for database connections."""
    
    def __init__(self):
        self.conn = None
        self.conn_id = None
    
    def __enter__(self):
        self.conn, self.conn_id = get_db_connection()
        return self.conn
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return_db_connection(self.conn, self.conn_id)
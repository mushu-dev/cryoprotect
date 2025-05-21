"""
Lightweight connection pool for Supabase REST API requests.
"""

import os
import time
import threading
import logging
import requests
from queue import Queue, Empty
from typing import Dict, Any, Optional, List, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ConnectionPool:
    """
    A simple connection pool that maintains a queue of reusable connection objects.
    Each connection is a tuple of (headers, created_at) where headers contain
    the necessary authentication information for Supabase REST API requests.
    """
    
    def __init__(self, min_size=2, max_size=10, timeout=1800):
        """
        Initialize the connection pool.
        
        Args:
            min_size: Minimum number of connections to maintain
            max_size: Maximum number of connections allowed
            timeout: Connection timeout in seconds (connections older than this will be recreated)
        """
        self.min_size = min_size
        self.max_size = max_size
        self.timeout = timeout
        self.pool = Queue(maxsize=max_size)
        self.active_connections = 0
        self.lock = threading.RLock()
        self.stats = {
            "created": 0,
            "reused": 0,
            "discarded": 0,
            "errors": 0
        }
        
        # Initialize the pool with minimum connections
        self._initialize_pool()
    
    def _initialize_pool(self):
        """Initialize the pool with minimum number of connections."""
        logger.info(f"Initializing connection pool with {self.min_size} connections")
        for _ in range(self.min_size):
            try:
                connection = self._create_connection()
                self.pool.put(connection)
            except Exception as e:
                logger.error(f"Error initializing connection: {str(e)}")
                self.stats["errors"] += 1
    
    def _create_connection(self) -> Dict[str, Any]:
        """
        Create a new connection (headers with authentication info).
        
        Returns:
            Dict with connection information
        """
        # Get API keys from environment
        api_key = os.environ.get("SUPABASE_KEY")
        service_key = os.environ.get("SUPABASE_SERVICE_KEY", api_key)
        
        if not api_key:
            raise ValueError("SUPABASE_KEY environment variable is not set")
        
        # Create headers for both anonymous and service role access
        anon_headers = {
            "apikey": api_key,
            "Authorization": f"Bearer {api_key}"
        }
        
        service_headers = {
            "apikey": service_key,
            "Authorization": f"Bearer {service_key}"
        }
        
        # Store both header sets in the connection object
        connection = {
            "anon": anon_headers,
            "service": service_headers,
            "created_at": time.time()
        }
        
        with self.lock:
            self.active_connections += 1
            self.stats["created"] += 1
        
        return connection
    
    def get_connection(self) -> Dict[str, Any]:
        """
        Get a connection from the pool or create a new one if needed.
        
        Returns:
            Connection object with headers
        """
        connection = None
        created_new = False
        
        # Try to get a connection from the pool
        try:
            connection = self.pool.get(block=False)
            
            # Check if the connection is still valid (not timed out)
            if time.time() - connection["created_at"] > self.timeout:
                # Connection is too old, discard it and create a new one
                with self.lock:
                    self.active_connections -= 1
                    self.stats["discarded"] += 1
                
                connection = None
            else:
                with self.lock:
                    self.stats["reused"] += 1
        except Empty:
            # Pool is empty, we'll create a new connection below
            pass
        
        # Create a new connection if needed
        if not connection:
            try:
                with self.lock:
                    if self.active_connections < self.max_size:
                        connection = self._create_connection()
                        created_new = True
                    else:
                        # We've reached max_size, wait for a connection to become available
                        connection = self.pool.get(block=True, timeout=10)
                        self.stats["reused"] += 1
            except Exception as e:
                logger.error(f"Error getting connection: {str(e)}")
                self.stats["errors"] += 1
                raise
        
        return connection
    
    def release_connection(self, connection: Dict[str, Any]):
        """
        Release a connection back to the pool.
        
        Args:
            connection: Connection object to release
        """
        try:
            # Check if the connection is still valid (not timed out)
            if time.time() - connection["created_at"] > self.timeout:
                # Connection is too old, discard it
                with self.lock:
                    self.active_connections -= 1
                    self.stats["discarded"] += 1
            else:
                # Put the connection back in the pool
                self.pool.put(connection)
        except Exception as e:
            logger.error(f"Error releasing connection: {str(e)}")
            self.stats["errors"] += 1
            
            # Ensure we don't leak connections
            with self.lock:
                self.active_connections -= 1
    
    def close_all(self):
        """Close all connections in the pool."""
        logger.info("Closing all connections in the pool")
        
        # Clear the pool and reset counters
        while not self.pool.empty():
            try:
                self.pool.get(block=False)
                with self.lock:
                    self.active_connections -= 1
            except Empty:
                break
        
        with self.lock:
            self.active_connections = 0
    
    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about the connection pool."""
        with self.lock:
            stats = self.stats.copy()
            stats["pool_size"] = self.pool.qsize()
            stats["active_connections"] = self.active_connections
            stats["min_size"] = self.min_size
            stats["max_size"] = self.max_size
            return stats

# Global connection pool instance
_pool_instance = None

def get_pool() -> ConnectionPool:
    """
    Get the global connection pool instance.
    
    Returns:
        ConnectionPool instance
    """
    global _pool_instance
    
    if _pool_instance is None:
        # Get pool configuration from environment variables
        min_size = int(os.environ.get("DB_POOL_MIN_SIZE", "2"))
        max_size = int(os.environ.get("DB_POOL_MAX_SIZE", "10"))
        timeout = int(os.environ.get("DB_POOL_TIMEOUT", "1800"))
        
        _pool_instance = ConnectionPool(
            min_size=min_size,
            max_size=max_size,
            timeout=timeout
        )
    
    return _pool_instance

# Convenience functions for database operations

def rest_request(
    method: str,
    path: str,
    data: Any = None,
    params: Dict[str, Any] = None,
    use_service_role: bool = False,
    timeout: int = 10
) -> requests.Response:
    """
    Make a request to Supabase REST API using the connection pool.
    
    Args:
        method: HTTP method (GET, POST, PUT, PATCH, DELETE)
        path: API path (without /rest/v1/ prefix)
        data: Request data (for POST, PUT, PATCH)
        params: Query parameters
        use_service_role: Whether to use service role key instead of anon key
        timeout: Request timeout in seconds
    
    Returns:
        requests.Response object
    """
    # Get Supabase URL from environment
    supabase_url = os.environ.get("SUPABASE_URL")
    if not supabase_url:
        raise ValueError("SUPABASE_URL environment variable is not set")
    
    # Build URL (remove leading slash if present)
    if path.startswith('/'):
        path = path[1:]
    
    url = f"{supabase_url}/rest/v1/{path}"
    
    # Get a connection from the pool
    pool = get_pool()
    connection = pool.get_connection()
    
    try:
        # Use service role or anonymous headers
        headers = connection["service"] if use_service_role else connection["anon"]
        
        # Make the request
        response = requests.request(
            method=method,
            url=url,
            headers=headers,
            json=data,
            params=params,
            timeout=timeout
        )
        
        return response
    finally:
        # Always release the connection back to the pool
        pool.release_connection(connection)

def get_table_data(
    table_name: str,
    filters: Dict[str, Any] = None,
    limit: int = 100,
    offset: int = 0,
    order_by: str = None,
    select: str = "*",
    use_service_role: bool = False
) -> List[Dict[str, Any]]:
    """
    Get data from a table using the connection pool.
    
    Args:
        table_name: Table name
        filters: Filter conditions
        limit: Maximum number of records to return
        offset: Offset for pagination
        order_by: Order by clause
        select: Columns to select
        use_service_role: Whether to use service role key
    
    Returns:
        List of records
    """
    params = {
        "limit": limit,
        "offset": offset
    }
    
    # Add select parameter if not default
    if select != "*":
        params["select"] = select
    
    # Add order parameter if provided
    if order_by:
        params["order"] = order_by
    
    # Add filter parameters
    if filters:
        params.update(filters)
    
    response = rest_request(
        method="GET",
        path=table_name,
        params=params,
        use_service_role=use_service_role
    )
    
    response.raise_for_status()
    return response.json()

def get_pool_stats() -> Dict[str, Any]:
    """
    Get statistics about the connection pool.
    
    Returns:
        Dict with pool statistics
    """
    pool = get_pool()
    return pool.get_stats()
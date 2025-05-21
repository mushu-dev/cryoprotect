# Simplified Connection Pooling for Fedora

This document provides a practical approach to implementing connection pooling in the Fedora environment for CryoProtect.

## Overview

Connection pooling allows us to:
1. Reuse database connections instead of creating new ones for each request
2. Reduce connection establishment overhead
3. Handle peak loads more efficiently
4. Set limits on maximum concurrent connections

In our Fedora environment, we'll implement a simple connection pool using a lightweight approach that works with our simplified application.

## Implementation

### 1. Create Connection Pool Module

Create a new file `db_pool.py` with a lightweight connection pool implementation:

```python
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
    
    def _create_connection(self) -> Tuple[Dict[str, str], float]:
        """
        Create a new connection (headers with authentication info).
        
        Returns:
            Tuple of (headers, created_at_timestamp)
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
```

### 2. Update Simplified App to Use Connection Pool

Modify `simplified_app.py` to use the connection pool:

```python
"""
Simplified CryoProtect app with connection pooling.
"""

import os
import sys
import logging
import json
from flask import Flask, jsonify, request
from dotenv import load_dotenv

# Import our connection pool
from db_pool import get_table_data, rest_request, get_pool_stats

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Create Flask app
app = Flask(__name__)

# Get Supabase credentials from environment
SUPABASE_URL = os.environ.get('SUPABASE_URL')
SUPABASE_KEY = os.environ.get('SUPABASE_KEY')

if not SUPABASE_URL or not SUPABASE_KEY:
    logger.error("Supabase URL or key not configured. Please set SUPABASE_URL and SUPABASE_KEY environment variables.")

@app.route('/')
def index():
    """Root endpoint with basic information"""
    return jsonify({
        'name': 'CryoProtect Simplified App with Connection Pooling',
        'status': 'running',
        'environment': os.environ.get('FLASK_ENV', 'development'),
        'python_version': sys.version,
        'supabase_configured': bool(SUPABASE_URL and SUPABASE_KEY)
    })

@app.route('/supabase/tables')
def list_tables():
    """List tables in Supabase using the connection pool"""
    try:
        # Get table information using connection pool
        response = rest_request('GET', '')
        response.raise_for_status()
        
        # Parse OpenAPI schema to get table names
        schema = response.json()
        paths = schema.get('paths', {})
        
        tables = []
        for path in paths:
            parts = path.strip('/').split('/')
            if len(parts) >= 3 and parts[0] == 'rest' and parts[1] == 'v1':
                table_name = parts[2]
                if table_name and '.' not in table_name and '{' not in table_name:
                    tables.append(table_name)
        
        return jsonify({
            'status': 'success',
            'tables': sorted(list(set(tables)))
        })
    except Exception as e:
        logger.error(f'Error listing tables: {str(e)}')
        return jsonify({
            'status': 'error',
            'message': f'Exception: {str(e)}'
        }), 500

@app.route('/api/molecules')
def list_molecules():
    """List molecules using the connection pool"""
    try:
        # Get query parameters
        limit = request.args.get('limit', default=10, type=int)
        offset = request.args.get('offset', default=0, type=int)
        name_filter = request.args.get('name', default=None)
        
        # Build filters
        filters = {}
        if name_filter:
            filters['name'] = f'ilike.{name_filter}%'
        
        # Get data using the connection pool
        molecules = get_table_data(
            table_name='molecules',
            filters=filters,
            limit=limit,
            offset=offset,
            order_by='name.asc'
        )
        
        return jsonify({
            'status': 'success',
            'count': len(molecules),
            'molecules': molecules
        })
    except Exception as e:
        logger.error(f'Error listing molecules: {str(e)}')
        return jsonify({
            'status': 'error',
            'message': f'Exception: {str(e)}'
        }), 500

@app.route('/pool/stats')
def pool_stats():
    """Show connection pool statistics"""
    try:
        stats = get_pool_stats()
        
        return jsonify({
            'status': 'success',
            'stats': stats
        })
    except Exception as e:
        logger.error(f'Error getting pool stats: {str(e)}')
        return jsonify({
            'status': 'error',
            'message': f'Exception: {str(e)}'
        }), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    host = os.environ.get('HOST', '0.0.0.0')
    debug = os.environ.get('FLASK_ENV', 'development') == 'development'
    
    print(f"Starting CryoProtect simplified app with connection pooling on {host}:{port} (debug={debug})")
    
    # Initialize connection pool
    from db_pool import get_pool
    pool = get_pool()
    print(f"Connection pool initialized (min={pool.min_size}, max={pool.max_size})")
    
    app.run(host=host, port=port, debug=debug)
```

## Testing Connection Pooling

You can test the connection pooling implementation with these steps:

1. Create the `db_pool.py` module
2. Update the simplified app to use the connection pool
3. Run the app and test the endpoints

```bash
# Start the app with connection pooling
python simplified_app.py

# In another terminal, test endpoints and pool
curl http://localhost:5000/
curl http://localhost:5000/supabase/tables
curl http://localhost:5000/api/molecules?limit=5
curl http://localhost:5000/pool/stats
```

## Performance Testing

To verify the performance benefits of connection pooling, you can create a simple load test:

```python
import requests
import time
import threading
from concurrent.futures import ThreadPoolExecutor

def make_request(i):
    """Make a request to the API"""
    start = time.time()
    response = requests.get('http://localhost:5000/api/molecules?limit=1')
    duration = time.time() - start
    
    return {
        'request_id': i,
        'status_code': response.status_code,
        'duration': duration
    }

def run_test(num_requests=100, max_workers=10):
    """Run a load test with multiple concurrent requests"""
    results = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(make_request, i) for i in range(num_requests)]
        
        for future in futures:
            results.append(future.result())
    
    # Calculate statistics
    durations = [r['duration'] for r in results]
    avg_duration = sum(durations) / len(durations)
    min_duration = min(durations)
    max_duration = max(durations)
    
    print(f"Results for {num_requests} requests with {max_workers} workers:")
    print(f"  Average: {avg_duration:.3f}s")
    print(f"  Min: {min_duration:.3f}s")
    print(f"  Max: {max_duration:.3f}s")
    
    # Get pool stats
    response = requests.get('http://localhost:5000/pool/stats')
    stats = response.json()['stats']
    print(f"Pool stats: {stats}")

if __name__ == '__main__':
    print("Running load test...")
    run_test(num_requests=100, max_workers=10)
```

Save this as `test_pool_performance.py` and run it after starting the app.

## Benefits

This simple connection pooling implementation provides significant benefits:

1. **Performance**: Reuses connections instead of creating new ones for each request
2. **Stability**: Limits the number of concurrent connections to prevent overloading
3. **Reliability**: Handles connection errors and timeouts gracefully
4. **Monitoring**: Provides statistics for monitoring and troubleshooting
5. **Simplicity**: Uses a straightforward approach without complex dependencies

## When to Use

Use connection pooling when:
1. Your application handles multiple concurrent requests
2. Connection establishment is expensive (e.g., SSL handshake, authentication)
3. You need to limit the number of concurrent connections
4. You want to monitor connection usage and performance

For very simple applications with few users, connection pooling might not be necessary, but it's still a good practice to implement.
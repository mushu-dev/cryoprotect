#!/usr/bin/env python3
"""
CryoProtect v2 - Implement Connection Pooling

This script implements connection pooling for the CryoProtect database with minimal downtime.
It configures connection pool settings, implements retry logic, includes monitoring functions,
and provides a graceful shutdown mechanism.

Usage:
    python implement_connection_pooling.py [--project-id PROJECT_ID] [--test]

Options:
    --project-id PROJECT_ID    Supabase project ID (default: tsdlmynydfuypiugmkev)
    --test                     Run a load test after implementing connection pooling
"""

import os
import sys
import time
import json
import logging
import argparse
import threading
import queue
import signal
import datetime
from typing import Dict, List, Any, Optional, Callable, Union
from pathlib import Path
from contextlib import contextmanager

# Import logging configuration
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from logging_config import setup_logging
except ImportError:
    # Fallback logging setup if the import fails
    def setup_logging():
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler("connection_pooling.log"),
                logging.StreamHandler()
            ]
        )

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import execute_sql_on_supabase, apply_migration
except ImportError:
    print("Error: supabase_mcp_tools.py not found. Make sure it's in the same directory.")
    sys.exit(1)

# Set up logging
setup_logging()
logger = logging.getLogger("connection_pooling")

# Connection pool settings
DEFAULT_MIN_CONNECTIONS = 2
DEFAULT_MAX_CONNECTIONS = 10
DEFAULT_CONNECTION_TIMEOUT = 30  # seconds
DEFAULT_CONNECTION_LIFETIME = 3600  # 1 hour
DEFAULT_IDLE_TIMEOUT = 300  # 5 minutes
DEFAULT_RETRY_ATTEMPTS = 3
DEFAULT_RETRY_DELAY = 1  # seconds

# Global connection pool
connection_pool = None
pool_lock = threading.RLock()
pool_stats = {
    "created": 0,
    "acquired": 0,
    "released": 0,
    "errors": 0,
    "retries": 0,
    "peak_usage": 0,
    "last_reset": datetime.datetime.now().isoformat()
}

class ConnectionPoolError(Exception):
    """Exception raised for connection pool errors."""
    pass

class Connection:
    """Represents a database connection in the pool."""
    
    def __init__(self, project_id: str):
        """
        Initialize a new connection.
        
        Args:
            project_id: The Supabase project ID
        """
        self.project_id = project_id
        self.created_at = datetime.datetime.now()
        self.last_used_at = self.created_at
        self.in_use = False
        self.error_count = 0
        self.query_count = 0
        self.id = f"conn-{id(self)}"
        logger.debug(f"Created new connection {self.id}")
        
    def execute(self, query: str) -> List[Dict[str, Any]]:
        """
        Execute a query using this connection.
        
        Args:
            query: The SQL query to execute
            
        Returns:
            The query results
            
        Raises:
            Exception: If the query execution fails
        """
        try:
            self.last_used_at = datetime.datetime.now()
            result = execute_sql_on_supabase(self.project_id, query)
            self.query_count += 1
            return result
        except Exception as e:
            self.error_count += 1
            raise e
            
    def is_expired(self, lifetime: int) -> bool:
        """
        Check if the connection has expired based on its lifetime.
        
        Args:
            lifetime: Maximum connection lifetime in seconds
            
        Returns:
            True if the connection has expired, False otherwise
        """
        age = (datetime.datetime.now() - self.created_at).total_seconds()
        return age > lifetime
        
    def is_idle(self, idle_timeout: int) -> bool:
        """
        Check if the connection is idle based on the idle timeout.
        
        Args:
            idle_timeout: Idle timeout in seconds
            
        Returns:
            True if the connection is idle, False otherwise
        """
        idle_time = (datetime.datetime.now() - self.last_used_at).total_seconds()
        return idle_time > idle_timeout and not self.in_use

class ConnectionPool:
    """A connection pool for managing database connections."""
    
    def __init__(self, project_id: str, min_connections: int = DEFAULT_MIN_CONNECTIONS, 
                 max_connections: int = DEFAULT_MAX_CONNECTIONS, 
                 connection_timeout: int = DEFAULT_CONNECTION_TIMEOUT,
                 connection_lifetime: int = DEFAULT_CONNECTION_LIFETIME,
                 idle_timeout: int = DEFAULT_IDLE_TIMEOUT):
        """
        Initialize the connection pool.
        
        Args:
            project_id: The Supabase project ID
            min_connections: Minimum number of connections to maintain
            max_connections: Maximum number of connections allowed
            connection_timeout: Timeout for acquiring a connection in seconds
            connection_lifetime: Maximum lifetime of a connection in seconds
            idle_timeout: Timeout for idle connections in seconds
        """
        self.project_id = project_id
        self.min_connections = min_connections
        self.max_connections = max_connections
        self.connection_timeout = connection_timeout
        self.connection_lifetime = connection_lifetime
        self.idle_timeout = idle_timeout
        
        self.connections = []
        self.available = queue.Queue()
        self.size = 0
        self.running = True
        
        # Initialize the minimum number of connections
        self._initialize_connections()
        
        # Start maintenance thread
        self.maintenance_thread = threading.Thread(target=self._maintenance_task, daemon=True)
        self.maintenance_thread.start()
        
        logger.info(f"Connection pool initialized with {min_connections} connections " +
                   f"(max: {max_connections}, timeout: {connection_timeout}s)")
    
    def _initialize_connections(self):
        """Initialize the minimum number of connections."""
        for _ in range(self.min_connections):
            conn = self._create_connection()
            self.available.put(conn)
    
    def _create_connection(self) -> Connection:
        """
        Create a new database connection.
        
        Returns:
            A new Connection object
        """
        conn = Connection(self.project_id)
        with pool_lock:
            self.connections.append(conn)
            self.size += 1
            pool_stats["created"] += 1
        return conn
    
    def _maintenance_task(self):
        """Background task to maintain the connection pool."""
        while self.running:
            try:
                # Sleep for a short interval
                time.sleep(10)
                
                if not self.running:
                    break
                
                with pool_lock:
                    # Check for expired or idle connections
                    current_time = datetime.datetime.now()
                    to_remove = []
                    
                    for conn in self.connections:
                        if not conn.in_use:
                            if (conn.is_expired(self.connection_lifetime) or 
                                conn.is_idle(self.idle_timeout)):
                                to_remove.append(conn)
                    
                    # Remove expired/idle connections
                    for conn in to_remove:
                        self.connections.remove(conn)
                        self.size -= 1
                        logger.debug(f"Removed expired/idle connection {conn.id}")
                    
                    # Ensure minimum connections
                    available_count = self.available.qsize()
                    if available_count < self.min_connections and self.size < self.max_connections:
                        needed = min(self.min_connections - available_count, 
                                    self.max_connections - self.size)
                        for _ in range(needed):
                            conn = self._create_connection()
                            self.available.put(conn)
                            logger.debug(f"Added new connection {conn.id} to maintain minimum")
                
                # Update peak usage statistic
                with pool_lock:
                    in_use = sum(1 for conn in self.connections if conn.in_use)
                    if in_use > pool_stats["peak_usage"]:
                        pool_stats["peak_usage"] = in_use
                
            except Exception as e:
                logger.error(f"Error in connection pool maintenance: {str(e)}")
    
    def acquire(self) -> Connection:
        """
        Acquire a connection from the pool.
        
        Returns:
            A Connection object
            
        Raises:
            ConnectionPoolError: If unable to acquire a connection
        """
        start_time = time.time()
        
        while time.time() - start_time < self.connection_timeout:
            # Try to get an available connection
            try:
                conn = self.available.get(block=False)
                
                # Check if the connection is expired
                if conn.is_expired(self.connection_lifetime):
                    with pool_lock:
                        self.connections.remove(conn)
                        self.size -= 1
                    # Create a new connection instead
                    conn = self._create_connection()
                
                conn.in_use = True
                conn.last_used_at = datetime.datetime.now()
                
                with pool_lock:
                    pool_stats["acquired"] += 1
                
                return conn
                
            except queue.Empty:
                # No available connections, create a new one if possible
                with pool_lock:
                    if self.size < self.max_connections:
                        conn = self._create_connection()
                        conn.in_use = True
                        conn.last_used_at = datetime.datetime.now()
                        pool_stats["acquired"] += 1
                        return conn
                
                # Wait a bit before trying again
                time.sleep(0.1)
        
        # Timeout reached
        with pool_lock:
            pool_stats["errors"] += 1
        raise ConnectionPoolError(f"Timeout acquiring connection after {self.connection_timeout}s")
    
    def release(self, conn: Connection):
        """
        Release a connection back to the pool.
        
        Args:
            conn: The Connection object to release
        """
        if conn not in self.connections:
            logger.warning(f"Attempted to release unknown connection {conn.id}")
            return
        
        conn.in_use = False
        self.available.put(conn)
        
        with pool_lock:
            pool_stats["released"] += 1
    
    def execute_with_retry(self, query: str, max_retries: int = DEFAULT_RETRY_ATTEMPTS, 
                          retry_delay: int = DEFAULT_RETRY_DELAY) -> List[Dict[str, Any]]:
        """
        Execute a query with retry logic.
        
        Args:
            query: The SQL query to execute
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries in seconds
            
        Returns:
            The query results
            
        Raises:
            Exception: If all retry attempts fail
        """
        last_error = None
        
        for attempt in range(max_retries + 1):
            conn = None
            try:
                conn = self.acquire()
                result = conn.execute(query)
                return result
                
            except Exception as e:
                last_error = e
                with pool_lock:
                    pool_stats["retries"] += 1
                
                if attempt < max_retries:
                    retry_wait = retry_delay * (2 ** attempt)  # Exponential backoff
                    logger.warning(f"Query failed, retrying in {retry_wait}s " +
                                  f"(attempt {attempt + 1}/{max_retries + 1}): {str(e)}")
                    time.sleep(retry_wait)
                
            finally:
                if conn:
                    self.release(conn)
        
        with pool_lock:
            pool_stats["errors"] += 1
        
        raise last_error or ConnectionPoolError("Query failed after retries")
    
    def shutdown(self, timeout: int = 30):
        """
        Shutdown the connection pool gracefully.
        
        Args:
            timeout: Maximum time to wait for connections to be released in seconds
        """
        logger.info("Shutting down connection pool...")
        self.running = False
        
        # Wait for maintenance thread to exit
        if self.maintenance_thread.is_alive():
            self.maintenance_thread.join(timeout=5)
        
        # Wait for all connections to be released
        start_time = time.time()
        while time.time() - start_time < timeout:
            with pool_lock:
                in_use = sum(1 for conn in self.connections if conn.in_use)
                if in_use == 0:
                    break
            time.sleep(0.5)
        
        # Log final stats
        with pool_lock:
            in_use = sum(1 for conn in self.connections if conn.in_use)
            logger.info(f"Connection pool shutdown complete. " +
                       f"{in_use} connections still in use out of {self.size} total.")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the connection pool.
        
        Returns:
            Dictionary with pool statistics
        """
        with pool_lock:
            stats = {
                "size": self.size,
                "available": self.available.qsize(),
                "in_use": sum(1 for conn in self.connections if conn.in_use),
                "min_connections": self.min_connections,
                "max_connections": self.max_connections,
                **pool_stats
            }
        return stats
    
    def reset_stats(self):
        """Reset the connection pool statistics."""
        with pool_lock:
            pool_stats["created"] = self.size
            pool_stats["acquired"] = 0
            pool_stats["released"] = 0
            pool_stats["errors"] = 0
            pool_stats["retries"] = 0
            pool_stats["peak_usage"] = sum(1 for conn in self.connections if conn.in_use)
            pool_stats["last_reset"] = datetime.datetime.now().isoformat()

@contextmanager
def get_connection():
    """
    Context manager for acquiring and releasing a connection.
    
    Yields:
        A Connection object
    """
    global connection_pool
    if connection_pool is None:
        raise ConnectionPoolError("Connection pool not initialized")
    
    conn = connection_pool.acquire()
    try:
        yield conn
    finally:
        connection_pool.release(conn)

def execute_query(query: str, max_retries: int = DEFAULT_RETRY_ATTEMPTS, 
                 retry_delay: int = DEFAULT_RETRY_DELAY) -> List[Dict[str, Any]]:
    """
    Execute a query using the connection pool with retry logic.
    
    Args:
        query: The SQL query to execute
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retries in seconds
        
    Returns:
        The query results
        
    Raises:
        Exception: If all retry attempts fail
    """
    global connection_pool
    if connection_pool is None:
        raise ConnectionPoolError("Connection pool not initialized")
    
    return connection_pool.execute_with_retry(query, max_retries, retry_delay)

def initialize_connection_pool(project_id: str, min_connections: int = DEFAULT_MIN_CONNECTIONS, 
                              max_connections: int = DEFAULT_MAX_CONNECTIONS, 
                              connection_timeout: int = DEFAULT_CONNECTION_TIMEOUT,
                              connection_lifetime: int = DEFAULT_CONNECTION_LIFETIME,
                              idle_timeout: int = DEFAULT_IDLE_TIMEOUT) -> ConnectionPool:
    """
    Initialize the global connection pool.
    
    Args:
        project_id: The Supabase project ID
        min_connections: Minimum number of connections to maintain
        max_connections: Maximum number of connections allowed
        connection_timeout: Timeout for acquiring a connection in seconds
        connection_lifetime: Maximum lifetime of a connection in seconds
        idle_timeout: Timeout for idle connections in seconds
        
    Returns:
        The initialized ConnectionPool instance
    """
    global connection_pool
    
    if connection_pool is not None:
        logger.warning("Connection pool already initialized, shutting down existing pool")
        connection_pool.shutdown()
    
    connection_pool = ConnectionPool(
        project_id=project_id,
        min_connections=min_connections,
        max_connections=max_connections,
        connection_timeout=connection_timeout,
        connection_lifetime=connection_lifetime,
        idle_timeout=idle_timeout
    )
    
    return connection_pool

def shutdown_connection_pool(timeout: int = 30):
    """
    Shutdown the global connection pool gracefully.
    
    Args:
        timeout: Maximum time to wait for connections to be released in seconds
    """
    global connection_pool
    if connection_pool is not None:
        connection_pool.shutdown(timeout)
        connection_pool = None

def get_connection_pool_stats() -> Dict[str, Any]:
    """
    Get statistics about the global connection pool.
    
    Returns:
        Dictionary with pool statistics
    """
    global connection_pool
    if connection_pool is None:
        return {"status": "not_initialized"}
    
    return connection_pool.get_stats()

def reset_connection_pool_stats():
    """Reset the global connection pool statistics."""
    global connection_pool
    if connection_pool is not None:
        connection_pool.reset_stats()

def configure_database_pooling(project_id: str) -> bool:
    """
    Configure connection pooling at the database level.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        True if successful, False otherwise
    """
    logger.info("Configuring database-level connection pooling...")
    
    # SQL for configuring connection pooling
    pooling_sql = """
    -- Configure connection pooling
    ALTER SYSTEM SET max_connections = '100';
    ALTER SYSTEM SET superuser_reserved_connections = '3';
    ALTER SYSTEM SET idle_in_transaction_session_timeout = '30000';  -- 30 seconds
    ALTER SYSTEM SET statement_timeout = '30000';  -- 30 seconds
    
    -- Configure connection pooler settings
    ALTER SYSTEM SET pgbouncer.max_client_conn = '1000';
    ALTER SYSTEM SET pgbouncer.default_pool_size = '50';
    ALTER SYSTEM SET pgbouncer.min_pool_size = '5';
    ALTER SYSTEM SET pgbouncer.reserve_pool_size = '10';
    ALTER SYSTEM SET pgbouncer.max_db_connections = '80';
    ALTER SYSTEM SET pgbouncer.max_user_connections = '50';
    
    -- Apply changes
    SELECT pg_reload_conf();
    """
    
    try:
        # Apply the configuration using a migration
        apply_migration(
            project_id=project_id,
            name="configure_connection_pooling",
            query=pooling_sql
        )
        logger.info("Successfully configured database-level connection pooling")
        return True
    except Exception as e:
        logger.error(f"Error configuring database-level connection pooling: {str(e)}")
        return False

def monitor_connection_pool(interval: int = 60, output_file: str = "connection_pool_stats.json"):
    """
    Monitor the connection pool and save statistics to a file.
    
    Args:
        interval: Monitoring interval in seconds
        output_file: File to save statistics to
    """
    def _monitor():
        while connection_pool and connection_pool.running:
            try:
                stats = get_connection_pool_stats()
                stats["timestamp"] = datetime.datetime.now().isoformat()
                
                # Save to file
                with open(output_file, "w") as f:
                    json.dump(stats, f, indent=2)
                
                logger.debug(f"Connection pool stats: {stats}")
                time.sleep(interval)
            except Exception as e:
                logger.error(f"Error monitoring connection pool: {str(e)}")
                time.sleep(interval)
    
    # Start monitoring in a background thread
    monitor_thread = threading.Thread(target=_monitor, daemon=True)
    monitor_thread.start()
    
    return monitor_thread

def run_load_test(project_id: str, num_queries: int = 100, 
                 concurrency: int = 5, query_type: str = "select") -> Dict[str, Any]:
    """
    Run a load test on the connection pool.
    
    Args:
        project_id: The Supabase project ID
        num_queries: Number of queries to execute
        concurrency: Number of concurrent queries
        query_type: Type of query to execute (select, insert, update)
        
    Returns:
        Dictionary with test results
    """
    logger.info(f"Running load test with {num_queries} queries, " +
               f"{concurrency} concurrent connections, query type: {query_type}")
    
    # Reset stats before the test
    reset_connection_pool_stats()
    
    # Prepare test queries
    if query_type == "select":
        query = "SELECT COUNT(*) FROM public.molecules;"
    elif query_type == "insert":
        query = """
        INSERT INTO public.test_connection_pool (name, created_at)
        VALUES ('test-' || floor(random() * 1000)::text, NOW())
        RETURNING id;
        """
    elif query_type == "update":
        query = """
        UPDATE public.test_connection_pool
        SET updated_at = NOW()
        WHERE id = (SELECT id FROM public.test_connection_pool ORDER BY RANDOM() LIMIT 1)
        RETURNING id;
        """
    else:
        raise ValueError(f"Invalid query type: {query_type}")
    
    # Create test table if needed
    if query_type in ("insert", "update"):
        try:
            create_table_query = """
            CREATE TABLE IF NOT EXISTS public.test_connection_pool (
                id SERIAL PRIMARY KEY,
                name TEXT NOT NULL,
                created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                updated_at TIMESTAMP WITH TIME ZONE
            );
            """
            execute_query(create_table_query)
            
            # Insert some initial data for update tests
            if query_type == "update":
                for i in range(10):
                    insert_query = f"""
                    INSERT INTO public.test_connection_pool (name, created_at)
                    VALUES ('initial-{i}', NOW())
                    """
                    execute_query(insert_query)
        except Exception as e:
            logger.error(f"Error creating test table: {str(e)}")
            return {"status": "error", "message": str(e)}
    
    # Run the load test
    results = []
    errors = []
    
    def worker(worker_id):
        for i in range(num_queries // concurrency):
            try:
                start_time = time.time()
                result = execute_query(query)
                end_time = time.time()
                
                results.append({
                    "worker_id": worker_id,
                    "query_id": i,
                    "duration": end_time - start_time,
                    "success": True
                })
            except Exception as e:
                errors.append({
                    "worker_id": worker_id,
                    "query_id": i,
                    "error": str(e)
                })
                logger.error(f"Worker {worker_id}, Query {i} failed: {str(e)}")
    
    # Start worker threads
    threads = []
    start_time = time.time()
    
    for i in range(concurrency):
        t = threading.Thread(target=worker, args=(i,))
        threads.append(t)
        t.start()
    
    # Wait for all threads to complete
    for t in threads:
        t.join()
    
    end_time = time.time()
    
    # Calculate statistics
    durations = [r["duration"] for r in results]
    
    stats = {
        "status": "success",
        "query_type": query_type,
        "num_queries": num_queries,
        "concurrency": concurrency,
        "total_time": end_time - start_time,
        "queries_per_second": len(results) / (end_time - start_time),
        "success_rate": len(results) / (len(results) + len(errors)) if (len(results) + len(errors)) > 0 else 0,
        "min_duration": min(durations) if durations else 0,
        "max_duration": max(durations) if durations else 0,
        "avg_duration": sum(durations) / len(durations) if durations else 0,
        "num_errors": len(errors),
        "pool_stats": get_connection_pool_stats()
    }
    
    logger.info(f"Load test completed: {stats['queries_per_second']:.2f} queries/sec, " +
               f"{stats['success_rate'] * 100:.2f}% success rate")
    
    return stats

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Implement connection pooling for CryoProtect database")
    parser.add_argument("--project-id", default="tsdlmynydfuypiugmkev",
                        help="Supabase project ID (default: tsdlmynydfuypiugmkev)")
    parser.add_argument("--min-connections", type=int, default=DEFAULT_MIN_CONNECTIONS,
                        help=f"Minimum number of connections (default: {DEFAULT_MIN_CONNECTIONS})")
    parser.add_argument("--max-connections", type=int, default=DEFAULT_MAX_CONNECTIONS,
                        help=f"Maximum number of connections (default: {DEFAULT_MAX_CONNECTIONS})")
    parser.add_argument("--connection-timeout", type=int, default=DEFAULT_CONNECTION_TIMEOUT,
                        help=f"Connection timeout in seconds (default: {DEFAULT_CONNECTION_TIMEOUT})")
    parser.add_argument("--test", action="store_true",
                        help="Run a load test after implementing connection pooling")
    parser.add_argument("--test-queries", type=int, default=100,
                        help="Number of queries to execute in the load test (default: 100)")
    parser.add_argument("--test-concurrency", type=int, default=5,
                        help="Number of concurrent connections in the load test (default: 5)")
    parser.add_argument("--monitor", action="store_true",
                        help="Enable connection pool monitoring")
    parser.add_argument("--monitor-interval", type=int, default=60,
                        help="Monitoring interval in seconds (default: 60)")
    return parser.parse_args()

def main():
    """Main function to implement connection pooling."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        
        # Configure database-level connection pooling
        logger.info(f"Implementing connection pooling for project {args.project_id}")
        db_pooling_success = configure_database_pooling(args.project_id)
        
        if not db_pooling_success:
            logger.warning("Database-level connection pooling configuration failed, " +
                          "continuing with application-level pooling only")
        
        # Initialize the connection pool
        pool = initialize_connection_pool(
            project_id=args.project_id,
            min_connections=args.min_connections,
            max_connections=args.max_connections,
            connection_timeout=args.connection_timeout
        )
        
        # Start monitoring if requested
        monitor_thread = None
        if args.monitor:
            logger.info(f"Starting connection pool monitoring (interval: {args.monitor_interval}s)")
            monitor_thread = monitor_connection_pool(interval=args.monitor_interval)
        
        # Run a load test if requested
        if args.test:
            logger.info("Running connection pool load test")
            test_results = run_load_test(
                project_id=args.project_id,
                num_queries=args.test_queries,
                concurrency=args.test_concurrency
            )
            
            # Save test results to file
            with open("connection_pool_test_results.json", "w") as f:
                json.dump(test_results, f, indent=2)
            
            logger.info(f"Load test results saved to connection_pool_test_results.json")
        
        # Keep the script running if monitoring is enabled
        if args.monitor:
            logger.info("Connection pooling implemented and monitoring active")
            logger.info("Press Ctrl+C to exit")
            
            # Set up signal handler for graceful shutdown
            def signal_handler(sig, frame):
                logger.info("Shutting down...")
                shutdown_connection_pool()
                sys.exit(0)
            
            signal.signal(signal.SIGINT, signal_handler)
            signal.signal(signal.SIGTERM, signal_handler)
            
            # Keep the main thread alive
            while True:
                time.sleep(1)
        else:
            # Display final stats
            stats = get_connection_pool_stats()
            logger.info(f"Connection pooling implemented successfully: {stats}")
            
            # Shutdown the pool
            shutdown_connection_pool()
            
            return 0
    
    except Exception as e:
        logger.error(f"Error implementing connection pooling: {str(e)}")
        return 1
    
    finally:
        # Ensure the connection pool is shut down
        if 'connection_pool' in globals() and connection_pool is not None:
            shutdown_connection_pool()

if __name__ == "__main__":
    sys.exit(main())
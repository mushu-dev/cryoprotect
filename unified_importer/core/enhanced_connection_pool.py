"""
Enhanced connection pool with health checks and dynamic sizing.

This module provides an improved connection pool implementation for the
unified molecular importer with advanced features like:
- Connection health checking
- Dynamic pool sizing based on workload
- Connection validation
- Detailed metrics collection
- Circuit breaker pattern for error handling
"""

import os
import time
import asyncio
import logging
import threading
import queue
import uuid
import json
from typing import Dict, List, Any, Optional, Tuple, Union, Set, Callable, AsyncIterator
from dataclasses import dataclass, field
import contextlib
from concurrent.futures import ThreadPoolExecutor
import random
from enum import Enum

# Optional imports based on available backends
try:
    import psycopg2
    import psycopg2.extras
    PSYCOPG2_AVAILABLE = True
except ImportError:
    PSYCOPG2_AVAILABLE = False

try:
    from supabase import create_client, Client
    SUPABASE_AVAILABLE = True
except ImportError:
    SUPABASE_AVAILABLE = False


class ConnectionState(Enum):
    """Enum for tracking connection state."""
    IDLE = "idle"
    IN_USE = "in_use"
    VALIDATING = "validating"
    UNHEALTHY = "unhealthy"
    CLOSED = "closed"


@dataclass
class EnhancedConnectionInfo:
    """Enhanced information about a database connection."""
    id: str = field(default_factory=lambda: str(uuid.uuid4()))
    host: str = "localhost"
    port: int = 5432
    database: str = "postgres"
    user: str = "postgres"
    password: str = ""
    connection_string: str = ""
    created_at: float = field(default_factory=time.time)
    last_used: float = field(default_factory=time.time)
    last_validated: float = field(default_factory=time.time)
    state: ConnectionState = ConnectionState.IDLE
    validation_failures: int = 0
    query_count: int = 0
    error_count: int = 0
    consecutive_errors: int = 0
    connection_type: str = "direct"
    extra: Dict[str, Any] = field(default_factory=dict)


class CircuitBreaker:
    """Circuit breaker pattern implementation for error handling."""
    
    def __init__(
        self,
        failure_threshold: int = 5,
        reset_timeout: float = 60.0,
        half_open_timeout: float = 5.0,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the circuit breaker.
        
        Args:
            failure_threshold: Number of failures before opening circuit
            reset_timeout: Time to wait before resetting to half-open state
            half_open_timeout: Time to wait in half-open state before closing
            logger: Logger instance
        """
        self.failure_threshold = failure_threshold
        self.reset_timeout = reset_timeout
        self.half_open_timeout = half_open_timeout
        self.logger = logger or logging.getLogger(__name__)
        
        # State variables
        self.failure_count = 0
        self.last_failure_time = 0
        self.state = "closed"  # closed, open, half-open
        self.lock = threading.RLock()
    
    def record_failure(self) -> None:
        """Record a failure and potentially open the circuit."""
        with self.lock:
            self.failure_count += 1
            self.last_failure_time = time.time()
            
            if self.state == "closed" and self.failure_count >= self.failure_threshold:
                self.state = "open"
                self.logger.warning(
                    f"Circuit opened after {self.failure_count} consecutive failures"
                )
            elif self.state == "half-open":
                self.state = "open"
                self.logger.warning("Circuit reopened after failure in half-open state")
    
    def record_success(self) -> None:
        """Record a success and potentially close the circuit."""
        with self.lock:
            if self.state == "half-open":
                self.state = "closed"
                self.failure_count = 0
                self.logger.info("Circuit closed after successful operation")
    
    def allow_request(self) -> bool:
        """
        Check if a request is allowed to proceed.
        
        Returns:
            True if the request is allowed, False otherwise
        """
        with self.lock:
            if self.state == "closed":
                return True
            
            if self.state == "open":
                time_since_failure = time.time() - self.last_failure_time
                
                if time_since_failure >= self.reset_timeout:
                    self.state = "half-open"
                    self.logger.info(
                        f"Circuit half-opened after {time_since_failure:.2f}s timeout"
                    )
                    return True
                
                return False
            
            # Half-open state - only allow occasional requests
            return True
    
    def get_state(self) -> str:
        """Get the current circuit breaker state."""
        with self.lock:
            return self.state


class EnhancedConnectionPool:
    """
    Enhanced connection pool with advanced features.
    
    This class extends the basic connection pool with:
    - Health checking of connections
    - Dynamic pool sizing based on workload
    - Connection validation
    - Detailed metrics
    - Circuit breaker pattern for handling persistent errors
    """
    
    def __init__(
        self,
        min_size: int = 2,
        max_size: int = 10,
        target_utilization: float = 0.7,
        timeout: float = 30.0,
        max_lifetime: float = 3600.0,
        validation_interval: float = 300.0,
        health_check_interval: float = 60.0,
        connection_type: str = "direct",
        connection_params: Optional[Dict[str, Any]] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the enhanced connection pool.
        
        Args:
            min_size: Minimum number of connections to maintain
            max_size: Maximum number of connections allowed
            target_utilization: Target pool utilization (0.0-1.0)
            timeout: Connection timeout in seconds
            max_lifetime: Maximum lifetime of a connection in seconds
            validation_interval: Interval between connection validations
            health_check_interval: Interval between pool health checks
            connection_type: Type of connection ("direct", "supabase")
            connection_params: Parameters for creating connections
            logger: Logger instance
        """
        self.min_size = min_size
        self.max_size = max_size
        self.target_utilization = target_utilization
        self.timeout = timeout
        self.max_lifetime = max_lifetime
        self.validation_interval = validation_interval
        self.health_check_interval = health_check_interval
        self.connection_type = connection_type
        self.connection_params = connection_params or {}
        self.logger = logger or logging.getLogger(__name__)
        
        # Pool data structures
        self.idle_connections = queue.Queue()
        self.all_connections: Dict[str, Tuple[Any, EnhancedConnectionInfo]] = {}
        self.active_count = 0
        self.lock = threading.RLock()
        
        # Async pool data structures
        self.async_idle_connections = asyncio.Queue()
        self.async_lock = asyncio.Lock()
        self.executor = ThreadPoolExecutor(max_workers=max_size)
        
        # Metrics and statistics
        self.stats = {
            "created": 0,
            "acquired": 0,
            "released": 0,
            "discarded": 0,
            "errors": 0,
            "timeouts": 0,
            "validations": 0,
            "validation_failures": 0,
            "health_checks": 0,
            "health_check_actions": 0,
            "peak_active": 0,
            "resize_events": 0,
            "last_error": None,
            "last_resize": 0,
            "last_health_check": 0
        }
        
        # Circuit breakers for error handling
        self.conn_circuit_breaker = CircuitBreaker(
            failure_threshold=5,
            reset_timeout=60.0,
            logger=self.logger
        )
        
        # Create validation query based on connection type
        if connection_type == "direct":
            self.validation_query = "SELECT 1"
        else:
            # Supabase validation is simpler
            self.validation_query = None
        
        # Initialize the pool
        self._initialize_pool()
        
        # Start background tasks
        self._start_background_tasks()
    
    def _initialize_pool(self) -> None:
        """Initialize the pool with minimum connections."""
        self.logger.info(
            f"Initializing enhanced connection pool with {self.min_size} connections"
        )
        
        initial_size = self.min_size
        created_count = 0
        
        for _ in range(initial_size):
            try:
                conn, conn_info = self._create_connection()
                
                # Add to pool data structures
                with self.lock:
                    self.all_connections[conn_info.id] = (conn, conn_info)
                    self.idle_connections.put(conn_info.id)
                    self.stats["created"] += 1
                    created_count += 1
                
                # Also initialize async pool
                self.async_idle_connections.put_nowait(conn_info.id)
            except Exception as e:
                self.logger.error(f"Error initializing connection: {str(e)}")
                self.stats["errors"] += 1
                self.conn_circuit_breaker.record_failure()
        
        self.logger.info(f"Successfully created {created_count}/{initial_size} initial connections")
    
    def _create_connection(self) -> Tuple[Any, EnhancedConnectionInfo]:
        """
        Create a new database connection.
        
        Returns:
            Tuple of (connection, connection_info)
        """
        start_time = time.time()
        
        if not self.conn_circuit_breaker.allow_request():
            self.logger.warning("Connection circuit breaker is open, using delayed retry")
            time.sleep(1.0)  # Brief delay to prevent hammering
            
            if not self.conn_circuit_breaker.allow_request():
                raise ConnectionError("Connection circuit breaker is open")
        
        with self.lock:
            self.active_count += 1
            
            if self.active_count > self.stats["peak_active"]:
                self.stats["peak_active"] = self.active_count
        
        try:
            if self.connection_type == "direct":
                # Direct PostgreSQL connection
                if not PSYCOPG2_AVAILABLE:
                    raise ImportError("psycopg2 is required for direct PostgreSQL connections")
                
                conn_params = {
                    'dbname': self.connection_params.get('database', 'postgres'),
                    'user': self.connection_params.get('user', 'postgres'),
                    'password': self.connection_params.get('password', ''),
                    'host': self.connection_params.get('host', 'localhost'),
                    'port': self.connection_params.get('port', 5432),
                    'connect_timeout': min(int(self.timeout), 60)  # Maximum allowed timeout
                }
                
                conn = psycopg2.connect(**conn_params)
                
                # Configure connection for optimal performance
                conn.set_session(autocommit=True)  # Non-transaction state by default
                
                # Create connection info
                conn_info = EnhancedConnectionInfo(
                    host=conn_params['host'],
                    port=conn_params['port'],
                    database=conn_params['dbname'],
                    user=conn_params['user'],
                    created_at=time.time(),
                    last_used=time.time(),
                    last_validated=time.time(),
                    state=ConnectionState.IDLE,
                    connection_type=self.connection_type
                )
                
                # Validate the connection immediately
                self._validate_connection(conn, conn_info)
                
                # Record success
                self.conn_circuit_breaker.record_success()
                
                # Log connection creation time
                create_time = time.time() - start_time
                self.logger.debug(f"Created PostgreSQL connection in {create_time:.3f}s")
                
                return conn, conn_info
            elif self.connection_type == "supabase":
                # Supabase connection
                if not SUPABASE_AVAILABLE:
                    raise ImportError("supabase-py is required for Supabase connections")
                
                supabase_url = self.connection_params.get('url') or os.environ.get('SUPABASE_URL')
                supabase_key = self.connection_params.get('key') or os.environ.get('SUPABASE_KEY')
                
                if not supabase_url or not supabase_key:
                    raise ValueError("Supabase URL and key are required")
                
                conn = create_client(supabase_url, supabase_key)
                
                # Create connection info
                conn_info = EnhancedConnectionInfo(
                    host=supabase_url,
                    created_at=time.time(),
                    last_used=time.time(),
                    last_validated=time.time(),
                    state=ConnectionState.IDLE,
                    connection_type=self.connection_type,
                    extra={"url": supabase_url}
                )
                
                # Record success
                self.conn_circuit_breaker.record_success()
                
                # Log connection creation time
                create_time = time.time() - start_time
                self.logger.debug(f"Created Supabase connection in {create_time:.3f}s")
                
                return conn, conn_info
            else:
                raise ValueError(f"Unsupported connection type: {self.connection_type}")
        except Exception as e:
            with self.lock:
                self.active_count -= 1
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)
            
            self.conn_circuit_breaker.record_failure()
            raise ConnectionError(f"Failed to create connection: {str(e)}")
    
    def _validate_connection(self, conn: Any, conn_info: EnhancedConnectionInfo) -> bool:
        """
        Validate a connection by executing a simple query.
        
        Args:
            conn: Database connection
            conn_info: Connection information
            
        Returns:
            True if connection is valid, False otherwise
        """
        if conn_info.state == ConnectionState.CLOSED:
            return False
        
        with self.lock:
            self.stats["validations"] += 1
        
        try:
            conn_info.state = ConnectionState.VALIDATING
            
            if self.connection_type == "direct":
                # For PostgreSQL, run a simple query
                with conn.cursor() as cursor:
                    cursor.execute(self.validation_query)
                    cursor.fetchone()
            else:
                # For Supabase, checking the connection state is enough
                pass
            
            # Connection is valid
            conn_info.state = ConnectionState.IDLE
            conn_info.last_validated = time.time()
            conn_info.validation_failures = 0
            conn_info.consecutive_errors = 0
            return True
        except Exception as e:
            # Connection validation failed
            with self.lock:
                self.stats["validation_failures"] += 1
                conn_info.validation_failures += 1
                conn_info.consecutive_errors += 1
                conn_info.state = ConnectionState.UNHEALTHY
            
            self.logger.warning(
                f"Connection validation failed for {conn_info.id}: {str(e)}"
            )
            return False
    
    def get_connection(self) -> Tuple[Any, str, bool]:
        """
        Get a connection from the pool.
        
        Returns:
            Tuple of (connection, connection_id, is_new)
        """
        conn = None
        conn_id = None
        is_new = False
        start_time = time.time()
        
        try:
            # Try to get a connection from the idle queue
            try:
                conn_id = self.idle_connections.get(block=False)
                conn, conn_info = self.all_connections[conn_id]
                
                # Check connection age and validity
                if time.time() - conn_info.created_at > self.max_lifetime:
                    # Connection is too old, discard it
                    self._close_connection(conn, conn_info)
                    with self.lock:
                        del self.all_connections[conn_id]
                        self.stats["discarded"] += 1
                    
                    conn = None
                    conn_id = None
                elif time.time() - conn_info.last_validated > self.validation_interval:
                    # Connection needs validation
                    if not self._validate_connection(conn, conn_info):
                        # Validation failed, discard and create a new one
                        self._close_connection(conn, conn_info)
                        with self.lock:
                            del self.all_connections[conn_id]
                            self.stats["discarded"] += 1
                        
                        conn = None
                        conn_id = None
                else:
                    # Connection is good to use
                    with self.lock:
                        conn_info.last_used = time.time()
                        conn_info.state = ConnectionState.IN_USE
                        self.stats["acquired"] += 1
            except queue.Empty:
                # No idle connections available
                pass
            
            # Create a new connection if needed
            if conn is None:
                with self.lock:
                    total_connections = len(self.all_connections)
                    
                    if total_connections < self.max_size:
                        # We can create a new connection
                        try:
                            conn, conn_info = self._create_connection()
                            conn_id = conn_info.id
                            is_new = True
                            
                            # Track the connection
                            self.all_connections[conn_id] = (conn, conn_info)
                            
                            # Mark as in use
                            conn_info.state = ConnectionState.IN_USE
                            self.stats["acquired"] += 1
                            
                            # Check if we need to increase the pool size
                            self._check_pool_size()
                        except Exception as e:
                            # Failed to create a new connection
                            self.logger.error(f"Failed to create new connection: {str(e)}")
                            
                            # Re-raise with more context
                            raise ConnectionError(f"Failed to get connection: {str(e)}")
                    else:
                        # We've reached max connections, wait for an idle one
                        try:
                            conn_id = self.idle_connections.get(block=True, timeout=self.timeout)
                            conn, conn_info = self.all_connections[conn_id]
                            
                            # Update usage information
                            conn_info.last_used = time.time()
                            conn_info.state = ConnectionState.IN_USE
                            
                            with self.lock:
                                self.stats["acquired"] += 1
                        except queue.Empty:
                            with self.lock:
                                self.stats["timeouts"] += 1
                                self.stats["errors"] += 1
                                self.stats["last_error"] = "Timed out waiting for connection"
                            
                            raise TimeoutError("Timed out waiting for database connection")
            
            # Log acquisition time
            acquisition_time = time.time() - start_time
            self.logger.debug(
                f"Acquired {'new' if is_new else 'existing'} connection in {acquisition_time:.3f}s"
            )
            
            return conn, conn_id, is_new
        except Exception as e:
            with self.lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)
            
            self.logger.error(f"Error getting connection: {str(e)}")
            raise
    
    def release_connection(self, conn: Any, conn_id: str) -> None:
        """
        Release a connection back to the pool.
        
        Args:
            conn: Database connection
            conn_id: Connection ID
        """
        try:
            # Get connection info
            if conn_id not in self.all_connections:
                self.logger.warning(f"Attempted to release unknown connection: {conn_id}")
                return
            
            _, conn_info = self.all_connections[conn_id]
            
            # Update usage information
            conn_info.last_used = time.time()
            
            # Check age and validity
            if time.time() - conn_info.created_at > self.max_lifetime:
                # Connection is too old, discard it
                self._close_connection(conn, conn_info)
                with self.lock:
                    del self.all_connections[conn_id]
                    self.stats["discarded"] += 1
            elif conn_info.consecutive_errors > 3:
                # Too many errors, discard and replace
                self.logger.warning(
                    f"Discarding connection with {conn_info.consecutive_errors} consecutive errors"
                )
                self._close_connection(conn, conn_info)
                with self.lock:
                    del self.all_connections[conn_id]
                    self.stats["discarded"] += 1
                
                # Create a replacement connection
                self._ensure_min_connections()
            else:
                # Return connection to the idle pool
                conn_info.state = ConnectionState.IDLE
                
                # Put back in the idle queue
                self.idle_connections.put(conn_id)
                
                with self.lock:
                    self.stats["released"] += 1
                
                # Check if we should resize the pool
                self._check_pool_size()
        except Exception as e:
            with self.lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)
            
            self.logger.error(f"Error releasing connection: {str(e)}")
    
    def _close_connection(self, conn: Any, conn_info: EnhancedConnectionInfo) -> None:
        """
        Close a database connection.
        
        Args:
            conn: Database connection
            conn_info: Connection information
        """
        try:
            conn_info.state = ConnectionState.CLOSED
            
            if self.connection_type == "direct":
                conn.close()
            # Supabase connections don't need explicit closing
            
            self.logger.debug(f"Closed connection {conn_info.id}")
        except Exception as e:
            self.logger.error(f"Error closing connection {conn_info.id}: {str(e)}")
    
    def _ensure_min_connections(self) -> None:
        """Ensure the pool has at least min_size connections."""
        with self.lock:
            current_size = len(self.all_connections)
            
            if current_size < self.min_size:
                connections_to_add = self.min_size - current_size
                
                self.logger.info(
                    f"Pool below minimum size ({current_size}/{self.min_size}), "
                    f"adding {connections_to_add} connections"
                )
                
                # Add connections
                for _ in range(connections_to_add):
                    try:
                        conn, conn_info = self._create_connection()
                        
                        # Add to pool data structures
                        self.all_connections[conn_info.id] = (conn, conn_info)
                        self.idle_connections.put(conn_info.id)
                        self.stats["created"] += 1
                        
                        # Also add to async pool
                        self.async_idle_connections.put_nowait(conn_info.id)
                    except Exception as e:
                        self.logger.error(f"Error creating connection: {str(e)}")
                        self.stats["errors"] += 1
    
    def _check_pool_size(self) -> None:
        """
        Check if the pool size should be adjusted based on utilization.
        
        This implements dynamic pool sizing based on current load.
        """
        with self.lock:
            # Only resize occasionally
            if time.time() - self.stats["last_resize"] < 60.0:
                return
            
            total_connections = len(self.all_connections)
            idle_count = self.idle_connections.qsize()
            active_count = total_connections - idle_count
            
            if total_connections == 0:
                return
            
            utilization = active_count / total_connections
            
            # Grow pool if utilization is high
            if utilization > self.target_utilization and total_connections < self.max_size:
                connections_to_add = min(2, self.max_size - total_connections)
                
                self.logger.info(
                    f"Pool utilization high ({utilization:.2f}), "
                    f"adding {connections_to_add} connections"
                )
                
                # Mark as resized
                self.stats["last_resize"] = time.time()
                self.stats["resize_events"] += 1
                
                # Add connections outside the lock
                for _ in range(connections_to_add):
                    try:
                        conn, conn_info = self._create_connection()
                        
                        # Add to pool data structures
                        with self.lock:
                            self.all_connections[conn_info.id] = (conn, conn_info)
                            self.idle_connections.put(conn_info.id)
                            self.stats["created"] += 1
                        
                        # Also add to async pool
                        self.async_idle_connections.put_nowait(conn_info.id)
                    except Exception as e:
                        self.logger.error(f"Error creating connection: {str(e)}")
                        self.stats["errors"] += 1
            
            # Shrink pool if utilization is low
            elif utilization < self.target_utilization / 2 and total_connections > self.min_size:
                idle_excess = idle_count - self.min_size
                
                if idle_excess > 2:
                    connections_to_remove = min(2, idle_excess)
                    
                    self.logger.info(
                        f"Pool utilization low ({utilization:.2f}), "
                        f"removing {connections_to_remove} idle connections"
                    )
                    
                    # Mark as resized
                    self.stats["last_resize"] = time.time()
                    self.stats["resize_events"] += 1
                    
                    # Remove connections
                    for _ in range(connections_to_remove):
                        try:
                            conn_id = self.idle_connections.get(block=False)
                            conn, conn_info = self.all_connections[conn_id]
                            
                            # Close and remove
                            self._close_connection(conn, conn_info)
                            del self.all_connections[conn_id]
                            self.stats["discarded"] += 1
                        except queue.Empty:
                            break
    
    def _run_health_check(self) -> None:
        """
        Run a health check on the connection pool.
        
        This verifies all connections and takes corrective actions as needed.
        """
        with self.lock:
            self.stats["health_checks"] += 1
            self.stats["last_health_check"] = time.time()
            
            # Snapshot current idle connections
            idle_connections = []
            while not self.idle_connections.empty():
                try:
                    conn_id = self.idle_connections.get(block=False)
                    idle_connections.append(conn_id)
                except queue.Empty:
                    break
            
            # Check each connection
            validated_count = 0
            failed_count = 0
            requeued_count = 0
            
            for conn_id in idle_connections:
                if conn_id in self.all_connections:
                    conn, conn_info = self.all_connections[conn_id]
                    
                    # Check if validation is needed
                    needs_validation = (
                        time.time() - conn_info.last_validated > self.validation_interval
                    )
                    
                    if needs_validation:
                        # Validate the connection
                        if self._validate_connection(conn, conn_info):
                            # Connection is valid, requeue it
                            self.idle_connections.put(conn_id)
                            validated_count += 1
                            requeued_count += 1
                        else:
                            # Connection is invalid, close and discard it
                            self._close_connection(conn, conn_info)
                            del self.all_connections[conn_id]
                            self.stats["discarded"] += 1
                            failed_count += 1
                            
                            # Record action
                            self.stats["health_check_actions"] += 1
                    else:
                        # No validation needed, requeue
                        self.idle_connections.put(conn_id)
                        requeued_count += 1
            
            # Ensure we maintain minimum pool size
            self._ensure_min_connections()
            
            # Log health check results
            self.logger.debug(
                f"Health check: validated {validated_count}, failed {failed_count}, "
                f"requeued {requeued_count}, pool size {len(self.all_connections)}"
            )
    
    def _start_background_tasks(self) -> None:
        """Start background tasks for pool maintenance."""
        # Health check thread
        self.health_check_thread = threading.Thread(
            target=self._health_check_loop,
            daemon=True,
            name="pool-health-check"
        )
        self.health_check_thread.start()
    
    def _health_check_loop(self) -> None:
        """Background loop for running periodic health checks."""
        while True:
            try:
                # Sleep first to allow initialization to complete
                time.sleep(self.health_check_interval)
                
                # Run the health check
                self._run_health_check()
            except Exception as e:
                self.logger.error(f"Error in health check loop: {str(e)}")
    
    def close_all(self) -> None:
        """Close all connections in the pool."""
        self.logger.info("Closing all connections in the pool")
        
        # Get all connections
        with self.lock:
            all_conn_ids = list(self.all_connections.keys())
        
        # Close each connection
        closed_count = 0
        error_count = 0
        
        for conn_id in all_conn_ids:
            if conn_id in self.all_connections:
                conn, conn_info = self.all_connections[conn_id]
                
                try:
                    self._close_connection(conn, conn_info)
                    closed_count += 1
                except Exception as e:
                    self.logger.error(f"Error closing connection {conn_id}: {str(e)}")
                    error_count += 1
                finally:
                    with self.lock:
                        if conn_id in self.all_connections:
                            del self.all_connections[conn_id]
        
        # Clear queues
        while not self.idle_connections.empty():
            try:
                self.idle_connections.get(block=False)
            except queue.Empty:
                break
        
        # Reset state
        with self.lock:
            self.active_count = 0
        
        # Shutdown executor
        self.executor.shutdown(wait=False)
        
        self.logger.info(f"Closed {closed_count} connections with {error_count} errors")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get detailed statistics about the connection pool.
        
        Returns:
            Dictionary of pool statistics
        """
        with self.lock:
            # Clone the stats to avoid race conditions
            stats = self.stats.copy()
            
            # Add current state
            stats.update({
                "pool_size": len(self.all_connections),
                "idle_connections": self.idle_connections.qsize(),
                "active_connections": self.active_count,
                "min_size": self.min_size,
                "max_size": self.max_size,
                "target_utilization": self.target_utilization,
                "circuit_breaker_state": self.conn_circuit_breaker.get_state(),
                "connection_age_stats": self._get_connection_age_stats(),
                "connection_usage_stats": self._get_connection_usage_stats(),
                "connection_error_stats": self._get_connection_error_stats()
            })
            
            return stats
    
    def _get_connection_age_stats(self) -> Dict[str, Any]:
        """Get statistics about connection ages."""
        now = time.time()
        
        # Extract connection info
        connection_infos = [info for _, info in self.all_connections.values()]
        
        if not connection_infos:
            return {
                "count": 0,
                "min_age": 0,
                "max_age": 0,
                "avg_age": 0
            }
        
        # Calculate ages
        ages = [now - info.created_at for info in connection_infos]
        
        return {
            "count": len(ages),
            "min_age": min(ages),
            "max_age": max(ages),
            "avg_age": sum(ages) / len(ages)
        }
    
    def _get_connection_usage_stats(self) -> Dict[str, Any]:
        """Get statistics about connection usage."""
        # Extract connection info
        connection_infos = [info for _, info in self.all_connections.values()]
        
        if not connection_infos:
            return {
                "count": 0,
                "min_queries": 0,
                "max_queries": 0,
                "avg_queries": 0,
                "total_queries": 0
            }
        
        # Calculate query counts
        query_counts = [info.query_count for info in connection_infos]
        
        return {
            "count": len(query_counts),
            "min_queries": min(query_counts) if query_counts else 0,
            "max_queries": max(query_counts) if query_counts else 0,
            "avg_queries": sum(query_counts) / len(query_counts) if query_counts else 0,
            "total_queries": sum(query_counts)
        }
    
    def _get_connection_error_stats(self) -> Dict[str, Any]:
        """Get statistics about connection errors."""
        # Extract connection info
        connection_infos = [info for _, info in self.all_connections.values()]
        
        if not connection_infos:
            return {
                "count": 0,
                "min_errors": 0,
                "max_errors": 0,
                "avg_errors": 0,
                "total_errors": 0
            }
        
        # Calculate error counts
        error_counts = [info.error_count for info in connection_infos]
        
        return {
            "count": len(error_counts),
            "min_errors": min(error_counts) if error_counts else 0,
            "max_errors": max(error_counts) if error_counts else 0,
            "avg_errors": sum(error_counts) / len(error_counts) if error_counts else 0,
            "total_errors": sum(error_counts)
        }
    
    # Async methods follow the same pattern as the sync methods
    
    async def get_connection_async(self) -> Tuple[Any, str, bool]:
        """
        Get a connection asynchronously from the pool.
        
        Returns:
            Tuple of (connection, connection_id, is_new)
        """
        conn = None
        conn_id = None
        is_new = False
        start_time = time.time()
        
        try:
            # Try to get a connection from the async pool
            try:
                conn_id = await self.async_idle_connections.get_nowait()
                conn, conn_info = self.all_connections[conn_id]
                
                # Check connection age and validity
                if time.time() - conn_info.created_at > self.max_lifetime:
                    # Connection is too old, discard it
                    await asyncio.to_thread(self._close_connection, conn, conn_info)
                    async with self.async_lock:
                        del self.all_connections[conn_id]
                        self.stats["discarded"] += 1
                    
                    conn = None
                    conn_id = None
                elif time.time() - conn_info.last_validated > self.validation_interval:
                    # Connection needs validation
                    is_valid = await asyncio.to_thread(self._validate_connection, conn, conn_info)
                    
                    if not is_valid:
                        # Validation failed, discard and create a new one
                        await asyncio.to_thread(self._close_connection, conn, conn_info)
                        async with self.async_lock:
                            del self.all_connections[conn_id]
                            self.stats["discarded"] += 1
                        
                        conn = None
                        conn_id = None
                else:
                    # Connection is good to use
                    async with self.async_lock:
                        conn_info.last_used = time.time()
                        conn_info.state = ConnectionState.IN_USE
                        self.stats["acquired"] += 1
            except asyncio.QueueEmpty:
                # No idle connections available
                pass
            
            # Create a new connection if needed
            if conn is None:
                async with self.async_lock:
                    total_connections = len(self.all_connections)
                    
                    if total_connections < self.max_size:
                        # We can create a new connection
                        try:
                            # Use thread to run blocking connection creation
                            conn, conn_info = await asyncio.to_thread(self._create_connection)
                            conn_id = conn_info.id
                            is_new = True
                            
                            # Track the connection
                            self.all_connections[conn_id] = (conn, conn_info)
                            
                            # Mark as in use
                            conn_info.state = ConnectionState.IN_USE
                            self.stats["acquired"] += 1
                            
                            # Check if we need to increase the pool size
                            await asyncio.to_thread(self._check_pool_size)
                        except Exception as e:
                            # Failed to create a new connection
                            self.logger.error(f"Failed to create new connection: {str(e)}")
                            
                            # Re-raise with more context
                            raise ConnectionError(f"Failed to get connection: {str(e)}")
                    else:
                        # We've reached max connections, wait for an idle one
                        try:
                            # Use asyncio.wait_for to implement timeout
                            conn_id = await asyncio.wait_for(
                                self.async_idle_connections.get(),
                                timeout=self.timeout
                            )
                            conn, conn_info = self.all_connections[conn_id]
                            
                            # Update usage information
                            conn_info.last_used = time.time()
                            conn_info.state = ConnectionState.IN_USE
                            
                            self.stats["acquired"] += 1
                        except asyncio.TimeoutError:
                            async with self.async_lock:
                                self.stats["timeouts"] += 1
                                self.stats["errors"] += 1
                                self.stats["last_error"] = "Timed out waiting for connection"
                            
                            raise TimeoutError("Timed out waiting for database connection")
            
            # Log acquisition time
            acquisition_time = time.time() - start_time
            self.logger.debug(
                f"Acquired {'new' if is_new else 'existing'} connection asynchronously "
                f"in {acquisition_time:.3f}s"
            )
            
            return conn, conn_id, is_new
        except Exception as e:
            async with self.async_lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)
            
            self.logger.error(f"Error getting connection asynchronously: {str(e)}")
            raise
    
    async def release_connection_async(self, conn: Any, conn_id: str) -> None:
        """
        Release a connection back to the async pool.
        
        Args:
            conn: Database connection
            conn_id: Connection ID
        """
        try:
            # Get connection info
            if conn_id not in self.all_connections:
                self.logger.warning(f"Attempted to release unknown connection: {conn_id}")
                return
            
            _, conn_info = self.all_connections[conn_id]
            
            # Update usage information
            conn_info.last_used = time.time()
            
            # Check age and validity
            if time.time() - conn_info.created_at > self.max_lifetime:
                # Connection is too old, discard it
                await asyncio.to_thread(self._close_connection, conn, conn_info)
                async with self.async_lock:
                    del self.all_connections[conn_id]
                    self.stats["discarded"] += 1
            elif conn_info.consecutive_errors > 3:
                # Too many errors, discard and replace
                self.logger.warning(
                    f"Discarding connection with {conn_info.consecutive_errors} consecutive errors"
                )
                await asyncio.to_thread(self._close_connection, conn, conn_info)
                async with self.async_lock:
                    del self.all_connections[conn_id]
                    self.stats["discarded"] += 1
                
                # Create a replacement connection asynchronously
                await asyncio.to_thread(self._ensure_min_connections)
            else:
                # Return connection to the idle pool
                conn_info.state = ConnectionState.IDLE
                
                # Put back in the idle queue
                await self.async_idle_connections.put(conn_id)
                self.idle_connections.put(conn_id)  # Also add to sync pool
                
                async with self.async_lock:
                    self.stats["released"] += 1
                
                # Check if we should resize the pool
                await asyncio.to_thread(self._check_pool_size)
        except Exception as e:
            async with self.async_lock:
                self.stats["errors"] += 1
                self.stats["last_error"] = str(e)
            
            self.logger.error(f"Error releasing connection asynchronously: {str(e)}")
    
    async def close_all_async(self) -> None:
        """Close all connections in the pool asynchronously."""
        self.logger.info("Closing all connections in the pool asynchronously")
        
        # Get all connections
        async with self.async_lock:
            all_conn_ids = list(self.all_connections.keys())
        
        # Close each connection
        closed_count = 0
        error_count = 0
        
        for conn_id in all_conn_ids:
            if conn_id in self.all_connections:
                conn, conn_info = self.all_connections[conn_id]
                
                try:
                    await asyncio.to_thread(self._close_connection, conn, conn_info)
                    closed_count += 1
                except Exception as e:
                    self.logger.error(f"Error closing connection asynchronously: {str(e)}")
                    error_count += 1
                finally:
                    async with self.async_lock:
                        if conn_id in self.all_connections:
                            del self.all_connections[conn_id]
        
        # Clear queues
        while not self.async_idle_connections.empty():
            try:
                await self.async_idle_connections.get_nowait()
            except asyncio.QueueEmpty:
                break
        
        while not self.idle_connections.empty():
            try:
                self.idle_connections.get(block=False)
            except queue.Empty:
                break
        
        # Reset state
        async with self.async_lock:
            self.active_count = 0
        
        # Shutdown executor
        self.executor.shutdown(wait=False)
        
        self.logger.info(f"Closed {closed_count} connections with {error_count} errors")
    
    async def get_stats_async(self) -> Dict[str, Any]:
        """
        Get detailed statistics about the connection pool asynchronously.
        
        Returns:
            Dictionary of pool statistics
        """
        async with self.async_lock:
            # Get stats from sync method
            stats = await asyncio.to_thread(self.get_stats)
            
            # Add async-specific stats
            stats.update({
                "async_idle_connections": self.async_idle_connections.qsize()
            })
            
            return stats
    
    async def run_health_check_async(self) -> Dict[str, Any]:
        """
        Run a health check asynchronously and return results.
        
        Returns:
            Dictionary with health check results
        """
        # Run health check in a thread
        await asyncio.to_thread(self._run_health_check)
        
        # Get stats
        stats = await self.get_stats_async()
        
        # Add health check specific info
        health_info = {
            "pool_size": stats["pool_size"],
            "idle_connections": stats["idle_connections"],
            "active_connections": stats["active_connections"],
            "health_checks": stats["health_checks"],
            "health_check_actions": stats["health_check_actions"],
            "last_health_check": stats["last_health_check"],
            "connection_age_stats": stats["connection_age_stats"],
            "connection_error_stats": stats["connection_error_stats"],
            "circuit_breaker_state": stats["circuit_breaker_state"]
        }
        
        return health_info
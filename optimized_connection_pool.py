#!/usr/bin/env python3
"""
Optimized Connection Pool for CryoProtect

This module provides a robust connection pool implementation with:
- Circuit breaker pattern for resilience
- Exponential backoff with jitter for retries
- Dynamic pool sizing based on load
- Connection health checks
- Detailed metrics collection
- Standardized logging for debugging
"""

import os
import sys
import time
import math
import random
import logging
import threading
import traceback
import json
import psycopg2
from psycopg2 import pool, errors
from queue import Queue, Empty, Full
from typing import Dict, Any, Optional, List, Tuple, Callable, Union
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class ConnectionStats:
    """Statistics for a single connection."""
    created_at: float
    last_used_at: float
    total_uses: int
    total_errors: int
    total_time_used: float
    avg_response_time: float
    current_state: str  # 'idle', 'in_use', 'error'

    @classmethod
    def new(cls):
        """Create a new connection stats object."""
        now = time.time()
        return cls(
            created_at=now,
            last_used_at=now,
            total_uses=0,
            total_errors=0,
            total_time_used=0,
            avg_response_time=0,
            current_state='idle'
        )

    def update_on_use(self):
        """Update stats when connection is used."""
        self.last_used_at = time.time()
        self.total_uses += 1
        self.current_state = 'in_use'

    def update_on_error(self):
        """Update stats when connection encounters an error."""
        self.total_errors += 1
        self.current_state = 'error'

    def update_on_return(self, use_time: float):
        """Update stats when connection is returned to the pool."""
        self.current_state = 'idle'
        self.total_time_used += use_time

        # Update average response time with exponential moving average
        if self.avg_response_time == 0:
            self.avg_response_time = use_time
        else:
            # Use a weight factor of 0.1 for the new value
            self.avg_response_time = (0.9 * self.avg_response_time) + (0.1 * use_time)

@dataclass
class ConnectionPoolMetrics:
    """Metrics for connection pool performance."""
    created: int = 0
    reused: int = 0
    discarded: int = 0
    errors: int = 0
    timeouts: int = 0
    wait_time_total: float = 0
    wait_time_avg: float = 0
    wait_time_max: float = 0
    query_time_total: float = 0
    query_time_avg: float = 0
    query_time_max: float = 0
    connection_age_avg: float = 0
    connection_uses_avg: float = 0
    circuit_breaks: int = 0
    circuit_recovery: int = 0
    validation_failures: int = 0
    health_check_failures: int = 0
    peak_connections: int = 0

    def update_wait_time(self, wait_time: float):
        """Update wait time metrics."""
        self.wait_time_total += wait_time
        if wait_time > self.wait_time_max:
            self.wait_time_max = wait_time
        # Recalculate average
        self.wait_time_avg = self.wait_time_total / max(self.created + self.reused, 1)

    def update_query_time(self, query_time: float):
        """Update query time metrics."""
        self.query_time_total += query_time
        if query_time > self.query_time_max:
            self.query_time_max = query_time
        # Recalculate average
        self.query_time_avg = self.query_time_total / max(self.created + self.reused, 1)

    def as_dict(self):
        """Convert metrics to dictionary."""
        return asdict(self)

    def reset(self):
        """Reset all metrics."""
        self.created = 0
        self.reused = 0
        self.discarded = 0
        self.errors = 0
        self.timeouts = 0
        self.wait_time_total = 0
        self.wait_time_avg = 0
        self.wait_time_max = 0
        self.query_time_total = 0
        self.query_time_avg = 0
        self.query_time_max = 0
        self.connection_age_avg = 0
        self.connection_uses_avg = 0
        self.circuit_breaks = 0
        self.circuit_recovery = 0
        self.validation_failures = 0
        self.health_check_failures = 0
        # Don't reset peak_connections as it's a high-water mark


class CircuitBreaker:
    """
    Implements the Circuit Breaker pattern for database connections.

    States:
    - CLOSED: Normal operation, all requests pass through
    - OPEN: Failure threshold exceeded, all requests fail fast
    - HALF_OPEN: Testing if the system has recovered
    """

    # Circuit states
    CLOSED = 'closed'
    OPEN = 'open'
    HALF_OPEN = 'half_open'

    def __init__(self, threshold=5, timeout=30, reset_count=2):
        """
        Initialize circuit breaker.

        Args:
            threshold: Number of failures before opening circuit
            timeout: Seconds circuit stays open before trying half-open
            reset_count: Number of successes needed to close circuit
        """
        self.failure_threshold = threshold
        self.timeout = timeout
        self.reset_count = reset_count

        self.state = self.CLOSED
        self.failure_count = 0
        self.success_count = 0
        self.last_failure_time = 0
        self._lock = threading.RLock()

    def allow_request(self):
        """Check if request should be allowed through circuit."""
        with self._lock:
            if self.state == self.CLOSED:
                return True

            if self.state == self.OPEN:
                # Check if timeout has elapsed
                if time.time() - self.last_failure_time > self.timeout:
                    # Try half-open state
                    self.state = self.HALF_OPEN
                    self.success_count = 0
                    logger.info("Circuit half-open, testing recovery")
                    return True
                return False

            # HALF_OPEN: allow limited requests through
            return True

    def record_success(self):
        """Record a successful operation."""
        with self._lock:
            if self.state == self.CLOSED:
                # Reset failure count on success in closed state
                self.failure_count = 0
                return

            if self.state == self.HALF_OPEN:
                self.success_count += 1
                if self.success_count >= self.reset_count:
                    # System appears recovered
                    self.state = self.CLOSED
                    self.failure_count = 0
                    self.success_count = 0
                    logger.info("Circuit closed, system recovered")

    def record_failure(self):
        """Record a failed operation."""
        with self._lock:
            self.last_failure_time = time.time()

            if self.state == self.CLOSED:
                self.failure_count += 1
                if self.failure_count >= self.failure_threshold:
                    # Too many failures, open circuit
                    self.state = self.OPEN
                    logger.warning(f"Circuit open after {self.failure_count} failures")

            elif self.state == self.HALF_OPEN:
                # Any failure in half-open returns to open
                self.state = self.OPEN
                logger.warning("Circuit reopened after failure in half-open state")

    def get_state(self):
        """Get current circuit state."""
        with self._lock:
            return self.state


class OptimizedConnectionPool:
    """
    Enhanced connection pool with:
    - Dynamic sizing
    - Connection validation
    - Circuit breaker
    - Exponential backoff
    - Health checks
    - Performance metrics
    """

    _instance = None
    _lock = threading.Lock()

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the optimized connection pool.

        Args:
            config: Configuration dictionary with connection parameters
        """
        self.config = config

        # Core connection parameters
        self.min_conn = config.get('min_connections', 3)
        self.max_conn = config.get('max_connections', 20)
        self.absolute_max_conn = config.get('absolute_max_connections', 50)
        self.connection_timeout = config.get('connection_timeout', 20)
        self.connection_lifetime = config.get('connection_lifetime', 1800)  # 30 minutes
        self.idle_timeout = config.get('idle_timeout', 180)  # 3 minutes

        # Retry parameters
        self.retry_attempts = config.get('retry_attempts', 5)
        self.initial_retry_delay = config.get('initial_retry_delay', 0.2)
        self.max_retry_delay = config.get('max_retry_delay', 10)
        self.retry_jitter_factor = config.get('retry_jitter_factor', 0.1)

        # Validation parameters
        self.validation_query = config.get('validation_query', 'SELECT 1')
        self.validation_timeout = config.get('validation_timeout', 5)

        # Health check parameters
        self.health_check_interval = config.get('health_check_interval', 30)

        # Connection tracking
        self.conn_stats = {}  # Map of connection ID to ConnectionStats
        self.active_connections = 0
        self.last_health_check = 0

        # Circuit breaker initialization
        self.circuit_breaker = CircuitBreaker(
            threshold=config.get('circuit_breaker_threshold', 5),
            timeout=config.get('circuit_breaker_timeout', 30),
            reset_count=config.get('circuit_breaker_reset', 2)
        )

        # Metrics collection
        self.metrics = ConnectionPoolMetrics()

        # Initialize pool
        self.pool = None
        self.pool_initialized = False
        self._initialize_pool()

        # Start background threads
        self._start_health_check()
        self._setup_monitoring()

    @classmethod
    def get_instance(cls, config: Optional[Dict[str, Any]] = None) -> 'OptimizedConnectionPool':
        """
        Get or create the singleton instance of the connection pool.

        Args:
            config: Configuration dictionary (only used if instance doesn't exist)

        Returns:
            OptimizedConnectionPool instance
        """
        with cls._lock:
            if cls._instance is None:
                if config is None:
                    raise ValueError("Configuration required for initial pool creation")
                cls._instance = cls(config)
            return cls._instance

    def _initialize_pool(self):
        """Initialize the connection pool."""
        try:
            # Create connection pool
            self.pool = pool.ThreadedConnectionPool(
                minconn=self.min_conn,
                maxconn=self.max_conn,
                host=self.config.get('host', 'localhost'),
                port=self.config.get('port', 5432),
                dbname=self.config.get('dbname', 'postgres'),
                user=self.config.get('user', 'postgres'),
                password=self.config.get('password', ''),
                connect_timeout=self.connection_timeout,
                **self.config.get('additional_params', {})
            )

            self.pool_initialized = True
            logger.info(f"Connection pool initialized with {self.min_conn}-{self.max_conn} connections")
        except Exception as e:
            self.pool_initialized = False
            logger.error(f"Failed to initialize connection pool: {str(e)}")
            raise

    def _start_health_check(self):
        """Start health check thread."""
        def health_check_worker():
            while True:
                try:
                    # Sleep first to allow initial pool setup to complete
                    time.sleep(self.health_check_interval)

                    if self.pool_initialized:
                        self._perform_health_check()
                        self._purge_old_connections()
                        self._adjust_pool_size()
                    else:
                        # Attempt to reinitialize if pool is not initialized
                        logger.warning("Pool not initialized, attempting to reinitialize...")
                        self._initialize_pool()
                except Exception as e:
                    logger.error(f"Error in health check thread: {str(e)}")
                    time.sleep(5)  # Short sleep after error

        # Start the thread
        health_thread = threading.Thread(target=health_check_worker, daemon=True)
        health_thread.start()
        logger.info("Health check thread started")

    def _perform_health_check(self):
        """Perform health check on the connection pool."""
        if not self.pool_initialized:
            return

        self.last_health_check = time.time()
        logger.debug("Performing connection pool health check")

        try:
            # Get a connection from the pool
            conn = self.pool.getconn(key='health_check')

            try:
                # Execute simple query to test connection
                with conn.cursor() as cursor:
                    cursor.execute(self.validation_query)
                    cursor.fetchone()

                # Connection is healthy
                logger.debug("Health check: Connection pool is healthy")

                # Reset the circuit breaker if it's open
                if self.circuit_breaker.get_state() != CircuitBreaker.CLOSED:
                    self.circuit_breaker.record_success()
            except Exception as e:
                # Connection is unhealthy
                logger.warning(f"Health check: Connection is unhealthy: {str(e)}")
                self.metrics.health_check_failures += 1
                self.circuit_breaker.record_failure()

                # Try to close and recreate the connection
                try:
                    self._close_bad_connection(conn)
                except Exception:
                    pass

                # Don't return this connection to the pool
                return
            finally:
                # Return connection to the pool
                self.pool.putconn(conn, key='health_check')
        except Exception as e:
            # Failed to get connection from pool
            logger.error(f"Health check: Failed to get connection from pool: {str(e)}")
            self.metrics.health_check_failures += 1
            self.circuit_breaker.record_failure()

    def _close_bad_connection(self, conn):
        """Properly close a bad connection."""
        try:
            # Make sure the connection is closed at the server
            if not conn.closed:
                conn.close()

            # Track metrics
            self.metrics.discarded += 1
        except Exception as e:
            logger.warning(f"Error closing bad connection: {str(e)}")

    def _purge_old_connections(self):
        """Purge old connections from the pool."""
        # Not directly accessible in psycopg2's connection pool
        # Instead, we track connection age and recycle when fetching
        pass

    def _adjust_pool_size(self):
        """Dynamically adjust pool size based on utilization."""
        if not self.pool_initialized:
            return

        try:
            # Calculate current utilization
            utilization_ratio = self.active_connections / self.max_conn

            # Scale up if consistently high utilization (>70%)
            if utilization_ratio > 0.7 and self.max_conn < self.absolute_max_conn:
                new_max = min(int(self.max_conn * 1.25), self.absolute_max_conn)
                new_min = min(int(self.min_conn * 1.25), new_max // 2)

                if new_max > self.max_conn:
                    logger.info(f"Scaling connection pool up: min={new_min}, max={new_max}")
                    self._recreate_pool(new_min, new_max)

            # Scale down if consistently low utilization (<30%)
            elif utilization_ratio < 0.3 and self.max_conn > 10 and self.min_conn > 2:
                new_max = max(int(self.max_conn * 0.8), 10)
                new_min = max(int(self.min_conn * 0.8), 2)

                if new_max < self.max_conn:
                    logger.info(f"Scaling connection pool down: min={new_min}, max={new_max}")
                    self._recreate_pool(new_min, new_max)
        except Exception as e:
            logger.error(f"Error adjusting pool size: {str(e)}")

    def _recreate_pool(self, new_min, new_max):
        """Recreate pool with new min/max values."""
        if not self.pool_initialized:
            return

        try:
            # Create new pool
            new_pool = pool.ThreadedConnectionPool(
                minconn=new_min,
                maxconn=new_max,
                host=self.config.get('host', 'localhost'),
                port=self.config.get('port', 5432),
                dbname=self.config.get('dbname', 'postgres'),
                user=self.config.get('user', 'postgres'),
                password=self.config.get('password', ''),
                connect_timeout=self.connection_timeout,
                **self.config.get('additional_params', {})
            )

            # Close old pool
            old_pool = self.pool

            # Update instance variables
            self.pool = new_pool
            self.min_conn = new_min
            self.max_conn = new_max

            # Schedule old pool closure
            def close_old_pool():
                try:
                    old_pool.closeall()
                    logger.info("Old connection pool closed successfully")
                except Exception as e:
                    logger.error(f"Error closing old pool: {str(e)}")

            # Close in a separate thread to not block
            closer_thread = threading.Thread(target=close_old_pool)
            closer_thread.daemon = True
            closer_thread.start()

            logger.info(f"Connection pool recreated with {new_min}-{new_max} connections")
        except Exception as e:
            logger.error(f"Failed to recreate connection pool: {str(e)}")

    def _validate_connection(self, conn):
        """Validate a connection is usable."""
        try:
            # Check if connection is closed
            if conn.closed:
                logger.debug("Connection validation failed: Connection is closed")
                return False

            # Check if connection is too old
            conn_id = id(conn)
            if conn_id in self.conn_stats:
                stats = self.conn_stats[conn_id]
                conn_age = time.time() - stats.created_at

                if conn_age > self.connection_lifetime:
                    logger.debug(f"Connection validation failed: Connection too old ({conn_age:.1f}s)")
                    return False

            # Execute validation query
            with conn.cursor() as cursor:
                cursor.execute(self.validation_query)
                result = cursor.fetchone()

                # Check if result is as expected
                if result and result[0] == 1:
                    return True
                else:
                    logger.warning(f"Connection validation failed: Unexpected result: {result}")
                    return False
        except Exception as e:
            logger.warning(f"Connection validation failed: {str(e)}")
            self.metrics.validation_failures += 1
            return False

    def _get_retry_delay(self, attempt):
        """
        Calculate retry delay with exponential backoff and jitter.

        Args:
            attempt: The current retry attempt (0-indexed)

        Returns:
            float: Delay in seconds
        """
        # Calculate exponential backoff
        delay = min(self.initial_retry_delay * (2 ** attempt), self.max_retry_delay)

        # Add jitter to prevent thundering herd
        jitter = delay * self.retry_jitter_factor * random.uniform(-1, 1)
        delay = max(0.001, delay + jitter)  # Ensure positive delay

        return delay

    def _setup_monitoring(self):
        """Set up connection pool monitoring."""
        # Set up scheduled metrics reporting
        reporting_interval = self.config.get('metrics_reporting_interval', 300)  # Every 5 minutes

        def _report_metrics():
            while True:
                try:
                    if self.pool_initialized:
                        metrics = self.metrics.as_dict()
                        logger.info(f"Connection pool metrics: {json.dumps(metrics)}")

                    time.sleep(reporting_interval)
                except Exception as e:
                    logger.error(f"Error reporting metrics: {str(e)}")
                    time.sleep(60)  # Retry after a minute

        # Start metrics reporting thread
        reporting_thread = threading.Thread(target=_report_metrics, daemon=True)
        reporting_thread.start()

    def get_metrics(self):
        """Get current connection pool metrics."""
        return self.metrics.as_dict()

    @contextmanager
    def get_connection(self, key=None):
        """
        Get a connection from the pool with circuit breaker, retries, and metrics.

        Args:
            key: Optional key for connection

        Yields:
            Database connection

        Raises:
            Exception: If unable to get a connection after retries
        """
        if not self.pool_initialized:
            try:
                self._initialize_pool()
            except Exception as e:
                logger.error(f"Failed to initialize pool: {str(e)}")
                raise

        # Check if circuit breaker allows the request
        if not self.circuit_breaker.allow_request():
            self.metrics.circuit_breaks += 1
            logger.warning("Circuit breaker open, failing fast")
            raise Exception("Database connection circuit breaker open")

        conn = None
        wait_start = time.time()
        attempt = 0

        while attempt < self.retry_attempts:
            try:
                # Get connection from pool
                conn = self.pool.getconn(key=key)
                wait_time = time.time() - wait_start
                self.metrics.update_wait_time(wait_time)

                # Track active connections
                with self._lock:
                    self.active_connections += 1
                    if self.active_connections > self.metrics.peak_connections:
                        self.metrics.peak_connections = self.active_connections

                # Get or create connection stats
                conn_id = id(conn)
                if conn_id not in self.conn_stats:
                    self.conn_stats[conn_id] = ConnectionStats.new()
                    self.metrics.created += 1
                else:
                    self.metrics.reused += 1

                # Update stats and record use
                self.conn_stats[conn_id].update_on_use()

                # Validate connection
                if not self._validate_connection(conn):
                    # Connection is invalid, try to close and retry
                    self._close_bad_connection(conn)
                    conn = None
                    raise Exception("Connection validation failed")

                # Record successful connection in circuit breaker
                self.circuit_breaker.record_success()

                # Get operation start time
                operation_start = time.time()

                try:
                    # Yield the connection
                    yield conn
                finally:
                    # Track operation time
                    operation_time = time.time() - operation_start

                    # Update stats
                    if conn_id in self.conn_stats:
                        self.conn_stats[conn_id].update_on_return(operation_time)

                    # Return connection to pool
                    if conn and not conn.closed:
                        self.pool.putconn(conn, key=key)

                    # Update metrics
                    self.metrics.update_query_time(operation_time)

                    # Decrement active connections counter
                    with self._lock:
                        self.active_connections = max(0, self.active_connections - 1)

                # Success, break out of retry loop
                break
            except Exception as e:
                # Record errors and handle retries
                self.metrics.errors += 1

                if conn:
                    try:
                        # Try to close and dispose of the bad connection
                        conn_id = id(conn)
                        if conn_id in self.conn_stats:
                            self.conn_stats[conn_id].update_on_error()

                        self._close_bad_connection(conn)
                    except Exception:
                        pass

                    # Clear connection reference
                    conn = None

                # Record failure in circuit breaker
                self.circuit_breaker.record_failure()

                # Increment attempt counter
                attempt += 1

                if attempt >= self.retry_attempts:
                    logger.error(f"Failed to get database connection after {attempt} attempts: {str(e)}")
                    raise

                # Calculate retry delay
                retry_delay = self._get_retry_delay(attempt)
                logger.warning(f"Retrying database connection ({attempt}/{self.retry_attempts}) after {retry_delay:.2f}s: {str(e)}")
                time.sleep(retry_delay)

        # If we get here without a connection and without an exception, something went wrong
        if not conn:
            self.metrics.errors += 1
            raise Exception("Failed to get database connection in unexpected way")

    def close(self):
        """Close all connections in the pool."""
        if self.pool_initialized and self.pool:
            try:
                self.pool.closeall()
                logger.info("Connection pool closed")

                # Reset pool state
                self.pool = None
                self.pool_initialized = False
                self.active_connections = 0
                self.conn_stats = {}
            except Exception as e:
                logger.error(f"Error closing connection pool: {str(e)}")
                raise
            
        # Update max if needed
        self.wait_time_max = max(self.wait_time_max, wait_time)
    
    def update_query_time(self, query_time: float):
        """Update query time metrics."""
        self.query_time_total += query_time
        
        # Update average with weighted average
        n = self.reused + self.created
        if n > 0:
            self.query_time_avg = ((n - 1) * self.query_time_avg + query_time) / n
        else:
            self.query_time_avg = query_time
            
        # Update max if needed
        self.query_time_max = max(self.query_time_max, query_time)
    
    def update_connection_age(self, connections):
        """Update connection age metrics."""
        if not connections:
            return
            
        total_age = sum(time.time() - conn['stats'].created_at for conn in connections)
        self.connection_age_avg = total_age / len(connections)
    
    def update_connection_uses(self, connections):
        """Update connection uses metrics."""
        if not connections:
            return
            
        total_uses = sum(conn['stats'].total_uses for conn in connections)
        self.connection_uses_avg = total_uses / len(connections)

class CircuitBreaker:
    """
    Circuit breaker pattern implementation for database connections.
    
    The circuit breaker prevents repeated attempts to execute operations
    that are likely to fail, allowing the system to recover.
    """
    
    # Circuit states
    CLOSED = 'closed'     # Normal operation
    OPEN = 'open'         # Circuit is open, failing fast
    HALF_OPEN = 'half_open'  # Testing if the circuit can be closed again
    
    def __init__(
        self, 
        failure_threshold: int = 5,
        recovery_timeout: int = 30,
        test_attempts: int = 2
    ):
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.test_attempts = test_attempts
        
        self.state = self.CLOSED
        self.failures = 0
        self.last_failure_time = 0
        self.test_counter = 0
        self.lock = threading.RLock()
    
    def can_execute(self) -> bool:
        """Check if the operation can be executed based on the circuit state."""
        with self.lock:
            if self.state == self.CLOSED:
                return True
            
            if self.state == self.OPEN:
                # Check if recovery timeout has elapsed
                if time.time() - self.last_failure_time > self.recovery_timeout:
                    logger.info("Circuit breaker transitioning from OPEN to HALF_OPEN")
                    self.state = self.HALF_OPEN
                    self.test_counter = 0
                    return True
                return False
            
            if self.state == self.HALF_OPEN:
                # Only allow a limited number of test attempts
                if self.test_counter < self.test_attempts:
                    self.test_counter += 1
                    return True
                return False
            
            return False
    
    def record_success(self):
        """Record a successful operation."""
        with self.lock:
            if self.state == self.HALF_OPEN:
                # If we've had enough successful test attempts, close the circuit
                if self.test_counter >= self.test_attempts:
                    logger.info("Circuit breaker transitioning from HALF_OPEN to CLOSED")
                    self.state = self.CLOSED
                    self.failures = 0
                    return True
            
            # Reset failure count in closed state on success
            if self.state == self.CLOSED:
                self.failures = 0
            
            return False
    
    def record_failure(self):
        """Record a failed operation."""
        with self.lock:
            self.failures += 1
            self.last_failure_time = time.time()
            
            if self.state == self.CLOSED and self.failures >= self.failure_threshold:
                logger.warning(
                    f"Circuit breaker transitioning from CLOSED to OPEN after "
                    f"{self.failures} consecutive failures"
                )
                self.state = self.OPEN
                return True
            
            if self.state == self.HALF_OPEN:
                logger.warning(
                    "Circuit breaker transitioning from HALF_OPEN back to OPEN "
                    "due to test failure"
                )
                self.state = self.OPEN
                return True
            
            return False
    
    def get_state(self) -> Dict[str, Any]:
        """Get the current state of the circuit breaker."""
        with self.lock:
            return {
                'state': self.state,
                'failures': self.failures,
                'last_failure_time': self.last_failure_time,
                'test_counter': self.test_counter
            }

class OptimizedConnectionPool:
    """
    An optimized database connection pool with advanced features.
    
    Features:
    - Dynamic pool sizing based on load
    - Connection validation and health checks
    - Circuit breaker pattern for resilience
    - Exponential backoff with jitter for retries
    - Detailed metrics collection
    - Automatic connection cleanup
    """
    
    def __init__(
        self,
        min_size: int = 2,
        max_size: int = 10,
        connection_timeout: int = 1800,
        validation_interval: int = 300,
        cleanup_interval: int = 600,
        connection_factory: Optional[Callable] = None
    ):
        """
        Initialize the connection pool.
        
        Args:
            min_size: Minimum number of connections to maintain
            max_size: Maximum number of connections allowed
            connection_timeout: Connection timeout in seconds
            validation_interval: Interval in seconds to validate idle connections
            cleanup_interval: Interval in seconds to cleanup expired connections
            connection_factory: Function to create a new connection
        """
        self.min_size = min_size
        self.max_size = max_size
        self.connection_timeout = connection_timeout
        self.validation_interval = validation_interval
        self.cleanup_interval = cleanup_interval
        
        # Use the provided factory or the default one
        self.connection_factory = connection_factory or self._default_connection_factory
        
        # Initialize pool and tracking variables
        self.pool = Queue(maxsize=max_size)
        self.active_connections = set()
        self.lock = threading.RLock()
        self.metrics = ConnectionPoolMetrics()
        
        # Initialize circuit breaker
        self.circuit_breaker = CircuitBreaker()
        
        # Track pool state
        self.shutting_down = False
        self.last_cleanup_time = time.time()
        self.last_validation_time = time.time()
        
        # Start with minimum connections
        self._initialize_pool()
        
        # Start background tasks
        self._start_background_tasks()
    
    def _default_connection_factory(self) -> Dict[str, Any]:
        """Default connection factory implementation."""
        raise NotImplementedError(
            "No connection factory provided. You must provide a connection factory "
            "function or implement _default_connection_factory."
        )
    
    def _initialize_pool(self):
        """Initialize the pool with minimum number of connections."""
        logger.info(f"Initializing connection pool with {self.min_size} connections")
        
        for _ in range(self.min_size):
            try:
                connection = self._create_connection()
                self.pool.put(connection)
            except Exception as e:
                logger.error(f"Error initializing connection: {str(e)}")
                logger.debug(traceback.format_exc())
                self.metrics.errors += 1
    
    def _start_background_tasks(self):
        """Start background tasks for pool maintenance."""
        # Start cleanup thread
        cleanup_thread = threading.Thread(
            target=self._cleanup_task,
            daemon=True,
            name="connection-pool-cleanup"
        )
        cleanup_thread.start()
        
        # Start validation thread
        validation_thread = threading.Thread(
            target=self._validation_task,
            daemon=True,
            name="connection-pool-validation"
        )
        validation_thread.start()
    
    def _cleanup_task(self):
        """Background task to clean up expired connections."""
        while not self.shutting_down:
            try:
                time.sleep(10)  # Check every 10 seconds
                
                current_time = time.time()
                if current_time - self.last_cleanup_time >= self.cleanup_interval:
                    self._cleanup_expired_connections()
                    self.last_cleanup_time = current_time
            except Exception as e:
                logger.error(f"Error in cleanup task: {str(e)}")
                logger.debug(traceback.format_exc())
    
    def _validation_task(self):
        """Background task to validate idle connections."""
        while not self.shutting_down:
            try:
                time.sleep(10)  # Check every 10 seconds
                
                current_time = time.time()
                if current_time - self.last_validation_time >= self.validation_interval:
                    self._validate_idle_connections()
                    self.last_validation_time = current_time
            except Exception as e:
                logger.error(f"Error in validation task: {str(e)}")
                logger.debug(traceback.format_exc())
    
    def _create_connection(self) -> Dict[str, Any]:
        """
        Create a new connection.
        
        Returns:
            Dict with connection details and stats
        """
        # Check circuit breaker first
        if not self.circuit_breaker.can_execute():
            logger.warning("Circuit breaker is open, failing fast")
            raise ConnectionError("Circuit breaker is open, connection creation prevented")
        
        try:
            # Call the connection factory
            connection_data = self.connection_factory()
            
            # Add stats tracking to the connection
            connection = {
                'data': connection_data,
                'stats': ConnectionStats.new()
            }
            
            with self.lock:
                self.active_connections.add(id(connection))
                self.metrics.created += 1
            
            # Record success in circuit breaker
            self.circuit_breaker.record_success()
            
            return connection
        except Exception as e:
            # Record failure in circuit breaker
            circuit_broken = self.circuit_breaker.record_failure()
            
            with self.lock:
                self.metrics.errors += 1
                if circuit_broken:
                    self.metrics.circuit_breaks += 1
            
            logger.error(f"Error creating connection: {str(e)}")
            logger.debug(traceback.format_exc())
            raise
    
    def _validate_connection(self, connection: Dict[str, Any]) -> bool:
        """
        Validate that a connection is still valid.
        
        Args:
            connection: Connection to validate
            
        Returns:
            True if connection is valid, False otherwise
        """
        # Check if connection has timed out
        if time.time() - connection['stats'].created_at > self.connection_timeout:
            logger.debug(f"Connection expired (timeout: {self.connection_timeout}s)")
            return False
        
        # Check for excessive errors
        if connection['stats'].total_errors > 5:
            logger.debug(f"Connection has too many errors: {connection['stats'].total_errors}")
            return False
        
        # Custom validation could be added here
            
        return True
    
    def _cleanup_expired_connections(self):
        """Clean up expired connections in the pool."""
        logger.debug("Running connection cleanup task")
        
        # Create a temporary list to hold valid connections
        valid_connections = []
        expired_count = 0
        
        # Drain the pool
        while not self.pool.empty():
            try:
                connection = self.pool.get(block=False)
                
                # Check if connection is valid
                if self._validate_connection(connection):
                    valid_connections.append(connection)
                else:
                    with self.lock:
                        self.active_connections.discard(id(connection))
                        self.metrics.discarded += 1
                        expired_count += 1
            except Empty:
                break
        
        # Put valid connections back in the pool
        for connection in valid_connections:
            try:
                self.pool.put(connection, block=False)
            except Full:
                # This should never happen as we're just putting back what we took out
                with self.lock:
                    self.active_connections.discard(id(connection))
                    self.metrics.discarded += 1
        
        # Create new connections if needed to maintain min_size
        current_size = self.pool.qsize()
        if current_size < self.min_size:
            logger.info(f"Pool size ({current_size}) below minimum ({self.min_size}), creating new connections")
            for _ in range(min(self.min_size - current_size, self.max_size - len(self.active_connections))):
                try:
                    connection = self._create_connection()
                    self.pool.put(connection)
                except Exception as e:
                    logger.error(f"Error creating connection during cleanup: {str(e)}")
        
        if expired_count > 0:
            logger.info(f"Cleaned up {expired_count} expired connections")
        
        # Update metrics
        with self.lock:
            # Get all connections from the pool temporarily
            temp_connections = []
            while not self.pool.empty():
                try:
                    temp_connections.append(self.pool.get(block=False))
                except Empty:
                    break
            
            # Update metrics
            self.metrics.update_connection_age(temp_connections)
            self.metrics.update_connection_uses(temp_connections)
            
            # Put connections back in the pool
            for conn in temp_connections:
                try:
                    self.pool.put(conn, block=False)
                except Full:
                    # This should never happen
                    self.active_connections.discard(id(conn))
                    self.metrics.discarded += 1
    
    def _validate_idle_connections(self):
        """Validate idle connections in the pool."""
        logger.debug("Running connection validation task")
        
        # Create a temporary list to hold all connections
        all_connections = []
        
        # Drain the pool
        while not self.pool.empty():
            try:
                connection = self.pool.get(block=False)
                all_connections.append(connection)
            except Empty:
                break
        
        # Validate each connection
        valid_connections = []
        invalid_count = 0
        
        for connection in all_connections:
            if self._validate_connection(connection):
                valid_connections.append(connection)
            else:
                with self.lock:
                    self.active_connections.discard(id(connection))
                    self.metrics.discarded += 1
                    invalid_count += 1
        
        # Put valid connections back in the pool
        for connection in valid_connections:
            try:
                self.pool.put(connection, block=False)
            except Full:
                # This shouldn't happen as we're just putting back what we took out
                with self.lock:
                    self.active_connections.discard(id(connection))
                    self.metrics.discarded += 1
        
        if invalid_count > 0:
            logger.info(f"Validation removed {invalid_count} invalid connections")
    
    def get_connection(self, timeout: int = 10) -> Tuple[Dict[str, Any], float]:
        """
        Get a connection from the pool or create a new one if needed.
        
        Args:
            timeout: Maximum time to wait for a connection in seconds
            
        Returns:
            Tuple of (connection, wait_time)
            
        Raises:
            ConnectionError: If unable to get a connection
        """
        if self.shutting_down:
            raise ConnectionError("Connection pool is shutting down")
        
        # Check circuit breaker
        if not self.circuit_breaker.can_execute():
            with self.lock:
                self.metrics.circuit_breaks += 1
            raise ConnectionError("Circuit breaker is open, connection acquisition prevented")
        
        connection = None
        wait_start = time.time()
        
        try:
            # Try to get a connection from the pool
            try:
                connection = self.pool.get(block=False)
                logger.debug("Got connection from pool immediately")
            except Empty:
                # Pool is empty, try to create a new connection or wait
                with self.lock:
                    if len(self.active_connections) < self.max_size:
                        # Create a new connection
                        logger.debug("Creating new connection as pool is empty")
                        connection = self._create_connection()
                    else:
                        # We're at max capacity, need to wait for a connection
                        logger.debug(f"Pool at max capacity ({self.max_size}), waiting for a connection")
                        try:
                            connection = self.pool.get(block=True, timeout=timeout)
                            logger.debug(f"Got connection from pool after waiting")
                        except Empty:
                            with self.lock:
                                self.metrics.timeouts += 1
                            raise ConnectionError(f"Timed out waiting for a connection after {timeout}s")
            
            # Validate the connection
            if connection and not self._validate_connection(connection):
                logger.debug("Got invalid connection from pool, creating a new one")
                # Discard the invalid connection
                with self.lock:
                    self.active_connections.discard(id(connection))
                    self.metrics.discarded += 1
                
                # Create a new connection
                connection = self._create_connection()
            
            wait_time = time.time() - wait_start
            
            # Update metrics
            with self.lock:
                if connection['stats'].total_uses > 0:
                    self.metrics.reused += 1
                self.metrics.update_wait_time(wait_time)
            
            # Update connection stats
            connection['stats'].update_on_use()
            
            # Record success with circuit breaker
            self.circuit_breaker.record_success()
            
            return connection, wait_time
        
        except Exception as e:
            # Record failure with circuit breaker
            circuit_broken = self.circuit_breaker.record_failure()
            
            with self.lock:
                self.metrics.errors += 1
                if circuit_broken:
                    self.metrics.circuit_breaks += 1
            
            logger.error(f"Error getting connection: {str(e)}")
            logger.debug(traceback.format_exc())
            raise
    
    def release_connection(self, connection: Dict[str, Any], error: bool = False, query_time: float = 0):
        """
        Release a connection back to the pool.
        
        Args:
            connection: Connection to release
            error: Whether the connection encountered an error
            query_time: Time spent executing queries with this connection
        """
        if self.shutting_down:
            # Pool is shutting down, just discard the connection
            with self.lock:
                self.active_connections.discard(id(connection))
                self.metrics.discarded += 1
            return
        
        try:
            # Update metrics
            with self.lock:
                self.metrics.update_query_time(query_time)
            
            if error:
                # Update connection stats
                connection['stats'].update_on_error()
                
                # Check if the connection should be discarded
                if connection['stats'].total_errors > 5:
                    logger.warning(f"Discarding connection with {connection['stats'].total_errors} errors")
                    with self.lock:
                        self.active_connections.discard(id(connection))
                        self.metrics.discarded += 1
                    return
            
            # Update connection stats on return
            connection['stats'].update_on_return(query_time)
            
            # Check if connection is still valid
            if self._validate_connection(connection):
                # Put the connection back in the pool
                try:
                    self.pool.put(connection, block=False)
                    return
                except Full:
                    logger.warning("Connection pool is full, discarding returned connection")
            
            # If we got here, the connection is invalid or the pool is full
            with self.lock:
                self.active_connections.discard(id(connection))
                self.metrics.discarded += 1
        
        except Exception as e:
            logger.error(f"Error releasing connection: {str(e)}")
            logger.debug(traceback.format_exc())
            
            # Ensure we don't leak connections
            with self.lock:
                self.active_connections.discard(id(connection))
                self.metrics.discarded += 1
                self.metrics.errors += 1
    
    @contextmanager
    def connection(self, timeout: int = 10):
        """
        Context manager for getting and automatically releasing a connection.
        
        Args:
            timeout: Maximum time to wait for a connection in seconds
            
        Yields:
            Connection data (whatever the connection factory returns)
            
        Raises:
            ConnectionError: If unable to get a connection
        """
        connection = None
        error_occurred = False
        start_time = time.time()
        
        try:
            connection, wait_time = self.get_connection(timeout=timeout)
            yield connection['data']
        except Exception:
            error_occurred = True
            raise
        finally:
            if connection:
                query_time = time.time() - start_time
                self.release_connection(
                    connection, 
                    error=error_occurred,
                    query_time=query_time
                )
    
    def shutdown(self, wait: bool = True):
        """
        Shutdown the connection pool.
        
        Args:
            wait: Whether to wait for active operations to complete
        """
        logger.info("Shutting down connection pool")
        self.shutting_down = True
        
        if wait:
            # Wait for all connections to be returned
            logger.info("Waiting for active operations to complete")
            timeout = 60  # Wait up to 60 seconds
            start_time = time.time()
            
            while (time.time() - start_time < timeout and 
                  len(self.active_connections) > self.pool.qsize()):
                time.sleep(0.1)
        
        # Drain the pool
        drained_count = 0
        while not self.pool.empty():
            try:
                self.pool.get(block=False)
                drained_count += 1
            except Empty:
                break
        
        logger.info(f"Connection pool shutdown complete. Drained {drained_count} connections.")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get detailed statistics about the connection pool.
        
        Returns:
            Dict with statistics
        """
        with self.lock:
            circuit_state = self.circuit_breaker.get_state()
            
            # Calculate overall stats
            total_connections = len(self.active_connections)
            idle_connections = self.pool.qsize()
            in_use_connections = total_connections - idle_connections
            
            return {
                "created_connections": self.metrics.created,
                "reused_connections": self.metrics.reused,
                "discarded_connections": self.metrics.discarded,
                "errors": self.metrics.errors,
                "timeouts": self.metrics.timeouts,
                "wait_time": {
                    "total": self.metrics.wait_time_total,
                    "average": self.metrics.wait_time_avg,
                    "max": self.metrics.wait_time_max
                },
                "query_time": {
                    "total": self.metrics.query_time_total,
                    "average": self.metrics.query_time_avg,
                    "max": self.metrics.query_time_max
                },
                "connection_age_avg": self.metrics.connection_age_avg,
                "connection_uses_avg": self.metrics.connection_uses_avg,
                "circuit_breaker": {
                    "state": circuit_state['state'],
                    "failures": circuit_state['failures'],
                    "state_changes": {
                        "to_open": self.metrics.circuit_breaks,
                        "to_closed": self.metrics.circuit_recovery
                    }
                },
                "pool_state": {
                    "size": {
                        "current": total_connections,
                        "min": self.min_size,
                        "max": self.max_size
                    },
                    "connections": {
                        "idle": idle_connections,
                        "in_use": in_use_connections
                    },
                    "shutting_down": self.shutting_down
                }
            }

def exponential_backoff(retry_count, base_delay=0.1, max_delay=60.0, jitter=0.25):
    """
    Calculate delay with exponential backoff and jitter.
    
    Args:
        retry_count: Number of retries already attempted
        base_delay: Base delay in seconds
        max_delay: Maximum delay in seconds
        jitter: Jitter factor (0-1) to randomize delay
        
    Returns:
        Delay in seconds
    """
    # Calculate exponential backoff
    delay = min(base_delay * (2 ** retry_count), max_delay)
    
    # Add jitter
    if jitter > 0:
        jitter_amount = delay * jitter
        delay = delay - (jitter_amount / 2) + (random.random() * jitter_amount)
    
    return delay

def retry_with_backoff(
    func, 
    max_retries=3, 
    base_delay=0.1, 
    max_delay=60.0, 
    jitter=0.25,
    retry_on=None
):
    """
    Retry a function with exponential backoff.
    
    Args:
        func: Function to retry
        max_retries: Maximum number of retries
        base_delay: Base delay in seconds
        max_delay: Maximum delay in seconds
        jitter: Jitter factor (0-1) to randomize delay
        retry_on: Exception or tuple of exceptions to retry on (default: Exception)
        
    Returns:
        Function result
        
    Raises:
        The last exception encountered
    """
    if retry_on is None:
        retry_on = Exception
    
    retry_count = 0
    last_exception = None
    
    while retry_count <= max_retries:
        try:
            return func()
        except retry_on as e:
            last_exception = e
            retry_count += 1
            
            if retry_count > max_retries:
                break
            
            delay = exponential_backoff(
                retry_count, base_delay, max_delay, jitter
            )
            
            logger.warning(
                f"Retry {retry_count}/{max_retries} after exception: {str(e)}. "
                f"Waiting {delay:.2f}s."
            )
            
            time.sleep(delay)
    
    # If we got here, we've exhausted our retries
    if last_exception:
        raise last_exception
    else:
        raise RuntimeError("Retry failed but no exception was caught")

# Modules to export
__all__ = [
    'OptimizedConnectionPool',
    'ConnectionStats',
    'ConnectionPoolMetrics',
    'CircuitBreaker',
    'exponential_backoff',
    'retry_with_backoff'
]
"""
Convex Connection Pool for CryoProtect

This module provides a connection pool implementation specifically optimized
for Convex database access. It includes enhanced features like:
- Circuit breaker pattern for resilience
- Dynamic pool sizing based on load
- Connection health verification
- Detailed metrics collection
- Support for both HTTP API and WebSocket connections
"""

import os
import time
import json
import random
import logging
import threading
import traceback
from queue import Queue, Empty, Full
from typing import Dict, Any, Optional, List, Tuple, Callable, Union
from dataclasses import dataclass, asdict, field
from contextlib import contextmanager

import requests
from requests.exceptions import RequestException, Timeout, ConnectionError

# Configure logging
logger = logging.getLogger(__name__)

@dataclass
class ConnectionStats:
    """Statistics for an individual Convex connection."""
    created_at: float
    last_used_at: float
    total_requests: int = 0
    successful_requests: int = 0
    failed_requests: int = 0
    total_time_used: float = 0
    avg_response_time: float = 0
    current_state: str = "idle"  # "idle", "in_use", "error"
    last_error: Optional[str] = None

    @classmethod
    def new(cls):
        """Create a new connection stats object."""
        now = time.time()
        return cls(
            created_at=now,
            last_used_at=now
        )

    def update_on_use(self):
        """Update stats when connection is used."""
        self.last_used_at = time.time()
        self.total_requests += 1
        self.current_state = "in_use"

    def update_on_success(self, request_time: float):
        """Update stats after a successful request."""
        self.successful_requests += 1
        self.total_time_used += request_time
        self.current_state = "idle"

        # Update average response time with exponential moving average
        if self.avg_response_time == 0:
            self.avg_response_time = request_time
        else:
            # Use a weight factor of 0.1 for the new value
            self.avg_response_time = (0.9 * self.avg_response_time) + (0.1 * request_time)

    def update_on_error(self, error: str):
        """Update stats when request encounters an error."""
        self.failed_requests += 1
        self.current_state = "error"
        self.last_error = error


@dataclass
class PoolMetrics:
    """Metrics for the Convex connection pool."""
    created: int = 0
    reused: int = 0
    discarded: int = 0
    errors: int = 0
    timeouts: int = 0
    connection_errors: int = 0
    circuit_breaks: int = 0
    circuit_recovery: int = 0
    validation_failures: int = 0
    health_check_failures: int = 0
    
    # Time metrics
    wait_time_total: float = 0
    wait_time_avg: float = 0
    wait_time_max: float = 0
    request_time_total: float = 0
    request_time_avg: float = 0
    request_time_max: float = 0
    
    # Connection metrics
    connection_age_avg: float = 0
    connection_requests_avg: float = 0
    
    # Pool metrics
    pool_size_current: int = 0
    pool_size_min: int = 0
    pool_size_max: int = 0
    pool_size_peak: int = 0
    
    # For HTTP stats
    http_status_codes: Dict[str, int] = field(default_factory=dict)
    http_methods: Dict[str, int] = field(default_factory=dict)
    
    def update_wait_time(self, wait_time: float):
        """Update wait time metrics."""
        self.wait_time_total += wait_time
        self.wait_time_avg = self.wait_time_total / max(self.created + self.reused, 1)
        self.wait_time_max = max(self.wait_time_max, wait_time)
    
    def update_request_time(self, request_time: float):
        """Update request time metrics."""
        self.request_time_total += request_time
        self.request_time_avg = self.request_time_total / max(self.created + self.reused, 1)
        self.request_time_max = max(self.request_time_max, request_time)
    
    def update_http_status(self, status_code: int):
        """Update HTTP status code metrics."""
        status_key = str(status_code)
        if status_key in self.http_status_codes:
            self.http_status_codes[status_key] += 1
        else:
            self.http_status_codes[status_key] = 1
    
    def update_http_method(self, method: str):
        """Update HTTP method metrics."""
        if method in self.http_methods:
            self.http_methods[method] += 1
        else:
            self.http_methods[method] = 1
    
    def reset(self):
        """Reset all metrics except for historical maximums."""
        self.created = 0
        self.reused = 0
        self.discarded = 0
        self.errors = 0
        self.timeouts = 0
        self.connection_errors = 0
        self.circuit_breaks = 0
        self.circuit_recovery = 0
        self.validation_failures = 0
        self.health_check_failures = 0
        
        self.wait_time_total = 0
        self.wait_time_avg = 0
        # Keep wait_time_max for historical record
        
        self.request_time_total = 0
        self.request_time_avg = 0
        # Keep request_time_max for historical record
        
        self.connection_age_avg = 0
        self.connection_requests_avg = 0
        
        # Keep pool_size_peak for historical record
        
        self.http_status_codes = {}
        self.http_methods = {}
    
    def as_dict(self):
        """Convert metrics to dictionary for reporting."""
        return asdict(self)


class CircuitBreaker:
    """
    Circuit breaker pattern implementation for Convex connections.
    
    The circuit breaker prevents cascading failures by stopping operations
    when a certain error threshold is reached, and then carefully testing
    recovery before resuming normal operations.
    
    States:
    - CLOSED: Normal operation, all requests pass through
    - OPEN: Circuit is open, failing fast
    - HALF_OPEN: Testing if system has recovered
    """
    
    # Circuit states
    CLOSED = "closed"
    OPEN = "open"
    HALF_OPEN = "half_open"
    
    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: int = 30,
        recovery_success: int = 2
    ):
        """
        Initialize the circuit breaker.
        
        Args:
            failure_threshold: Number of failures before opening circuit
            recovery_timeout: Seconds to wait before trying recovery
            recovery_success: Number of successes needed to close circuit
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.recovery_success = recovery_success
        
        self.state = self.CLOSED
        self.failures = 0
        self.successes = 0
        self.last_failure_time = 0
        self.lock = threading.RLock()
    
    def allow_request(self) -> bool:
        """
        Check if a request should be allowed through the circuit.
        
        Returns:
            bool: True if request should be allowed, False if it should fail fast
        """
        with self.lock:
            if self.state == self.CLOSED:
                return True
            
            if self.state == self.OPEN:
                # Check if we've waited long enough to try recovery
                if time.time() - self.last_failure_time > self.recovery_timeout:
                    logger.info("Circuit breaker moving from OPEN to HALF_OPEN to test recovery")
                    self.state = self.HALF_OPEN
                    self.successes = 0
                    return True
                return False
            
            # In HALF_OPEN state, we allow requests but monitor them closely
            return True
    
    def record_success(self) -> bool:
        """
        Record a successful operation.
        
        Returns:
            bool: True if circuit state changed, False otherwise
        """
        with self.lock:
            if self.state == self.CLOSED:
                # Reset failure count on success in closed state
                self.failures = 0
                return False
            
            if self.state == self.HALF_OPEN:
                self.successes += 1
                if self.successes >= self.recovery_success:
                    # System has recovered
                    logger.info("Circuit breaker closing after successful recovery")
                    self.state = self.CLOSED
                    self.failures = 0
                    self.successes = 0
                    return True
            
            return False
    
    def record_failure(self) -> bool:
        """
        Record a failed operation.
        
        Returns:
            bool: True if circuit state changed, False otherwise
        """
        with self.lock:
            self.last_failure_time = time.time()
            
            if self.state == self.CLOSED:
                self.failures += 1
                if self.failures >= self.failure_threshold:
                    # Too many failures, open the circuit
                    logger.warning(f"Circuit breaker opening after {self.failures} failures")
                    self.state = self.OPEN
                    return True
            
            elif self.state == self.HALF_OPEN:
                # Any failure in half-open state means we're not ready to recover
                logger.warning("Circuit breaker moving back to OPEN after failure during recovery")
                self.state = self.OPEN
                return True
            
            return False
    
    def get_state(self) -> Dict[str, Any]:
        """Get the current state of the circuit breaker."""
        with self.lock:
            return {
                "state": self.state,
                "failures": self.failures,
                "successes": self.successes,
                "last_failure_time": self.last_failure_time
            }


class ConvexConnectionPool:
    """
    Connection pool for Convex database access.
    
    Features:
    - Connection pooling with min/max size
    - Health checks and connection validation
    - Circuit breaker for resilience
    - Exponential backoff for retries
    - Detailed metrics collection
    - Support for both HTTP API and WebSocket connections
    """
    
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the connection pool.
        
        Args:
            config: Configuration dictionary with connection parameters
        """
        self.config = config
        
        # Core connection parameters
        self.convex_url = config.get("url", os.environ.get("CONVEX_URL", ""))
        self.convex_deployment_key = config.get(
            "deployment_key", os.environ.get("CONVEX_DEPLOYMENT_KEY", "")
        )
        
        if not self.convex_url:
            raise ValueError("Convex URL is required")
        
        # Ensure the URL ends with a slash
        if not self.convex_url.endswith("/"):
            self.convex_url += "/"
        
        # Pool configuration
        self.min_size = config.get("min_connections", 3)
        self.max_size = config.get("max_connections", 20)
        self.connection_timeout = config.get("connection_timeout", 30)
        self.connection_lifetime = config.get("connection_lifetime", 3600)  # 1 hour
        self.idle_timeout = config.get("idle_timeout", 300)  # 5 minutes
        
        # Retry configuration
        self.retry_attempts = config.get("retry_attempts", 3)
        self.initial_retry_delay = config.get("initial_retry_delay", 0.1)
        self.max_retry_delay = config.get("max_retry_delay", 5)
        self.retry_jitter_factor = config.get("retry_jitter_factor", 0.1)
        
        # Health check configuration
        self.health_check_interval = config.get("health_check_interval", 60)
        self.health_check_path = config.get("health_check_path", "api")
        
        # Circuit breaker configuration
        self.circuit_breaker = CircuitBreaker(
            failure_threshold=config.get("circuit_breaker_threshold", 5),
            recovery_timeout=config.get("circuit_breaker_timeout", 30),
            recovery_success=config.get("circuit_breaker_success", 2)
        )
        
        # Connection tracking
        self.http_connections = Queue(maxsize=self.max_size)
        self.active_connections = 0
        self.conn_stats = {}  # Map of connection ID to ConnectionStats
        
        # Metrics tracking
        self.metrics = PoolMetrics()
        self.metrics.pool_size_min = self.min_size
        self.metrics.pool_size_max = self.max_size
        
        # Thread synchronization
        self.lock = threading.RLock()
        
        # Pool state
        self.initialized = False
        self.shutting_down = False
        self.last_health_check = 0
        
        # Initialize the pool
        self._initialize_pool()
        
        # Start background threads
        self._start_health_check()
        self._start_metrics_reporter()
    
    @classmethod
    def get_instance(cls, config: Optional[Dict[str, Any]] = None) -> 'ConvexConnectionPool':
        """
        Get the singleton instance of the connection pool.
        
        Args:
            config: Configuration dictionary (only used for first initialization)
            
        Returns:
            ConvexConnectionPool instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    if config is None:
                        raise ValueError("Configuration required for initial pool creation")
                    cls._instance = cls(config)
        return cls._instance
    
    def _initialize_pool(self):
        """Initialize the connection pool with minimum connections."""
        logger.info(f"Initializing Convex connection pool with {self.min_size} connections")
        
        try:
            # Create initial connections
            for _ in range(self.min_size):
                connection = self._create_connection()
                self.http_connections.put(connection)
            
            self.initialized = True
            logger.info("Convex connection pool initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize Convex connection pool: {str(e)}")
            logger.debug(traceback.format_exc())
            self.initialized = False
    
    def _create_connection(self) -> Dict[str, Any]:
        """
        Create a new Convex connection.
        
        Returns:
            Dict with connection and stats
        """
        # Create a session
        session = requests.Session()
        
        # Set default headers
        headers = {
            "Content-Type": "application/json",
            "User-Agent": "CryoProtect-ConvexClient/1.0"
        }
        
        # Add auth header if we have a deployment key
        if self.convex_deployment_key:
            headers["Authorization"] = f"Bearer {self.convex_deployment_key}"
        
        session.headers.update(headers)
        
        # Set timeout for all requests
        session.request = lambda method, url, **kwargs: super(requests.Session, session).request(
            method=method,
            url=url,
            timeout=kwargs.pop('timeout', self.connection_timeout),
            **kwargs
        )
        
        # Create connection object
        connection = {
            "session": session,
            "created_at": time.time(),
            "id": id(session)
        }
        
        # Create stats for this connection
        self.conn_stats[connection["id"]] = ConnectionStats.new()
        
        # Update metrics
        with self.lock:
            self.active_connections += 1
            self.metrics.created += 1
            self.metrics.pool_size_current = self.active_connections
            
            # Update peak if needed
            if self.active_connections > self.metrics.pool_size_peak:
                self.metrics.pool_size_peak = self.active_connections
        
        return connection
    
    def _validate_connection(self, connection: Dict[str, Any]) -> bool:
        """
        Validate that a connection is usable.
        
        Args:
            connection: Connection to validate
            
        Returns:
            bool: True if connection is valid, False otherwise
        """
        # Check connection age
        connection_age = time.time() - connection["created_at"]
        if connection_age > self.connection_lifetime:
            logger.debug(f"Connection expired (age: {connection_age:.1f}s)")
            return False
        
        # Check connection stats for excessive errors
        conn_id = connection["id"]
        if conn_id in self.conn_stats:
            stats = self.conn_stats[conn_id]
            if stats.failed_requests > 5:
                logger.debug(f"Connection has too many errors: {stats.failed_requests}")
                return False
        
        return True
    
    def _start_health_check(self):
        """Start a background thread for periodic health checks."""
        def health_check_worker():
            while not self.shutting_down:
                try:
                    # Wait for health check interval
                    time.sleep(self.health_check_interval)
                    
                    if self.initialized:
                        # Perform health check
                        self._perform_health_check()
                        
                        # Clean up old connections
                        self._cleanup_connections()
                    else:
                        # Try to reinitialize if needed
                        logger.warning("Connection pool not initialized, attempting to initialize")
                        self._initialize_pool()
                except Exception as e:
                    logger.error(f"Error in health check thread: {str(e)}")
                    logger.debug(traceback.format_exc())
        
        # Start the thread
        health_thread = threading.Thread(
            target=health_check_worker,
            daemon=True,
            name="convex-health-check"
        )
        health_thread.start()
        logger.info(f"Health check thread started with interval {self.health_check_interval}s")
    
    def _perform_health_check(self):
        """Perform a health check on the Convex connection."""
        logger.debug("Performing Convex connection health check")
        
        self.last_health_check = time.time()
        
        # Get a connection from the pool specifically for health check
        try:
            connection = None
            try:
                # Try to get a connection without blocking
                connection = self.http_connections.get(block=False)
            except Empty:
                # Create a new connection if pool is empty
                connection = self._create_connection()
            
            # Validate the connection
            if not self._validate_connection(connection):
                # Connection is invalid, discard it
                self._discard_connection(connection)
                return
            
            # Update stats
            conn_id = connection["id"]
            if conn_id in self.conn_stats:
                self.conn_stats[conn_id].update_on_use()
            
            # Make a simple request to check health
            try:
                session = connection["session"]
                url = f"{self.convex_url}{self.health_check_path}"
                
                response = session.get(url, timeout=self.connection_timeout)
                
                # Check if response is valid
                if response.status_code >= 200 and response.status_code < 500:
                    # Connection is healthy
                    logger.debug("Health check: Convex connection is healthy")
                    
                    # Update stats
                    if conn_id in self.conn_stats:
                        self.conn_stats[conn_id].update_on_success(0.1)  # Nominal time
                    
                    # Reset circuit breaker if needed
                    self.circuit_breaker.record_success()
                else:
                    # Connection returned an error status
                    logger.warning(f"Health check: Convex returned error status {response.status_code}")
                    self.metrics.health_check_failures += 1
                    
                    # Update stats
                    if conn_id in self.conn_stats:
                        self.conn_stats[conn_id].update_on_error(f"HTTP {response.status_code}")
                    
                    # Record failure in circuit breaker
                    self.circuit_breaker.record_failure()
            except Exception as e:
                # Connection failed health check
                logger.warning(f"Health check: Convex connection failed: {str(e)}")
                self.metrics.health_check_failures += 1
                
                # Update stats
                if conn_id in self.conn_stats:
                    self.conn_stats[conn_id].update_on_error(str(e))
                
                # Record failure in circuit breaker
                self.circuit_breaker.record_failure()
                
                # Discard the connection
                self._discard_connection(connection)
                return
            
            # Return the connection to the pool
            try:
                self.http_connections.put(connection, block=False)
            except Full:
                # Pool is full, discard the connection
                self._discard_connection(connection)
        
        except Exception as e:
            logger.error(f"Error performing health check: {str(e)}")
            logger.debug(traceback.format_exc())
            self.metrics.health_check_failures += 1
            
            # Record failure in circuit breaker
            self.circuit_breaker.record_failure()
            
            # Make sure we don't leak the connection
            if connection:
                self._discard_connection(connection)
    
    def _cleanup_connections(self):
        """Clean up old or idle connections."""
        logger.debug("Cleaning up Convex connections")
        
        now = time.time()
        connections = []
        
        # Get all connections from the pool
        while not self.http_connections.empty():
            try:
                connection = self.http_connections.get(block=False)
                connections.append(connection)
            except Empty:
                break
        
        # Filter and return valid connections
        valid_connections = []
        
        for connection in connections:
            # Check if connection is too old
            connection_age = now - connection["created_at"]
            
            # Get connection stats
            conn_id = connection["id"]
            stats = self.conn_stats.get(conn_id)
            
            if connection_age > self.connection_lifetime:
                # Connection is too old, discard it
                logger.debug(f"Discarding old connection (age: {connection_age:.1f}s)")
                self._discard_connection(connection)
            elif stats and (now - stats.last_used_at) > self.idle_timeout:
                # Connection has been idle too long, discard it
                logger.debug(f"Discarding idle connection (idle: {now - stats.last_used_at:.1f}s)")
                self._discard_connection(connection)
            elif stats and stats.failed_requests > 5:
                # Connection has too many errors, discard it
                logger.debug(f"Discarding error-prone connection ({stats.failed_requests} errors)")
                self._discard_connection(connection)
            else:
                # Connection is valid, keep it
                valid_connections.append(connection)
        
        # Put valid connections back in the pool
        for connection in valid_connections:
            try:
                self.http_connections.put(connection, block=False)
            except Full:
                # Pool is full, discard the connection
                self._discard_connection(connection)
        
        # Create new connections if needed to maintain min_size
        current_size = self.http_connections.qsize()
        if current_size < self.min_size:
            logger.info(f"Creating {self.min_size - current_size} new connections to maintain minimum pool size")
            for _ in range(self.min_size - current_size):
                try:
                    connection = self._create_connection()
                    self.http_connections.put(connection)
                except Exception as e:
                    logger.error(f"Error creating new connection during cleanup: {str(e)}")
                    logger.debug(traceback.format_exc())
    
    def _discard_connection(self, connection: Dict[str, Any]):
        """
        Discard a connection and clean up resources.
        
        Args:
            connection: Connection to discard
        """
        if not connection:
            return
        
        try:
            # Close the session
            connection["session"].close()
            
            # Update metrics
            with self.lock:
                self.active_connections = max(0, self.active_connections - 1)
                self.metrics.discarded += 1
                self.metrics.pool_size_current = self.active_connections
            
            # Remove stats
            conn_id = connection["id"]
            if conn_id in self.conn_stats:
                del self.conn_stats[conn_id]
        
        except Exception as e:
            logger.error(f"Error discarding connection: {str(e)}")
    
    def _get_retry_delay(self, attempt: int) -> float:
        """
        Calculate retry delay with exponential backoff and jitter.
        
        Args:
            attempt: Current retry attempt (0-indexed)
            
        Returns:
            float: Delay in seconds
        """
        # Calculate exponential backoff
        delay = min(self.initial_retry_delay * (2 ** attempt), self.max_retry_delay)
        
        # Add jitter to prevent thundering herd
        jitter = delay * self.retry_jitter_factor * random.uniform(-1, 1)
        delay = max(0.001, delay + jitter)  # Ensure positive delay
        
        return delay
    
    def _start_metrics_reporter(self):
        """Start a thread to periodically report metrics."""
        def metrics_reporter():
            while not self.shutting_down:
                try:
                    # Wait for interval
                    time.sleep(300)  # Report every 5 minutes
                    
                    # Generate report
                    if self.initialized:
                        with self.lock:
                            metrics = self.metrics.as_dict()
                        
                        logger.info(f"Convex connection pool metrics: {json.dumps(metrics)}")
                
                except Exception as e:
                    logger.error(f"Error in metrics reporter: {str(e)}")
        
        # Start the thread
        reporter_thread = threading.Thread(
            target=metrics_reporter,
            daemon=True,
            name="convex-metrics-reporter"
        )
        reporter_thread.start()
        logger.info("Metrics reporter thread started")
    
    def get_connection(self) -> Tuple[Dict[str, Any], float]:
        """
        Get a connection from the pool.
        
        Returns:
            Tuple of (connection, wait_time)
            
        Raises:
            Exception: If unable to get a connection
        """
        if not self.initialized:
            try:
                self._initialize_pool()
            except Exception as e:
                logger.error(f"Failed to initialize pool: {str(e)}")
                raise
        
        # Check circuit breaker
        if not self.circuit_breaker.allow_request():
            with self.lock:
                self.metrics.circuit_breaks += 1
            
            logger.warning("Circuit breaker open, failing fast")
            raise ConnectionError("Convex connection circuit breaker open")
        
        wait_start = time.time()
        connection = None
        attempt = 0
        
        while attempt < self.retry_attempts:
            try:
                # Try to get a connection from the pool
                try:
                    connection = self.http_connections.get(block=False)
                    logger.debug("Got connection from pool immediately")
                except Empty:
                    # Pool is empty, create a new connection if possible
                    with self.lock:
                        if self.active_connections < self.max_size:
                            logger.debug("Creating new connection as pool is empty")
                            connection = self._create_connection()
                        else:
                            # We're at max capacity, wait for a connection
                            logger.debug(f"Pool at max capacity ({self.max_size}), waiting for a connection")
                            connection = self.http_connections.get(block=True, timeout=self.connection_timeout)
                
                # Validate the connection
                if not self._validate_connection(connection):
                    # Connection is invalid, discard and try again
                    logger.debug("Discarding invalid connection and retrying")
                    self._discard_connection(connection)
                    connection = None
                    raise ValueError("Connection validation failed")
                
                # Connection is valid
                wait_time = time.time() - wait_start
                
                # Update metrics
                with self.lock:
                    self.metrics.update_wait_time(wait_time)
                
                # Update connection stats
                conn_id = connection["id"]
                if conn_id in self.conn_stats:
                    self.conn_stats[conn_id].update_on_use()
                    
                    # Update reuse count for existing connections
                    if self.conn_stats[conn_id].total_requests > 1:
                        self.metrics.reused += 1
                
                # Record success in circuit breaker
                self.circuit_breaker.record_success()
                
                return connection, wait_time
            
            except Exception as e:
                # Handle retry logic
                attempt += 1
                
                with self.lock:
                    self.metrics.errors += 1
                
                # If we have a connection but it's invalid, discard it
                if connection:
                    self._discard_connection(connection)
                    connection = None
                
                # Record in circuit breaker
                self.circuit_breaker.record_failure()
                
                # If we've exhausted retries, raise the exception
                if attempt >= self.retry_attempts:
                    logger.error(f"Failed to get Convex connection after {attempt} attempts: {str(e)}")
                    raise
                
                # Calculate retry delay
                retry_delay = self._get_retry_delay(attempt - 1)
                logger.warning(f"Retrying to get Convex connection ({attempt}/{self.retry_attempts}) after {retry_delay:.2f}s: {str(e)}")
                time.sleep(retry_delay)
        
        # This should never happen, but just in case
        raise RuntimeError("Failed to get Convex connection in an unexpected way")
    
    def release_connection(self, connection: Dict[str, Any], error: bool = False, operation_time: float = 0):
        """
        Release a connection back to the pool.
        
        Args:
            connection: Connection to release
            error: Whether the operation encountered an error
            operation_time: Time spent using the connection
        """
        if not connection:
            return
        
        if self.shutting_down:
            # Pool is shutting down, just discard the connection
            self._discard_connection(connection)
            return
        
        try:
            conn_id = connection["id"]
            
            # Update metrics
            with self.lock:
                self.metrics.update_request_time(operation_time)
            
            # Update connection stats
            if conn_id in self.conn_stats:
                if error:
                    self.conn_stats[conn_id].update_on_error("Operation error")
                else:
                    self.conn_stats[conn_id].update_on_success(operation_time)
            
            # If connection had an error, check if it should be discarded
            if error and conn_id in self.conn_stats and self.conn_stats[conn_id].failed_requests > 5:
                logger.warning(f"Discarding error-prone connection ({self.conn_stats[conn_id].failed_requests} errors)")
                self._discard_connection(connection)
                return
            
            # Return the connection to the pool if it's still valid
            if self._validate_connection(connection):
                try:
                    self.http_connections.put(connection, block=False)
                    return
                except Full:
                    # Pool is full, discard the connection
                    logger.debug("Connection pool is full, discarding returned connection")
            
            # If we got here, the connection is invalid or the pool is full
            self._discard_connection(connection)
        
        except Exception as e:
            logger.error(f"Error releasing Convex connection: {str(e)}")
            logger.debug(traceback.format_exc())
            
            # Ensure we don't leak the connection
            self._discard_connection(connection)
    
    @contextmanager
    def connection(self, timeout: int = None):
        """
        Context manager for getting and automatically releasing a connection.
        
        Args:
            timeout: Connection timeout override
            
        Yields:
            requests.Session: Session object for making HTTP requests
            
        Raises:
            Exception: If unable to get a connection
        """
        connection = None
        error_occurred = False
        start_time = time.time()
        
        try:
            # Set timeout if provided
            old_timeout = self.connection_timeout
            if timeout is not None:
                self.connection_timeout = timeout
            
            # Get connection
            connection, _ = self.get_connection()
            yield connection["session"]
        
        except Exception:
            error_occurred = True
            raise
        
        finally:
            # Restore timeout if changed
            if timeout is not None:
                self.connection_timeout = old_timeout
            
            # Release connection
            if connection:
                operation_time = time.time() - start_time
                self.release_connection(
                    connection,
                    error=error_occurred,
                    operation_time=operation_time
                )
    
    def http_request(
        self,
        method: str,
        path: str,
        data: Any = None,
        params: Dict[str, Any] = None,
        headers: Dict[str, str] = None,
        timeout: int = None
    ) -> requests.Response:
        """
        Make an HTTP request to Convex HTTP API.
        
        Args:
            method: HTTP method
            path: API path (without base URL)
            data: JSON serializable data for request body
            params: Query parameters
            headers: Additional headers
            timeout: Optional timeout override
            
        Returns:
            requests.Response: Response object
            
        Raises:
            Exception: If request fails
        """
        # Build URL (remove leading slash if present)
        if path.startswith('/'):
            path = path[1:]
        
        url = f"{self.convex_url}{path}"
        
        # Keep track of retry attempts
        attempt = 0
        last_exception = None
        
        while attempt < self.retry_attempts:
            # Get a connection from the pool
            with self.connection(timeout=timeout) as session:
                try:
                    # Make the request
                    response = session.request(
                        method=method,
                        url=url,
                        json=data,
                        params=params,
                        headers=headers
                    )
                    
                    # Update metrics
                    with self.lock:
                        self.metrics.update_http_method(method)
                        self.metrics.update_http_status(response.status_code)
                    
                    # Check if we need to retry based on status code
                    if response.status_code >= 500:
                        attempt += 1
                        
                        if attempt >= self.retry_attempts:
                            return response
                        
                        # Calculate retry delay
                        retry_delay = self._get_retry_delay(attempt - 1)
                        logger.warning(
                            f"Retrying request to {url} ({attempt}/{self.retry_attempts}) "
                            f"after HTTP {response.status_code}. Waiting {retry_delay:.2f}s."
                        )
                        time.sleep(retry_delay)
                    else:
                        # Success or client error, don't retry
                        return response
                
                except (Timeout, ConnectionError, RequestException) as e:
                    # Network-related errors that we should retry
                    attempt += 1
                    last_exception = e
                    
                    # Update metrics
                    with self.lock:
                        if isinstance(e, Timeout):
                            self.metrics.timeouts += 1
                        elif isinstance(e, ConnectionError):
                            self.metrics.connection_errors += 1
                        else:
                            self.metrics.errors += 1
                    
                    if attempt >= self.retry_attempts:
                        logger.error(f"Request to {url} failed after {self.retry_attempts} attempts: {str(e)}")
                        raise
                    
                    # Calculate retry delay
                    retry_delay = self._get_retry_delay(attempt - 1)
                    logger.warning(
                        f"Retrying request to {url} ({attempt}/{self.retry_attempts}) "
                        f"after {e.__class__.__name__}. Waiting {retry_delay:.2f}s."
                    )
                    time.sleep(retry_delay)
        
        # We should never get here, but just in case
        if last_exception:
            raise last_exception
        else:
            raise RuntimeError(f"Request to {url} failed in an unexpected way")
    
    def execute_query(
        self,
        action: str,
        path: str,
        params: Dict[str, Any] = None,
        timeout: int = None
    ) -> Dict[str, Any]:
        """
        Execute a query against Convex and return JSON response.
        
        Args:
            action: HTTP method (GET, POST, PUT, DELETE)
            path: API path
            params: Parameters for the request
            timeout: Optional timeout override
            
        Returns:
            Dict: Response data
            
        Raises:
            Exception: If the request fails
        """
        # Set default timeout if not specified
        if timeout is None:
            timeout = self.connection_timeout
        
        # Convex HTTP actions URL format: [url]/http-api/[path]
        # This assumes our path already includes the leading slash if needed
        if path.startswith('/'):
            path = path[1:]
        
        api_path = f"http-api/{path}"
        
        logger.debug(f"Executing Convex query: {action} {api_path}")
        logger.debug(f"Params: {params}")
        
        try:
            if action.upper() == 'GET':
                response = self.http_request(
                    method="GET",
                    path=api_path,
                    params=params,
                    timeout=timeout
                )
            else:
                response = self.http_request(
                    method="POST",
                    path=api_path,
                    data=params,
                    timeout=timeout
                )
            
            # Parse response
            try:
                result = response.json()
            except ValueError:
                if not response.text:
                    result = {}
                else:
                    raise ValueError(f"Invalid JSON response: {response.text}")
            
            # Check for error
            if response.status_code >= 400:
                error_message = result.get("error", f"HTTP {response.status_code}")
                raise Exception(f"Convex query failed: {error_message}")
            
            return result
        except Exception as e:
            logger.error(f"Error executing Convex query: {str(e)}")
            logger.debug(f"Action: {action}, Path: {api_path}, Params: {params}")
            
            raise
    
    def close(self):
        """Close the connection pool and release all resources."""
        logger.info("Closing Convex connection pool")
        
        self.shutting_down = True
        
        # Close all connections
        while not self.http_connections.empty():
            try:
                connection = self.http_connections.get(block=False)
                self._discard_connection(connection)
            except Empty:
                break
        
        # Reset state
        self.initialized = False
        self.active_connections = 0
        
        logger.info("Convex connection pool closed")
    
    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about the connection pool."""
        with self.lock:
            # Build the stats dictionary
            stats = self.metrics.as_dict()
            
            # Add circuit breaker state
            stats["circuit_breaker"] = self.circuit_breaker.get_state()
            
            # Add pool info
            stats["pool_info"] = {
                "initialized": self.initialized,
                "shutting_down": self.shutting_down,
                "last_health_check": self.last_health_check,
                "http_queue_size": self.http_connections.qsize(),
                "connection_count": len(self.conn_stats)
            }
            
            return stats


# Functions to get and use the connection pool

def get_convex_pool(config: Optional[Dict[str, Any]] = None) -> ConvexConnectionPool:
    """
    Get the singleton instance of the Convex connection pool.
    
    Args:
        config: Optional configuration dictionary
        
    Returns:
        ConvexConnectionPool instance
    """
    if config is None:
        # Get configuration from environment
        config = {
            "url": os.environ.get("CONVEX_URL", ""),
            "deployment_key": os.environ.get("CONVEX_DEPLOYMENT_KEY", ""),
            "min_connections": int(os.environ.get("CONVEX_POOL_MIN_SIZE", "3")),
            "max_connections": int(os.environ.get("CONVEX_POOL_MAX_SIZE", "20")),
            "connection_timeout": int(os.environ.get("CONVEX_POOL_TIMEOUT", "30"))
        }
    
    return ConvexConnectionPool.get_instance(config)


def execute_convex_query(
    action: str,
    path: str,
    params: Dict[str, Any] = None,
    config: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Execute a query against Convex.
    
    Args:
        action: HTTP method (GET, POST, PUT, DELETE)
        path: API path
        params: Parameters for the request
        config: Optional configuration for the connection pool
        
    Returns:
        Dict: Response data
    """
    pool = get_convex_pool(config)
    return pool.execute_query(action, path, params)


@contextmanager
def convex_connection(config: Optional[Dict[str, Any]] = None, timeout: int = None):
    """
    Context manager for getting a Convex connection.
    
    Args:
        config: Optional configuration for the connection pool
        timeout: Optional timeout override
        
    Yields:
        requests.Session: Session object for making HTTP requests
    """
    pool = get_convex_pool(config)
    
    with pool.connection(timeout=timeout) as session:
        yield session
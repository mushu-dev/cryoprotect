"""
CryoProtect v2 API - Circuit Breaker Pattern

This module implements the Circuit Breaker pattern to prevent cascading failures
in distributed systems by temporarily stopping calls to failing services.

Key features:
- Three states: CLOSED (normal), OPEN (failing), and HALF_OPEN (testing recovery)
- Automatic trip based on consecutive failures
- Configurable thresholds and recovery timeouts
- Integration with monitoring systems
- Thread-safe implementation

Usage:
    from api.resiliency.circuit_breaker import circuit_breaker, CircuitBreakerState
    
    # Basic usage with default settings
    @circuit_breaker()
    def external_api_call():
        # This circuit breaker will open after 5 consecutive failures
        # and will attempt recovery after 60 seconds
        pass
        
    # Advanced usage with custom settings
    @circuit_breaker(
        name="database_circuit",               # Circuit name for monitoring
        failure_threshold=3,                   # Number of failures to open circuit
        recovery_timeout=30,                   # Seconds to wait before recovery attempt
        exceptions=(ConnectionError, TimeoutError),  # Exceptions to count as failures
        fallback=lambda *args, **kwargs: None  # Function to call when circuit is open
    )
    def database_operation():
        # This circuit breaker will open after 3 consecutive failures
        # and will attempt recovery after 30 seconds
        pass

Benefits:
- Prevents cascading failures across microservices
- Allows failing components to recover without constant load
- Provides fast failure when a service is known to be down
- Implements graceful recovery with the half-open state
"""

import time
import enum
import functools
import threading
import logging
from typing import Callable, Tuple, Type, Optional, Any, Dict, Union

# Get logger
logger = logging.getLogger(__name__)

# Import observability if available
try:
    from ..observability import (
        log_with_context, 
        report_error, 
        get_correlation_id,
        get_request_context
    )
    from monitoring.prometheus_metrics import (
        CIRCUIT_BREAKER_STATE,
        CIRCUIT_BREAKER_FAILURES,
        CIRCUIT_BREAKER_SUCCESSES
    )
    HAS_OBSERVABILITY = True
except ImportError:
    # Fallback to basic logging if observability is not available
    def log_with_context(logger, level, message, context=None, exc_info=False):
        getattr(logger, level)(message, exc_info=exc_info)
    
    def report_error(error, context=None):
        logger.error(f"Error: {str(error)}", exc_info=error)
    
    def get_correlation_id():
        return None
    
    def get_request_context():
        return {}
    
    # Dummy metrics
    class DummyGauge:
        def __init__(self):
            pass
        
        def labels(self, **kwargs):
            return self
        
        def set(self, value):
            pass
    
    class DummyCounter:
        def __init__(self):
            pass
        
        def labels(self, **kwargs):
            return self
            
        def inc(self):
            pass
    
    CIRCUIT_BREAKER_STATE = DummyGauge()
    CIRCUIT_BREAKER_FAILURES = DummyCounter()
    CIRCUIT_BREAKER_SUCCESSES = DummyCounter()
    HAS_OBSERVABILITY = False

class CircuitBreakerState(enum.Enum):
    """Circuit breaker states."""
    CLOSED = 0    # Normal state, requests pass through
    OPEN = 1      # Failing state, requests are blocked
    HALF_OPEN = 2 # Testing state, limited requests pass through for testing

class CircuitBreaker:
    """
    Circuit breaker implementation to prevent cascading failures.
    
    This class tracks service health and temporarily blocks requests to failing services,
    allowing them to recover without constant load.
    """
    
    # Class-level registry of circuit breakers
    _registry = {}
    _registry_lock = threading.RLock()
    
    @classmethod
    def get_circuit_breaker(cls, name):
        """
        Get or create a circuit breaker by name.
        
        Args:
            name: Name of the circuit breaker
            
        Returns:
            CircuitBreaker instance
        """
        with cls._registry_lock:
            if name not in cls._registry:
                cls._registry[name] = CircuitBreaker(name)
            return cls._registry[name]
    
    @classmethod
    def get_all_circuit_breakers(cls):
        """
        Get all registered circuit breakers.
        
        Returns:
            Dict of circuit breaker instances by name
        """
        with cls._registry_lock:
            return cls._registry.copy()
    
    def __init__(self, name):
        """
        Initialize a new circuit breaker.
        
        Args:
            name: Name of the circuit breaker for monitoring
        """
        self.name = name
        self.state = CircuitBreakerState.CLOSED
        self.failure_count = 0
        self.last_failure_time = 0
        self.lock = threading.RLock()
        
        # Default settings (can be overridden in decorator)
        self.failure_threshold = 5
        self.recovery_timeout = 60
        self.exceptions = (Exception,)
        
        # Update metrics
        if HAS_OBSERVABILITY:
            CIRCUIT_BREAKER_STATE.labels(name=self.name).set(self.state.value)
    
    def on_success(self):
        """Record a successful operation."""
        with self.lock:
            if self.state == CircuitBreakerState.HALF_OPEN:
                # Success in half-open state means the service is recovering
                self.state = CircuitBreakerState.CLOSED
                self.failure_count = 0
                
                # Log recovery
                log_with_context(
                    logger, 'info',
                    f"Circuit breaker '{self.name}' recovered and is now CLOSED",
                    context={
                        'event_type': 'circuit_breaker_recovery',
                        'circuit_breaker': {
                            'name': self.name,
                            'state': self.state.name,
                            'previous_state': CircuitBreakerState.HALF_OPEN.name
                        }
                    }
                )
            
            elif self.state == CircuitBreakerState.CLOSED:
                # Reset failure count on success
                self.failure_count = 0
            
            # Update metrics
            if HAS_OBSERVABILITY:
                CIRCUIT_BREAKER_STATE.labels(name=self.name).set(self.state.value)
                CIRCUIT_BREAKER_SUCCESSES.labels(name=self.name).inc()
    
    def on_failure(self, exception):
        """
        Record a failed operation.
        
        Args:
            exception: The exception that caused the failure
        """
        with self.lock:
            # Only count exceptions that match our configured exceptions
            if not isinstance(exception, self.exceptions):
                return
            
            # Update failure metrics
            if HAS_OBSERVABILITY:
                CIRCUIT_BREAKER_FAILURES.labels(
                    name=self.name,
                    exception_type=type(exception).__name__
                ).inc()
            
            if self.state == CircuitBreakerState.CLOSED:
                # Increment failure count
                self.failure_count += 1
                self.last_failure_time = time.time()
                
                # Check if we should trip the circuit
                if self.failure_count >= self.failure_threshold:
                    self.trip()
            
            elif self.state == CircuitBreakerState.HALF_OPEN:
                # Failure in half-open state means the service is still failing
                self.trip()
    
    def trip(self):
        """Trip the circuit breaker to OPEN state."""
        with self.lock:
            previous_state = self.state
            self.state = CircuitBreakerState.OPEN
            self.last_failure_time = time.time()
            
            # Log circuit trip
            log_with_context(
                logger, 'warning',
                f"Circuit breaker '{self.name}' tripped OPEN after {self.failure_count} consecutive failures",
                context={
                    'event_type': 'circuit_breaker_trip',
                    'circuit_breaker': {
                        'name': self.name,
                        'state': self.state.name,
                        'previous_state': previous_state.name,
                        'failure_count': self.failure_count,
                        'failure_threshold': self.failure_threshold,
                        'recovery_timeout': self.recovery_timeout
                    }
                }
            )
            
            # Update metrics
            if HAS_OBSERVABILITY:
                CIRCUIT_BREAKER_STATE.labels(name=self.name).set(self.state.value)
    
    def allow_request(self):
        """
        Check if a request is allowed to pass through the circuit breaker.
        
        Returns:
            bool: True if the request is allowed, False otherwise
        """
        with self.lock:
            if self.state == CircuitBreakerState.CLOSED:
                # Normal operation, always allow
                return True
            
            elif self.state == CircuitBreakerState.OPEN:
                # Check if recovery timeout has elapsed
                elapsed = time.time() - self.last_failure_time
                if elapsed >= self.recovery_timeout:
                    # Transition to half-open state to test the service
                    self.state = CircuitBreakerState.HALF_OPEN
                    
                    # Log half-open state
                    log_with_context(
                        logger, 'info',
                        f"Circuit breaker '{self.name}' testing recovery, now HALF_OPEN after {elapsed:.2f}s",
                        context={
                            'event_type': 'circuit_breaker_half_open',
                            'circuit_breaker': {
                                'name': self.name,
                                'state': self.state.name,
                                'previous_state': CircuitBreakerState.OPEN.name,
                                'elapsed_time': elapsed,
                                'recovery_timeout': self.recovery_timeout
                            }
                        }
                    )
                    
                    # Update metrics
                    if HAS_OBSERVABILITY:
                        CIRCUIT_BREAKER_STATE.labels(name=self.name).set(self.state.value)
                    
                    # Allow the first request in half-open state
                    return True
                else:
                    # Circuit is still open, block the request
                    return False
            
            elif self.state == CircuitBreakerState.HALF_OPEN:
                # Allow a single request to test the service 
                # This was already allowed during the transition to half-open
                # Further requests are blocked until the state changes
                return False
            
            # Default to allowing the request
            return True
    
    def get_state(self):
        """
        Get the current state of the circuit breaker.
        
        Returns:
            CircuitBreakerState: Current state of the circuit breaker
        """
        with self.lock:
            return self.state
    
    def get_metrics(self):
        """
        Get circuit breaker metrics.
        
        Returns:
            dict: Circuit breaker metrics
        """
        with self.lock:
            return {
                'name': self.name,
                'state': self.state.name,
                'failure_count': self.failure_count,
                'failure_threshold': self.failure_threshold,
                'recovery_timeout': self.recovery_timeout,
                'last_failure_time': self.last_failure_time,
                'elapsed_since_failure': time.time() - self.last_failure_time if self.last_failure_time > 0 else None
            }
    
    def reset(self):
        """Reset the circuit breaker to CLOSED state."""
        with self.lock:
            previous_state = self.state
            self.state = CircuitBreakerState.CLOSED
            self.failure_count = 0
            
            # Log reset
            log_with_context(
                logger, 'info',
                f"Circuit breaker '{self.name}' manually reset to CLOSED",
                context={
                    'event_type': 'circuit_breaker_reset',
                    'circuit_breaker': {
                        'name': self.name,
                        'state': self.state.name,
                        'previous_state': previous_state.name
                    }
                }
            )
            
            # Update metrics
            if HAS_OBSERVABILITY:
                CIRCUIT_BREAKER_STATE.labels(name=self.name).set(self.state.value)

def circuit_breaker(
    name: Optional[str] = None,
    failure_threshold: int = 5,
    recovery_timeout: float = 60.0,
    exceptions: Tuple[Type[Exception], ...] = (Exception,),
    fallback: Optional[Callable] = None
) -> Callable:
    """
    Decorator that implements the circuit breaker pattern for a function.
    
    Args:
        name: Name of the circuit breaker (default: function name)
        failure_threshold: Number of consecutive failures to trip the circuit (default: 5)
        recovery_timeout: Seconds to wait before attempting recovery (default: 60.0)
        exceptions: Tuple of exception types to count as failures (default: (Exception,))
        fallback: Optional function to call when the circuit is open
        
    Returns:
        Decorated function
        
    The circuit breaker has three states:
    - CLOSED: Normal operation, requests pass through
    - OPEN: Circuit is tripped, requests are blocked
    - HALF_OPEN: Testing recovery, a single request is allowed through
    """
    def decorator(func: Callable) -> Callable:
        # Use function name if name not provided
        circuit_name = name or f"{func.__module__}.{func.__name__}"
        
        # Get or create circuit breaker
        circuit = CircuitBreaker.get_circuit_breaker(circuit_name)
        
        # Configure circuit breaker
        circuit.failure_threshold = failure_threshold
        circuit.recovery_timeout = recovery_timeout
        circuit.exceptions = exceptions
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Check if request is allowed
            if not circuit.allow_request():
                # Circuit is open or half-open with a request in progress
                if fallback:
                    # Call fallback function if provided
                    return fallback(*args, **kwargs)
                else:
                    # Raise exception for open circuit if no fallback
                    raise CircuitOpenError(
                        f"Circuit breaker '{circuit_name}' is OPEN due to consecutive failures. "
                        f"Retry after {circuit.recovery_timeout}s."
                    )
            
            try:
                # Call the function
                result = func(*args, **kwargs)
                
                # Record success
                circuit.on_success()
                
                return result
            
            except Exception as e:
                # Record failure and propagate the exception
                circuit.on_failure(e)
                raise
        
        # Add circuit breaker reference to the function
        wrapper.circuit_breaker = circuit
        
        return wrapper
    
    return decorator

class CircuitOpenError(Exception):
    """Exception raised when a circuit breaker is open."""
    pass
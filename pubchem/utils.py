"""
Utility functions for the PubChem client.

This module provides helper functions for the PubChem client, including
circuit breaker implementation and request handling utilities.
"""

import time
import logging
import functools
import random
from typing import Callable, Any, Dict, Optional, TypeVar, cast
from enum import Enum

logger = logging.getLogger(__name__)

# Type variable for generic function
T = TypeVar('T')

class CircuitState(Enum):
    """Circuit breaker states."""
    CLOSED = 'closed'  # Normal operation, requests allowed
    OPEN = 'open'      # Circuit is open, requests are blocked
    HALF_OPEN = 'half_open'  # Testing if service is back online


class CircuitBreaker:
    """
    Circuit breaker implementation to prevent repeated failures.
    
    Features:
    - Three states: CLOSED, OPEN, HALF_OPEN
    - Configurable failure threshold, timeout, and retry timeout
    - Exponential backoff for retries
    """
    
    def __init__(
        self,
        failure_threshold: int = 5,
        recovery_timeout: int = 30,
        expected_exceptions: tuple = (Exception,),
        name: str = "default"
    ):
        """
        Initialize the circuit breaker.
        
        Args:
            failure_threshold: Number of failures before opening the circuit
            recovery_timeout: Time in seconds to wait before trying again
            expected_exceptions: Tuple of exceptions that count as failures
            name: Name of this circuit breaker for logging
        """
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.expected_exceptions = expected_exceptions
        self.name = name
        
        self.state = CircuitState.CLOSED
        self.failure_count = 0
        self.last_failure_time = 0
        
        # Stats
        self.stats = {
            "success_count": 0,
            "failure_count": 0,
            "rejected_count": 0,
            "state_changes": []
        }
    
    def __call__(self, func: Callable[..., T]) -> Callable[..., T]:
        """
        Decorator to apply circuit breaker to a function.
        
        Args:
            func: Function to wrap with circuit breaker
            
        Returns:
            Wrapped function
        """
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            return self.call(func, *args, **kwargs)
        return wrapper
    
    def call(self, func: Callable[..., T], *args: Any, **kwargs: Any) -> T:
        """
        Call the function with circuit breaker protection.
        
        Args:
            func: Function to call
            *args: Positional arguments for the function
            **kwargs: Keyword arguments for the function
            
        Returns:
            Result of the function call
            
        Raises:
            Exception: If the circuit is open or the function call fails
        """
        if self.state == CircuitState.OPEN:
            if time.time() - self.last_failure_time >= self.recovery_timeout:
                # Try to recover
                self._change_state(CircuitState.HALF_OPEN)
            else:
                # Circuit is open, reject the request
                self.stats["rejected_count"] += 1
                recovery_time = self.recovery_timeout - (time.time() - self.last_failure_time)
                logger.warning(
                    f"Circuit {self.name} is OPEN. Request rejected. "
                    f"Retry after {recovery_time:.1f} seconds."
                )
                raise CircuitBreakerError(
                    f"Circuit {self.name} is OPEN. Request rejected. "
                    f"Retry after {recovery_time:.1f} seconds."
                )
        
        try:
            result = func(*args, **kwargs)
            self._on_success()
            return result
        except self.expected_exceptions as e:
            self._on_failure(e)
            raise
    
    def _change_state(self, new_state: CircuitState) -> None:
        """
        Change the state of the circuit breaker.
        
        Args:
            new_state: New state for the circuit breaker
        """
        if self.state != new_state:
            logger.info(f"Circuit {self.name} state change: {self.state.value} -> {new_state.value}")
            self.state = new_state
            self.stats["state_changes"].append({
                "time": time.time(),
                "from": self.state.value,
                "to": new_state.value
            })
    
    def _on_success(self) -> None:
        """Handle successful function call."""
        self.stats["success_count"] += 1
        
        if self.state == CircuitState.HALF_OPEN:
            # If we're testing the waters and it succeeded, close the circuit
            self._change_state(CircuitState.CLOSED)
            self.failure_count = 0
    
    def _on_failure(self, exception: Exception) -> None:
        """
        Handle failed function call.
        
        Args:
            exception: Exception that was raised
        """
        self.stats["failure_count"] += 1
        self.last_failure_time = time.time()
        
        if self.state == CircuitState.CLOSED:
            self.failure_count += 1
            
            if self.failure_count >= self.failure_threshold:
                self._change_state(CircuitState.OPEN)
        
        elif self.state == CircuitState.HALF_OPEN:
            # If we're testing the waters and it failed, open the circuit again
            self._change_state(CircuitState.OPEN)
    
    def reset(self) -> None:
        """Reset the circuit breaker to its initial state."""
        self._change_state(CircuitState.CLOSED)
        self.failure_count = 0
        self.last_failure_time = 0
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get circuit breaker statistics.
        
        Returns:
            Dictionary with circuit breaker statistics
        """
        return {
            "name": self.name,
            "state": self.state.value,
            "failure_count": self.failure_count,
            "last_failure_time": self.last_failure_time,
            "stats": self.stats
        }


class CircuitBreakerError(Exception):
    """Exception raised when a circuit breaker is open."""
    pass


def retry_with_backoff(
    max_retries: int = 5,
    base_delay: float = 1.0,
    max_delay: float = 60.0,
    backoff_factor: float = 2.0,
    jitter: bool = True,
    expected_exceptions: tuple = (Exception,)
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator for retrying a function with exponential backoff.
    
    Args:
        max_retries: Maximum number of retries
        base_delay: Initial delay between retries in seconds
        max_delay: Maximum delay between retries in seconds
        backoff_factor: Factor to increase delay for each retry
        jitter: Whether to add random jitter to delay
        expected_exceptions: Tuple of exceptions to catch and retry
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            retry_count = 0
            delay = base_delay
            
            while True:
                try:
                    return func(*args, **kwargs)
                except expected_exceptions as e:
                    retry_count += 1
                    
                    if retry_count > max_retries:
                        logger.error(
                            f"Maximum retries ({max_retries}) exceeded for {func.__name__}. "
                            f"Last error: {str(e)}"
                        )
                        raise
                    
                    # Calculate delay with optional jitter
                    if jitter:
                        actual_delay = delay * (0.5 + random.random())
                    else:
                        actual_delay = delay
                    
                    logger.warning(
                        f"Retry {retry_count}/{max_retries} for {func.__name__} "
                        f"after {actual_delay:.2f}s. Error: {str(e)}"
                    )
                    
                    time.sleep(actual_delay)
                    
                    # Increase delay for next retry, but cap at max_delay
                    delay = min(delay * backoff_factor, max_delay)
        
        return wrapper
    
    return decorator


def create_request_key(url: str, params: Optional[Dict[str, Any]] = None) -> str:
    """
    Create a unique key for caching based on request URL and parameters.
    
    Args:
        url: Request URL
        params: Request parameters
        
    Returns:
        String key for caching
    """
    if params:
        # Sort params to ensure consistent keys
        param_str = "&".join(f"{k}={v}" for k, v in sorted(params.items()))
        return f"{url}?{param_str}"
    return url
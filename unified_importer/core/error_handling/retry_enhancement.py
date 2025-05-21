"""
Enhanced retry mechanism for the unified molecular importer.

This module provides advanced retry capabilities including:
- Exponential backoff with jitter for API rate limiting
- Circuit breaker pattern for persistent errors
- Custom retry policies configurable through settings
"""

import time
import random
import logging
import asyncio
import json
import os
from enum import Enum, auto
from typing import Callable, TypeVar, Generic, Dict, List, Optional, Set, Tuple, Any, Union
from datetime import datetime, timedelta
from functools import wraps
from dataclasses import dataclass, field

from .error_classification import ErrorCategory, ErrorSeverity, RecoveryStrategy, ErrorContext, ClassifiedError, ErrorClassifier

# Type variables for generic functions
T = TypeVar('T')
U = TypeVar('U')

class CircuitState(Enum):
    """Represents the state of a circuit breaker."""
    CLOSED = auto()      # Normal operation - requests flow through
    OPEN = auto()        # Failure threshold reached - requests fail fast
    HALF_OPEN = auto()   # Testing if system recovered - limited requests allowed
    
@dataclass
class CircuitBreakerConfig:
    """Configuration for a circuit breaker."""
    failure_threshold: int = 5            # Number of failures before opening the circuit
    reset_timeout: float = 60.0           # Seconds before attempting half-open state
    half_open_max_calls: int = 1          # Maximum number of calls in half-open state
    success_threshold: int = 1            # Successes needed to close circuit from half-open
    open_timeout_increment_factor: float = 1.5  # Factor to increase timeout on repeated failures
    max_open_timeout: float = 300.0       # Maximum timeout in seconds (5 minutes)

@dataclass
class RetryConfig:
    """Configuration for retry behavior."""
    max_attempts: int = 3                  # Maximum number of retry attempts
    initial_delay: float = 0.1             # Initial delay before first retry (seconds)
    max_delay: float = 60.0                # Maximum delay between retries (seconds)
    backoff_factor: float = 2.0            # Exponential backoff multiplier
    jitter_factor: float = 0.1             # Random jitter factor (0.0-1.0) to add to delay
    non_retryable_errors: Set[type] = field(default_factory=set)  # Exception types that shouldn't be retried
    retryable_error_categories: Set[ErrorCategory] = field(default_factory=set)  # Error categories eligible for retry

@dataclass
class RetryState:
    """Tracks state for a retry operation."""
    attempts: int = 0                      # Current attempt number
    last_delay: float = 0.0                # Last delay used
    start_time: float = field(default_factory=time.time)  # Start time of first attempt
    total_delay: float = 0.0               # Total time spent in delay
    failures: List[Exception] = field(default_factory=list)  # List of failures encountered
    classified_errors: List[ClassifiedError] = field(default_factory=list)  # Classified errors

@dataclass
class CircuitBreakerState:
    """Tracks the state of a circuit breaker."""
    state: CircuitState = CircuitState.CLOSED    # Current circuit state
    failure_count: int = 0                       # Count of consecutive failures
    success_count: int = 0                       # Count of consecutive successes (for half-open)
    last_failure_time: float = 0.0               # Time of last failure
    open_until: float = 0.0                      # Time when circuit will try half-open
    current_timeout: float = 0.0                 # Current timeout duration
    half_open_calls: int = 0                     # Number of calls in half-open state
    service_identifier: str = ""                 # Identifier for the service
    operation_identifier: str = ""               # Identifier for the operation
    
    def should_execute(self) -> bool:
        """Determine if execution should be allowed based on circuit state."""
        now = time.time()
        
        if self.state == CircuitState.CLOSED:
            return True
        elif self.state == CircuitState.OPEN:
            # Check if we should move to half-open
            if now >= self.open_until:
                self.state = CircuitState.HALF_OPEN
                self.half_open_calls = 0
                return True
            else:
                return False
        elif self.state == CircuitState.HALF_OPEN:
            # Allow limited calls in half-open state
            if self.half_open_calls < self.current_half_open_max:
                self.half_open_calls += 1
                return True
            else:
                return False
        
        return False
    
    @property
    def current_half_open_max(self) -> int:
        """Maximum number of calls allowed in half-open state."""
        return getattr(self, 'half_open_max_calls', 1)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the state to a dictionary for serialization."""
        return {
            'state': self.state.name,
            'failure_count': self.failure_count,
            'success_count': self.success_count,
            'last_failure_time': self.last_failure_time,
            'open_until': self.open_until,
            'current_timeout': self.current_timeout,
            'half_open_calls': self.half_open_calls,
            'service_identifier': self.service_identifier,
            'operation_identifier': self.operation_identifier
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CircuitBreakerState':
        """Create a state object from a dictionary."""
        state = cls(
            state=CircuitState[data['state']],
            failure_count=data['failure_count'],
            success_count=data['success_count'],
            last_failure_time=data['last_failure_time'],
            open_until=data['open_until'],
            current_timeout=data['current_timeout'],
            half_open_calls=data['half_open_calls'],
            service_identifier=data['service_identifier'],
            operation_identifier=data['operation_identifier']
        )
        return state


class CircuitBreakerRegistry:
    """Registry for circuit breakers to maintain state across calls."""
    
    def __init__(self, state_file_path: Optional[str] = None, logger: Optional[logging.Logger] = None):
        """Initialize the circuit breaker registry.
        
        Args:
            state_file_path: Optional path to persist circuit breaker states
            logger: Optional logger instance
        """
        self._breakers: Dict[str, CircuitBreakerState] = {}
        self._state_file_path = state_file_path
        self._logger = logger or logging.getLogger(__name__)
        
        if state_file_path:
            self._load_states()
    
    def get_breaker(self, service: str, operation: str) -> CircuitBreakerState:
        """Get or create a circuit breaker for the specified service and operation.
        
        Args:
            service: Service identifier (e.g., "ChEMBL API")
            operation: Operation identifier (e.g., "fetch_molecules")
            
        Returns:
            Circuit breaker state
        """
        key = f"{service}:{operation}"
        if key not in self._breakers:
            self._breakers[key] = CircuitBreakerState(
                service_identifier=service,
                operation_identifier=operation
            )
        return self._breakers[key]
    
    def report_success(self, service: str, operation: str) -> None:
        """Report a successful operation to the circuit breaker.
        
        Args:
            service: Service identifier
            operation: Operation identifier
        """
        breaker = self.get_breaker(service, operation)
        
        if breaker.state == CircuitState.CLOSED:
            # Reset failure count if we're in closed state
            breaker.failure_count = 0
        elif breaker.state == CircuitState.HALF_OPEN:
            # In half-open, track successes and potentially close circuit
            breaker.success_count += 1
            if breaker.success_count >= getattr(breaker, 'success_threshold', 1):
                # Reset to closed state
                breaker.state = CircuitState.CLOSED
                breaker.failure_count = 0
                breaker.success_count = 0
                breaker.current_timeout = 0.0
                self._logger.info(
                    f"Circuit breaker for {service}:{operation} closed after successful recovery"
                )
        
        # Save state changes
        if self._state_file_path:
            self._save_states()
    
    def report_failure(self, service: str, operation: str, config: CircuitBreakerConfig) -> None:
        """Report a failed operation to the circuit breaker.
        
        Args:
            service: Service identifier
            operation: Operation identifier
            config: Circuit breaker configuration
        """
        breaker = self.get_breaker(service, operation)
        now = time.time()
        
        # Increment failure count regardless of state
        breaker.failure_count += 1
        breaker.last_failure_time = now
        
        if breaker.state == CircuitState.CLOSED:
            # Check if we should open the circuit
            if breaker.failure_count >= config.failure_threshold:
                # Calculate timeout (possibly increasing if this isn't the first time)
                if breaker.current_timeout == 0.0:
                    timeout = config.reset_timeout
                else:
                    # Increase timeout by the configured factor, up to the maximum
                    timeout = min(
                        breaker.current_timeout * config.open_timeout_increment_factor,
                        config.max_open_timeout
                    )
                
                breaker.state = CircuitState.OPEN
                breaker.current_timeout = timeout
                breaker.open_until = now + timeout
                breaker.success_count = 0
                
                self._logger.warning(
                    f"Circuit breaker for {service}:{operation} opened for {timeout:.2f}s "
                    f"after {breaker.failure_count} failures"
                )
                
        elif breaker.state == CircuitState.HALF_OPEN:
            # If we fail during half-open, go back to open with increased timeout
            timeout = min(
                breaker.current_timeout * config.open_timeout_increment_factor,
                config.max_open_timeout
            )
            
            breaker.state = CircuitState.OPEN
            breaker.current_timeout = timeout
            breaker.open_until = now + timeout
            breaker.success_count = 0
            
            self._logger.warning(
                f"Circuit breaker for {service}:{operation} re-opened for {timeout:.2f}s "
                f"after failure in half-open state"
            )
        
        # Save state changes
        if self._state_file_path:
            self._save_states()
    
    def _save_states(self) -> None:
        """Save circuit breaker states to disk."""
        try:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(self._state_file_path), exist_ok=True)
            
            # Convert breaker states to serializable format
            data = {
                key: breaker.to_dict() 
                for key, breaker in self._breakers.items()
            }
            
            with open(self._state_file_path, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            self._logger.error(f"Failed to save circuit breaker states: {e}")
    
    def _load_states(self) -> None:
        """Load circuit breaker states from disk."""
        if not os.path.exists(self._state_file_path):
            return
        
        try:
            with open(self._state_file_path, 'r') as f:
                data = json.load(f)
            
            for key, state_data in data.items():
                try:
                    # Only load if the circuit is still open
                    if state_data['state'] == 'OPEN' and state_data['open_until'] > time.time():
                        service, operation = key.split(':', 1)
                        self._breakers[key] = CircuitBreakerState.from_dict(state_data)
                    elif state_data['state'] == 'HALF_OPEN':
                        # Also load half-open circuits
                        service, operation = key.split(':', 1)
                        self._breakers[key] = CircuitBreakerState.from_dict(state_data)
                except Exception as e:
                    self._logger.warning(f"Failed to load circuit breaker state for {key}: {e}")
        except Exception as e:
            self._logger.error(f"Failed to load circuit breaker states: {e}")


class EnhancedRetryManager:
    """Enhanced retry manager with circuit breaker pattern and configurable policies."""
    
    def __init__(
        self, 
        default_config: Optional[RetryConfig] = None,
        circuit_breaker_config: Optional[CircuitBreakerConfig] = None,
        circuit_breaker_registry: Optional[CircuitBreakerRegistry] = None,
        error_classifier: Optional[ErrorClassifier] = None,
        logger: Optional[logging.Logger] = None
    ):
        """Initialize the enhanced retry manager.
        
        Args:
            default_config: Default retry configuration
            circuit_breaker_config: Circuit breaker configuration
            circuit_breaker_registry: Registry for circuit breakers
            error_classifier: Error classifier for categorizing exceptions
            logger: Logger instance
        """
        self.default_config = default_config or RetryConfig()
        self.circuit_breaker_config = circuit_breaker_config or CircuitBreakerConfig()
        self.circuit_breaker_registry = circuit_breaker_registry or CircuitBreakerRegistry()
        self.error_classifier = error_classifier or ErrorClassifier()
        self.logger = logger or logging.getLogger(__name__)
        
        # Set up default retryable error categories if not specified
        if not self.default_config.retryable_error_categories:
            self.default_config.retryable_error_categories = {
                ErrorCategory.NETWORK,
                ErrorCategory.TIMEOUT,
                ErrorCategory.DATABASE,
                ErrorCategory.EXTERNAL_SERVICE
            }
        
        # Set up default non-retryable errors if not specified
        if not self.default_config.non_retryable_errors:
            self.default_config.non_retryable_errors = {
                ValueError,
                TypeError,
                KeyError,
                AttributeError,
                NotImplementedError,
                PermissionError
            }
    
    def execute(
        self, 
        func: Callable[..., T], 
        component: str, 
        operation: str, 
        args: Optional[Tuple] = None, 
        kwargs: Optional[Dict[str, Any]] = None,
        config: Optional[RetryConfig] = None,
        use_circuit_breaker: bool = True,
        circuit_breaker_config: Optional[CircuitBreakerConfig] = None
    ) -> T:
        """Execute a function with retry and circuit breaker logic.
        
        Args:
            func: Function to execute
            component: Component name for error context
            operation: Operation name for error context
            args: Positional arguments for the function
            kwargs: Keyword arguments for the function
            config: Custom retry configuration
            use_circuit_breaker: Whether to use circuit breaker pattern
            circuit_breaker_config: Custom circuit breaker configuration
            
        Returns:
            Result from the function
            
        Raises:
            Exception: If max retries are exceeded or circuit is open
        """
        args = args or ()
        kwargs = kwargs or {}
        config = config or self.default_config
        circuit_breaker_config = circuit_breaker_config or self.circuit_breaker_config
        
        retry_state = RetryState()
        
        # Check circuit breaker before attempting execution
        if use_circuit_breaker:
            circuit = self.circuit_breaker_registry.get_breaker(component, operation)
            if not circuit.should_execute():
                until = datetime.fromtimestamp(circuit.open_until).strftime('%H:%M:%S')
                raise Exception(
                    f"Circuit breaker open for {component}:{operation} until {until}"
                )
        
        while retry_state.attempts <= config.max_attempts:
            try:
                retry_state.attempts += 1
                
                # Execute the function
                result = func(*args, **kwargs)
                
                # Report success to circuit breaker if used
                if use_circuit_breaker:
                    self.circuit_breaker_registry.report_success(component, operation)
                
                return result
            
            except Exception as e:
                # Create error context
                context = ErrorContext(
                    component=component,
                    operation=operation,
                    data={"attempt": retry_state.attempts, "args": str(args), "kwargs": str(kwargs)},
                    metadata={"start_time": retry_state.start_time}
                )
                
                # Classify the error
                classified_error = self.error_classifier.classify(e, context)
                retry_state.classified_errors.append(classified_error)
                retry_state.failures.append(e)
                
                # Check if we've reached max attempts
                if retry_state.attempts >= config.max_attempts:
                    # Report failure to circuit breaker if used
                    if use_circuit_breaker:
                        self.circuit_breaker_registry.report_failure(
                            component, operation, circuit_breaker_config
                        )
                    
                    # Construct detailed error message
                    attempts_str = f"Max retry attempts ({config.max_attempts}) exceeded for {component}:{operation}. "
                    errors_str = ", ".join([f"{type(err).__name__}: {str(err)}" for err in retry_state.failures])
                    
                    raise Exception(f"{attempts_str}Errors: {errors_str}") from e
                
                # Check if this error is retryable
                if (
                    type(e) in config.non_retryable_errors or
                    classified_error.category not in config.retryable_error_categories or
                    classified_error.recovery_strategy not in {RecoveryStrategy.RETRY, RecoveryStrategy.DELAYED_RETRY, RecoveryStrategy.CIRCUIT_BREAKER}
                ):
                    # Report failure to circuit breaker if used
                    if use_circuit_breaker:
                        self.circuit_breaker_registry.report_failure(
                            component, operation, circuit_breaker_config
                        )
                    
                    # Re-raise non-retryable errors
                    raise
                
                # Calculate delay with exponential backoff and jitter
                delay = min(
                    config.initial_delay * (config.backoff_factor ** (retry_state.attempts - 1)),
                    config.max_delay
                )
                
                # Add jitter to prevent retry storms in distributed systems
                jitter = random.uniform(-config.jitter_factor, config.jitter_factor) * delay
                delay = max(0, delay + jitter)
                
                retry_state.last_delay = delay
                retry_state.total_delay += delay
                
                self.logger.warning(
                    f"Retry {retry_state.attempts}/{config.max_attempts} for {component}:{operation} "
                    f"after error: {type(e).__name__}: {str(e)}. "
                    f"Retrying in {delay:.2f}s"
                )
                
                # Wait before retry
                time.sleep(delay)
    
    async def execute_async(
        self, 
        func: Callable[..., T], 
        component: str, 
        operation: str, 
        args: Optional[Tuple] = None, 
        kwargs: Optional[Dict[str, Any]] = None,
        config: Optional[RetryConfig] = None,
        use_circuit_breaker: bool = True,
        circuit_breaker_config: Optional[CircuitBreakerConfig] = None
    ) -> T:
        """Execute an async function with retry and circuit breaker logic.
        
        Args:
            func: Async function to execute
            component: Component name for error context
            operation: Operation name for error context
            args: Positional arguments for the function
            kwargs: Keyword arguments for the function
            config: Custom retry configuration
            use_circuit_breaker: Whether to use circuit breaker pattern
            circuit_breaker_config: Custom circuit breaker configuration
            
        Returns:
            Result from the function
            
        Raises:
            Exception: If max retries are exceeded or circuit is open
        """
        args = args or ()
        kwargs = kwargs or {}
        config = config or self.default_config
        circuit_breaker_config = circuit_breaker_config or self.circuit_breaker_config
        
        retry_state = RetryState()
        
        # Check circuit breaker before attempting execution
        if use_circuit_breaker:
            circuit = self.circuit_breaker_registry.get_breaker(component, operation)
            if not circuit.should_execute():
                until = datetime.fromtimestamp(circuit.open_until).strftime('%H:%M:%S')
                raise Exception(
                    f"Circuit breaker open for {component}:{operation} until {until}"
                )
        
        while retry_state.attempts <= config.max_attempts:
            try:
                retry_state.attempts += 1
                
                # Execute the async function
                result = await func(*args, **kwargs)
                
                # Report success to circuit breaker if used
                if use_circuit_breaker:
                    self.circuit_breaker_registry.report_success(component, operation)
                
                return result
            
            except Exception as e:
                # Create error context
                context = ErrorContext(
                    component=component,
                    operation=operation,
                    data={"attempt": retry_state.attempts, "args": str(args), "kwargs": str(kwargs)},
                    metadata={"start_time": retry_state.start_time}
                )
                
                # Classify the error
                classified_error = self.error_classifier.classify(e, context)
                retry_state.classified_errors.append(classified_error)
                retry_state.failures.append(e)
                
                # Check if we've reached max attempts
                if retry_state.attempts >= config.max_attempts:
                    # Report failure to circuit breaker if used
                    if use_circuit_breaker:
                        self.circuit_breaker_registry.report_failure(
                            component, operation, circuit_breaker_config
                        )
                    
                    # Construct detailed error message
                    attempts_str = f"Max retry attempts ({config.max_attempts}) exceeded for {component}:{operation}. "
                    errors_str = ", ".join([f"{type(err).__name__}: {str(err)}" for err in retry_state.failures])
                    
                    raise Exception(f"{attempts_str}Errors: {errors_str}") from e
                
                # Check if this error is retryable
                if (
                    type(e) in config.non_retryable_errors or
                    classified_error.category not in config.retryable_error_categories or
                    classified_error.recovery_strategy not in {RecoveryStrategy.RETRY, RecoveryStrategy.DELAYED_RETRY, RecoveryStrategy.CIRCUIT_BREAKER}
                ):
                    # Report failure to circuit breaker if used
                    if use_circuit_breaker:
                        self.circuit_breaker_registry.report_failure(
                            component, operation, circuit_breaker_config
                        )
                    
                    # Re-raise non-retryable errors
                    raise
                
                # Calculate delay with exponential backoff and jitter
                delay = min(
                    config.initial_delay * (config.backoff_factor ** (retry_state.attempts - 1)),
                    config.max_delay
                )
                
                # Add jitter to prevent retry storms in distributed systems
                jitter = random.uniform(-config.jitter_factor, config.jitter_factor) * delay
                delay = max(0, delay + jitter)
                
                retry_state.last_delay = delay
                retry_state.total_delay += delay
                
                self.logger.warning(
                    f"Retry {retry_state.attempts}/{config.max_attempts} for {component}:{operation} "
                    f"after error: {type(e).__name__}: {str(e)}. "
                    f"Retrying in {delay:.2f}s"
                )
                
                # Wait before retry
                await asyncio.sleep(delay)

    def retry(
        self, 
        func: Union[Callable[..., T], Callable[..., Awaitable[T]]],
        component: str, 
        operation: str,
        *args,
        max_attempts: Optional[int] = None,
        initial_delay: Optional[float] = None,
        max_delay: Optional[float] = None,
        backoff_factor: Optional[float] = None,
        jitter_factor: Optional[float] = None,
        use_circuit_breaker: bool = True,
        **kwargs
    ) -> Union[T, Awaitable[T]]:
        """Decorator-style function for retrying operations.
        
        Args:
            func: Function to retry
            component: Component name for context
            operation: Operation name for context
            *args: Arguments to pass to the function
            max_attempts: Maximum retry attempts (overrides default)
            initial_delay: Initial delay before retry (overrides default)
            max_delay: Maximum delay between retries (overrides default)
            backoff_factor: Backoff multiplier (overrides default)
            jitter_factor: Jitter factor (overrides default)
            use_circuit_breaker: Whether to use circuit breaker pattern
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            Result from the function or awaitable
        """
        # Create custom config if any parameters were specified
        config = None
        if any(param is not None for param in [max_attempts, initial_delay, max_delay, backoff_factor, jitter_factor]):
            config = RetryConfig(
                max_attempts=max_attempts or self.default_config.max_attempts,
                initial_delay=initial_delay or self.default_config.initial_delay,
                max_delay=max_delay or self.default_config.max_delay,
                backoff_factor=backoff_factor or self.default_config.backoff_factor,
                jitter_factor=jitter_factor or self.default_config.jitter_factor,
                non_retryable_errors=self.default_config.non_retryable_errors,
                retryable_error_categories=self.default_config.retryable_error_categories
            )
        
        # Check if this is an async function
        if asyncio.iscoroutinefunction(func):
            return self.execute_async(
                func, component, operation, args=args, kwargs=kwargs,
                config=config, use_circuit_breaker=use_circuit_breaker
            )
        else:
            return self.execute(
                func, component, operation, args=args, kwargs=kwargs,
                config=config, use_circuit_breaker=use_circuit_breaker
            )
    
    def retry_decorator(
        self, 
        component: str, 
        operation: str,
        max_attempts: Optional[int] = None,
        initial_delay: Optional[float] = None,
        max_delay: Optional[float] = None,
        backoff_factor: Optional[float] = None,
        jitter_factor: Optional[float] = None,
        use_circuit_breaker: bool = True
    ) -> Callable[[Callable[..., T]], Callable[..., T]]:
        """Create a decorator for retrying functions.
        
        Args:
            component: Component name for context
            operation: Operation name for context
            max_attempts: Maximum retry attempts (overrides default)
            initial_delay: Initial delay before retry (overrides default)
            max_delay: Maximum delay between retries (overrides default)
            backoff_factor: Backoff multiplier (overrides default)
            jitter_factor: Jitter factor (overrides default)
            use_circuit_breaker: Whether to use circuit breaker pattern
            
        Returns:
            Decorator function
        """
        # Create custom config if any parameters were specified
        config = None
        if any(param is not None for param in [max_attempts, initial_delay, max_delay, backoff_factor, jitter_factor]):
            config = RetryConfig(
                max_attempts=max_attempts or self.default_config.max_attempts,
                initial_delay=initial_delay or self.default_config.initial_delay,
                max_delay=max_delay or self.default_config.max_delay,
                backoff_factor=backoff_factor or self.default_config.backoff_factor,
                jitter_factor=jitter_factor or self.default_config.jitter_factor,
                non_retryable_errors=self.default_config.non_retryable_errors,
                retryable_error_categories=self.default_config.retryable_error_categories
            )
        
        def decorator(func: Callable[..., T]) -> Callable[..., T]:
            if asyncio.iscoroutinefunction(func):
                @wraps(func)
                async def async_wrapper(*args, **kwargs):
                    return await self.execute_async(
                        func, component, operation, args=args, kwargs=kwargs,
                        config=config, use_circuit_breaker=use_circuit_breaker
                    )
                return async_wrapper
            else:
                @wraps(func)
                def wrapper(*args, **kwargs):
                    return self.execute(
                        func, component, operation, args=args, kwargs=kwargs,
                        config=config, use_circuit_breaker=use_circuit_breaker
                    )
                return wrapper
        
        return decorator
    
    def configure_from_dict(self, config_dict: Dict[str, Any]) -> None:
        """Configure retry settings from a dictionary.
        
        Args:
            config_dict: Dictionary with configuration parameters
        """
        if 'retry' in config_dict:
            retry_config = config_dict['retry']
            
            # Update default retry config
            if 'max_attempts' in retry_config:
                self.default_config.max_attempts = retry_config['max_attempts']
            if 'initial_delay' in retry_config:
                self.default_config.initial_delay = retry_config['initial_delay']
            if 'max_delay' in retry_config:
                self.default_config.max_delay = retry_config['max_delay']
            if 'backoff_factor' in retry_config:
                self.default_config.backoff_factor = retry_config['backoff_factor']
            if 'jitter_factor' in retry_config:
                self.default_config.jitter_factor = retry_config['jitter_factor']
        
        if 'circuit_breaker' in config_dict:
            cb_config = config_dict['circuit_breaker']
            
            # Update circuit breaker config
            if 'failure_threshold' in cb_config:
                self.circuit_breaker_config.failure_threshold = cb_config['failure_threshold']
            if 'reset_timeout' in cb_config:
                self.circuit_breaker_config.reset_timeout = cb_config['reset_timeout']
            if 'half_open_max_calls' in cb_config:
                self.circuit_breaker_config.half_open_max_calls = cb_config['half_open_max_calls']
            if 'success_threshold' in cb_config:
                self.circuit_breaker_config.success_threshold = cb_config['success_threshold']
            if 'open_timeout_increment_factor' in cb_config:
                self.circuit_breaker_config.open_timeout_increment_factor = cb_config['open_timeout_increment_factor']
            if 'max_open_timeout' in cb_config:
                self.circuit_breaker_config.max_open_timeout = cb_config['max_open_timeout']
    
    def get_circuit_breaker_status(self) -> Dict[str, Dict[str, Any]]:
        """Get the status of all circuit breakers.
        
        Returns:
            Dictionary of circuit breaker states
        """
        return {
            key: breaker.to_dict() 
            for key, breaker in self.circuit_breaker_registry._breakers.items()
        }
    
    def reset_circuit_breaker(self, service: str, operation: str) -> None:
        """Reset a circuit breaker to closed state.
        
        Args:
            service: Service identifier
            operation: Operation identifier
        """
        circuit = self.circuit_breaker_registry.get_breaker(service, operation)
        circuit.state = CircuitState.CLOSED
        circuit.failure_count = 0
        circuit.success_count = 0
        circuit.current_timeout = 0.0
        circuit.half_open_calls = 0
        
        self.logger.info(f"Circuit breaker for {service}:{operation} manually reset to CLOSED state")
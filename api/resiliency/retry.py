"""
CryoProtect v2 API - Retry Mechanism

This module provides a retry mechanism with exponential backoff for handling
transient failures in API calls and other operations.

Key features:
- Decorator-based approach for easy integration
- Exponential backoff with jitter for optimal retry timing
- Configurable retry conditions based on exception types
- Automatic logging and observability integration
- Prometheus metrics for tracking retry attempts

Usage:
    from api.resiliency.retry import retry_with_backoff
    
    # Basic usage with default settings
    @retry_with_backoff()
    def external_api_call():
        # This function will retry up to 3 times with exponential backoff
        # if it raises any exception
        pass
        
    # Advanced usage with custom settings
    @retry_with_backoff(
        max_retries=5,                              # Maximum number of retry attempts
        exceptions=(ConnectionError, TimeoutError), # Specific exceptions to retry on
        base_delay=0.5,                             # Initial delay in seconds
        max_delay=30,                               # Maximum delay in seconds
        backoff_factor=2,                           # Multiplicative factor for backoff
        jitter=0.1                                  # Random jitter factor (0.0-1.0)
    )
    def database_operation():
        # This function will retry up to 5 times with exponential backoff
        # if it raises ConnectionError or TimeoutError
        pass

Benefits:
- Improves system stability by handling transient failures
- Reduces cascading failures in distributed systems
- Provides detailed logging and metrics for debugging
- Implements industry best practices for backoff algorithms
"""

import time
import random
import functools
import logging
from typing import Callable, Tuple, Type, Optional, Any, Dict, Union, List
import traceback

# Get logger
logger = logging.getLogger(__name__)

# Import observability if available
try:
    from ..observability import (
        log_with_context, 
        report_error, 
        get_correlation_id, 
        get_request_context,
        trace_request
    )
    from monitoring.prometheus_metrics import RETRY_ATTEMPTS
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
    
    def trace_request(func):
        return func
    
    # Dummy metrics
    class DummyCounter:
        def __init__(self):
            pass
        
        def labels(self, **kwargs):
            return self
            
        def inc(self):
            pass
    
    RETRY_ATTEMPTS = DummyCounter()
    HAS_OBSERVABILITY = False

def retry_with_backoff(
    max_retries: int = 3,
    exceptions: Tuple[Type[Exception], ...] = (Exception,),
    base_delay: float = 0.1,
    max_delay: float = 10.0,
    backoff_factor: float = 2.0,
    jitter: float = 0.1,
    on_retry: Optional[Callable[[Exception, int, float], None]] = None,
    giveup_after: Optional[float] = None,
    retry_condition: Optional[Callable[[Exception], bool]] = None
) -> Callable:
    """
    Decorator that retries a function with exponential backoff when it raises specified exceptions.
    
    Args:
        max_retries: Maximum number of retry attempts (default: 3)
        exceptions: Tuple of exception types to retry on (default: (Exception,))
        base_delay: Initial delay in seconds (default: 0.1)
        max_delay: Maximum delay in seconds (default: 10.0)
        backoff_factor: Multiplicative factor for backoff (default: 2.0)
        jitter: Random jitter factor to add (0.0-1.0) (default: 0.1)
        on_retry: Optional callback function called before each retry attempt
        giveup_after: Optional time limit in seconds after which to stop retrying
        retry_condition: Optional function to determine whether to retry based on the exception
        
    Returns:
        Decorated function
        
    The retry delay is calculated as follows:
        delay = min(base_delay * (backoff_factor ** retry_attempt), max_delay)
        jitter_amount = random.uniform(-jitter * delay, jitter * delay)
        final_delay = delay + jitter_amount
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Initialize retry attempt counter
            retry_attempt = 0
            start_time = time.time()
            
            # Keep track of exceptions for better error reporting
            exceptions_history = []
            
            # Generate operation ID for logging
            operation_id = f"{func.__module__}.{func.__name__}"
            
            # Get request context if available
            context = get_request_context() if HAS_OBSERVABILITY else {}
            
            # Add function information to context
            context.update({
                'function': {
                    'name': func.__name__,
                    'module': func.__module__,
                    'max_retries': max_retries,
                    'base_delay': base_delay,
                    'max_delay': max_delay,
                    'backoff_factor': backoff_factor,
                    'jitter': jitter
                }
            })
            
            while True:
                try:
                    # Attempt to call the function
                    return func(*args, **kwargs)
                
                except exceptions as e:
                    # Check if we've exceeded the maximum retry attempts
                    if retry_attempt >= max_retries:
                        # Log failure after maximum retries
                        log_with_context(
                            logger, 'error',
                            f"Operation {operation_id} failed after {retry_attempt} retries: {str(e)}",
                            context={
                                **context,
                                'event_type': 'retry_exhausted',
                                'retry': {
                                    'attempt': retry_attempt,
                                    'max_retries': max_retries,
                                    'exceptions_history': exceptions_history
                                },
                                'error': {
                                    'type': type(e).__name__,
                                    'message': str(e),
                                    'traceback': traceback.format_exc().split('\n')
                                }
                            },
                            exc_info=True
                        )
                        
                        # Re-raise the last exception
                        raise
                    
                    # Check if the retry condition is met
                    if retry_condition and not retry_condition(e):
                        # Log non-retryable exception
                        log_with_context(
                            logger, 'error',
                            f"Operation {operation_id} failed with non-retryable exception: {str(e)}",
                            context={
                                **context,
                                'event_type': 'retry_condition_failed',
                                'retry': {
                                    'attempt': retry_attempt,
                                    'max_retries': max_retries
                                },
                                'error': {
                                    'type': type(e).__name__,
                                    'message': str(e),
                                    'traceback': traceback.format_exc().split('\n')
                                }
                            },
                            exc_info=True
                        )
                        
                        # Re-raise the exception
                        raise
                    
                    # Check if we've exceeded the time limit
                    elapsed_time = time.time() - start_time
                    if giveup_after is not None and elapsed_time > giveup_after:
                        # Log failure after time limit
                        log_with_context(
                            logger, 'error',
                            f"Operation {operation_id} failed after {elapsed_time:.2f}s (time limit: {giveup_after}s): {str(e)}",
                            context={
                                **context,
                                'event_type': 'retry_timeout',
                                'retry': {
                                    'attempt': retry_attempt,
                                    'max_retries': max_retries,
                                    'elapsed_time': elapsed_time,
                                    'giveup_after': giveup_after,
                                    'exceptions_history': exceptions_history
                                },
                                'error': {
                                    'type': type(e).__name__,
                                    'message': str(e),
                                    'traceback': traceback.format_exc().split('\n')
                                }
                            },
                            exc_info=True
                        )
                        
                        # Re-raise the exception
                        raise
                    
                    # Calculate delay with exponential backoff
                    delay = min(base_delay * (backoff_factor ** retry_attempt), max_delay)
                    
                    # Add jitter to prevent synchronized retries in distributed systems
                    jitter_amount = random.uniform(-jitter * delay, jitter * delay)
                    delay = max(0, delay + jitter_amount)  # Ensure delay is not negative
                    
                    # Record the exception
                    exception_info = {
                        'type': type(e).__name__,
                        'message': str(e),
                        'attempt': retry_attempt,
                        'delay': delay
                    }
                    exceptions_history.append(exception_info)
                    
                    # Update Prometheus metrics if available
                    if HAS_OBSERVABILITY:
                        try:
                            RETRY_ATTEMPTS.labels(
                                function=operation_id,
                                exception_type=type(e).__name__
                            ).inc()
                        except Exception:
                            # Ignore metrics errors
                            pass
                    
                    # Log retry attempt
                    log_with_context(
                        logger, 'warning',
                        f"Retrying operation {operation_id} after error: {str(e)}. "
                        f"Attempt {retry_attempt + 1}/{max_retries} in {delay:.2f}s",
                        context={
                            **context,
                            'event_type': 'retry_attempt',
                            'retry': {
                                'attempt': retry_attempt,
                                'max_retries': max_retries,
                                'delay': delay,
                                'elapsed_time': elapsed_time
                            },
                            'error': {
                                'type': type(e).__name__,
                                'message': str(e)
                            }
                        }
                    )
                    
                    # Call custom retry callback if provided
                    if on_retry:
                        try:
                            on_retry(e, retry_attempt, delay)
                        except Exception as callback_error:
                            # Log callback error but continue with retry
                            logger.warning(
                                f"Error in retry callback for {operation_id}: {str(callback_error)}"
                            )
                    
                    # Increment retry attempt counter
                    retry_attempt += 1
                    
                    # Wait before retrying
                    time.sleep(delay)
                
                except Exception as e:
                    # For exceptions not in the retry list, log and re-raise
                    log_with_context(
                        logger, 'error',
                        f"Operation {operation_id} failed with unhandled exception: {str(e)}",
                        context={
                            **context,
                            'event_type': 'unhandled_exception',
                            'error': {
                                'type': type(e).__name__,
                                'message': str(e),
                                'traceback': traceback.format_exc().split('\n')
                            }
                        },
                        exc_info=True
                    )
                    
                    # Re-raise the exception
                    raise
        
        return wrapper
    
    return decorator

# Define commonly used retry configurations
def retry_network_operations(max_retries: int = 3, **kwargs):
    """
    Retry decorator specifically for network operations.
    
    Args:
        max_retries: Maximum number of retry attempts (default: 3)
        **kwargs: Additional arguments to pass to retry_with_backoff
        
    Returns:
        Decorated function
    """
    return retry_with_backoff(
        max_retries=max_retries,
        exceptions=(
            ConnectionError,
            TimeoutError,
            ConnectionRefusedError,
            ConnectionResetError,
            # Add any HTTP client library exceptions here
        ),
        base_delay=0.5,
        max_delay=10.0,
        **kwargs
    )

def retry_database_operations(max_retries: int = 3, **kwargs):
    """
    Retry decorator specifically for database operations.
    
    Args:
        max_retries: Maximum number of retry attempts (default: 3)
        **kwargs: Additional arguments to pass to retry_with_backoff
        
    Returns:
        Decorated function
    """
    # Import database-specific exceptions here to avoid circular imports
    try:
        # PostgreSQL exceptions
        from psycopg2 import OperationalError, InterfaceError
        db_exceptions = (OperationalError, InterfaceError)
    except ImportError:
        # Fallback to generic exceptions
        db_exceptions = (ConnectionError, TimeoutError)
    
    return retry_with_backoff(
        max_retries=max_retries,
        exceptions=db_exceptions,
        base_delay=0.2,
        max_delay=5.0,
        **kwargs
    )
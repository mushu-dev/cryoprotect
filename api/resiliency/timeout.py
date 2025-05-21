"""
CryoProtect v2 API - Timeout Management

This module provides a mechanism to enforce timeouts on operations to prevent
resource blockage in case of slow or hanging external systems.

Key features:
- Decorator-based approach for easy integration
- Thread-safe implementation
- Integration with logging and observability
- Configurable timeout thresholds

Usage:
    from api.resiliency.timeout import with_timeout
    
    # Basic usage with default settings
    @with_timeout()
    def external_api_call():
        # This function will timeout after 30 seconds
        pass
        
    # Advanced usage with custom settings
    @with_timeout(
        seconds=5,                     # Timeout in seconds
        exception_class=CustomTimeout  # Custom exception class to raise
    )
    def database_operation():
        # This function will timeout after 5 seconds
        pass

Benefits:
- Prevents operations from hanging indefinitely
- Frees up resources that would otherwise be blocked
- Provides clear error messages when operations time out
- Can be combined with retry mechanisms for robust error handling
"""

import time
import signal
import functools
import threading
import logging
from typing import Callable, Type, Optional, Any

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
    from monitoring.prometheus_metrics import TIMEOUT_ERRORS
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
    class DummyCounter:
        def __init__(self):
            pass
        
        def labels(self, **kwargs):
            return self
            
        def inc(self):
            pass
    
    TIMEOUT_ERRORS = DummyCounter()
    HAS_OBSERVABILITY = False

class TimeoutError(Exception):
    """Exception raised when an operation times out."""
    pass

class ThreadingTimeout:
    """
    Timeout implementation using threading for Python applications.
    
    This implementation is safe for multi-threaded applications and doesn't
    use signals, which can interfere with other signal handlers.
    """
    
    def __init__(self, seconds, exception_class=TimeoutError):
        """
        Initialize a new threading timeout.
        
        Args:
            seconds: Timeout in seconds
            exception_class: Exception class to raise when timeout occurs
        """
        self.seconds = seconds
        self.exception_class = exception_class
        self.timer = None
        self.timed_out = False
        self.original_thread_id = None
        self.done = threading.Event()
    
    def _timeout_handler(self):
        """
        Handler called when timeout occurs.
        
        This handler sets the timed_out flag and tries to interrupt the main thread.
        """
        if not self.done.is_set():
            self.timed_out = True
            
            # Log timeout
            operation_id = f"thread_{self.original_thread_id}"
            log_with_context(
                logger, 'warning',
                f"Operation {operation_id} timed out after {self.seconds}s",
                context={
                    'event_type': 'operation_timeout',
                    'timeout': {
                        'seconds': self.seconds,
                        'thread_id': self.original_thread_id
                    }
                }
            )
            
            # Update metrics
            if HAS_OBSERVABILITY:
                TIMEOUT_ERRORS.labels(operation=operation_id).inc()
    
    def __enter__(self):
        """
        Start the timeout timer.
        
        Returns:
            self
        """
        self.original_thread_id = threading.current_thread().ident
        self.timed_out = False
        self.done.clear()
        
        # Create and start the timer
        self.timer = threading.Timer(self.seconds, self._timeout_handler)
        self.timer.daemon = True
        self.timer.start()
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Cancel the timeout timer.
        
        Args:
            exc_type: Exception type
            exc_val: Exception value
            exc_tb: Exception traceback
            
        Returns:
            bool: True if exception should be suppressed, False otherwise
        """
        # Mark as done to prevent timeout handler from firing
        self.done.set()
        
        # Cancel the timer
        if self.timer:
            self.timer.cancel()
        
        # Check if we timed out
        if self.timed_out:
            # Raise timeout exception
            raise self.exception_class(f"Operation timed out after {self.seconds} seconds")
        
        # Don't suppress other exceptions
        return False

def with_timeout(
    seconds: float = 30.0,
    exception_class: Type[Exception] = TimeoutError
) -> Callable:
    """
    Decorator that enforces a timeout on a function.
    
    Args:
        seconds: Timeout in seconds (default: 30.0)
        exception_class: Exception class to raise when timeout occurs (default: TimeoutError)
        
    Returns:
        Decorated function
        
    The decorated function will raise the specified exception if it doesn't
    complete within the specified timeout.
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Get operation ID for logging
            operation_id = f"{func.__module__}.{func.__name__}"
            
            # Get request context if available
            context = get_request_context() if HAS_OBSERVABILITY else {}
            
            # Add function information to context
            context.update({
                'function': {
                    'name': func.__name__,
                    'module': func.__module__,
                    'timeout': seconds
                }
            })
            
            try:
                # Use threading timeout
                with ThreadingTimeout(seconds, exception_class):
                    return func(*args, **kwargs)
            
            except exception_class as e:
                # Log timeout
                log_with_context(
                    logger, 'warning',
                    f"Operation {operation_id} timed out after {seconds}s: {str(e)}",
                    context={
                        **context,
                        'event_type': 'function_timeout',
                        'timeout': {
                            'seconds': seconds,
                            'exception': str(e)
                        }
                    }
                )
                
                # Update metrics
                if HAS_OBSERVABILITY:
                    TIMEOUT_ERRORS.labels(function=operation_id).inc()
                
                # Re-raise the exception
                raise
        
        return wrapper
    
    return decorator
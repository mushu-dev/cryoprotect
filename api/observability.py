#!/usr/bin/env python3
"""
CryoProtect v2 - API Observability Module

This module provides comprehensive observability features for the CryoProtect API:
- Request tracing with correlation IDs
- Performance timing middleware
- Detailed error reporting with full context
- Integration with structured logging system

Usage:
    from api.observability import (
        init_observability,
        trace_request,
        track_performance,
        report_error,
        get_correlation_id,
        ObservabilityMiddleware
    )
"""

import time
import uuid
import functools
import traceback
import logging
import json
from typing import Dict, Any, Optional, Callable, Union, Tuple
from datetime import datetime
from flask import Flask, request, g, Response, current_app
from werkzeug.exceptions import HTTPException

from logging_enhanced import log_with_context, get_logger
from monitoring.prometheus_metrics import (
    REQUEST_LATENCY, ERROR_COUNT, REQUEST_COUNT
)

# Set up logger
logger = get_logger(__name__)

# Constants
CORRELATION_ID_HEADER = 'X-Correlation-ID'
REQUEST_ID_HEADER = 'X-Request-ID'
TRACE_ID_HEADER = 'X-Trace-ID'
PARENT_SPAN_ID_HEADER = 'X-Parent-Span-ID'
SPAN_ID_HEADER = 'X-Span-ID'
TIMING_HEADER = 'X-Response-Time-Ms'

# Performance thresholds (in seconds)
PERFORMANCE_THRESHOLDS = {
    'warning': 1.0,  # Log warning if request takes more than 1 second
    'critical': 3.0,  # Log critical if request takes more than 3 seconds
}

class ObservabilityMiddleware:
    """
    Middleware for API observability.
    
    This middleware adds request tracing, performance timing, and error reporting
    to all API requests.
    """
    
    def __init__(self, app: Optional[Flask] = None):
        """
        Initialize the middleware.
        
        Args:
            app: Optional Flask application instance
        """
        self.app = app
        if app is not None:
            self.init_app(app)
    
    def init_app(self, app: Flask) -> None:
        """
        Initialize the middleware with a Flask application.
        
        Args:
            app: Flask application instance
        """
        # Register before_request handler
        @app.before_request
        def before_request():
            # Start request timing (using datetime for consistency with logging_enhanced.py)
            from datetime import datetime
            g.start_time = datetime.utcnow()
            
            # Generate request ID if not already present
            g.request_id = request.headers.get(REQUEST_ID_HEADER) or str(uuid.uuid4())
            
            # Use existing correlation ID from header or generate a new one
            g.correlation_id = request.headers.get(CORRELATION_ID_HEADER) or str(uuid.uuid4())
            
            # Generate trace ID if not already present
            g.trace_id = request.headers.get(TRACE_ID_HEADER) or str(uuid.uuid4())
            
            # Generate span ID for this request
            g.span_id = str(uuid.uuid4())
            
            # Store parent span ID if present
            g.parent_span_id = request.headers.get(PARENT_SPAN_ID_HEADER)
            
            # Extract user ID if authenticated
            g.user_id = getattr(g, 'user_id', 'anonymous')
            
            # Log request with tracing information
            log_with_context(
                logger, 'info',
                f"Request started: {request.method} {request.path}",
                context={
                    'event_type': 'request_started',
                    'trace': {
                        'correlation_id': g.correlation_id,
                        'request_id': g.request_id,
                        'trace_id': g.trace_id,
                        'span_id': g.span_id,
                        'parent_span_id': g.parent_span_id
                    },
                    'request': {
                        'method': request.method,
                        'path': request.path,
                        'remote_addr': request.remote_addr,
                        'user_agent': request.headers.get('User-Agent', 'N/A'),
                        'content_length': request.content_length,
                        'content_type': request.content_type,
                        'query_string': request.query_string.decode('utf-8', errors='replace') if request.query_string else '',
                        'headers': {k: v for k, v in request.headers.items() if k.lower() not in ('authorization', 'cookie')}
                    },
                    'user_id': g.user_id
                }
            )
        
        # Register after_request handler
        @app.after_request
        def after_request(response):
            # Calculate request duration
            duration = None
            if hasattr(g, 'start_time'):
                # Use datetime for consistency with logging_enhanced.py
                from datetime import datetime
                duration = (datetime.utcnow() - g.start_time).total_seconds()
                
                # Add timing header to response
                if duration is not None:
                    response.headers[TIMING_HEADER] = str(int(duration * 1000))
                
                # Log performance warning if request is slow
                if duration > PERFORMANCE_THRESHOLDS['critical']:
                    log_level = 'critical'
                elif duration > PERFORMANCE_THRESHOLDS['warning']:
                    log_level = 'warning'
                else:
                    log_level = 'info'
                
                if log_level != 'info':
                    log_with_context(
                        logger, log_level,
                        f"Slow request: {request.method} {request.path} took {duration:.2f}s",
                        context={
                            'event_type': 'slow_request',
                            'trace': {
                                'correlation_id': getattr(g, 'correlation_id', 'N/A'),
                                'request_id': getattr(g, 'request_id', 'N/A'),
                                'trace_id': getattr(g, 'trace_id', 'N/A'),
                                'span_id': getattr(g, 'span_id', 'N/A')
                            },
                            'request': {
                                'method': request.method,
                                'path': request.path,
                                'duration': duration
                            },
                            'performance': {
                                'threshold_exceeded': log_level,
                                'duration': duration,
                                'threshold': PERFORMANCE_THRESHOLDS[log_level]
                            }
                        }
                    )
            
            # Add tracing headers to response
            if hasattr(g, 'correlation_id'):
                response.headers[CORRELATION_ID_HEADER] = g.correlation_id
            if hasattr(g, 'request_id'):
                response.headers[REQUEST_ID_HEADER] = g.request_id
            if hasattr(g, 'trace_id'):
                response.headers[TRACE_ID_HEADER] = g.trace_id
            if hasattr(g, 'span_id'):
                response.headers[SPAN_ID_HEADER] = g.span_id
            
            # Log response
            status_code = response.status_code
            log_level = 'warning' if status_code >= 400 else 'info'
            
            log_with_context(
                logger, log_level,
                f"Request completed: {request.method} {request.path} {status_code}",
                context={
                    'event_type': 'request_completed',
                    'trace': {
                        'correlation_id': getattr(g, 'correlation_id', 'N/A'),
                        'request_id': getattr(g, 'request_id', 'N/A'),
                        'trace_id': getattr(g, 'trace_id', 'N/A'),
                        'span_id': getattr(g, 'span_id', 'N/A')
                    },
                    'request': {
                        'method': request.method,
                        'path': request.path,
                        'remote_addr': request.remote_addr
                    },
                    'response': {
                        'status_code': status_code,
                        'content_length': response.content_length,
                        'content_type': response.content_type,
                        'duration': duration
                    },
                    'user_id': getattr(g, 'user_id', 'anonymous')
                }
            )
            
            # Update Prometheus metrics
            endpoint = request.endpoint or 'unknown'
            method = request.method
            
            if duration is not None:
                REQUEST_LATENCY.labels(method=method, endpoint=endpoint).observe(duration)
            
            REQUEST_COUNT.labels(method=method, endpoint=endpoint, status=status_code).inc()
            
            return response
        
        # Register error handler
        @app.errorhandler(Exception)
        def handle_exception(error):
            # Get error details
            error_type = type(error).__name__
            error_message = str(error)
            status_code = 500
            
            # Get status code for HTTP exceptions
            if isinstance(error, HTTPException):
                status_code = error.code
            
            # Get traceback
            tb = traceback.format_exc()
            
            # Log error with context
            log_with_context(
                logger, 'error',
                f"Request failed: {request.method} {request.path} - {error_type}: {error_message}",
                context={
                    'event_type': 'request_failed',
                    'trace': {
                        'correlation_id': getattr(g, 'correlation_id', 'N/A'),
                        'request_id': getattr(g, 'request_id', 'N/A'),
                        'trace_id': getattr(g, 'trace_id', 'N/A'),
                        'span_id': getattr(g, 'span_id', 'N/A')
                    },
                    'request': {
                        'method': request.method,
                        'path': request.path,
                        'remote_addr': request.remote_addr,
                        'query_string': request.query_string.decode('utf-8', errors='replace') if request.query_string else '',
                        'headers': {k: v for k, v in request.headers.items() if k.lower() not in ('authorization', 'cookie')}
                    },
                    'error': {
                        'type': error_type,
                        'message': error_message,
                        'traceback': tb.split('\n')
                    },
                    'user_id': getattr(g, 'user_id', 'anonymous')
                },
                exc_info=True
            )
            
            # Update Prometheus metrics
            endpoint = request.endpoint or 'unknown'
            ERROR_COUNT.labels(type=error_type, endpoint=endpoint).inc()
            
            # Let the default error handler handle the exception
            raise error
        
        # Log initialization
        logger.info("API observability middleware initialized")

def trace_request(func: Callable) -> Callable:
    """
    Decorator to trace API requests.
    
    This decorator adds tracing information to API requests, including
    correlation ID, request ID, and span ID.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Generate span ID for this function call
        function_span_id = str(uuid.uuid4())
        
        # Get parent span ID (current request span ID)
        parent_span_id = getattr(g, 'span_id', None)
        
        # Get correlation ID and trace ID
        correlation_id = getattr(g, 'correlation_id', None)
        trace_id = getattr(g, 'trace_id', None)
        
        # Log function entry with tracing information
        log_with_context(
            logger, 'debug',
            f"Function {func.__name__} started",
            context={
                'event_type': 'function_started',
                'trace': {
                    'correlation_id': correlation_id,
                    'trace_id': trace_id,
                    'span_id': function_span_id,
                    'parent_span_id': parent_span_id
                },
                'function': {
                    'name': func.__name__,
                    'module': func.__module__,
                    'args': [str(arg) for arg in args],
                    'kwargs': {k: str(v) for k, v in kwargs.items()}
                }
            }
        )
        
        try:
            # Call the function
            result = func(*args, **kwargs)
            
            # Log function exit
            log_with_context(
                logger, 'debug',
                f"Function {func.__name__} completed",
                context={
                    'event_type': 'function_completed',
                    'trace': {
                        'correlation_id': correlation_id,
                        'trace_id': trace_id,
                        'span_id': function_span_id,
                        'parent_span_id': parent_span_id
                    },
                    'function': {
                        'name': func.__name__,
                        'module': func.__module__
                    }
                }
            )
            
            return result
        except Exception as e:
            # Log function error
            log_with_context(
                logger, 'error',
                f"Function {func.__name__} failed: {str(e)}",
                context={
                    'event_type': 'function_failed',
                    'trace': {
                        'correlation_id': correlation_id,
                        'trace_id': trace_id,
                        'span_id': function_span_id,
                        'parent_span_id': parent_span_id
                    },
                    'function': {
                        'name': func.__name__,
                        'module': func.__module__
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
    
    return wrapper

def track_performance(threshold: float = None) -> Callable:
    """
    Decorator to track function performance.
    
    This decorator measures the execution time of a function and logs a warning
    if it exceeds the specified threshold.
    
    Args:
        threshold: Threshold in seconds (default: PERFORMANCE_THRESHOLDS['warning'])
        
    Returns:
        Decorator function
    """
    if threshold is None:
        threshold = PERFORMANCE_THRESHOLDS['warning']
    
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Start timing
            start_time = time.time()
            
            try:
                # Call the function
                result = func(*args, **kwargs)
                
                # Calculate duration
                duration = time.time() - start_time
                
                # Log performance warning if function is slow
                if duration > threshold:
                    log_with_context(
                        logger, 'warning',
                        f"Slow function: {func.__name__} took {duration:.2f}s (threshold: {threshold:.2f}s)",
                        context={
                            'event_type': 'slow_function',
                            'trace': {
                                'correlation_id': getattr(g, 'correlation_id', 'N/A'),
                                'request_id': getattr(g, 'request_id', 'N/A')
                            },
                            'function': {
                                'name': func.__name__,
                                'module': func.__module__
                            },
                            'performance': {
                                'duration': duration,
                                'threshold': threshold
                            }
                        }
                    )
                
                return result
            except Exception as e:
                # Calculate duration even if function fails
                duration = time.time() - start_time
                
                # Re-raise the exception
                raise
        
        return wrapper
    
    return decorator

def report_error(error: Exception, context: Dict[str, Any] = None) -> None:
    """
    Report an error with detailed context.
    
    Args:
        error: Exception to report
        context: Additional context to include in the error report
    """
    if context is None:
        context = {}
    
    # Get error details
    error_type = type(error).__name__
    error_message = str(error)
    
    # Get traceback
    tb = traceback.format_exc()
    
    # Build error context
    error_context = {
        'event_type': 'error',
        'error': {
            'type': error_type,
            'message': error_message,
            'traceback': tb.split('\n')
        }
    }
    
    # Add tracing information if available
    if hasattr(g, 'correlation_id') or hasattr(g, 'request_id'):
        error_context['trace'] = {
            'correlation_id': getattr(g, 'correlation_id', 'N/A'),
            'request_id': getattr(g, 'request_id', 'N/A'),
            'trace_id': getattr(g, 'trace_id', 'N/A'),
            'span_id': getattr(g, 'span_id', 'N/A')
        }
    
    # Add request information if available
    if request:
        error_context['request'] = {
            'method': request.method,
            'path': request.path,
            'remote_addr': request.remote_addr,
            'user_agent': request.headers.get('User-Agent', 'N/A')
        }
    
    # Add user information if available
    if hasattr(g, 'user_id'):
        error_context['user_id'] = g.user_id
    
    # Add additional context
    error_context.update(context)
    
    # Log error with context
    log_with_context(
        logger, 'error',
        f"Error: {error_type}: {error_message}",
        context=error_context,
        exc_info=True
    )
    
    # Update Prometheus metrics if in a request context
    if request and hasattr(request, 'endpoint'):
        endpoint = request.endpoint or 'unknown'
        ERROR_COUNT.labels(type=error_type, endpoint=endpoint).inc()

def get_correlation_id() -> str:
    """
    Get the current correlation ID.
    
    Returns:
        str: Current correlation ID or a new one if not available
    """
    if hasattr(g, 'correlation_id'):
        return g.correlation_id
    
    # Generate a new correlation ID if not available
    correlation_id = str(uuid.uuid4())
    g.correlation_id = correlation_id
    
    return correlation_id

def get_request_context() -> Dict[str, Any]:
    """
    Get the current request context for logging.
    
    Returns:
        Dict[str, Any]: Request context
    """
    context = {
        'timestamp': datetime.utcnow().isoformat()
    }
    
    # Add tracing information if available
    if hasattr(g, 'correlation_id') or hasattr(g, 'request_id'):
        context['trace'] = {
            'correlation_id': getattr(g, 'correlation_id', 'N/A'),
            'request_id': getattr(g, 'request_id', 'N/A'),
            'trace_id': getattr(g, 'trace_id', 'N/A'),
            'span_id': getattr(g, 'span_id', 'N/A')
        }
    
    # Add request information if available
    if request:
        context['request'] = {
            'method': request.method,
            'path': request.path,
            'remote_addr': request.remote_addr,
            'user_agent': request.headers.get('User-Agent', 'N/A')
        }
    
    # Add user information if available
    if hasattr(g, 'user_id'):
        context['user_id'] = g.user_id
    
    return context

def init_observability(app: Flask) -> None:
    """
    Initialize API observability for a Flask application.
    
    Args:
        app: Flask application instance
    """
    # Initialize observability middleware
    ObservabilityMiddleware(app)
    
    # Log initialization
    logger.info("API observability initialized")
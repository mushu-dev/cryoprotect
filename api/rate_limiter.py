"""
CryoProtect Analyzer - Rate Limiting Module

This module provides rate limiting functionality for the CryoProtect Analyzer API.
It supports user-based, IP-based, and endpoint-based rate limiting.
"""

import time
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Tuple, Optional, Union, Callable
from functools import wraps
from flask import request, g, jsonify, current_app, Response
import redis
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

# Setup logging
logger = logging.getLogger(__name__)

# Initialize rate limiter
limiter = None

def get_user_identifier():
    """
    Get the user identifier for rate limiting.
    Uses user ID if authenticated, otherwise falls back to IP address.
    Based on the RATE_LIMIT_BY configuration.
    
    Returns:
        str: User ID or IP address
    """
    # Import here to avoid circular imports
    from .utils import get_user_id
    from flask import current_app
    
    # Get rate limit by configuration
    rate_limit_by = current_app.config.get('RATE_LIMIT_BY', 'hybrid')
    
    user_id = get_user_id()
    ip_address = get_remote_address()
    
    if rate_limit_by == 'user' and user_id:
        return f"user:{user_id}"
    elif rate_limit_by == 'ip':
        return ip_address
    else:  # hybrid - default
        if user_id:
            return f"user:{user_id}"
        return ip_address

def configure_rate_limiter(app):
    """
    Configure the rate limiter for the application.
    
    Args:
        app: Flask application
    """
    global limiter
    
    # Default configuration
    default_limits = [
        "200 per day",
        "50 per hour",
        "10 per minute"
    ]
    
    # Load configuration from app config if available
    if app.config.get('RATE_LIMIT_ENABLED', True):
        default_limits = app.config.get('RATE_LIMIT_DEFAULT', default_limits)
    else:
        # If rate limiting is disabled, set very high limits
        default_limits = ["1000000 per day"]
    
    # Initialize limiter
    limiter = Limiter(
        app=app,
        key_func=get_user_identifier,
        default_limits=default_limits,
        headers_enabled=app.config.get('RATE_LIMIT_HEADERS_ENABLED', True),
        strategy=app.config.get('RATE_LIMIT_STRATEGY', "fixed-window")
    )
    
    # Configure storage
    if app.config.get('RATE_LIMIT_STORAGE_URL'):
        limiter.storage_uri = app.config.get('RATE_LIMIT_STORAGE_URL')
    
    # Apply endpoint-specific rate limits
    endpoint_limits = app.config.get('RATE_LIMIT_ENDPOINTS', {})
    for endpoint, limits in endpoint_limits.items():
        view_function = app.view_functions.get(endpoint.replace('/', '_').strip('_'))
        if view_function is not None:
            limiter.limit(limits)(view_function)
        else:
            app.logger.warning(f"Rate limit endpoint not found: {endpoint}")
    
    # Apply exempt endpoints
    exempt_endpoints = app.config.get('RATE_LIMIT_EXEMPT', [])
    for endpoint in exempt_endpoints:
        view_function = app.view_functions.get(endpoint.replace('/', '_').strip('_'))
        if view_function is not None:
            limiter.exempt(view_function)
        else:
            app.logger.warning(f"Rate limit exempt endpoint not found: {endpoint}")
    
    # Register error handler
    @app.errorhandler(429)
    def ratelimit_handler(e):
        retry_after = app.config.get('RATE_LIMIT_RETRY_AFTER', 60)
        
        # Try to get the actual retry-after value
        if hasattr(e, 'description') and 'retry after' in e.description:
            try:
                retry_after = int(e.description.split('retry after ')[1].split(' ')[0])
            except (IndexError, ValueError):
                pass
        
        response = jsonify({
            'status': 'error',
            'message': 'Rate limit exceeded',
            'details': str(e.description),
            'retry_after': retry_after
        })
        response.status_code = 429
        response.headers['Retry-After'] = str(retry_after)
        return response
    
    # Apply role-based rate limits
    @app.before_request
    def check_role_limits():
        # Import here to avoid circular imports
        from .jwt_auth import get_current_user
        
        # Get role-based limits
        role_limits = app.config.get('RATE_LIMIT_ROLES', {})
        if not role_limits:
            return
        
        # Get current user
        user = get_current_user()
        if not user:
            return
        
        # Get user role
        role = user.get('role', 'basic')
        
        # Apply role-based limits if defined
        if role in role_limits:
            # Store the role limits in g for potential use in other parts of the request
            g.role_rate_limits = role_limits[role]
    
    logger.info("Rate limiter configured with endpoint-specific and role-based limits")
    return limiter

def get_limiter():
    """
    Get the configured rate limiter instance.
    
    Returns:
        Limiter: Configured rate limiter
    """
    global limiter
    if limiter is None:
        raise RuntimeError("Rate limiter not configured. Call configure_rate_limiter first.")
    return limiter

def endpoint_rate_limit(limit_string: str):
    """
    Decorator to apply rate limiting to a specific endpoint.
    
    Args:
        limit_string: Rate limit string (e.g., "5 per minute")
        
    Returns:
        Decorated function
    """
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            # Get limiter
            limiter = get_limiter()
            
            # Apply rate limit
            limiter.limit(limit_string)(f)(*args, **kwargs)
            
            return f(*args, **kwargs)
        return decorated_function
    return decorator

def exempt_from_rate_limit(f):
    """
    Decorator to exempt an endpoint from rate limiting.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Get limiter
        limiter = get_limiter()
        
        # Exempt from rate limiting
        limiter.exempt(f)
        
        return f(*args, **kwargs)
    return decorated_function

def add_rate_limit_headers(response):
    """
    Add rate limit headers to the response.
    
    Args:
        response: Flask response
        
    Returns:
        Flask response with rate limit headers
    """
    # Check if headers are enabled
    if not current_app.config.get('RATE_LIMIT_HEADERS_ENABLED', True):
        return response
    
    # Get limiter
    try:
        limiter = get_limiter()
        
        # Get current limits
        current_limit = getattr(g, 'view_rate_limit', None)
        if current_limit is not None:
            window_stats = getattr(g, 'view_rate_limit_window', None)
            reset = 0
            if window_stats:
                reset = window_stats.get('reset', 0)
            
            # Standard rate limit headers
            response.headers.add('X-RateLimit-Limit', str(current_limit.limit))
            response.headers.add('X-RateLimit-Remaining', str(current_limit.remaining))
            response.headers.add('X-RateLimit-Reset', str(reset))
            
            # Add RFC 6585 compliant header if we're near the limit
            if current_limit.remaining < (current_limit.limit * 0.1):  # Less than 10% remaining
                retry_after = reset - int(time.time())
                if retry_after < 0:
                    retry_after = 0
                response.headers.add('Retry-After', str(retry_after))
            
            # Add additional headers for more context
            response.headers.add('X-RateLimit-Policy', str(current_limit.key_for_scope))
            
            # Add role-based limit info if available
            role_limits = getattr(g, 'role_rate_limits', None)
            if role_limits:
                response.headers.add('X-RateLimit-Role', str(role_limits))
    except Exception as e:
        logger.warning(f"Error adding rate limit headers: {str(e)}")
    
    return response

def check_rate_limit(key, limit=None, period=None):
    """
    Check if a request is allowed based on rate limits.
    
    Args:
        key: Key to use for rate limiting
        limit: Rate limit to apply (default: from config)
        period: Period in seconds for the rate limit (default: 60)
        
    Returns:
        Tuple[bool, int, datetime]: (allowed, current count, reset time)
    """
    # Get limiter
    limiter = get_limiter()
    
    # Default values
    if period is None:
        period = 60  # Default to 60 seconds
    
    if limit is None:
        limit = current_app.config.get('RATE_LIMIT_DEFAULT', [20])[0]
        if isinstance(limit, str) and 'per' in limit:
            limit = int(limit.split(' ')[0])
    
    # Check if we're using flask-limiter
    if hasattr(limiter, 'limiter'):
        # Get current usage
        current = limiter.limiter.get(key)
        if current is None:
            current = 0
        
        # Calculate reset time
        reset_time = datetime.now().replace(microsecond=0) + timedelta(seconds=period)
        
        # Check if allowed
        allowed = current < limit
        
        return allowed, current, reset_time
    
    # Fallback to simple check
    return True, 0, datetime.now() + timedelta(seconds=period)

def update_rate_limit(key, period=None):
    """
    Update the rate limit counter for a key.
    
    Args:
        key: Key to update
        period: Period in seconds for the rate limit (default: 60)
    """
    # Get limiter
    limiter = get_limiter()
    
    # Update counter if using flask-limiter
    if hasattr(limiter, 'limiter'):
        limiter.limiter.hit(key)
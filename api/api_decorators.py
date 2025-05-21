"""
CryoProtect Analyzer API - Decorators

This module provides decorators for API endpoints to ensure standardized responses,
error handling, and proper HTTP status code usage.
"""

import functools
import logging
from typing import Callable, Any, Dict, Tuple
from flask import request, current_app, abort
from marshmallow import ValidationError
from datetime import datetime

from .api_standards import (
    create_standard_response,
    create_error_response,
    create_success_response,
    jsonify_standard_response,
    HTTP_STATUS_CODES
)

# Set up logging
logger = logging.getLogger(__name__)

def standardize_response(f: Callable) -> Callable:
    """
    Decorator to standardize API responses.
    
    This decorator wraps API endpoint functions to ensure they return
    standardized responses with proper status codes. It also handles
    exceptions and converts them to standardized error responses.
    
    Args:
        f: The function to decorate
        
    Returns:
        Decorated function that returns standardized responses
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        try:
            # Call the original function
            result = f(*args, **kwargs)
            
            # Handle different return types
            if isinstance(result, tuple) and len(result) >= 2:
                # Function returned (data, status_code, [headers])
                data = result[0]
                status_code = result[1]
                headers = result[2] if len(result) > 2 else None
                
                # Check if data is already in standard format
                if isinstance(data, dict) and "status" in data and "timestamp" in data:
                    # Already standardized, just return it
                    response = jsonify_standard_response(data, status_code)
                else:
                    # Create a standardized success response
                    response = jsonify_standard_response(
                        *create_success_response(
                            data=data,
                            status_code=status_code
                        )
                    )
                
                # Add headers if provided
                if headers:
                    response[0].headers.extend(headers)
                    
                return response
            else:
                # Function returned just data
                # Create a standardized success response
                return jsonify_standard_response(
                    *create_success_response(data=result)
                )
                
        except Exception as e:
            # Log the exception
            logger.exception(f"Error in {f.__name__}: {str(e)}")
            
            # Get the endpoint name for context
            endpoint = request.endpoint or f.__name__
            
            # Handle validation errors specifically
            if isinstance(e, ValidationError):
                return jsonify_standard_response(
                    *create_error_response(
                        error=e,
                        status_code=400,
                        context=f"Validation error in {endpoint}",
                        details={"validation_errors": e.messages}
                    )
                )
            
            # Handle other exceptions
            return jsonify_standard_response(
                *create_error_response(
                    error=e,
                    context=f"Error in {endpoint}"
                )
            )
    
    return decorated

def validate_request_schema(schema_class):
    """
    Decorator to validate request data against a schema.
    
    Args:
        schema_class: Marshmallow schema class to validate against
        
    Returns:
        Decorator function
    """
    def decorator(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            try:
                # Create schema instance
                schema = schema_class()
                
                # Get request data based on content type
                if request.is_json:
                    data = request.get_json()
                elif request.form:
                    data = request.form.to_dict()
                else:
                    data = {}
                
                # Validate data against schema
                validated_data = schema.load(data)
                
                # Add validated data to kwargs
                kwargs.update(validated_data)
                
                # Call the original function
                return f(*args, **kwargs)
                
            except ValidationError as e:
                # Log the validation error
                logger.warning(f"Validation error in {f.__name__}: {e.messages}")
                
                # Return standardized error response
                return jsonify_standard_response(
                    *create_error_response(
                        error="Request validation failed",
                        status_code=400,
                        context=f"Validation error in {f.__name__}",
                        details={"validation_errors": e.messages}
                    )
                )
                
            except Exception as e:
                # Log the exception
                logger.exception(f"Error in {f.__name__}: {str(e)}")
                
                # Return standardized error response
                return jsonify_standard_response(
                    *create_error_response(
                        error=e,
                        context=f"Error in {f.__name__}"
                    )
                )
                
        return decorated
    
    return decorator

def rate_limit(limit_key=None, limit=None, period=None):
    """
    Decorator to apply rate limiting to an endpoint.
    
    Args:
        limit_key: Key to use for rate limiting (default: endpoint name)
        limit: Rate limit to apply (default: from config)
        period: Period in seconds for the rate limit (default: 60)
        
    Returns:
        Decorator function
    """
    def decorator(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            # Import here to avoid circular imports
            from api.rate_limiter import check_rate_limit, update_rate_limit
            
            # Get endpoint name if limit_key not provided
            key = limit_key or request.endpoint or f.__name__
            
            # Check rate limit
            allowed, current, reset_time = check_rate_limit(key, limit, period)
            
            if not allowed:
                # Return rate limit exceeded error
                response = jsonify_standard_response(
                    *create_error_response(
                        error="Rate limit exceeded",
                        status_code=429,
                        context=f"Rate limit exceeded for {key}",
                        details={
                            "limit": limit,
                            "current": current,
                            "reset": reset_time.isoformat() if reset_time else None
                        }
                    )
                )
                
                # Add rate limit headers
                response[0].headers.update({
                    "X-RateLimit-Limit": str(limit),
                    "X-RateLimit-Remaining": "0",
                    "X-RateLimit-Reset": str(int(reset_time.timestamp())) if reset_time else "",
                    "Retry-After": str(int((reset_time - datetime.now()).total_seconds())) if reset_time else "60"
                })
                
                return response
            
            # Update rate limit counter
            update_rate_limit(key, period)
            
            # Call the original function
            result = f(*args, **kwargs)
            
            # Add rate limit headers to response
            if isinstance(result, tuple) and len(result) >= 2:
                response = result[0]
                if hasattr(response, "headers"):
                    response.headers.update({
                        "X-RateLimit-Limit": str(limit),
                        "X-RateLimit-Remaining": str(limit - current - 1),
                        "X-RateLimit-Reset": str(int(reset_time.timestamp())) if reset_time else ""
                    })
            
            return result
        
        return decorated
    
    return decorator

def csrf_protect():
    """
    Decorator to protect API endpoints from CSRF attacks.
    
    This decorator should be applied to all state-changing API endpoints
    (POST, PUT, PATCH, DELETE) to ensure they require and validate CSRF tokens.
    
    Returns:
        Decorator function
    """
    def decorator(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            # Skip CSRF check for GET, HEAD, OPTIONS requests
            if request.method in ('GET', 'HEAD', 'OPTIONS'):
                return f(*args, **kwargs)
            
            # Skip CSRF check if explicitly disabled for testing
            if current_app.config.get('TESTING') and current_app.config.get('CSRF_DISABLED'):
                return f(*args, **kwargs)
            
            # Import here to avoid circular imports
            from api.csrf import validate_csrf_token, CSRF_HEADER_NAME
            
            # Get token from header or form data
            token = request.headers.get(CSRF_HEADER_NAME)
            if not token and request.form:
                token = request.form.get('csrf_token')
            if not token and request.json:
                token = request.json.get('csrf_token')
            
            if not validate_csrf_token(token):
                # Log the CSRF failure
                logger.warning(
                    f"CSRF validation failed for {request.path}",
                    extra={
                        "event_type": "security",
                        "security_event": "csrf_failure",
                        "path": request.path,
                        "method": request.method,
                        "remote_addr": request.remote_addr
                    }
                )
                
                # Return CSRF validation error
                return jsonify_standard_response(
                    *create_error_response(
                        error="CSRF token validation failed",
                        status_code=403,
                        context="Security validation failed"
                    )
                )
            
            # Call the original function
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator

def csrf_exempt():
    """
    Decorator to exempt an API endpoint from CSRF protection.
    
    This should only be used for endpoints that don't change state or for
    external API endpoints that can't include CSRF tokens.
    
    Returns:
        Decorator function
    """
    def decorator(f):
        @functools.wraps(f)
        def decorated(*args, **kwargs):
            # Mark the view as exempt from CSRF protection
            setattr(f, '_csrf_exempt', True)
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator

def require_auth(f):
    """
    Decorator to require authentication for an endpoint.
    
    This is a temporary implementation for testing. In production, this would
    validate JWT tokens or other authentication mechanisms.
    
    Returns:
        Decorator function
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        # For testing, we'll always allow access
        # In production, this would check for valid authentication
        request.user = type('User', (), {'id': 'test-user-id'})
        return f(*args, **kwargs)
    
    return decorated

def require_service_role(f):
    """
    Decorator to require service role authentication for an endpoint.
    
    This is a temporary implementation for testing. In production, this would
    validate service role tokens or other authorization mechanisms.
    
    Returns:
        Decorator function
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        # For testing, we'll always allow access with service role
        # In production, this would check for valid service role authorization
        request.user = type('ServiceUser', (), {'id': 'service-role-user-id'})
        return f(*args, **kwargs)
    
    return decorated
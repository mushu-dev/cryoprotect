"""
CryoProtect Analyzer - Utility Functions

This module provides utility functions for the CryoProtect Analyzer API.

Rate Limiting:
------------
The CryoProtect API implements a comprehensive rate limiting system to ensure fair usage
and protect against abuse. Rate limits are applied at multiple levels:

1. Global Rate Limits:
   - Default limits applied to all endpoints (e.g., 1000 per day, 100 per hour, 20 per minute)
   - Configurable via environment variables or configuration files
   - Applied to all users regardless of authentication status

2. Endpoint-Specific Rate Limits:
   - Different limits for different endpoints based on resource intensity
   - Higher limits for common operations (e.g., GET /molecules)
   - Lower limits for computationally intensive operations (e.g., /rdkit/similarity)
   - Very low limits for batch operations to prevent abuse

3. Role-Based Rate Limits:
   - Different limits based on user roles (admin, premium, basic)
   - Higher limits for premium users and administrators
   - Configurable via environment variables or directly in the database

4. Identification Strategy:
   - User-based: Limits applied per user ID for authenticated requests
   - IP-based: Limits applied per IP address for unauthenticated requests
   - Hybrid: User ID if authenticated, otherwise IP address (default)

5. Rate Limit Response:
   - 429 Too Many Requests status code when limit is exceeded
   - Retry-After header indicating when to retry
   - JSON response with details about the limit and retry time
   - Rate limit headers in all responses (X-RateLimit-*)

6. Exemptions:
   - Critical endpoints exempt from rate limiting (e.g., health checks)
   - Documentation endpoints exempt to ensure accessibility
   - Configurable exemptions via configuration files

For detailed information on rate limiting, see README_RATE_LIMITING.md
"""

import os
import re
import json
import uuid
import logging
import traceback
from typing import Dict, Any, Tuple, Optional, List, Union
from datetime import datetime, timedelta
from flask import current_app, g, request, jsonify
from marshmallow import ValidationError
from functools import wraps

# Import Supabase client
from .supabase_client import create_client

def get_supabase_client():
    """
    Get the Supabase client.
    
    Returns:
        Supabase client
    """
    if 'supabase_client' not in g:
        g.supabase_client = create_client()
    return g.supabase_client

def release_supabase_connection():
    """
    Release the Supabase connection.
    This function is called when the request context is torn down.
    """
    client = g.pop('supabase_client', None)
    if client:
        # Close any open connections
        pass  # Supabase client doesn't have an explicit close method

def authenticate_user(email=None, password=None):
    """
    Authenticate a user with email and password or from request token.
    
    Args:
        email: User email (optional if authenticating from token)
        password: User password (optional if authenticating from token)
        
    Returns:
        User object if authentication is successful, None otherwise
    """
    # Import here to avoid circular imports
    from .jwt_auth import get_token_from_request, extract_user_from_token
    
    # Check if we have a token in the request
    if email is None and password is None:
        token = get_token_from_request()
        if token:
            try:
                user_data, user_id = extract_user_from_token(token)
                # Create a user object with the extracted data
                class User:
                    def __init__(self, data):
                        self.id = data.get('id')
                        self.email = data.get('email')
                        self.user_metadata = data.get('user_metadata', {})
                        self.app_metadata = data.get('app_metadata', {})
                        self.role = data.get('role', 'user')
                
                return User(user_data)
            except Exception as e:
                logging.error(f"Token authentication error: {str(e)}")
                return None
    
    # If email and password are provided, authenticate with Supabase
    if email and password:
        try:
            client = get_supabase_client()
            response = client.auth.sign_in_with_password({
                'email': email,
                'password': password
            })
            
            if response.user:
                return response.user
            
            return None
        except Exception as e:
            logging.error(f"Password authentication error: {str(e)}")
            return None
    
    return None

def token_required(f):
    """
    Decorator to require a valid token for API endpoints.
    This is a wrapper around jwt_required for backward compatibility.
    
    This decorator enforces authentication via JWT tokens. Protected endpoints
    will return a 401 Unauthorized response if no valid token is provided.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
        
    Authentication Flow:
        1. Client obtains a JWT token by authenticating with credentials
        2. Client includes token in Authorization header: "Bearer <token>"
        3. Server validates token and extracts user information
        4. If token is valid, the endpoint function is executed
        5. If token is invalid or missing, a 401 response is returned
    
    Example Usage:
        @token_required
        def protected_endpoint():
            # This function will only execute if a valid token is provided
            return {"message": "This is protected data"}
    """
    # Import here to avoid circular imports
    from .jwt_auth import jwt_required
    
    # Use the jwt_required decorator
    decorated = jwt_required(f)
    
    # Add OpenAPI documentation if flask-apispec is available
    try:
        decorated = document_auth_requirements(decorated, [{'Bearer': []}])
    except (NameError, AttributeError):
        # document_auth_requirements might not be available yet
        pass
    
    return decorated

def _handle_json_serialization(obj):
    """
    Handle JSON serialization for custom types.
    
    Args:
        obj: Object to serialize
        
    Returns:
        JSON serializable representation
    """
    if isinstance(obj, datetime):
        return obj.isoformat()
    elif isinstance(obj, uuid.UUID):
        return str(obj)
    raise TypeError(f"Type {type(obj)} not serializable")

def get_user_id():
    """
    Get the current user ID from the request context.
    
    Returns:
        User ID or None if not authenticated
    """
    # Import here to avoid circular imports
    from .jwt_auth import get_current_user
    
    # First check if user_id is already in the request context
    if hasattr(g, 'user_id') and g.user_id:
        return g.user_id
    
    # Try to get the current user
    user = get_current_user()
    if user:
        return user.get('id')
    
    return None

def handle_supabase_error(response):
    """
    Handle Supabase response errors.
    
    Args:
        response: Supabase response
        
    Returns:
        Tuple of (error_message, status_code) or (None, None) if no error
    """
    if response.error:
        error_data = response.error
        error_message = error_data.get('message', 'Unknown error')
        error_code = error_data.get('code', 500)
        
        # Map Supabase error codes to HTTP status codes
        status_code = 500  # Default to internal server error
        if error_code == 'PGRST301':  # Foreign key violation
            status_code = 400
        elif error_code == 'PGRST302':  # Unique violation
            status_code = 409
        elif error_code == 'PGRST401':  # Not authenticated
            status_code = 401
        elif error_code == 'PGRST403':  # Permission denied
            status_code = 403
        elif error_code == 'PGRST404':  # Not found
            status_code = 404
        
        return error_message, status_code
    
    return None, None

def generate_id():
    """
    Generate a unique ID.
    
    Returns:
        Unique ID string
    """
    return str(uuid.uuid4())

def format_timestamp(timestamp=None):
    """
    Format a timestamp in ISO format.
    
    Args:
        timestamp: Timestamp to format (default: current time)
        
    Returns:
        Formatted timestamp string
    """
    if timestamp is None:
        timestamp = datetime.now()
    return timestamp.isoformat()

def parse_timestamp(timestamp_str):
    """
    Parse a timestamp string.
    
    Args:
        timestamp_str: Timestamp string to parse
        
    Returns:
        Datetime object
    """
    return datetime.fromisoformat(timestamp_str)

def sanitize_input(input_str):
    """
    Sanitize input string to prevent injection attacks.
    
    Args:
        input_str: Input string to sanitize
        
    Returns:
        Sanitized string
    """
    if not input_str:
        return input_str
    
    # Remove any potentially dangerous characters
    sanitized = re.sub(r'[;<>&\'"\\]', '', input_str)
    return sanitized

def validate_smiles(smiles):
    """
    Validate a SMILES string.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    # Basic validation - check for valid characters
    if not smiles:
        return False
    
    # Check for balanced brackets
    bracket_count = 0
    for char in smiles:
        if char == '(':
            bracket_count += 1
        elif char == ')':
            bracket_count -= 1
        if bracket_count < 0:
            return False
    
    return bracket_count == 0

def validate_molecular_formula(formula):
    """
    Validate a molecular formula.
    
    Args:
        formula: Molecular formula to validate
        
    Returns:
        True if valid, False otherwise
    """
    # Basic validation - check for valid characters
    if not formula:
        return False
    
    # Check for valid format (e.g., C6H12O6)
    pattern = r'^([A-Z][a-z]?\d*)+$'
    return bool(re.match(pattern, formula))

def format_response(data=None, message=None, status="success", meta=None):
    """
    Format a standardized API response.
    
    Args:
        data: Response data
        message: Response message
        status: Response status ("success" or "error")
        meta: Metadata (e.g., pagination info)
        
    Returns:
        Formatted response dictionary
    """
    response = {
        "status": status,
        "timestamp": format_timestamp()
    }
    
    if data is not None:
        response["data"] = data
    
    if message:
        response["message"] = message
    
    if meta:
        response["meta"] = meta
    
    return response

def handle_error(error, context=None, log_level='error', return_response=False, status_code=None):
    """
    Standardized error handling utility for consistent logging and error reporting.
    
    Args:
        error: The exception or error message
        context: Additional context information (e.g., function name, operation being performed)
        log_level: Logging level ('error', 'warning', 'info', 'debug')
        return_response: Whether to return a structured response for API endpoints
        status_code: HTTP status code for API responses (default: inferred from error type)
    
    Returns:
        If return_response is True, returns a tuple (response_dict, status_code)
        Otherwise, returns None (just logs the error)
    
    Raises:
        The original exception if return_response is False and no status_code is provided
    """
    # Get logger
    logger = current_app.logger if 'current_app' in globals() and hasattr(current_app, 'logger') else logging.getLogger(__name__)
    
    # Format error message
    if isinstance(error, Exception):
        error_message = str(error)
        error_type = error.__class__.__name__
        exception = error
    else:
        error_message = str(error)
        error_type = "Error"
        exception = None
    
    # Add context to message if provided
    if context:
        log_message = f"{context}: {error_message}"
    else:
        log_message = error_message
    
    # Determine appropriate status code if not provided
    if status_code is None:
        if isinstance(error, ValueError) or isinstance(error, ValidationError):
            status_code = 400
        elif isinstance(error, PermissionError) or isinstance(error, ConnectionRefusedError):
            status_code = 403
        elif isinstance(error, FileNotFoundError) or isinstance(error, LookupError):
            status_code = 404
        elif isinstance(error, TimeoutError):
            status_code = 408
        elif isinstance(error, ConnectionError):
            status_code = 503
        else:
            status_code = 500
    
    # Log the error with appropriate level and include stack trace for errors
    if log_level == 'error':
        logger.error(log_message, exc_info=exception is not None)
    elif log_level == 'warning':
        logger.warning(log_message, exc_info=exception is not None)
    elif log_level == 'info':
        logger.info(log_message)
    elif log_level == 'debug':
        logger.debug(log_message)
    
    # Return structured response if requested
    if return_response:
        response = format_response(
            message=error_message,
            status="error",
            meta={
                "error_type": error_type,
                "context": context
            }
        )
        return response, status_code
    
    # Re-raise the exception if not returning a response and no status code provided
    if exception and not status_code:
        raise exception
    
    return None

def get_pagination_params(request):
    """
    Extract pagination parameters from request.
    
    Args:
        request: Flask request object
        
    Returns:
        Dictionary with page and per_page values
    """
    try:
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 20))
        
        # Validate and constrain values
        page = max(1, page)  # Minimum page is 1
        per_page = max(1, min(100, per_page))  # Between 1 and 100
        
        return {
            'page': page,
            'per_page': per_page,
            'offset': (page - 1) * per_page,
            'limit': per_page
        }
    except (ValueError, TypeError):
        # Default values if parsing fails
        return {
            'page': 1,
            'per_page': 20,
            'offset': 0,
            'limit': 20
        }

def get_sort_params(request, allowed_fields=None):
    """
    Extract sorting parameters from request.
    
    Args:
        request: Flask request object
        allowed_fields: List of fields that are allowed to be sorted
        
    Returns:
        List of (field, direction) tuples
    """
    sort_param = request.args.get('sort', '')
    if not sort_param:
        return []
    
    sort_fields = []
    for field in sort_param.split(','):
        field = field.strip()
        if field.startswith('-'):
            direction = 'desc'
            field = field[1:]
        else:
            direction = 'asc'
        
        # Validate field if allowed_fields is provided
        if allowed_fields and field not in allowed_fields:
            continue
        
        sort_fields.append((field, direction))
    
    return sort_fields

def get_filter_params(request, allowed_fields=None):
    """
    Extract filter parameters from request.
    
    Args:
        request: Flask request object
        allowed_fields: List of fields that are allowed to be filtered
        
    Returns:
        Dictionary of filter conditions
    """
    filters = {}
    for key, value in request.args.items():
        # Skip pagination and sorting params
        if key in ['page', 'per_page', 'sort']:
            continue
        
        # Validate field if allowed_fields is provided
        if allowed_fields and key not in allowed_fields:
            continue
        
        filters[key] = value
    
    return filters

def serialize_datetime(obj):
    """
    JSON serializer for datetime objects.
    
    Args:
        obj: Object to serialize
        
    Returns:
        Serialized object
    """
    if isinstance(obj, datetime):
        return obj.isoformat()
    raise TypeError(f"Type {type(obj)} not serializable")

def to_json(data):
    """
    Convert data to JSON string.
    
    Args:
        data: Data to convert
        
    Returns:
        JSON string
    """
    return json.dumps(data, default=serialize_datetime)

def from_json(json_str):
    """
    Convert JSON string to data.
    
    Args:
        json_str: JSON string to convert
        
    Returns:
        Parsed data
    """
    return json.loads(json_str)

def marshal_with(fields):
    """
    Decorator to marshal response with the given fields.
    
    Args:
        fields: Field definitions for marshalling
        
    Returns:
        Decorated function
    """
    from flask_restful import marshal
    
    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            resp = f(*args, **kwargs)
            
            if isinstance(resp, tuple):
                data, code = resp
                return marshal(data, fields), code
            else:
                return marshal(resp, fields)
        return wrapper
    return decorator

def get_file_extension(filename):
    """
    Get the file extension from a filename.
    
    Args:
        filename: Filename to extract extension from
        
    Returns:
        File extension (lowercase, without the dot)
    """
    if not filename:
        return ''
    
    parts = filename.rsplit('.', 1)
    if len(parts) > 1:
        return parts[1].lower()
    return ''

def is_valid_uuid(uuid_str):
    """
    Check if a string is a valid UUID.
    
    Args:
        uuid_str: String to check
        
    Returns:
        True if valid UUID, False otherwise
    """
    try:
        uuid_obj = uuid.UUID(uuid_str)
        return str(uuid_obj) == uuid_str
    except (ValueError, AttributeError, TypeError):
        return False

def get_env_var(name, default=None):
    """
    Get an environment variable.
    
    Args:
        name: Name of the environment variable
        default: Default value if not found
        
    Returns:
        Value of the environment variable or default
    """
    return os.environ.get(name, default)

def get_config_value(name, default=None):
    """
    Get a configuration value.
    
    Args:
        name: Name of the configuration value
        default: Default value if not found
        
    Returns:
        Configuration value or default
    """
    if 'current_app' in globals() and hasattr(current_app, 'config'):
        return current_app.config.get(name, default)
    return default

# === API Documentation Utilities ===

try:
    from apispec import APISpec
    from apispec.ext.marshmallow import MarshmallowPlugin
    from flask_apispec.extension import FlaskApiSpec
except ImportError:
    APISpec = None
    MarshmallowPlugin = None
    FlaskApiSpec = None

def get_apispec() -> "APISpec":
    """
    Create and return a configured APISpec instance for OpenAPI documentation.

    Returns:
        APISpec: Configured APISpec object (OpenAPI 3.0.2, Marshmallow/Flask plugins)
        
    This function configures the OpenAPI specification with:
    - Basic API information (title, version, description)
    - Security schemes (JWT Bearer token, API key)
    - Marshmallow plugin for schema conversion
    """
    if APISpec is None or MarshmallowPlugin is None:
        raise ImportError("apispec and flask-apispec must be installed to use API documentation utilities.")

    spec = APISpec(
        title="CryoProtect v2 API",
        version="1.0.0",
        openapi_version="3.0.2",
        plugins=[MarshmallowPlugin()],
        info={
            "description": "OpenAPI documentation for the CryoProtect v2 API.",
            "contact": {"email": "support@cryoprotect.com"},
            "license": {"name": "MIT License"}
        }
    )
    
    # Register authentication schemes
    spec.components.security_scheme('Bearer', {
        'type': 'http',
        'scheme': 'bearer',
        'bearerFormat': 'JWT',
        'description': 'JWT token authentication. Provide the token in the format: Bearer <token>'
    })
    
    spec.components.security_scheme('ApiKey', {
        'type': 'apiKey',
        'in': 'header',
        'name': 'X-API-Key',
        'description': 'API key authentication (legacy support)'
    })
    
    return spec

# Placeholder for future documentation helpers (e.g., resource registration, schema helpers)

# === Authentication Documentation ===

def register_auth_documentation(spec):
    """
    Register authentication documentation with the APISpec instance.
    
    This function adds security schemes and requirements to the OpenAPI specification.
    
    Args:
        spec (APISpec): The APISpec instance to register with
    """
    # Define Bearer token security scheme
    spec.components.security_scheme('Bearer', {
        'type': 'http',
        'scheme': 'bearer',
        'bearerFormat': 'JWT',
        'description': 'JWT token authentication. Provide the token in the format: Bearer <token>'
    })
    
    # Define API key security scheme (for legacy support)
    spec.components.security_scheme('ApiKey', {
        'type': 'apiKey',
        'in': 'header',
        'name': 'X-API-Key',
        'description': 'API key authentication (legacy support)'
    })

def document_auth_requirements(view_func, security_requirements):
    """
    Decorator to document authentication requirements for an endpoint.
    
    This function adds OpenAPI security requirement documentation to Flask view functions.
    It's used to specify which authentication mechanisms are required for each endpoint,
    making this information available in the generated API documentation.
    
    The security requirements are added to the view function's __apispec__ attribute,
    which is used by flask-apispec when generating the OpenAPI specification.
    
    Args:
        view_func: The view function to decorate
        security_requirements: List of security requirements in OpenAPI format.
            Example: [{'Bearer': []}] for JWT authentication or
                    [{'ApiKey': []}] for API key authentication
    
    Returns:
        Decorated function with OpenAPI security documentation
        
    Example Usage:
        @app.route('/api/protected')
        @document_auth_requirements([{'Bearer': []}])
        def protected_endpoint():
            # This endpoint will be documented as requiring Bearer authentication
            return {"message": "Protected data"}
            
    Note:
        This function does not enforce authentication itself; it only documents
        the requirements. Use with token_required or jwt_required to enforce
        authentication.
    """
    if hasattr(view_func, '__apispec__'):
        if 'security' not in view_func.__apispec__:
            view_func.__apispec__['security'] = security_requirements
    else:
        view_func.__apispec__ = {'security': security_requirements}
    return view_func

def jwt_auth_required_with_doc(view_func):
    """
    Decorator that combines token_required with documentation.
    
    This is a convenience decorator that both enforces JWT authentication
    and documents the authentication requirements in the OpenAPI specification.
    It combines the functionality of token_required and document_auth_requirements
    in a single decorator, making it easier to create protected endpoints
    that are properly documented.
    
    When applied to a Flask route function, this decorator will:
    1. Enforce JWT token authentication (reject requests without valid tokens)
    2. Document the endpoint as requiring Bearer authentication in the OpenAPI spec
    3. Make the endpoint visible in the Swagger UI with the authentication requirements
    
    Args:
        view_func: The view function to decorate
        
    Returns:
        Decorated function that enforces and documents JWT authentication
        
    Example Usage:
        @app.route('/api/protected-resource')
        @jwt_auth_required_with_doc
        def get_protected_resource():
            # This endpoint requires a valid JWT token
            # and is documented as requiring Bearer authentication
            return jsonify({"data": "Protected resource"})
    """
    # Apply token_required decorator
    decorated_func = token_required(view_func)
    
    # Add documentation
    return document_auth_requirements(decorated_func, [{'Bearer': []}])
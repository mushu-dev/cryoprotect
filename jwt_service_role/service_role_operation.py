"""
CryoProtect - Service Role Operation Decorator

This module provides a decorator for API functions that require service role
authentication, with automatic scope verification and audit logging.
"""

import time
import logging
import functools
import threading
from typing import Dict, List, Optional, Union, Tuple, Any, Callable

from flask import request, g, jsonify, current_app

from .token_manager import ServiceRoleTokenManager

# Configure logging
logger = logging.getLogger(__name__)

# Thread-local storage for request tracking
request_local = threading.local()

def get_token_from_request() -> Optional[str]:
    """
    Extract the JWT token from the request.
    Checks Authorization header.
    
    Returns:
        The token or None if not found
    """
    # Check Authorization header
    auth_header = request.headers.get('Authorization', '')
    if auth_header.startswith('Bearer '):
        return auth_header.split(' ')[1]
    
    return None

def service_role_operation(required_scope: Union[str, List[str]], audit_operation: str = None):
    """
    Decorator for API functions that require service role authentication.
    
    This decorator:
    1. Extracts the service role token from the request
    2. Validates the token
    3. Verifies the token has the required scope
    4. Records audit logs for the operation
    5. Provides the authenticated client ID to the decorated function
    
    Usage:
        @app.route('/api/admin/users', methods=['GET'])
        @service_role_operation(required_scope='users:admin')
        def get_all_users(client_id):
            # This function will only execute if:
            # 1. A valid service role token is provided
            # 2. The token has the 'users:admin' scope
            # client_id is passed automatically from the token
            return {"users": [...]}
    
    Args:
        required_scope: Scope or list of scopes required for the operation
        audit_operation: Name of the operation for audit logs (defaults to function name)
        
    Returns:
        Decorated function
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Start timing for performance logging
            start_time = time.time()
            
            # Extract token
            token = get_token_from_request()
            if not token:
                logger.warning("Service role token missing in request")
                return jsonify({'message': 'Service role token required', 'error': 'unauthorized'}), 401
            
            # Create token manager
            token_manager = ServiceRoleTokenManager()
            
            # Validate token
            is_valid, token_data, error_message = token_manager.validate_token(token)
            if not is_valid:
                logger.warning(f"Invalid service role token: {error_message}")
                return jsonify({'message': f'Invalid service role token: {error_message}', 'error': 'unauthorized'}), 401
            
            # Verify token type
            if token_data.get('type') != 'service_role':
                logger.warning("Token is not a service role token")
                return jsonify({'message': 'Service role token required', 'error': 'unauthorized'}), 401
            
            # Get client ID from token
            client_id = token_data.get('sub')
            if not client_id:
                logger.warning("Token missing client ID")
                return jsonify({'message': 'Invalid service role token: missing client ID', 'error': 'unauthorized'}), 401
            
            # Convert single scope to list
            required_scopes = required_scope if isinstance(required_scope, list) else [required_scope]
            
            # Verify scopes
            has_required_scope = False
            for scope in required_scopes:
                if token_manager.check_scope(token_data, scope):
                    has_required_scope = True
                    break
            
            if not has_required_scope:
                logger.warning(f"Client {client_id} lacks required scope: {required_scope}")
                return jsonify({
                    'message': f'Insufficient permissions. Required scope: {required_scope}',
                    'error': 'forbidden'
                }), 403
            
            # Track request in thread-local storage for audit logging
            request_local.client_id = client_id
            request_local.token_data = token_data
            request_local.start_time = start_time
            
            # Store client ID in flask.g for potential use in the view function
            g.service_role_client_id = client_id
            g.service_role_scopes = token_data.get('scopes', [])
            
            try:
                # Call the decorated function with client_id as first argument
                result = func(client_id, *args, **kwargs)
                
                # Log the operation for audit purposes
                operation = audit_operation or func.__name__
                elapsed_time = time.time() - start_time
                logger.info(
                    f"Service role operation: {operation} | "
                    f"Client: {client_id} | "
                    f"Status: success | "
                    f"Duration: {elapsed_time:.3f}s"
                )
                
                return result
            except Exception as e:
                # Log the operation failure
                operation = audit_operation or func.__name__
                elapsed_time = time.time() - start_time
                logger.error(
                    f"Service role operation: {operation} | "
                    f"Client: {client_id} | "
                    f"Status: error | "
                    f"Duration: {elapsed_time:.3f}s | "
                    f"Error: {str(e)}"
                )
                raise
            finally:
                # Clean up thread-local storage
                if hasattr(request_local, 'client_id'):
                    del request_local.client_id
                if hasattr(request_local, 'token_data'):
                    del request_local.token_data
                if hasattr(request_local, 'start_time'):
                    del request_local.start_time
        
        return wrapper
    
    return decorator

def get_current_service_role_client():
    """
    Get the current service role client ID from the request context.
    
    Returns:
        Client ID or None if not authenticated as service role
    """
    return getattr(g, 'service_role_client_id', None)

def get_current_service_role_scopes():
    """
    Get the current service role scopes from the request context.
    
    Returns:
        List of scopes or empty list if not authenticated as service role
    """
    return getattr(g, 'service_role_scopes', [])
#!/usr/bin/env python3
"""
Service Role Operation decorator for API endpoints.

This module provides a decorator for securing API endpoints
with service role authentication using JWT tokens.
"""

import os
import time
import json
import logging
import functools
from typing import Dict, Any, Optional, Callable, List, Union
from jwt.exceptions import PyJWTError

from flask import request, jsonify, g, current_app
from werkzeug.exceptions import Unauthorized, Forbidden

from .token_manager import validate_service_token

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def service_role_required(
    scopes: Optional[Union[str, List[str]]] = None,
    allow_services: Optional[List[str]] = None
):
    """
    Decorator to require service role authentication for an API endpoint.
    
    Args:
        scopes: Required scope or list of scopes (any match is sufficient)
        allow_services: List of allowed service names
        
    Returns:
        Decorated function
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            
            try:
                # Get Bearer token from Authorization header
                auth_header = request.headers.get('Authorization')
                
                if not auth_header or not auth_header.startswith('Bearer '):
                    raise Unauthorized("Missing or invalid Authorization header")
                
                token = auth_header.split(' ')[1]
                
                # Validate token
                try:
                    claims = validate_service_token(token)
                except PyJWTError as e:
                    logger.warning(f"Invalid service token: {str(e)}")
                    raise Unauthorized(f"Invalid service token: {str(e)}")
                
                # Check service name if restricted
                if allow_services:
                    service_name = claims.get('client')
                    if not service_name or service_name not in allow_services:
                        logger.warning(f"Unauthorized service: {service_name}")
                        raise Forbidden(f"Service '{service_name}' is not authorized")
                
                # Check scopes if required
                if scopes:
                    # Convert to list for consistent handling
                    required_scopes = [scopes] if isinstance(scopes, str) else scopes
                    token_scopes = claims.get('scopes', [])
                    
                    # Check if user has any of the required scopes
                    # or if they have the wildcard scope
                    if (not any(scope in token_scopes for scope in required_scopes) and
                        '*' not in token_scopes):
                        logger.warning(
                            f"Missing required scopes. Required: {required_scopes}, "
                            f"Token scopes: {token_scopes}"
                        )
                        raise Forbidden("Missing required scopes")
                
                # Store claims in Flask g object for later use
                g.service_claims = claims
                
                # Call the decorated function
                result = func(*args, **kwargs)
                
                # Log access with timing
                elapsed = time.time() - start_time
                logger.info(
                    f"Service role access: {claims.get('client')} - "
                    f"{request.method} {request.path} - {elapsed:.3f}s"
                )
                
                return result
            
            except (Unauthorized, Forbidden) as e:
                # Log the error
                logger.warning(
                    f"Service role access denied: {request.method} {request.path} - {str(e)}"
                )
                
                # Return appropriate status code and error message
                status_code = e.code if hasattr(e, 'code') else 401
                
                return jsonify({
                    'error': str(e),
                    'status': 'error',
                    'code': status_code
                }), status_code
            
            except Exception as e:
                # Log unexpected errors
                logger.error(
                    f"Unexpected error in service role authentication: {str(e)}"
                )
                
                # Return 500 for unexpected errors
                return jsonify({
                    'error': 'Internal server error',
                    'status': 'error',
                    'code': 500
                }), 500
        
        # Add authentication_required attribute for introspection
        wrapper.service_role_required = True
        
        return wrapper
    
    return decorator

def get_current_service() -> Optional[Dict[str, Any]]:
    """
    Get the current authenticated service from Flask g object.
    
    Returns:
        Dict with service claims or None if not authenticated
    """
    return getattr(g, 'service_claims', None)

def require_service_scope(scope: str) -> bool:
    """
    Check if the current service has a required scope.
    
    Args:
        scope: Required scope
        
    Returns:
        True if the service has the scope, False otherwise
    """
    claims = get_current_service()
    
    if not claims:
        return False
    
    scopes = claims.get('scopes', [])
    return scope in scopes or '*' in scopes

# Convenience decorator for common operations
def admin_service_required(func: Callable) -> Callable:
    """
    Decorator to require admin service role authentication.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    return service_role_required(scopes='admin')(func)

def data_service_required(func: Callable) -> Callable:
    """
    Decorator to require data service role authentication.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    return service_role_required(scopes='data')(func)

def analytics_service_required(func: Callable) -> Callable:
    """
    Decorator to require analytics service role authentication.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    return service_role_required(scopes='analytics')(func)

# Flask extension class for easier integration
class ServiceRoleAuth:
    """
    Flask extension for service role authentication.
    
    Usage:
        app = Flask(__name__)
        service_auth = ServiceRoleAuth(app)
        
        # Or with factory pattern
        service_auth = ServiceRoleAuth()
        service_auth.init_app(app)
    """
    
    def __init__(self, app=None):
        """
        Initialize the extension.
        
        Args:
            app: Optional Flask application instance
        """
        self.app = app
        if app is not None:
            self.init_app(app)
    
    def init_app(self, app):
        """
        Initialize the extension with a Flask application.
        
        Args:
            app: Flask application instance
        """
        # Register extension with app
        app.extensions = getattr(app, 'extensions', {})
        app.extensions['service_role_auth'] = self
        
        # Add convenience methods to app context
        app.service_role_required = service_role_required
        app.admin_service_required = admin_service_required
        app.data_service_required = data_service_required
        app.analytics_service_required = analytics_service_required
        app.get_current_service = get_current_service
        app.require_service_scope = require_service_scope
        
        # Optional configuration
        self.config = app.config.get('SERVICE_ROLE_AUTH', {})
        
        # Log initialization
        logger.info("ServiceRoleAuth extension initialized")

# Example usage in a Flask application
"""
from flask import Flask, jsonify
from auth.service_role_operation import ServiceRoleAuth, service_role_required

app = Flask(__name__)
service_auth = ServiceRoleAuth(app)

@app.route('/api/admin/config', methods=['GET'])
@service_role_required(scopes=['admin'])
def get_admin_config():
    # This endpoint is protected and only accessible with a valid service token
    # that has the 'admin' scope
    return jsonify({
        'status': 'success',
        'message': 'Admin configuration retrieved'
    })

@app.route('/api/data/query', methods=['POST'])
@service_role_required(scopes=['data', 'analytics'])
def query_data():
    # This endpoint is accessible with 'data' OR 'analytics' scope
    return jsonify({
        'status': 'success',
        'message': 'Data query executed'
    })

@app.route('/api/health', methods=['GET'])
def health_check():
    # Public endpoint, no authentication required
    return jsonify({
        'status': 'ok'
    })

if __name__ == '__main__':
    app.run(debug=True)
"""
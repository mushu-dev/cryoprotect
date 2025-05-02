"""
CSRF Protection Module for CryoProtect API

This module provides CSRF protection for all state-changing API endpoints
(POST, PUT, PATCH, DELETE) as required by the security audit remediation plan.
"""

import os
import secrets
from functools import wraps
from flask import request, abort, current_app, g, session, jsonify
from datetime import datetime, timedelta

# Constants
CSRF_HEADER_NAME = 'X-CSRF-Token'
CSRF_COOKIE_NAME = 'csrf_token'
CSRF_TOKEN_LENGTH = 32  # 256 bits
CSRF_TOKEN_EXPIRY = 3600  # 1 hour in seconds

def generate_csrf_token():
    """
    Generate a secure CSRF token using the secrets module.
    
    Returns:
        str: A secure random token
    """
    return secrets.token_hex(CSRF_TOKEN_LENGTH // 2)  # token_hex returns 2 hex chars per byte

def get_csrf_token():
    """
    Get the current CSRF token from the session or generate a new one.
    
    Returns:
        str: The CSRF token
    """
    if 'csrf_token' not in session:
        session['csrf_token'] = generate_csrf_token()
        session['csrf_token_expiry'] = (datetime.utcnow() + timedelta(seconds=CSRF_TOKEN_EXPIRY)).timestamp()
    
    # Check if token is expired
    if datetime.utcnow().timestamp() > session.get('csrf_token_expiry', 0):
        session['csrf_token'] = generate_csrf_token()
        session['csrf_token_expiry'] = (datetime.utcnow() + timedelta(seconds=CSRF_TOKEN_EXPIRY)).timestamp()
    
    return session['csrf_token']

def validate_csrf_token(token):
    """
    Validate the provided CSRF token against the one in the session.
    
    Args:
        token (str): The CSRF token to validate
        
    Returns:
        bool: True if the token is valid, False otherwise
    """
    if not token or not session.get('csrf_token'):
        return False
    
    # Check if token is expired
    if datetime.utcnow().timestamp() > session.get('csrf_token_expiry', 0):
        return False
    
    # Use constant-time comparison to prevent timing attacks
    return secrets.compare_digest(token, session.get('csrf_token', ''))

def csrf_protect():
    """
    Decorator to protect routes from CSRF attacks.
    
    This decorator should be applied to all state-changing API endpoints
    (POST, PUT, PATCH, DELETE).
    
    Returns:
        function: The decorated function
    """
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            # Skip CSRF check for GET, HEAD, OPTIONS requests
            if request.method in ('GET', 'HEAD', 'OPTIONS'):
                return f(*args, **kwargs)
            
            # Skip CSRF check if explicitly disabled for testing
            if current_app.config.get('TESTING') and current_app.config.get('CSRF_DISABLED'):
                return f(*args, **kwargs)
            
            # Get token from header or form data
            token = request.headers.get(CSRF_HEADER_NAME)
            if not token and request.form:
                token = request.form.get('csrf_token')
            if not token and request.json:
                token = request.json.get('csrf_token')
            
            if not validate_csrf_token(token):
                current_app.logger.warning(
                    f"CSRF validation failed for {request.path}",
                    extra={
                        "event_type": "security",
                        "security_event": "csrf_failure",
                        "path": request.path,
                        "method": request.method,
                        "remote_addr": request.remote_addr
                    }
                )
                abort(403, description="CSRF token validation failed")
            
            return f(*args, **kwargs)
        return decorated_function
    return decorator

def csrf_exempt():
    """
    Decorator to exempt a route from CSRF protection.
    
    This should only be used for routes that don't change state or for
    external API endpoints that can't include CSRF tokens.
    
    Returns:
        function: The decorated function
    """
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            # Mark the view as exempt from CSRF protection
            setattr(f, '_csrf_exempt', True)
            return f(*args, **kwargs)
        return decorated_function
    return decorator

def is_csrf_exempt(view_func):
    """
    Check if a view function is exempt from CSRF protection.
    
    Args:
        view_func (function): The view function to check
        
    Returns:
        bool: True if the view is exempt, False otherwise
    """
    return getattr(view_func, '_csrf_exempt', False)

def set_csrf_cookie(response):
    """
    Set the CSRF token in a cookie for JavaScript access.
    
    Args:
        response (Response): The Flask response object
        
    Returns:
        Response: The modified response with the CSRF cookie
    """
    token = get_csrf_token()
    
    # Import secure cookie function if available
    try:
        from api.session_security import set_secure_cookie
        return set_secure_cookie(
            response,
            CSRF_COOKIE_NAME,
            token,
            max_age=CSRF_TOKEN_EXPIRY,
            secure=current_app.config.get('SECURE_COOKIES', True),
            httponly=False,  # Must be accessible to JavaScript
            samesite=current_app.config.get('SAME_SITE_COOKIES', 'Lax')
        )
    except ImportError:
        # Fallback to standard cookie setting
        response.set_cookie(
            CSRF_COOKIE_NAME,
            token,
            max_age=CSRF_TOKEN_EXPIRY,
            secure=current_app.config.get('SECURE_COOKIES', True),
            httponly=False,  # Must be accessible to JavaScript
            samesite=current_app.config.get('SAME_SITE_COOKIES', 'Lax')
        )
        return response

def init_csrf(app):
    """
    Initialize CSRF protection for a Flask application.
    
    This function sets up the necessary hooks to apply CSRF protection
    to all state-changing API endpoints.
    
    Args:
        app (Flask): The Flask application
    """
    # Set up CSRF configuration
    app.config.setdefault('CSRF_ENABLED', True)
    app.config.setdefault('CSRF_TOKEN_LENGTH', CSRF_TOKEN_LENGTH)
    app.config.setdefault('CSRF_TOKEN_EXPIRY', CSRF_TOKEN_EXPIRY)
    
    # Add CSRF token to all responses
    @app.after_request
    def add_csrf_token(response):
        # Skip for static files
        if request.path.startswith('/static/'):
            return response
        
        # Skip for non-HTML/JSON responses
        if not response.mimetype in ('text/html', 'application/json'):
            return response
        
        # Set CSRF token in cookie
        return set_csrf_cookie(response)
    
    # Add CSRF protection to all state-changing routes
    @app.before_request
    def csrf_protect_request():
        # Skip CSRF check for GET, HEAD, OPTIONS requests
        if request.method in ('GET', 'HEAD', 'OPTIONS'):
            return
        
        # Skip CSRF check if explicitly disabled for testing
        if app.config.get('TESTING') and app.config.get('CSRF_DISABLED'):
            return
        
        # Skip for static files
        if request.path.startswith('/static/'):
            return
        
        # Skip for CSRF exempt views
        if is_csrf_exempt(app.view_functions.get(request.endpoint)):
            return
        
        # Get token from header or form data
        token = request.headers.get(CSRF_HEADER_NAME)
        if not token and request.form:
            token = request.form.get('csrf_token')
        if not token and request.json:
            token = request.json.get('csrf_token')
        
        if not validate_csrf_token(token):
            app.logger.warning(
                f"CSRF validation failed for {request.path}",
                extra={
                    "event_type": "security",
                    "security_event": "csrf_failure",
                    "path": request.path,
                    "method": request.method,
                    "remote_addr": request.remote_addr
                }
            )
            abort(403, description="CSRF token validation failed")
    
    # Add route to get CSRF token
    @app.route('/api/v1/csrf-token', methods=['GET'])
    def get_csrf_token_route():
        """API endpoint to get a new CSRF token."""
        token = get_csrf_token()
        response = jsonify({'csrf_token': token})
        return set_csrf_cookie(response)
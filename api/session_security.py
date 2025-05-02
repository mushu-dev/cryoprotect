"""
Session Security Enhancement Module for CryoProtect API

This module provides enhanced security for session management, including:
- Secure cookie attributes (Secure, HttpOnly, SameSite)
- Session cookie rotation on security-sensitive events
- Session expiration management

Implements requirements from the security audit remediation plan.
"""

import os
import secrets
from datetime import datetime, timedelta
from flask import session, request, current_app, g
from functools import wraps

# Import authentication configuration
from auth_config import (
    SESSION_TIMEOUT, REFRESH_TOKEN_ROTATION,
    SECURE_COOKIES, HTTP_ONLY_COOKIES, SAME_SITE_COOKIES
)

def rotate_session():
    """
    Rotate the current session to prevent session fixation attacks.
    
    This creates a new session with the same data but a new session ID.
    """
    if not session:
        return
    
    # Store current session data
    session_data = dict(session)
    
    # Clear current session
    session.clear()
    
    # Create new session with same data
    for key, value in session_data.items():
        session[key] = value
    
    # Update session creation time
    session['created_at'] = datetime.utcnow().timestamp()
    
    # Log session rotation
    if hasattr(current_app, 'logger'):
        current_app.logger.info(
            "Session rotated",
            extra={
                "event_type": "security",
                "security_event": "session_rotation",
                "user_id": session.get('user_id', 'anonymous'),
                "remote_addr": request.remote_addr
            }
        )

def set_secure_cookie(response, key, value, max_age=None, expires=None, path='/', 
                     domain=None, secure=None, httponly=None, samesite=None):
    """
    Set a cookie with secure attributes.
    
    Args:
        response: Flask response object
        key: Cookie name
        value: Cookie value
        max_age: Cookie max age in seconds
        expires: Cookie expiration time
        path: Cookie path
        domain: Cookie domain
        secure: Whether the cookie should only be sent over HTTPS
        httponly: Whether the cookie should be accessible only via HTTP(S)
        samesite: SameSite attribute (Lax, Strict, None)
        
    Returns:
        Flask response with cookie set
    """
    # Use configured defaults if not specified
    if secure is None:
        secure = SECURE_COOKIES
    if httponly is None:
        httponly = HTTP_ONLY_COOKIES
    if samesite is None:
        samesite = SAME_SITE_COOKIES
    
    # Set the cookie with secure attributes
    response.set_cookie(
        key,
        value,
        max_age=max_age,
        expires=expires,
        path=path,
        domain=domain,
        secure=secure,
        httponly=httponly,
        samesite=samesite
    )
    
    return response

def session_expiry_required(f):
    """
    Decorator to check if the session has expired and clear it if so.
    
    This should be applied to routes that require an active session.
    """
    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Check if session exists and has a creation time
        if session and 'created_at' in session:
            created_at = session.get('created_at')
            now = datetime.utcnow().timestamp()
            
            # Check if session has expired
            if now - created_at > SESSION_TIMEOUT:
                # Clear session
                session.clear()
                
                # Log session expiry
                if hasattr(current_app, 'logger'):
                    current_app.logger.info(
                        "Session expired",
                        extra={
                            "event_type": "security",
                            "security_event": "session_expiry",
                            "remote_addr": request.remote_addr
                        }
                    )
        
        return f(*args, **kwargs)
    return decorated_function

def init_session_security(app):
    """
    Initialize session security for a Flask application.
    
    This sets up secure session configuration and hooks for session management.
    
    Args:
        app: Flask application
    """
    # Configure session security
    app.config.update({
        'SESSION_COOKIE_SECURE': SECURE_COOKIES,
        'SESSION_COOKIE_HTTPONLY': HTTP_ONLY_COOKIES,
        'SESSION_COOKIE_SAMESITE': SAME_SITE_COOKIES,
        'PERMANENT_SESSION_LIFETIME': timedelta(seconds=SESSION_TIMEOUT),
        'SESSION_REFRESH_EACH_REQUEST': True
    })
    
    # Set session to permanent but respect the configured timeout
    @app.before_request
    def make_session_permanent():
        session.permanent = True
        
        # Initialize session creation time if not set
        if 'created_at' not in session:
            session['created_at'] = datetime.utcnow().timestamp()
    
    # Check for privilege escalation that would require session rotation
    @app.before_request
    def check_privilege_changes():
        if not session or not g.get('user'):
            return
            
        # Check if user role has changed since last request
        current_role = g.user.get('role')
        session_role = session.get('user_role')
        
        if session_role and current_role and current_role != session_role:
            # Role has changed, rotate session
            rotate_session()
            
            # Update session role
            session['user_role'] = current_role
            
            # Log privilege change
            if hasattr(current_app, 'logger'):
                current_app.logger.warning(
                    f"User role changed from {session_role} to {current_role}",
                    extra={
                        "event_type": "security",
                        "security_event": "privilege_change",
                        "user_id": g.user.get('id', 'unknown'),
                        "old_role": session_role,
                        "new_role": current_role,
                        "remote_addr": request.remote_addr
                    }
                )
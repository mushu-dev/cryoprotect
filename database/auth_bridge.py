"""
Authentication Bridge for Convex and Flask

This module provides utilities to synchronize authentication between
our Flask API and Convex database by generating and validating JWTs
that both systems can use for authentication.

Enhanced with improved Convex compatibility, bidirectional auth support,
and comprehensive role-based access control.
"""

import os
import time
import json
import jwt
import logging
from datetime import datetime, timedelta
from functools import wraps
from flask import request, jsonify, g, current_app
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
logger = logging.getLogger(__name__)

# JWT configuration
JWT_SECRET = os.environ.get('JWT_SECRET', 'replace-with-your-secret-key')
JWT_ALGORITHM = 'HS256'
JWT_EXPIRY_SECONDS = 24 * 60 * 60  # 24 hours

# Convex configuration
CONVEX_URL = os.environ.get('CONVEX_URL', 'https://upbeat-parrot-866.convex.cloud')
CONVEX_DEPLOYMENT_KEY = os.environ.get('CONVEX_DEPLOYMENT_KEY', '')

class AuthBridge:
    """
    Authentication bridge for Convex integration.
    Handles JWT generation and validation compatible with both systems.
    """
    
    def __init__(self, app=None, secret_key=None, algorithm=None, expiry=None):
        """
        Initialize the auth bridge.
        
        Args:
            app: Optional Flask app to initialize with
            secret_key: Secret key for JWT encoding/decoding
            algorithm: Algorithm for JWT encoding/decoding
            expiry: Token expiry time in seconds
        """
        self.app = app
        self.secret_key = secret_key or JWT_SECRET
        self.algorithm = algorithm or JWT_ALGORITHM
        self.expiry_seconds = expiry or JWT_EXPIRY_SECONDS
        self.convex_url = CONVEX_URL
        self.convex_deployment_key = CONVEX_DEPLOYMENT_KEY
        
        if app is not None:
            self.init_app(app)
            
        logger.info("AuthBridge initialized with algorithm %s and expiry %s seconds", 
                   self.algorithm, self.expiry_seconds)
    
    def init_app(self, app):
        """
        Initialize with a Flask app.
        
        Args:
            app: Flask app instance
        """
        self.app = app
        app.auth_bridge = self
        
        # Register middleware
        app.before_request(self.jwt_middleware())
    
    def generate_token(self, user_id, user_data=None):
        """
        Generate a JWT token for a user that is compatible with Convex.
        
        Args:
            user_id: User ID to include in the token
            user_data: Additional user data to include
            
        Returns:
            str: JWT token
        """
        now = int(time.time())
        
        # Create payload with claims required by Convex
        payload = {
            'sub': str(user_id),
            'iat': now,
            'exp': now + self.expiry_seconds,
            # Convex-specific fields
            'iss': os.environ.get('JWT_ISSUER', 'cryoprotect-api'),
            'aud': self.convex_url
        }
        
        # Add user data if provided
        if user_data:
            # Include only safe fields
            safe_fields = ['name', 'email', 'role', 'permissions']
            for field in safe_fields:
                if field in user_data:
                    payload[field] = user_data[field]
        
        # Add Convex-specific structure for authentication
        payload['convex'] = {
            'identity': str(user_id),
            'tokenIdentifier': user_data.get('email') if user_data else str(user_id),
            'traits': {
                'role': user_data.get('role', 'user') if user_data else 'user',
                'name': user_data.get('name', '') if user_data else '',
            }
        }
        
        # Generate token
        token = jwt.encode(payload, self.secret_key, algorithm=self.algorithm)
        
        logger.debug("Generated token for user %s with expiry %s", user_id, now + self.expiry_seconds)
        return token
    
    def validate_token(self, token):
        """
        Validate a JWT token.
        
        Args:
            token: JWT token to validate
            
        Returns:
            dict: Decoded payload if valid, None otherwise
        """
        if not token:
            return None
            
        try:
            payload = jwt.decode(
                token,
                self.secret_key,
                algorithms=[self.algorithm],
                options={'verify_signature': True}
            )
            
            # Check if token has expired
            exp = payload.get('exp', 0)
            now = int(time.time())
            if exp < now:
                logger.warning("Token has expired for user %s", payload.get('sub', 'unknown'))
                return None
                
            return payload
        except jwt.PyJWTError as e:
            logger.error("JWT validation error: %s", str(e))
            return None
    
    def create_convex_identity(self, user_id, user_data=None):
        """
        Create a Convex-compatible identity object.
        
        Args:
            user_id: User ID
            user_data: Additional user data
            
        Returns:
            dict: Identity object for Convex
        """
        token = self.generate_token(user_id, user_data)
        
        # Extract payload for identity
        payload = self.validate_token(token)
        if not payload:
            logger.error("Failed to validate newly created token for user %s", user_id)
            return None
            
        convex_data = payload.get('convex', {})
        
        # Create identity object for Convex
        identity = {
            'identity': str(user_id),
            'tokenIdentifier': convex_data.get('tokenIdentifier', str(user_id)),
            'token': token,
            'authType': 'jwt',
            'expiry': payload.get('exp'),
            'traits': convex_data.get('traits', {})
        }
        
        return identity
        
    def jwt_middleware(self):
        """
        Create a Flask middleware function for JWT authentication.
        
        Returns:
            function: Middleware function to process JWTs
        """
        def process_jwt():
            auth_header = request.headers.get('Authorization')
            if auth_header and auth_header.startswith('Bearer '):
                token = auth_header[7:]  # Remove 'Bearer ' prefix
                
                # Validate token
                payload = self.validate_token(token)
                
                if payload:
                    # Store user info in Flask g object
                    g.user_id = payload.get('sub')
                    g.user_data = {
                        'email': payload.get('email'),
                        'role': payload.get('role', 'user'),
                        'name': payload.get('name', ''),
                        'permissions': payload.get('permissions', [])
                    }
                    
                    # Store convex-specific data
                    g.convex_identity = payload.get('convex')
                    g.convex_token = token
                    
                    logger.debug("Authenticated user %s with role %s", 
                               g.user_id, g.user_data.get('role'))
                else:
                    g.user_id = None
                    g.user_data = None
                    g.convex_identity = None
                    g.convex_token = None
            else:
                g.user_id = None
                g.user_data = None
                g.convex_identity = None
                g.convex_token = None
            
            return None
        
        return process_jwt
    
    def get_admin_token(self):
        """
        Generate an admin token for Convex operations.
        
        Returns:
            str: Admin token for Convex
        """
        # Create admin user data
        admin_data = {
            'email': 'admin@cryoprotect.internal',
            'role': 'admin',
            'name': 'System Admin',
            'permissions': ['*']
        }
        
        # Generate token with admin privileges
        return self.generate_token('admin', admin_data)

# Enhanced decorators for authentication and authorization

def jwt_required(f):
    """
    Decorator that requires a valid JWT for the endpoint.
    
    Args:
        f: Function to decorate
        
    Returns:
        function: Decorated function
    """
    @wraps(f)
    def decorated(*args, **kwargs):
        # Check if user_id is set in g by the middleware
        if not g.get('user_id'):
            return jsonify({'message': 'Authentication required', 'error': 'unauthorized'}), 401
        
        return f(*args, **kwargs)
    
    return decorated

def role_required(role):
    """
    Decorator that requires a specific role for the endpoint.
    
    Args:
        role: Required role or list of roles
        
    Returns:
        function: Decorator function
    """
    def decorator(f):
        @wraps(f)
        @jwt_required
        def decorated(*args, **kwargs):
            user_role = g.user_data.get('role')
            
            # Admin role has access to everything
            if user_role == 'admin':
                return f(*args, **kwargs)
            
            # Check if role matches required role
            if isinstance(role, list):
                if user_role not in role:
                    return jsonify({
                        'message': f'Access denied. Required role: {", ".join(role)}',
                        'error': 'forbidden'
                    }), 403
            else:
                if user_role != role:
                    return jsonify({
                        'message': f'Access denied. Required role: {role}',
                        'error': 'forbidden'
                    }), 403
            
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator

def convex_identity_required(f):
    """
    Decorator that requires a valid Convex identity for the endpoint.
    
    Args:
        f: Function to decorate
        
    Returns:
        function: Decorated function
    """
    @wraps(f)
    @jwt_required
    def decorated(*args, **kwargs):
        # Check if convex_identity is set in g by the middleware
        if not g.get('convex_identity'):
            return jsonify({
                'message': 'Convex authentication required',
                'error': 'unauthorized_convex'
            }), 401
        
        return f(*args, **kwargs)
    
    return decorated

# Integration function for Flask app
def init_auth_bridge(app, secret_key=None, algorithm=None, expiry=None):
    """
    Initialize the AuthBridge for a Flask app.
    
    Args:
        app: Flask app instance
        secret_key: Secret key for JWT encoding/decoding
        algorithm: Algorithm for JWT encoding/decoding
        expiry: Token expiry time in seconds
        
    Returns:
        AuthBridge: Configured auth bridge instance
    """
    # Create and initialize auth bridge
    auth_bridge = AuthBridge(app, secret_key, algorithm, expiry)
    
    # Register error handlers for auth errors
    @app.errorhandler(401)
    def unauthorized(error):
        return jsonify({
            'message': 'Authentication required',
            'error': 'unauthorized'
        }), 401
    
    @app.errorhandler(403)
    def forbidden(error):
        return jsonify({
            'message': 'Access denied',
            'error': 'forbidden'
        }), 403
    
    # Add convex-specific helpers to the app
    app.get_convex_admin_token = auth_bridge.get_admin_token
    app.create_convex_identity = auth_bridge.create_convex_identity
    
    logger.info("Auth bridge initialized for Flask app")
    return auth_bridge
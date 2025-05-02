"""
CryoProtect Analyzer - JWT Authentication

This module provides functions for JWT token generation, validation, and user extraction.
"""

import os
import jwt
import json
import logging
import requests
from datetime import datetime, timedelta
from functools import wraps
from flask import request, g, jsonify, current_app
from typing import Dict, Any, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# Cache for JWKs to avoid fetching them on every request
jwks_cache = {
    'keys': None,
    'last_updated': None,
    'cache_duration': timedelta(hours=24)  # Cache JWKs for 24 hours
}

def get_jwks() -> Dict:
    """
    Fetch JSON Web Key Set (JWKS) from Supabase.
    Uses caching to avoid fetching on every request.
    
    Returns:
        Dict: JWKS containing public keys
    """
    global jwks_cache
    
    now = datetime.now()
    
    # Return cached keys if they exist and are not expired
    if (jwks_cache['keys'] is not None and 
        jwks_cache['last_updated'] is not None and 
        now - jwks_cache['last_updated'] < jwks_cache['cache_duration']):
        return jwks_cache['keys']
    
    # Fetch new keys
    try:
        supabase_url = os.environ.get('SUPABASE_URL') or current_app.config.get('SUPABASE_URL')
        supabase_key = os.environ.get('SUPABASE_KEY') or current_app.config.get('SUPABASE_KEY')
        
        if not supabase_url:
            raise ValueError("SUPABASE_URL is not configured")
        if not supabase_key:
            raise ValueError("SUPABASE_KEY is not configured")
        
        # Construct the JWKS URL
        jwks_url = f"{supabase_url}/auth/v1/jwks"
        
        # Include the Supabase API key in the request headers
        headers = {
            'apikey': supabase_key,
            'Authorization': f'Bearer {supabase_key}'
        }
        
        response = requests.get(jwks_url, headers=headers)
        response.raise_for_status()
        
        jwks = response.json()
        
        # Update cache
        jwks_cache['keys'] = jwks
        jwks_cache['last_updated'] = now
        
        return jwks
    except Exception as e:
        logger.error(f"Failed to fetch JWKS: {str(e)}")
        # Return cached keys if available, even if expired
        if jwks_cache['keys'] is not None:
            logger.warning("Using expired JWKS cache due to fetch failure")
            return jwks_cache['keys']
        raise

def find_signing_key(token_header: Dict, jwks: Dict) -> Optional[Dict]:
    """
    Find the signing key in the JWKS that matches the key ID in the token header.
    
    Args:
        token_header: The decoded token header
        jwks: The JSON Web Key Set
        
    Returns:
        The signing key or None if not found
    """
    kid = token_header.get('kid')
    if not kid:
        return None
    
    for key in jwks.get('keys', []):
        if key.get('kid') == kid:
            return key
    
    return None

def decode_token(token: str) -> Dict:
    """
    Decode and validate a JWT token.
    
    Args:
        token: The JWT token to decode
        
    Returns:
        The decoded token payload
        
    Raises:
        jwt.InvalidTokenError: If the token is invalid
    """
    # Decode the token header without verification to get the key ID
    header = jwt.get_unverified_header(token)
    
    # Get the JWKS
    jwks = get_jwks()
    
    # Find the signing key
    signing_key = find_signing_key(header, jwks)
    if not signing_key:
        raise jwt.InvalidTokenError(f"Signing key not found for kid: {header.get('kid')}")
    
    # Construct the public key from the JWKS
    public_key = jwt.algorithms.RSAAlgorithm.from_jwk(json.dumps(signing_key))
    
    # Decode and verify the token
    supabase_url = os.environ.get('SUPABASE_URL') or current_app.config.get('SUPABASE_URL')
    audience = f"{supabase_url}/auth/v1"
    
    decoded = jwt.decode(
        token,
        public_key,
        algorithms=['RS256'],
        options={
            'verify_signature': True,
            'verify_exp': True,
            'verify_nbf': True,
            'verify_iat': True,
            'verify_aud': True,
            'verify_iss': True,
            'require_exp': True,
            'require_iat': True,
            'require_nbf': False
        },
        audience=audience,
        issuer=supabase_url
    )
    
    return decoded

def extract_user_from_token(token: str) -> Tuple[Dict, str]:
    """
    Extract user information from a JWT token.
    
    Args:
        token: The JWT token
        
    Returns:
        Tuple of (user_data, user_id)
        
    Raises:
        jwt.InvalidTokenError: If the token is invalid
    """
    payload = decode_token(token)
    
    # Extract user ID from the 'sub' claim
    user_id = payload.get('sub')
    if not user_id:
        raise jwt.InvalidTokenError("Token does not contain a subject (user ID)")
    
    # Extract user data
    user_data = {
        'id': user_id,
        'email': payload.get('email', ''),
        'app_metadata': payload.get('app_metadata', {}),
        'user_metadata': payload.get('user_metadata', {}),
        'role': payload.get('role', 'user'),  # Default role is 'user'
        'roles': payload.get('roles', []),    # All roles assigned to the user
        'permissions': payload.get('permissions', []),  # All permissions granted to the user
        'aud': payload.get('aud', ''),
        'exp': payload.get('exp', 0)
    }
    
    # If roles is not in the token, try to fetch from database
    if not user_data['roles'] and user_id:
        try:
            from api.rbac import UserRoleManager
            user_data['roles'] = UserRoleManager.get_user_roles(user_id)
        except Exception as e:
            logger.warning(f"Failed to fetch user roles: {str(e)}")
    
    # If permissions is not in the token, try to fetch from database
    if not user_data['permissions'] and user_id:
        try:
            from api.rbac import UserRoleManager
            user_data['permissions'] = UserRoleManager.get_user_permissions(user_id)
        except Exception as e:
            logger.warning(f"Failed to fetch user permissions: {str(e)}")
    
    return user_data, user_id

def get_token_from_request() -> Optional[str]:
    """
    Extract the JWT token from the request.
    Checks Authorization header and cookies.
    
    Returns:
        The token or None if not found
    """
    # Check Authorization header
    auth_header = request.headers.get('Authorization', '')
    if auth_header.startswith('Bearer '):
        return auth_header.split(' ')[1]
    
    # Check cookies
    token = request.cookies.get('access_token')
    if token:
        return token
    
    return None

def jwt_required(f):
    """
    Decorator to require a valid JWT token for API endpoints.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @wraps(f)
    def decorated(*args, **kwargs):
        token = get_token_from_request()
        
        if not token:
            return jsonify({'message': 'Authentication token is missing'}), 401
        
        try:
            # Decode and validate token
            user_data, user_id = extract_user_from_token(token)
            
            # Store user data in request context
            g.user = user_data
            g.user_id = user_id
            
            return f(*args, **kwargs)
        except jwt.ExpiredSignatureError:
            return jsonify({'message': 'Authentication token has expired'}), 401
        except jwt.InvalidTokenError as e:
            logger.warning(f"Invalid token: {str(e)}")
            return jsonify({'message': 'Invalid authentication token'}), 401
        except Exception as e:
            logger.error(f"Token validation error: {str(e)}")
            return jsonify({'message': 'Authentication error'}), 401
    
    return decorated

def get_current_user() -> Optional[Dict]:
    """
    Get the current authenticated user from the request context.
    
    Returns:
        User data or None if not authenticated
    """
    if hasattr(g, 'user'):
        return g.user
    
    # Try to authenticate from token
    token = get_token_from_request()
    if token:
        try:
            user_data, user_id = extract_user_from_token(token)
            g.user = user_data
            g.user_id = user_id
            return user_data
        except Exception as e:
            logger.warning(f"Failed to extract user from token: {str(e)}")
    
    return None

def has_role(required_role: Union[str, list]) -> bool:
    """
    Check if the current user has the required role.
    
    Args:
        required_role: Role or list of roles required
        
    Returns:
        True if the user has the required role, False otherwise
    """
    user = get_current_user()
    if not user:
        return False
    
    user_id = user.get('id')
    
    # Try to use the RBAC system first
    try:
        from api.rbac import UserRoleManager
        
        # Convert single role to list
        roles = required_role if isinstance(required_role, list) else [required_role]
        
        # Check if user has any of the required roles
        for role in roles:
            if UserRoleManager.has_role(user_id, role):
                return True
        
        return False
    except Exception as e:
        logger.warning(f"Failed to use RBAC for role check: {str(e)}")
        
        # Fall back to legacy role check
        user_role = user.get('role', 'user')
        
        if isinstance(required_role, list):
            return user_role in required_role
        
        return user_role == required_role

def role_required(required_role: Union[str, list]):
    """
    Decorator to require a specific role for API endpoints.
    
    Args:
        required_role: Role or list of roles required
        
    Returns:
        Decorated function
    """
    def decorator(f):
        @wraps(f)
        @jwt_required
        def decorated(*args, **kwargs):
            if not has_role(required_role):
                return jsonify({'message': 'Insufficient permissions'}), 403
            
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator
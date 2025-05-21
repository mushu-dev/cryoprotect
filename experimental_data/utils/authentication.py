#!/usr/bin/env python3
"""
Authentication utilities for CryoProtect Enhanced Experimental Data System.

This module provides authentication utilities for making API requests.
"""

from typing import Dict, Any, Optional
import os
import logging
import jwt
import time
from datetime import datetime, timedelta

# Configure logging
logger = logging.getLogger(__name__)

# JWT token cache
_token_cache = {
    "token": None,
    "expires_at": 0
}

async def get_auth_headers() -> Dict[str, str]:
    """
    Get authentication headers for API requests.
    
    Returns:
        Dict with authentication headers
    """
    token = await get_auth_token()
    return {"Authorization": f"Bearer {token}"}

async def get_auth_token() -> str:
    """
    Get an authentication token for API requests.
    If a valid token is cached, it will be returned.
    Otherwise, a new token will be generated.
    
    Returns:
        Authentication token
    """
    # Check if we have a cached token that's still valid
    now = time.time()
    if _token_cache["token"] and _token_cache["expires_at"] > now + 60:  # Add 60s buffer
        return _token_cache["token"]
    
    # Generate a new token
    token = _generate_token()
    
    # Cache the token
    _token_cache["token"] = token
    # Set expiration to 1 hour from now
    _token_cache["expires_at"] = now + 3600
    
    return token

def _generate_token() -> str:
    """
    Generate a JWT token for authentication.
    
    Returns:
        JWT token
    
    Raises:
        ValueError: If required environment variables are not set
    """
    try:
        # Get required environment variables
        secret_key = os.environ.get("JWT_SECRET_KEY")
        service_role = os.environ.get("SERVICE_ROLE")
        
        if not secret_key:
            raise ValueError("JWT_SECRET_KEY environment variable not set")
        
        if not service_role:
            raise ValueError("SERVICE_ROLE environment variable not set")
        
        # Create token payload
        now = datetime.utcnow()
        expires = now + timedelta(hours=1)
        
        payload = {
            "sub": "service",
            "role": service_role,
            "iat": now.timestamp(),
            "exp": expires.timestamp()
        }
        
        # Generate token
        token = jwt.encode(payload, secret_key, algorithm="HS256")
        
        return token
    
    except Exception as e:
        logger.error(f"Error generating JWT token: {str(e)}")
        raise ValueError(f"Failed to generate auth token: {str(e)}")

def validate_token(token: str) -> Dict[str, Any]:
    """
    Validate a JWT token.
    
    Args:
        token: JWT token to validate
    
    Returns:
        Dict with decoded token payload
    
    Raises:
        ValueError: If the token is invalid
    """
    try:
        # Get secret key
        secret_key = os.environ.get("JWT_SECRET_KEY")
        
        if not secret_key:
            raise ValueError("JWT_SECRET_KEY environment variable not set")
        
        # Decode token
        payload = jwt.decode(token, secret_key, algorithms=["HS256"])
        
        return payload
    
    except jwt.ExpiredSignatureError:
        raise ValueError("Token has expired")
    
    except jwt.InvalidTokenError as e:
        raise ValueError(f"Invalid token: {str(e)}")
    
    except Exception as e:
        logger.error(f"Error validating JWT token: {str(e)}")
        raise ValueError(f"Failed to validate token: {str(e)}")
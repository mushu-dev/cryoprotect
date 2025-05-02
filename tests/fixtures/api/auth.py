"""
Authentication fixtures for API testing.

This module provides fixtures for authentication testing, including
tokens for different user roles and authentication scenarios.
"""

import pytest
import jwt
import uuid
from datetime import datetime, timedelta
from typing import Dict, Any

# Secret key for JWT token generation (for testing only)
JWT_SECRET = "test-secret-key"

# User roles
ROLE_ADMIN = "admin"
ROLE_USER = "user"
ROLE_SCIENTIST = "scientist"

# Sample user IDs
ADMIN_USER_ID = str(uuid.uuid4())
REGULAR_USER_ID = str(uuid.uuid4())
SCIENTIST_USER_ID = str(uuid.uuid4())


def generate_token(user_id: str, role: str, expires_in: int = 3600) -> str:
    """
    Generate a JWT token for testing.
    
    Args:
        user_id: User ID
        role: User role
        expires_in: Token expiration time in seconds
        
    Returns:
        JWT token
    """
    now = datetime.utcnow()
    payload = {
        'sub': user_id,
        'role': role,
        'iat': now,
        'exp': now + timedelta(seconds=expires_in),
        'iss': 'cryoprotect-test'
    }
    
    return jwt.encode(payload, JWT_SECRET, algorithm='HS256')


@pytest.fixture
def auth_token() -> str:
    """
    Generate a regular user authentication token.
    
    Returns:
        JWT token
    """
    return generate_token(REGULAR_USER_ID, ROLE_USER)


@pytest.fixture
def admin_token() -> str:
    """
    Generate an admin user authentication token.
    
    Returns:
        JWT token
    """
    return generate_token(ADMIN_USER_ID, ROLE_ADMIN)


@pytest.fixture
def user_token() -> str:
    """
    Generate a regular user authentication token.
    
    Returns:
        JWT token
    """
    return generate_token(REGULAR_USER_ID, ROLE_USER)


@pytest.fixture
def scientist_token() -> str:
    """
    Generate a scientist user authentication token.
    
    Returns:
        JWT token
    """
    return generate_token(SCIENTIST_USER_ID, ROLE_SCIENTIST)


@pytest.fixture
def expired_token() -> str:
    """
    Generate an expired authentication token.
    
    Returns:
        Expired JWT token
    """
    return generate_token(REGULAR_USER_ID, ROLE_USER, expires_in=-3600)


@pytest.fixture
def invalid_token() -> str:
    """
    Generate an invalid authentication token.
    
    Returns:
        Invalid token
    """
    return "invalid-token"


@pytest.fixture
def mock_auth_middleware():
    """
    Mock authentication middleware for testing.
    
    This fixture patches the authentication middleware to allow
    testing without actual authentication.
    
    Returns:
        Function to set the current user for testing
    """
    import flask
    from unittest.mock import patch
    from functools import wraps
    
    # Store the original token_required decorator
    from api.utils import token_required
    original_token_required = token_required
    
    # Mock user for testing
    mock_user = None
    
    # Mock token_required decorator
    def mock_token_required(f):
        @wraps(f)
        def decorated(*args, **kwargs):
            # Set the user in the Flask g object
            flask.g.user_id = mock_user['id'] if mock_user else None
            flask.g.user_role = mock_user['role'] if mock_user else None
            
            return f(*args, **kwargs)
        return decorated
    
    # Patch the token_required decorator
    with patch('api.utils.token_required', mock_token_required):
        # Function to set the current user for testing
        def set_current_user(user=None):
            nonlocal mock_user
            mock_user = user
        
        yield set_current_user
    
    # Restore the original token_required decorator
    token_required = original_token_required


@pytest.fixture
def mock_current_user(mock_auth_middleware):
    """
    Set a mock current user for testing.
    
    Args:
        mock_auth_middleware: Mock authentication middleware
        
    Returns:
        Function to set the current user for testing
    """
    return mock_auth_middleware


@pytest.fixture
def mock_admin_user(mock_current_user):
    """
    Set a mock admin user for testing.
    
    Args:
        mock_current_user: Function to set the current user
        
    Returns:
        None
    """
    mock_current_user({
        'id': ADMIN_USER_ID,
        'role': ROLE_ADMIN,
        'email': 'admin@example.com'
    })


@pytest.fixture
def mock_regular_user(mock_current_user):
    """
    Set a mock regular user for testing.
    
    Args:
        mock_current_user: Function to set the current user
        
    Returns:
        None
    """
    mock_current_user({
        'id': REGULAR_USER_ID,
        'role': ROLE_USER,
        'email': 'user@example.com'
    })


@pytest.fixture
def mock_scientist_user(mock_current_user):
    """
    Set a mock scientist user for testing.
    
    Args:
        mock_current_user: Function to set the current user
        
    Returns:
        None
    """
    mock_current_user({
        'id': SCIENTIST_USER_ID,
        'role': ROLE_SCIENTIST,
        'email': 'scientist@example.com'
    })


@pytest.fixture
def mock_unauthenticated_user(mock_current_user):
    """
    Set no current user for testing unauthenticated scenarios.
    
    Args:
        mock_current_user: Function to set the current user
        
    Returns:
        None
    """
    mock_current_user(None)
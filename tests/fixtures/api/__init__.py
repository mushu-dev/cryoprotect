"""
API fixtures for testing.

This module provides fixtures for API testing, including
client setup, authentication, and request/response handling.
"""

from tests.fixtures.api.client import (
    api_client,
    mock_api_client,
    authenticated_client,
    admin_client,
    scientist_client
)

from tests.fixtures.api.auth import (
    auth_token,
    admin_token,
    user_token,
    scientist_token,
    expired_token,
    invalid_token
)

__all__ = [
    'api_client',
    'mock_api_client',
    'authenticated_client',
    'admin_client',
    'scientist_client',
    'auth_token',
    'admin_token',
    'user_token',
    'scientist_token',
    'expired_token',
    'invalid_token'
]
"""
JWT-based Service Role Authentication for CryoProtect.

This package provides components for secure service-to-service communication
using JWT tokens with role-based authentication and authorization.

Components:
- TokenManager: Manages JWT token creation, validation, and revocation
- ServiceRoleClient: Client for making authenticated API requests
- service_role_required: Decorator for securing API endpoints
"""

# Import main components for easier access
from .token_manager import (
    get_token_manager, 
    create_service_token, 
    validate_service_token, 
    revoke_service_token
)

from .service_role_client import (
    ServiceRoleClient,
    create_client
)

from .service_role_operation import (
    service_role_required,
    admin_service_required,
    data_service_required,
    analytics_service_required,
    get_current_service,
    require_service_scope,
    ServiceRoleAuth
)

# Package metadata
__version__ = '1.0.0'
__all__ = [
    # Token Manager
    'get_token_manager',
    'create_service_token',
    'validate_service_token',
    'revoke_service_token',
    
    # Service Role Client
    'ServiceRoleClient',
    'create_client',
    
    # Service Role Operation
    'service_role_required',
    'admin_service_required',
    'data_service_required',
    'analytics_service_required',
    'get_current_service',
    'require_service_scope',
    'ServiceRoleAuth'
]
"""
CryoProtect - JWT Service Role Authentication

This package provides enhanced JWT-based service role authentication for the
CryoProtect application. It includes token management, client implementation,
and a decorator for service role operations.

The service role authentication system allows for secure API access with
fine-grained permission control through scopes, token validation, and
comprehensive audit logging.

Components:
- TokenManager: Manages service role JWT tokens
- ServiceRoleClient: Client for making authenticated service role API requests
- service_role_operation: Decorator for API functions requiring service role auth

Usage:
    from jwt_service_role import ServiceRoleTokenManager, service_role_operation
    
    # In API endpoint
    @app.route('/api/admin/operation')
    @service_role_operation(required_scope='admin:operation')
    def admin_operation(client_id):
        # Authenticated operation
        return {"status": "success"}
    
    # In client code
    from jwt_service_role import ServiceRoleClient
    
    client = ServiceRoleClient(
        base_url='https://api.example.com',
        client_id='service-client-1',
        client_scopes=['admin:operation']
    )
    
    response = client.get('/api/admin/operation')
"""

from .token_manager import ServiceRoleTokenManager, VALID_SCOPES
from .service_role_client import ServiceRoleClient
from .service_role_operation import (
    service_role_operation,
    get_current_service_role_client,
    get_current_service_role_scopes
)

__all__ = [
    'ServiceRoleTokenManager',
    'ServiceRoleClient',
    'service_role_operation',
    'get_current_service_role_client',
    'get_current_service_role_scopes',
    'VALID_SCOPES'
]
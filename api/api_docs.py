"""
CryoProtect Analyzer API - Documentation Utilities

This module provides utilities for generating OpenAPI documentation for API endpoints.
It ensures consistent documentation across all endpoints and provides helpers for
documenting common patterns.
"""

import logging
from typing import Dict, Any, List, Optional, Union, Callable
from flask import current_app
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from flask_apispec.extension import FlaskApiSpec
from marshmallow import Schema, fields

# Set up logging
logger = logging.getLogger(__name__)

# Common response schemas
class ErrorResponseSchema(Schema):
    """Schema for error responses."""
    status = fields.String(required=True, example="error")
    timestamp = fields.String(required=True, example="2025-04-21T17:45:00.000Z")
    code = fields.Integer(required=True, example=400)
    message = fields.String(required=True, example="Bad request")
    errors = fields.List(fields.Dict(), required=True, example=[
        {
            "type": "ValidationError",
            "message": "Invalid input data",
            "context": "Validating request data"
        }
    ])

class SuccessResponseSchema(Schema):
    """Schema for success responses."""
    status = fields.String(required=True, example="success")
    timestamp = fields.String(required=True, example="2025-04-21T17:45:00.000Z")
    code = fields.Integer(required=True, example=200)
    message = fields.String(required=True, example="Request succeeded")
    data = fields.Raw(required=False)
    meta = fields.Dict(required=False)
    pagination = fields.Dict(required=False)

class PaginationSchema(Schema):
    """Schema for pagination information."""
    page = fields.Integer(required=True, example=1)
    per_page = fields.Integer(required=True, example=20)
    total_items = fields.Integer(required=True, example=100)
    total_pages = fields.Integer(required=True, example=5)
    has_next = fields.Boolean(required=True, example=True)
    has_prev = fields.Boolean(required=True, example=False)
    links = fields.Dict(required=False, example={
        "self": "/api/v1/molecules?page=1&per_page=20",
        "first": "/api/v1/molecules?page=1&per_page=20",
        "last": "/api/v1/molecules?page=5&per_page=20",
        "next": "/api/v1/molecules?page=2&per_page=20"
    })

# Common response examples
COMMON_RESPONSES = {
    200: {
        "description": "Successful response",
        "content": {
            "application/json": {
                "schema": SuccessResponseSchema
            }
        }
    },
    400: {
        "description": "Bad request",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    },
    401: {
        "description": "Unauthorized",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    },
    403: {
        "description": "Forbidden",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    },
    404: {
        "description": "Not found",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    },
    429: {
        "description": "Too many requests",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    },
    500: {
        "description": "Internal server error",
        "content": {
            "application/json": {
                "schema": ErrorResponseSchema
            }
        }
    }
}

def init_docs(app):
    """
    Initialize API documentation.
    
    Args:
        app: Flask application
        
    Returns:
        FlaskApiSpec instance
    """
    app.config.update({
        'APISPEC_SPEC': APISpec(
            title=app.config.get('API_TITLE', 'CryoProtect API'),
            version=app.config.get('API_VERSION', '1.0.0'),
            openapi_version=app.config.get('OPENAPI_VERSION', '3.0.2'),
            plugins=[MarshmallowPlugin()],
            info={
                'description': 'API for CryoProtect Analyzer',
                'contact': {
                    'name': app.config.get('API_CONTACT_NAME', 'CryoProtect API Support'),
                    'email': app.config.get('API_CONTACT_EMAIL', 'support@cryoprotect.com'),
                    'url': app.config.get('API_CONTACT_URL', 'https://cryoprotect.com/support')
                },
                'license': {
                    'name': app.config.get('API_LICENSE_NAME', 'MIT'),
                    'url': app.config.get('API_LICENSE_URL', 'https://opensource.org/licenses/MIT')
                }
            }
        ),
        'APISPEC_SWAGGER_URL': '/api/v1/swagger/',
        'APISPEC_SWAGGER_UI_URL': '/swagger-ui/'
    })
    
    # Register security schemes
    app.config['APISPEC_SPEC'].components.security_scheme(
        'bearerAuth', 
        {
            'type': 'http',
            'scheme': 'bearer',
            'bearerFormat': 'JWT',
            'description': 'Enter JWT token in the format: Bearer <token>'
        }
    )
    
    # Initialize FlaskApiSpec
    docs = FlaskApiSpec(app)
    
    return docs

def register_resource(docs, resource_class, endpoint_name):
    """
    Register a resource with the API documentation.
    
    Args:
        docs: FlaskApiSpec instance
        resource_class: Resource class to register
        endpoint_name: Endpoint name
    """
    try:
        docs.register(resource_class, endpoint=endpoint_name)
        logger.info(f"Registered resource {resource_class.__name__} with documentation")
    except Exception as e:
        logger.warning(f"Failed to register resource {resource_class.__name__} with documentation: {str(e)}")

def document_endpoint(
    summary: str,
    description: str = None,
    tags: List[str] = None,
    params: Dict[str, Dict[str, Any]] = None,
    request_body: Dict[str, Any] = None,
    responses: Dict[int, Dict[str, Any]] = None,
    security: List[Dict[str, List[str]]] = None,
    deprecated: bool = False
):
    """
    Decorator to document an API endpoint.
    
    Args:
        summary: Short summary of the endpoint
        description: Detailed description of the endpoint
        tags: List of tags for categorizing the endpoint
        params: Dictionary of parameters
        request_body: Request body schema
        responses: Dictionary of responses
        security: Security requirements
        deprecated: Whether the endpoint is deprecated
        
    Returns:
        Decorator function
    """
    def decorator(f):
        # Store documentation in function attributes
        f.__apidoc__ = getattr(f, '__apidoc__', {})
        f.__apidoc__['summary'] = summary
        
        if description:
            f.__apidoc__['description'] = description
            
        if tags:
            f.__apidoc__['tags'] = tags
            
        if params:
            f.__apidoc__['params'] = params
            
        if request_body:
            f.__apidoc__['requestBody'] = request_body
            
        # Merge common responses with provided responses
        endpoint_responses = {}
        for status_code, response in COMMON_RESPONSES.items():
            endpoint_responses[status_code] = response
            
        if responses:
            for status_code, response in responses.items():
                endpoint_responses[status_code] = response
                
        f.__apidoc__['responses'] = endpoint_responses
        
        if security:
            f.__apidoc__['security'] = security
            
        if deprecated:
            f.__apidoc__['deprecated'] = deprecated
            
        return f
    
    return decorator

def auth_required(security_scheme='bearerAuth'):
    """
    Decorator to indicate that an endpoint requires authentication.
    
    Args:
        security_scheme: Security scheme to use
        
    Returns:
        Decorator function
    """
    def decorator(f):
        f.__apidoc__ = getattr(f, '__apidoc__', {})
        f.__apidoc__['security'] = [{security_scheme: []}]
        return f
    
    return decorator

def paginated_response(data_schema):
    """
    Create a paginated response schema.
    
    Args:
        data_schema: Schema for the data items
        
    Returns:
        Response schema with pagination
    """
    class PaginatedResponseSchema(SuccessResponseSchema):
        data = fields.List(fields.Nested(data_schema))
        pagination = fields.Nested(PaginationSchema)
        
    return PaginatedResponseSchema

def generate_openapi_spec():
    """
    Generate OpenAPI specification for the API.
    
    Returns:
        OpenAPI specification as a dictionary
    """
    if 'current_app' in globals() and hasattr(current_app, 'config'):
        spec = current_app.config.get('APISPEC_SPEC')
        if spec:
            return spec.to_dict()
    
    # Create a new spec if not available in current_app
    spec = APISpec(
        title='CryoProtect API',
        version='1.0.0',
        openapi_version='3.0.2',
        plugins=[MarshmallowPlugin()],
        info={
            'description': 'API for CryoProtect Analyzer',
            'contact': {
                'name': 'CryoProtect API Support',
                'email': 'support@cryoprotect.com',
                'url': 'https://cryoprotect.com/support'
            },
            'license': {
                'name': 'MIT',
                'url': 'https://opensource.org/licenses/MIT'
            }
        }
    )
    
    # Register security schemes
    spec.components.security_scheme(
        'bearerAuth', 
        {
            'type': 'http',
            'scheme': 'bearer',
            'bearerFormat': 'JWT',
            'description': 'Enter JWT token in the format: Bearer <token>'
        }
    )
    
    return spec.to_dict()
"""
CryoProtect API Documentation Utilities

This module provides utilities for generating and serving API documentation
using flask-apispec and apispec with OpenAPI/Swagger standards.
"""

import os
from flask import Blueprint, jsonify, render_template, current_app, send_from_directory
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from flask_apispec.extension import FlaskApiSpec
from flask_apispec import use_kwargs, marshal_with
from marshmallow import Schema, fields
from api.schemas import (
    ErrorResponseSchema, SuccessResponseSchema, MetadataSchema,
    MoleculeSchema, MoleculeListSchema, MoleculeCreateSchema,
    MixtureSchema, MixtureListSchema, MixtureCreateSchema,
    ExperimentSchema, ExperimentListSchema, ExperimentCreateSchema
)

# Create a blueprint for API documentation
docs_bp = Blueprint('api_docs', __name__, url_prefix='/api/docs')

def init_docs(app):
    """
    Initialize API documentation for the Flask application.
    
    Args:
        app: Flask application instance
    
    Returns:
        FlaskApiSpec: Configured FlaskApiSpec instance
    """
    # Configure API documentation
    spec = APISpec(
        title=app.config.get('API_TITLE', 'CryoProtect API'),
        version=app.config.get('API_VERSION', '1.0.0'),
        openapi_version=app.config.get('OPENAPI_VERSION', '3.0.2'),
        plugins=[MarshmallowPlugin()],
        info={
            'description': 'CryoProtect API Documentation',
            'termsOfService': app.config.get('API_TERMS_URL', ''),
            'contact': {
                'name': app.config.get('API_CONTACT_NAME', 'API Support'),
                'email': app.config.get('API_CONTACT_EMAIL', 'support@cryoprotect.com'),
                'url': app.config.get('API_CONTACT_URL', '')
            },
            'license': {
                'name': app.config.get('API_LICENSE_NAME', 'MIT'),
                'url': app.config.get('API_LICENSE_URL', '')
            }
        }
    )
    
    # Add security schemes to the API specification
    add_security_schemes(spec)
    
    # Manually add paths to the API specification
    # Molecules endpoints
    spec.path(
        path="/api/v1/molecules",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get a list of molecules',
                'description': 'Returns a paginated list of molecules',
                'responses': {'200': {'description': 'List of molecules'}}
            },
            'post': {
                'summary': 'Create a new molecule',
                'description': 'Import a molecule from PubChem',
                'responses': {'201': {'description': 'Molecule created'}}
            }
        })
    )
    
    spec.path(
        path="/api/v1/molecules/{molecule_id}",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get a molecule',
                'description': 'Returns a molecule by ID',
                'parameters': [{'name': 'molecule_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Molecule details'}}
            },
            'put': {
                'summary': 'Update a molecule',
                'description': 'Update a molecule by ID',
                'parameters': [{'name': 'molecule_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Updated molecule'}}
            },
            'delete': {
                'summary': 'Delete a molecule',
                'description': 'Delete a molecule by ID',
                'parameters': [{'name': 'molecule_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Molecule deleted'}}
            }
        })
    )
    
    # Molecule score endpoint
    spec.path(
        path="/api/v1/molecules/{molecule_id}/score",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate cryoprotection score for a molecule',
                'description': 'Calculate and optionally store a cryoprotection score for a molecule',
                'parameters': [{'name': 'molecule_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Molecule score'}}
            }
        })
    )
    
    # Molecule property calculation endpoint
    spec.path(
        path="/api/v1/molecules/{molecule_id}/calculate-properties",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate properties for a molecule',
                'description': 'Calculate and store properties for a molecule using RDKit',
                'parameters': [{'name': 'molecule_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Properties calculated'}}
            }
        })
    )
    
    # Mixture endpoints
    spec.path(
        path="/api/v1/mixtures",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get a list of mixtures',
                'description': 'Returns a paginated list of mixtures',
                'parameters': [
                    {'name': 'page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Page number'},
                    {'name': 'per_page', 'in': 'query', 'schema': {'type': 'integer'}, 'description': 'Items per page'},
                    {'name': 'search', 'in': 'query', 'schema': {'type': 'string'}, 'description': 'Search term'}
                ],
                'responses': {
                    '200': {
                        'description': 'List of mixtures',
                        'content': {
                            'application/json': {
                                'schema': MixtureListSchema
                            }
                        }
                    }
                }
            },
            'post': {
                'summary': 'Create a new mixture',
                'description': 'Creates a new mixture from the provided data',
                'requestBody': {
                    'required': True,
                    'content': {
                        'application/json': {
                            'schema': MixtureCreateSchema
                        }
                    }
                },
                'responses': {
                    '201': {
                        'description': 'Mixture created successfully',
                        'content': {
                            'application/json': {
                                'schema': MixtureSchema
                            }
                        }
                    },
                    '400': {
                        'description': 'Invalid request data'
                    }
                }
            }
        })
    )
    
    spec.path(
        path="/api/v1/mixtures/{mixture_id}",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get a mixture',
                'description': 'Returns a mixture by ID with its components',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Mixture details'}}
            },
            'put': {
                'summary': 'Update a mixture',
                'description': 'Update a mixture by ID',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Updated mixture'}}
            },
            'delete': {
                'summary': 'Delete a mixture',
                'description': 'Delete a mixture by ID',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Mixture deleted'}}
            }
        })
    )
    
    # Mixture score endpoint
    spec.path(
        path="/api/v1/mixtures/{mixture_id}/score",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate cryoprotection score for a mixture',
                'description': 'Calculate and optionally store a cryoprotection score for a mixture',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Mixture score'}}
            }
        })
    )
    
    # Predictions endpoints
    spec.path(
        path="/api/v1/mixtures/{mixture_id}/predictions",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get predictions for a mixture',
                'description': 'Returns a list of predictions for a mixture',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'List of predictions'}}
            },
            'post': {
                'summary': 'Add a prediction for a mixture',
                'description': 'Add a new prediction for a mixture',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'201': {'description': 'Prediction created'}}
            }
        })
    )
    
    # Experiment endpoints
    spec.path(
        path="/api/v1/mixtures/{mixture_id}/experiments",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get experiments for a mixture',
                'description': 'Returns a list of experiments for a mixture',
                'parameters': [
                    {'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string', 'format': 'uuid'}}
                ],
                'responses': {
                    '200': {
                        'description': 'List of experiments',
                        'content': {
                            'application/json': {
                                'schema': ExperimentListSchema
                            }
                        }
                    }
                }
            },
            'post': {
                'summary': 'Record an experiment for a mixture',
                'description': 'Record a new experiment for a mixture',
                'parameters': [
                    {'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string', 'format': 'uuid'}}
                ],
                'requestBody': {
                    'required': True,
                    'content': {
                        'application/json': {
                            'schema': ExperimentCreateSchema
                        }
                    }
                },
                'responses': {
                    '201': {
                        'description': 'Experiment recorded successfully',
                        'content': {
                            'application/json': {
                                'schema': ExperimentSchema
                            }
                        }
                    },
                    '400': {
                        'description': 'Invalid request data'
                    }
                }
            }
        })
    )
    
    # Comparison endpoint
    spec.path(
        path="/api/v1/mixtures/{mixture_id}/compare",
        operations=add_standard_responses({
            'get': {
                'summary': 'Compare prediction with experiment',
                'description': 'Compare prediction with experiment for a mixture',
                'parameters': [{'name': 'mixture_id', 'in': 'path', 'required': True, 'schema': {'type': 'string'}}],
                'responses': {'200': {'description': 'Comparison results'}}
            }
        })
    )
    
    # Property comparison endpoint
    spec.path(
        path="/api/v1/compare-properties",
        operations=add_standard_responses({
            'post': {
                'summary': 'Compare properties of multiple entities',
                'description': 'Compare properties of molecules and mixtures side-by-side',
                'responses': {'200': {'description': 'Property comparison results'}}
            }
        })
    )
    
    # RDKit endpoints
    spec.path(
        path="/api/v1/rdkit/properties",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate properties for a molecule',
                'description': 'Calculate physicochemical properties for a molecule using RDKit',
                'responses': {'200': {'description': 'Calculated properties'}}
            }
        })
    )
    
    spec.path(
        path="/api/v1/rdkit/visualization",
        operations=add_standard_responses({
            'post': {
                'summary': 'Generate visualization for a molecule',
                'description': 'Generate SVG visualization for a molecule',
                'responses': {'200': {'description': 'Molecule visualization'}}
            }
        })
    )
    
    spec.path(
        path="/api/v1/rdkit/substructure",
        operations=add_standard_responses({
            'post': {
                'summary': 'Perform substructure search',
                'description': 'Search for a substructure within a molecule',
                'responses': {'200': {'description': 'Substructure search results'}}
            }
        })
    )
    
    spec.path(
        path="/api/v1/rdkit/similarity",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate similarity between molecules',
                'description': 'Calculate similarity between two molecules using various fingerprint methods',
                'responses': {'200': {'description': 'Similarity results'}}
            }
        })
    )
    
    # Batch operations
    spec.path(
        path="/api/v1/batch",
        operations=add_standard_responses({
            'post': {
                'summary': 'Perform batch operation',
                'description': 'Perform a batch operation on multiple entities',
                'responses': {'200': {'description': 'Batch operation results'}}
            }
        })
    )
    
    # Scoring endpoints
    spec.path(
        path="/api/v1/scoring/molecules",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate cryoprotection score for a molecule',
                'description': 'Calculate cryoprotection score for a molecule based on its properties',
                'responses': {'200': {'description': 'Molecule score'}}
            }
        })
    )
    
    spec.path(
        path="/api/v1/scoring/batch",
        operations=add_standard_responses({
            'post': {
                'summary': 'Calculate scores for multiple entities',
                'description': 'Calculate cryoprotection scores for multiple molecules or mixtures',
                'responses': {'200': {'description': 'Batch scoring results'}}
            }
        })
    )
    
    # User profile endpoint
    spec.path(
        path="/api/v1/user_profile",
        operations=add_standard_responses({
            'get': {
                'summary': 'Get user profile',
                'description': 'Get the current user\'s profile',
                'responses': {'200': {'description': 'User profile'}}
            },
            'post': {
                'summary': 'Create or update user profile',
                'description': 'Create or update the current user\'s profile',
                'responses': {'201': {'description': 'User profile created'}}
            },
            'put': {
                'summary': 'Update user profile',
                'description': 'Update the current user\'s profile',
                'responses': {'200': {'description': 'User profile updated'}}
            }
        })
    )
    
    app.config.update({
        'APISPEC_SPEC': spec,
        'APISPEC_SWAGGER_URL': '/api/docs/swagger/',
        'APISPEC_SWAGGER_UI_URL': '/api/docs/swagger-ui/',
        'APISPEC_SWAGGER_UI_BUNDLE_JS_PATH': app.config.get(
            'APISPEC_SWAGGER_UI_BUNDLE_JS_PATH', 
            'https://cdn.jsdelivr.net/npm/swagger-ui-dist@3/swagger-ui-bundle.js'
        ),
        'APISPEC_SWAGGER_UI_CSS_PATH': app.config.get(
            'APISPEC_SWAGGER_UI_CSS_PATH', 
            'https://cdn.jsdelivr.net/npm/swagger-ui-dist@3/swagger-ui.css'
        )
    })
    
    # Add routes for documentation
    register_doc_routes(app)
    
    # Initialize FlaskApiSpec
    docs = FlaskApiSpec(app)
    
    # Register the documentation blueprint
    app.register_blueprint(docs_bp)
    
    return docs

def register_doc_routes(app):
    """
    Register routes for API documentation.
    
    Args:
        app: Flask application instance
    """
    # Define the routes but don't register them yet
    def api_docs_index():
        """Render the API documentation index page."""
        return render_template('api_docs/index.html',
                              title=app.config.get('API_TITLE', 'CryoProtect API'),
                              version=app.config.get('API_VERSION', '1.0.0'))
    
    def get_openapi_json():
        """Return the OpenAPI specification as JSON."""
        return jsonify(app.config.get('APISPEC_SPEC').to_dict())
    
    def get_openapi_yaml():
        """Return the OpenAPI specification as YAML."""
        from apispec.yaml_utils import dict_to_yaml
        spec_dict = app.config.get('APISPEC_SPEC').to_dict()
        yaml_content = dict_to_yaml(spec_dict)
        return yaml_content, 200, {'Content-Type': 'text/yaml'}
    
    # Add the routes to the blueprint
    docs_bp.route('/')(api_docs_index)
    docs_bp.route('/openapi.json')(get_openapi_json)
    docs_bp.route('/openapi.yaml')(get_openapi_yaml)

def register_resource(docs, resource, endpoint=None, blueprint=None, **kwargs):
    """
    Register a resource with the API documentation.
    
    Args:
        docs: FlaskApiSpec instance
        resource: Resource class to document
        endpoint: Endpoint name (optional)
        blueprint: Blueprint name (optional)
        **kwargs: Additional arguments to pass to docs.register
    """
    import logging
    logger = logging.getLogger(__name__)
    
    try:
        # Register the resource directly without specifying an endpoint
        # This will use the resource class itself for documentation
        docs.register(resource, **kwargs)
        logger.info(f"Successfully registered resource {resource.__name__}")
    except Exception as e:
        # Log the error and continue
        logger.error(f"Error registering resource {resource.__name__}: {str(e)}")

def add_standard_responses(operations):
    """
    Add standard responses to API operations.
    
    Args:
        operations: Dictionary of operations to update
    
    Returns:
        Updated operations dictionary with standard responses
    """
    for operation in operations.values():
        if 'responses' not in operation:
            operation['responses'] = {}
            
        # Add 400 Bad Request if not present
        if '400' not in operation['responses']:
            operation['responses']['400'] = {
                'description': 'Bad Request - Invalid input data',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
            
        # Add 401 Unauthorized if not present
        if '401' not in operation['responses']:
            operation['responses']['401'] = {
                'description': 'Unauthorized - Authentication required',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
            
        # Add 404 Not Found for single resource operations
        if '{' in operation.get('path', ''):
            if '404' not in operation['responses']:
                operation['responses']['404'] = {
                    'description': 'Not Found - Resource does not exist',
                    'content': {
                        'application/json': {
                            'schema': ErrorResponseSchema
                        }
                    }
                }
                
        # Add 500 Server Error if not present
        if '500' not in operation['responses']:
            operation['responses']['500'] = {
                'description': 'Internal Server Error',
                'content': {
                    'application/json': {
                        'schema': ErrorResponseSchema
                    }
                }
            }
    
    return operations


# Security scheme definitions
def add_security_schemes(spec):
    """
    Add security scheme definitions to the API specification.
    
    Args:
        spec: APISpec instance
    """
    # JWT Bearer Authentication
    spec.components.security_scheme('bearerAuth', {
        'type': 'http',
        'scheme': 'bearer',
        'bearerFormat': 'JWT',
        'description': 'Enter JWT Bearer token in the format: Bearer {token}'
    })
    
    # API Key Authentication
    spec.components.security_scheme('apiKeyAuth', {
        'type': 'apiKey',
        'in': 'header',
        'name': 'X-API-Key',
        'description': 'API key authentication'
    })
    
    # Cookie Authentication
    spec.components.security_scheme('cookieAuth', {
        'type': 'apiKey',
        'in': 'cookie',
        'name': 'access_token',
        'description': 'Cookie-based authentication'
    })

# Decorator for documenting authentication requirements
def auth_required(security_type='bearerAuth'):
    """
    Decorator to document authentication requirements for an endpoint.
    
    Args:
        security_type: Type of security scheme to use (default: 'bearerAuth')
    
    Returns:
        Decorator function
    """
    def decorator(func):
        if not hasattr(func, '_apidoc'):
            func._apidoc = {}
        
        if 'security' not in func._apidoc:
            func._apidoc['security'] = []
        
        func._apidoc['security'].append({security_type: []})
        
        return func
    
    return decorator

# Helper function to generate OpenAPI YAML file
def generate_openapi_yaml(app, output_path='docs/api/openapi.yaml'):
    """
    Generate OpenAPI YAML file from the API specification.
    
    Args:
        app: Flask application instance
        output_path: Path to save the YAML file (default: 'docs/api/openapi.yaml')
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        from apispec.yaml_utils import dict_to_yaml
        
        # Get the API specification
        spec = app.config.get('APISPEC_SPEC')
        if not spec:
            return False
        
        # Convert to YAML
        yaml_content = dict_to_yaml(spec.to_dict())
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Write to file
        with open(output_path, 'w') as f:
            f.write(yaml_content)
        
        return True
    except Exception as e:
        current_app.logger.error(f"Error generating OpenAPI YAML: {str(e)}")
        return False
"""
CryoProtect Analyzer API - OpenAPI Documentation Generator

This module generates OpenAPI documentation for all API endpoints.
It provides a comprehensive documentation of the API, including:
- Endpoint paths and methods
- Request parameters and schemas
- Response schemas and examples
- Authentication requirements
- Error responses

The generated documentation follows the OpenAPI 3.0.2 specification.
"""

import os
import json
import yaml
import logging
from typing import Dict, Any, Optional
from flask import Flask, jsonify, current_app, Blueprint, send_from_directory
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin

from api.api_docs import generate_openapi_spec

# Set up logging
logger = logging.getLogger(__name__)

def create_openapi_blueprint(app: Flask) -> Blueprint:
    """
    Create a blueprint for OpenAPI documentation.
    
    Args:
        app: Flask application
        
    Returns:
        Blueprint for OpenAPI documentation
    """
    openapi_bp = Blueprint('openapi', __name__, url_prefix='/api/v1/docs')
    
    @openapi_bp.route('/openapi.json')
    def get_openapi_json():
        """Get OpenAPI specification in JSON format."""
        spec = generate_openapi_spec()
        return jsonify(spec)
    
    @openapi_bp.route('/openapi.yaml')
    def get_openapi_yaml():
        """Get OpenAPI specification in YAML format."""
        spec = generate_openapi_spec()
        yaml_spec = yaml.dump(spec, default_flow_style=False)
        return yaml_spec, 200, {'Content-Type': 'text/yaml'}
    
    @openapi_bp.route('/')
    def get_redoc():
        """Get ReDoc UI for API documentation."""
        return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>CryoProtect API Documentation</title>
            <meta charset="utf-8"/>
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <link href="https://fonts.googleapis.com/css?family=Montserrat:300,400,700|Roboto:300,400,700" rel="stylesheet">
            <style>
                body {{
                    margin: 0;
                    padding: 0;
                }}
            </style>
        </head>
        <body>
            <redoc spec-url='/api/v1/docs/openapi.json'></redoc>
            <script src="https://cdn.jsdelivr.net/npm/redoc@next/bundles/redoc.standalone.js"></script>
        </body>
        </html>
        """
    
    @openapi_bp.route('/swagger')
    def get_swagger():
        """Get Swagger UI for API documentation."""
        return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>CryoProtect API Documentation</title>
            <meta charset="utf-8"/>
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <link rel="stylesheet" type="text/css" href="https://unpkg.com/swagger-ui-dist@3/swagger-ui.css">
            <style>
                html {{
                    box-sizing: border-box;
                    overflow: -moz-scrollbars-vertical;
                    overflow-y: scroll;
                }}
                
                *,
                *:before,
                *:after {{
                    box-sizing: inherit;
                }}
                
                body {{
                    margin: 0;
                    background: #fafafa;
                }}
            </style>
        </head>
        <body>
            <div id="swagger-ui"></div>
            <script src="https://unpkg.com/swagger-ui-dist@3/swagger-ui-bundle.js"></script>
            <script src="https://unpkg.com/swagger-ui-dist@3/swagger-ui-standalone-preset.js"></script>
            <script>
                window.onload = function() {{
                    const ui = SwaggerUIBundle({{
                        url: '/api/v1/docs/openapi.json',
                        dom_id: '#swagger-ui',
                        deepLinking: true,
                        presets: [
                            SwaggerUIBundle.presets.apis,
                            SwaggerUIStandalonePreset
                        ],
                        plugins: [
                            SwaggerUIBundle.plugins.DownloadUrl
                        ],
                        layout: 'StandaloneLayout'
                    }});
                    window.ui = ui;
                }};
            </script>
        </body>
        </html>
        """
    
    return openapi_bp

def generate_openapi_file(output_path: str, format: str = 'json') -> None:
    """
    Generate OpenAPI specification file.
    
    Args:
        output_path: Path to save the file
        format: Format of the file ('json' or 'yaml')
    """
    spec = generate_openapi_spec()
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write to file
    if format.lower() == 'json':
        with open(output_path, 'w') as f:
            json.dump(spec, f, indent=2)
    elif format.lower() == 'yaml':
        with open(output_path, 'w') as f:
            yaml.dump(spec, f, default_flow_style=False)
    else:
        raise ValueError(f"Unsupported format: {format}")
    
    logger.info(f"Generated OpenAPI specification file: {output_path}")

def register_openapi_blueprint(app: Flask) -> None:
    """
    Register OpenAPI blueprint with the Flask app.
    
    Args:
        app: Flask application
    """
    openapi_bp = create_openapi_blueprint(app)
    app.register_blueprint(openapi_bp)
    logger.info("Registered OpenAPI blueprint")

def generate_openapi_files() -> None:
    """
    Generate OpenAPI specification files in both JSON and YAML formats.
    """
    # Generate JSON file
    generate_openapi_file('docs/openapi.json', 'json')
    
    # Generate YAML file
    generate_openapi_file('docs/openapi.yaml', 'yaml')
    
    logger.info("Generated OpenAPI specification files")
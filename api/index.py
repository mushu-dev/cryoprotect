#!/usr/bin/env python3
"""
CryoProtect v2 API Router for Vercel Serverless

This file provides a consolidated API router that handles all API requests on Vercel.
It maps routes to the appropriate handlers, reducing the number of serverless functions.
"""

import os
import sys
import json
import re
from http import HTTPStatus
from urllib.parse import parse_qs
from flask import Flask, request, Response, jsonify

# Add project directory to path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(current_dir)
sys.path.append(project_dir)

from app import create_app, app as flask_app
from api.api_standards import (
    create_standard_response,
    create_error_response,
    jsonify_standard_response
)
from api.utils import get_supabase_client, handle_supabase_error

# Initialize Flask app
app = create_app()
app.testing = True  # Required for proper request context handling

# Map URL patterns to handler functions
url_patterns = {
    r'^/api/v1/molecules(/.*)?$': 'handle_molecules',
    r'^/api/v1/mixtures(/.*)?$': 'handle_mixtures',
    r'^/api/v1/batch(/.*)?$': 'handle_batch',
    r'^/api/v1/rdkit(/.*)?$': 'handle_rdkit',
    r'^/api/v1/scoring(/.*)?$': 'handle_scoring',
    r'^/api/v1/compare-properties(/.*)?$': 'handle_compare_properties',
    r'^/api/v1/user_profile(/.*)?$': 'handle_user_profile',
    r'^/api/v1/protocol-designer(/.*)?$': 'handle_protocol_designer',
    r'^/api/v1/predictive-models(/.*)?$': 'handle_predictive_models',
    r'^/api/v1/verification(/.*)?$': 'handle_verification',
    r'^/api/v1/auth(/.*)?$': 'handle_auth',
    r'^/api/v1/teams(/.*)?$': 'handle_teams',
    r'^/api/v1/rbac(/.*)?$': 'handle_rbac',
    r'^/api/v1/docs(/.*)?$': 'handle_docs',
    r'^/api/v1/health(/.*)?$': 'handle_health',
}

def handler(event, context):
    """
    Main serverless handler function that processes all API requests.
    """
    try:
        # Extract path and HTTP method from the event
        path = event.get('path', '')
        http_method = event.get('httpMethod', 'GET')
        
        # Find matching handler for the path
        handler_name = None
        for pattern, handler_func in url_patterns.items():
            if re.match(pattern, path):
                handler_name = handler_func
                break
        
        # If no handler matches, return 404
        if not handler_name:
            return create_response(
                *create_error_response(
                    error="Endpoint not found",
                    status_code=HTTPStatus.NOT_FOUND,
                    context="The requested API endpoint does not exist"
                )
            )
        
        # Call appropriate handler based on the URL
        handler_func = globals().get(handler_name)
        if handler_func and callable(handler_func):
            return handler_func(event, context)
        else:
            return create_response(
                *create_error_response(
                    error="Internal server error",
                    status_code=HTTPStatus.INTERNAL_SERVER_ERROR,
                    context="Handler function not implemented"
                )
            )
    
    except Exception as e:
        # Log the error for debugging
        print(f"Error processing request: {str(e)}")
        
        # Return a standardized error response
        return create_response(
            *create_error_response(
                error=str(e),
                status_code=HTTPStatus.INTERNAL_SERVER_ERROR,
                context="Unexpected error processing request"
            )
        )

def create_response(response_data, status_code=200):
    """
    Create a standardized response object for Vercel Serverless Functions.
    
    Args:
        response_data: Response data (dict or JSON-serializable object)
        status_code: HTTP status code
        
    Returns:
        dict: Response object compatible with Vercel serverless functions
    """
    # Convert response data to JSON if it's not already a string
    if not isinstance(response_data, str):
        response_data = json.dumps(response_data)
    
    # Return response in format expected by Vercel
    return {
        'statusCode': status_code,
        'headers': {
            'Content-Type': 'application/json',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'Content-Type, Authorization',
            'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS'
        },
        'body': response_data
    }

def handle_molecules(event, context):
    """Handle requests to /api/v1/molecules endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_mixtures(event, context):
    """Handle requests to /api/v1/mixtures endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_batch(event, context):
    """Handle requests to /api/v1/batch endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_rdkit(event, context):
    """Handle requests to /api/v1/rdkit endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_scoring(event, context):
    """Handle requests to /api/v1/scoring endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_compare_properties(event, context):
    """Handle requests to /api/v1/compare-properties endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_user_profile(event, context):
    """Handle requests to /api/v1/user_profile endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_protocol_designer(event, context):
    """Handle requests to /api/v1/protocol-designer endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_predictive_models(event, context):
    """Handle requests to /api/v1/predictive-models endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_verification(event, context):
    """Handle requests to /api/v1/verification endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_auth(event, context):
    """Handle requests to /api/v1/auth endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_teams(event, context):
    """Handle requests to /api/v1/teams endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_rbac(event, context):
    """Handle requests to /api/v1/rbac endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_docs(event, context):
    """Handle requests to /api/v1/docs endpoints."""
    with app.test_request_context(event.get('path'), method=event.get('httpMethod')):
        # Set up the request context with the original request data
        request.args = parse_qs(event.get('queryStringParameters', {}))
        request.form = parse_qs(event.get('body', {}))
        
        # Handle the request using Flask routes
        response = flask_app.dispatch_request()
        
        # Return the response
        return create_response(response.get_data(as_text=True), response.status_code)

def handle_health(event, context):
    """Handle requests to /api/v1/health endpoints."""
    # Simple health check endpoint for monitoring
    try:
        # Check database connection
        supabase = get_supabase_client()
        response = supabase.table('property_types').select('*').limit(1).execute()
        db_status = "connected" if hasattr(response, 'data') else "unknown"
        
        # Return health status
        return create_response({
            'status': 'ok' if db_status == 'connected' else 'degraded',
            'version': os.environ.get('API_VERSION', '1.0.0'),
            'timestamp': datetime.now().isoformat(),
            'services': {
                'database': db_status
            }
        })
    except Exception as e:
        return create_response({
            'status': 'error',
            'version': os.environ.get('API_VERSION', '1.0.0'),
            'error': str(e),
            'timestamp': datetime.now().isoformat()
        }, 500)
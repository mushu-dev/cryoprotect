#!/usr/bin/env python3
"""
Example script demonstrating JWT-based service role authentication.

This script shows how to:
1. Create and validate service tokens
2. Make authenticated API requests using the service client
3. Secure API endpoints with service role authentication
"""

import os
import sys
import json
import time
import logging
import argparse
from typing import Dict, Any

from flask import Flask, jsonify, request

# Import authentication components
from auth import (
    create_service_token,
    validate_service_token,
    revoke_service_token,
    ServiceRoleClient,
    service_role_required,
    ServiceRoleAuth
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_example_app():
    """Create a Flask application with protected endpoints."""
    app = Flask(__name__)

    # Initialize service role authentication
    service_auth = ServiceRoleAuth(app)

    @app.route('/api/public', methods=['GET'])
    def public_endpoint():
        """Public endpoint that doesn't require authentication."""
        return jsonify({
            'status': 'success',
            'message': 'This is a public endpoint'
        })

    @app.route('/api/service/data', methods=['GET'])
    @service_role_required(scopes=['data'])
    def service_data_endpoint():
        """Protected endpoint that requires data scope."""
        return jsonify({
            'status': 'success',
            'message': 'This is a protected data endpoint',
            'service': request.g.service_claims['client']
        })

    @app.route('/api/service/admin', methods=['GET'])
    @service_role_required(scopes=['admin'])
    def service_admin_endpoint():
        """Protected endpoint that requires admin scope."""
        return jsonify({
            'status': 'success',
            'message': 'This is a protected admin endpoint',
            'service': request.g.service_claims['client']
        })

    @app.route('/api/service/analytics', methods=['GET', 'POST'])
    @service_role_required(scopes=['analytics'])
    def service_analytics_endpoint():
        """Protected endpoint that requires analytics scope."""
        if request.method == 'GET':
            return jsonify({
                'status': 'success',
                'message': 'This is a protected analytics GET endpoint',
                'service': request.g.service_claims['client']
            })
        else:
            data = request.json or {}
            return jsonify({
                'status': 'success',
                'message': 'This is a protected analytics POST endpoint',
                'service': request.g.service_claims['client'],
                'received_data': data
            })

    @app.route('/api/service/multi-role', methods=['GET'])
    @service_role_required(scopes=['data', 'admin', 'analytics'])
    def service_multi_role_endpoint():
        """Protected endpoint that allows multiple scopes."""
        return jsonify({
            'status': 'success',
            'message': 'This is a protected multi-role endpoint',
            'service': request.g.service_claims['client']
        })

    return app

def token_examples():
    """Run token creation and validation examples."""
    logger.info("=== JWT Token Examples ===")

    # Create tokens for different services and scopes
    data_token = create_service_token(
        service_name="data-service",
        scopes=["data"],
        expiration=3600  # 1 hour
    )

    admin_token = create_service_token(
        service_name="admin-service",
        scopes=["admin"],
        expiration=1800  # 30 minutes
    )

    analytics_token = create_service_token(
        service_name="analytics-service",
        scopes=["analytics"],
        expiration=7200  # 2 hours
    )

    # Print tokens
    logger.info("Generated Service Tokens:")
    logger.info(f"Data Service: {data_token}")
    logger.info(f"Admin Service: {admin_token}")
    logger.info(f"Analytics Service: {analytics_token}")

    # Validate tokens
    logger.info("\nValidating Tokens:")

    # Validate data token
    try:
        data_claims = validate_service_token(data_token)
        logger.info(f"Data Token Valid: {json.dumps(data_claims, indent=2)}")
    except Exception as e:
        logger.error(f"Data Token Validation Error: {str(e)}")

    # Validate admin token
    try:
        admin_claims = validate_service_token(admin_token)
        logger.info(f"Admin Token Valid: {json.dumps(admin_claims, indent=2)}")
    except Exception as e:
        logger.error(f"Admin Token Validation Error: {str(e)}")

    # Validate analytics token
    try:
        analytics_claims = validate_service_token(analytics_token)
        logger.info(f"Analytics Token Valid: {json.dumps(analytics_claims, indent=2)}")
    except Exception as e:
        logger.error(f"Analytics Token Validation Error: {str(e)}")

    # Revoke a token
    logger.info("\nRevoking Admin Token:")
    try:
        revoke_service_token(admin_token)
        logger.info("Admin Token Revoked")

        # Try to validate the revoked token
        validate_service_token(admin_token)
        logger.error("Error: Revoked token validated successfully!")
    except Exception as e:
        logger.info(f"Expected Error: {str(e)}")

    return {
        "data_token": data_token,
        "admin_token": admin_token,
        "analytics_token": analytics_token
    }

def client_examples(tokens: Dict[str, str], base_url: str = "http://localhost:5000"):
    """Run service client examples."""
    logger.info("\n=== Service Client Examples ===")

    # Initialize clients for different services
    data_client = ServiceRoleClient(
        base_url=base_url,
        service_name="data-service",
        scopes=["data"]
    )

    admin_client = ServiceRoleClient(
        base_url=base_url,
        service_name="admin-service",
        scopes=["admin"]
    )

    analytics_client = ServiceRoleClient(
        base_url=base_url,
        service_name="analytics-service",
        scopes=["analytics"]
    )

    # Make requests to different endpoints
    logger.info("\nMaking Requests to Public Endpoint:")

    try:
        response = data_client.get("/api/public")
        logger.info(f"Data Client - Public: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Data Client Error: {str(e)}")

    logger.info("\nMaking Requests to Data Endpoint:")

    try:
        response = data_client.get("/api/service/data")
        logger.info(f"Data Client - Data Endpoint: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Data Client Error: {str(e)}")

    try:
        response = admin_client.get("/api/service/data")
        logger.info(f"Admin Client - Data Endpoint: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Admin Client Error: {str(e)}")

    logger.info("\nMaking Requests to Admin Endpoint:")

    try:
        response = admin_client.get("/api/service/admin")
        logger.info(f"Admin Client - Admin Endpoint: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Admin Client Error: {str(e)}")

    try:
        response = data_client.get("/api/service/admin")
        logger.info(f"Data Client - Admin Endpoint: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Data Client Error: {str(e)}")

    logger.info("\nMaking Requests to Analytics Endpoint:")

    try:
        response = analytics_client.get("/api/service/analytics")
        logger.info(f"Analytics Client - Analytics GET: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Analytics Client Error: {str(e)}")

    try:
        data = {"query": "sample analytics data", "time_range": "last_week"}
        response = analytics_client.post("/api/service/analytics", json_data=data)
        logger.info(f"Analytics Client - Analytics POST: {response.status_code}")
        logger.info(f"Response: {json.dumps(response.json(), indent=2)}")
    except Exception as e:
        logger.error(f"Analytics Client Error: {str(e)}")

def run_app(host="localhost", port=5000):
    """Run the example Flask application."""
    app = create_example_app()
    logger.info(f"\nStarting example server at http://{host}:{port}")
    app.run(host=host, port=port, debug=True)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='JWT Service Auth Example')
    parser.add_argument('--tokens', action='store_true', help='Run token examples')
    parser.add_argument('--client', action='store_true', help='Run client examples')
    parser.add_argument('--server', action='store_true', help='Run example server')
    parser.add_argument('--all', action='store_true', help='Run all examples')
    parser.add_argument('--host', default='localhost', help='Server host')
    parser.add_argument('--port', type=int, default=5000, help='Server port')

    args = parser.parse_args()

    # Default to all if no arguments specified
    run_all = args.all or not (args.tokens or args.client or args.server)

    if args.tokens or run_all:
        tokens = token_examples()
    else:
        tokens = {}

    if args.client or run_all:
        base_url = f"http://{args.host}:{args.port}"
        client_examples(tokens, base_url)

    if args.server or run_all:
        run_app(args.host, args.port)

    return 0

if __name__ == "__main__":
    main()
"""
Tests for API fixtures.

This module contains tests that demonstrate how to use the API fixtures
for testing API endpoints with different authentication scenarios.
"""

import pytest
import json
from flask import Flask, jsonify, request, g
from flask_restful import Api, Resource

from tests.fixtures.api.client import (
    app, client, api_client, authenticated_client,
    admin_client, scientist_client, mock_api_client, APIClient, MockResponse
)
from tests.fixtures.api.auth import (
    auth_token, admin_token, user_token, scientist_token,
    expired_token, invalid_token, mock_current_user,
    mock_admin_user, mock_regular_user, mock_scientist_user,
    mock_unauthenticated_user
)


# Test API endpoints for demonstration
class TestAPI:
    """Test API for demonstrating fixtures."""
    
    @pytest.fixture
    def test_app(self, scientist_token):
        """Create a test Flask app with API endpoints.

        NOTE: scientist_token is passed in to ensure the endpoint and test use the same token instance.
        """
        app = Flask(__name__)
        api = Api(app)
        
        # Public endpoint
        @app.route('/api/public', methods=['GET'])
        def public_endpoint():
            return jsonify({'message': 'This is a public endpoint'})
        
        # Protected endpoint requiring authentication
        @app.route('/api/protected', methods=['GET'])
        def protected_endpoint():
            # Check for Authorization header
            auth_header = request.headers.get('Authorization')
            if not auth_header or not auth_header.startswith('Bearer '):
                return jsonify({'error': 'Authentication required'}), 401
            
            return jsonify({'message': 'This is a protected endpoint'})
        
        # Admin-only endpoint
        @app.route('/api/admin', methods=['GET'])
        def admin_endpoint():
            # Check for Authorization header
            auth_header = request.headers.get('Authorization')
            if not auth_header or not auth_header.startswith('Bearer '):
                return jsonify({'error': 'Authentication required'}), 401
            
            # In a real app, we would verify the token and check the role
            # For this test, we'll check if it's the admin token
            token = auth_header.split(' ')[1]
            if token != admin_token:
                return jsonify({'error': 'Admin access required'}), 403
            
            return jsonify({'message': 'This is an admin endpoint'})
        
        # Scientist-only endpoint
        @app.route('/api/scientist', methods=['GET'])
        def scientist_endpoint():
            # Check for Authorization header
            auth_header = request.headers.get('Authorization')
            if not auth_header or not auth_header.startswith('Bearer '):
                return jsonify({'error': 'Authentication required'}), 401
            
            # In a real app, we would verify the token and check the role
            # For this test, we'll check if it's the scientist token
            token = auth_header.split(' ')[1]
            if token != scientist_token:
                return jsonify({'error': 'Scientist access required'}), 403
            
            return jsonify({'message': 'This is a scientist endpoint'})
        
        # JSON request/response endpoint
        @app.route('/api/data', methods=['POST'])
        def data_endpoint():
            # Get JSON data from request
            data = request.get_json()
            if not data:
                return jsonify({'error': 'JSON data required'}), 400
            
            # Process data
            result = {
                'received': data,
                'processed': True
            }
            
            return jsonify(result)
        
        return app
    
    @pytest.fixture
    def test_client(self, test_app):
        """Create a test client for the test app."""
        return test_app.test_client()
    
    @pytest.fixture
    def local_api_client(self, test_client):
        """APIClient for the local test app."""
        from tests.fixtures.api.client import APIClient
        return APIClient(test_client)

    def test_public_endpoint_unauthenticated(self, local_api_client):
        """Test accessing a public endpoint without authentication."""
        response = local_api_client.get('/api/public')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['message'] == 'This is a public endpoint'

    def test_protected_endpoint_unauthenticated(self, local_api_client):
        """Test accessing a protected endpoint without authentication."""
        response = local_api_client.get('/api/protected')
        assert response.status_code == 401
        data = json.loads(response.data)
        assert 'error' in data

    # For authenticated tests, you would need to extend local_api_client to accept tokens, or mock the authentication as appropriate.

    def test_admin_endpoint_as_regular_user(self, local_regular_api_client):
        """Test accessing an admin endpoint as a regular user."""
        response = local_regular_api_client.get('/api/admin')
        assert response.status_code == 403
        data = json.loads(response.data)
        assert 'error' in data
    
    @pytest.fixture
    def local_scientist_api_client(self, test_client, scientist_token):
        """Scientist APIClient for the local test app."""
        from tests.fixtures.api.client import APIClient
        return APIClient(test_client, scientist_token)

    def test_scientist_endpoint_as_scientist(self, local_scientist_api_client):
        """Test accessing a scientist endpoint as a scientist user."""
        response = local_scientist_api_client.get('/api/scientist')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['message'] == 'This is a scientist endpoint'

    @pytest.fixture
    def local_regular_api_client(self, test_client, user_token):
        """Regular user APIClient for the local test app."""
        from tests.fixtures.api.client import APIClient
        return APIClient(test_client, user_token)

    def test_scientist_endpoint_as_regular_user(self, local_regular_api_client):
        """Test accessing a scientist endpoint as a regular user."""
        response = local_regular_api_client.get('/api/scientist')
        assert response.status_code == 403
        data = json.loads(response.data)
        assert 'error' in data
    
    @pytest.fixture
    def local_authenticated_api_client(self, test_client, auth_token):
        """Authenticated APIClient for the local test app."""
        from tests.fixtures.api.client import APIClient
        return APIClient(test_client, auth_token)

    def test_json_request_response(self, local_authenticated_api_client):
        """Test JSON request/response handling."""
        test_data = {
            'name': 'Test User',
            'email': 'test@example.com',
            'data': [1, 2, 3]
        }
        response = local_authenticated_api_client.post('/api/data', json_data=test_data)
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['received'] == test_data
        assert data['processed'] is True

    def test_mock_api_client(self, mock_api_client):
        """Test using the mock API client."""
        mock_api_client.responses = {
            '/api/test': {'message': 'Test response'},
            '/api/data': lambda data, **kwargs: MockResponse({'echo': data})
        }
        response = mock_api_client.get('/api/test')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['message'] == 'Test response'
        
        # Make request to data endpoint with dynamic response
        test_data = {'test': 'data'}
        # Use mock_api_client for POST as well
        response = mock_api_client.post('/api/data', json_data=test_data)
        
        # Verify response
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['echo'] == test_data
        
        # Verify request was recorded
        assert len(mock_api_client.requests) == 2
        assert mock_api_client.requests[0]['method'] == 'GET'
        assert mock_api_client.requests[0]['endpoint'] == '/api/test'
        assert mock_api_client.requests[1]['method'] == 'POST'
        assert mock_api_client.requests[1]['endpoint'] == '/api/data'
        assert mock_api_client.requests[1]['json_data'] == test_data


# Integration with Flask app tests
class TestFlaskAppIntegration:
    """Test integration with Flask app."""

    @pytest.mark.skip(reason="Health endpoint requires valid Supabase API key; skip in unit tests.")
    def test_api_client_with_flask_app(self, app, client):
        """Test using the API client with a Flask app."""
        # Create an API client with the Flask test client
        api_client_instance = APIClient(client)
        
        # Make request to health endpoint
        response = api_client_instance.get('/health')
        
        # Verify response
        assert response.status_code == 200
        data = json.loads(response.data)
        assert 'status' in data

    def test_authenticated_client_with_flask_app(self, app, client, auth_token):
        """Test using the authenticated API client with a Flask app."""
        # Create an authenticated API client with the Flask test client
        authenticated_client_instance = APIClient(client, auth_token)
        
        # Add a test endpoint to the app
        @app.route('/api/test-auth', methods=['GET'])
        def test_auth_endpoint():
            # Check for Authorization header
            auth_header = request.headers.get('Authorization')
            if not auth_header or not auth_header.startswith('Bearer '):
                return jsonify({'error': 'Authentication required'}), 401
            
            return jsonify({'message': 'Authenticated'})
        
        # Make request to test endpoint
        response = authenticated_client_instance.get('/api/test-auth')
        
        # Verify response
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['message'] == 'Authenticated'
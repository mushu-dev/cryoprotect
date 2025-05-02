"""
API client fixtures for testing.

This module provides fixtures for API client testing, including
authenticated and unauthenticated clients for different user roles.
"""

import pytest
import json
from flask import Flask
from flask.testing import FlaskClient
from typing import Dict, Any, Generator, Optional

from app import create_app
from tests.fixtures.api.auth import (
    auth_token, admin_token, user_token, scientist_token
)


class APIClient:
    """Wrapper around Flask test client with helper methods for API testing."""
    
    def __init__(self, client: FlaskClient, token: Optional[str] = None):
        """
        Initialize the API client.
        
        Args:
            client: Flask test client
            token: Optional authentication token
        """
        self.client = client
        self.token = token
        self.headers = {}
        
        if token:
            self.headers['Authorization'] = f'Bearer {token}'
    
    def get(self, endpoint: str, params: Dict[str, Any] = None, **kwargs) -> Any:
        """
        Make a GET request to the API.
        
        Args:
            endpoint: API endpoint
            params: Query parameters
            **kwargs: Additional arguments to pass to the client
            
        Returns:
            Response object
        """
        headers = kwargs.pop('headers', {})
        headers.update(self.headers)
        
        return self.client.get(
            endpoint,
            query_string=params,
            headers=headers,
            **kwargs
        )
    
    def post(self, endpoint: str, data: Dict[str, Any] = None, json_data: Dict[str, Any] = None, **kwargs) -> Any:
        """
        Make a POST request to the API.
        
        Args:
            endpoint: API endpoint
            data: Form data
            json_data: JSON data
            **kwargs: Additional arguments to pass to the client
            
        Returns:
            Response object
        """
        headers = kwargs.pop('headers', {})
        headers.update(self.headers)
        
        if json_data is not None:
            return self.client.post(
                endpoint,
                json=json_data,
                headers=headers,
                **kwargs
            )
        else:
            return self.client.post(
                endpoint,
                data=data,
                headers=headers,
                **kwargs
            )
    
    def put(self, endpoint: str, data: Dict[str, Any] = None, json_data: Dict[str, Any] = None, **kwargs) -> Any:
        """
        Make a PUT request to the API.
        
        Args:
            endpoint: API endpoint
            data: Form data
            json_data: JSON data
            **kwargs: Additional arguments to pass to the client
            
        Returns:
            Response object
        """
        headers = kwargs.pop('headers', {})
        headers.update(self.headers)
        
        if json_data is not None:
            return self.client.put(
                endpoint,
                json=json_data,
                headers=headers,
                **kwargs
            )
        else:
            return self.client.put(
                endpoint,
                data=data,
                headers=headers,
                **kwargs
            )
    
    def delete(self, endpoint: str, **kwargs) -> Any:
        """
        Make a DELETE request to the API.
        
        Args:
            endpoint: API endpoint
            **kwargs: Additional arguments to pass to the client
            
        Returns:
            Response object
        """
        headers = kwargs.pop('headers', {})
        headers.update(self.headers)
        
        return self.client.delete(
            endpoint,
            headers=headers,
            **kwargs
        )
    
    def parse_json(self, response) -> Dict[str, Any]:
        """
        Parse JSON response.
        
        Args:
            response: Response object
            
        Returns:
            Parsed JSON data
        """
        return json.loads(response.data)


@pytest.fixture
def app() -> Flask:
    """
    Create a Flask app for testing.
    
    Returns:
        Flask app
    """
    app = create_app(testing=True)
    return app


@pytest.fixture
def client(app) -> FlaskClient:
    """
    Create a Flask test client.
    
    Args:
        app: Flask app
        
    Returns:
        Flask test client
    """
    return app.test_client()


@pytest.fixture
def api_client(client) -> APIClient:
    """
    Create an unauthenticated API client.
    
    Args:
        client: Flask test client
        
    Returns:
        API client
    """
    return APIClient(client)


@pytest.fixture
def authenticated_client(client, auth_token) -> APIClient:
    """
    Create an authenticated API client with a regular user token.
    
    Args:
        client: Flask test client
        auth_token: Authentication token
        
    Returns:
        Authenticated API client
    """
    return APIClient(client, auth_token)


@pytest.fixture
def admin_client(client, admin_token) -> APIClient:
    """
    Create an authenticated API client with an admin token.
    
    Args:
        client: Flask test client
        admin_token: Admin authentication token
        
    Returns:
        Admin API client
    """
    return APIClient(client, admin_token)


@pytest.fixture
def scientist_client(client, scientist_token) -> APIClient:
    """
    Create an authenticated API client with a scientist token.
    
    Args:
        client: Flask test client
        scientist_token: Scientist authentication token
        
    Returns:
        Scientist API client
    """
    return APIClient(client, scientist_token)


class MockResponse:
    """Mock response object for testing."""
    
    def __init__(self, data: Dict[str, Any], status_code: int = 200):
        """
        Initialize the mock response.
        
        Args:
            data: Response data
            status_code: HTTP status code
        """
        self.data = json.dumps(data).encode('utf-8')
        self.status_code = status_code
        self.headers = {}
    
    def get_json(self) -> Dict[str, Any]:
        """
        Get JSON data from the response.
        
        Returns:
            Parsed JSON data
        """
        return json.loads(self.data)


class MockAPIClient:
    """Mock API client for testing."""
    
    def __init__(self, responses: Dict[str, Any] = None):
        """
        Initialize the mock API client.
        
        Args:
            responses: Dictionary of endpoint to response mappings
        """
        self.responses = responses or {}
        self.requests = []
    
    def get(self, endpoint: str, params: Dict[str, Any] = None, **kwargs) -> MockResponse:
        """
        Mock GET request.
        
        Args:
            endpoint: API endpoint
            params: Query parameters
            **kwargs: Additional arguments
            
        Returns:
            Mock response
        """
        self.requests.append({
            'method': 'GET',
            'endpoint': endpoint,
            'params': params,
            'kwargs': kwargs
        })
        
        if endpoint in self.responses:
            response_data = self.responses[endpoint]
            if callable(response_data):
                return response_data(params, **kwargs)
            else:
                return MockResponse(response_data)
        
        return MockResponse({'error': 'Not found'}, 404)
    
    def post(self, endpoint: str, data: Dict[str, Any] = None, json_data: Dict[str, Any] = None, **kwargs) -> MockResponse:
        """
        Mock POST request.
        
        Args:
            endpoint: API endpoint
            data: Form data
            json_data: JSON data
            **kwargs: Additional arguments
            
        Returns:
            Mock response
        """
        self.requests.append({
            'method': 'POST',
            'endpoint': endpoint,
            'data': data,
            'json_data': json_data,
            'kwargs': kwargs
        })
        
        if endpoint in self.responses:
            response_data = self.responses[endpoint]
            if callable(response_data):
                return response_data(data or json_data, **kwargs)
            else:
                return MockResponse(response_data)
        
        return MockResponse({'error': 'Not found'}, 404)
    
    def put(self, endpoint: str, data: Dict[str, Any] = None, json_data: Dict[str, Any] = None, **kwargs) -> MockResponse:
        """
        Mock PUT request.
        
        Args:
            endpoint: API endpoint
            data: Form data
            json_data: JSON data
            **kwargs: Additional arguments
            
        Returns:
            Mock response
        """
        self.requests.append({
            'method': 'PUT',
            'endpoint': endpoint,
            'data': data,
            'json_data': json_data,
            'kwargs': kwargs
        })
        
        if endpoint in self.responses:
            response_data = self.responses[endpoint]
            if callable(response_data):
                return response_data(data or json_data, **kwargs)
            else:
                return MockResponse(response_data)
        
        return MockResponse({'error': 'Not found'}, 404)
    
    def delete(self, endpoint: str, **kwargs) -> MockResponse:
        """
        Mock DELETE request.
        
        Args:
            endpoint: API endpoint
            **kwargs: Additional arguments
            
        Returns:
            Mock response
        """
        self.requests.append({
            'method': 'DELETE',
            'endpoint': endpoint,
            'kwargs': kwargs
        })
        
        if endpoint in self.responses:
            response_data = self.responses[endpoint]
            if callable(response_data):
                return response_data(**kwargs)
            else:
                return MockResponse(response_data)
        
        return MockResponse({'error': 'Not found'}, 404)
    
    def parse_json(self, response) -> Dict[str, Any]:
        """
        Parse JSON response.
        
        Args:
            response: Response object
            
        Returns:
            Parsed JSON data
        """
        return json.loads(response.data)


@pytest.fixture
def mock_api_client() -> MockAPIClient:
    """
    Create a mock API client.
    
    Returns:
        Mock API client
    """
    return MockAPIClient()
"""
CSRF Test Fixtures

This module provides fixtures for testing CSRF protection in the CryoProtect API.
"""

import pytest
import json
from flask import session

@pytest.fixture
def csrf_token(client):
    """
    Get a CSRF token from the API.
    
    Args:
        client: Flask test client
        
    Returns:
        str: CSRF token
    """
    response = client.get('/api/v1/csrf-token')
    data = json.loads(response.data)
    return data['csrf_token']

@pytest.fixture
def csrf_headers(csrf_token):
    """
    Create headers with CSRF token.
    
    Args:
        csrf_token: CSRF token from the csrf_token fixture
        
    Returns:
        dict: Headers with CSRF token
    """
    return {'X-CSRF-Token': csrf_token}

@pytest.fixture
def csrf_client(client, csrf_token):
    """
    Create a client with CSRF token in headers.
    
    This fixture creates a wrapper around the Flask test client that
    automatically includes the CSRF token in all requests.
    
    Args:
        client: Flask test client
        csrf_token: CSRF token from the csrf_token fixture
        
    Returns:
        object: Client wrapper with CSRF token
    """
    class CSRFClientWrapper:
        def __init__(self, client, token):
            self.client = client
            self.token = token
            
        def get(self, *args, **kwargs):
            return self.client.get(*args, **kwargs)
            
        def post(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.post(*args, **kwargs)
            
        def put(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.put(*args, **kwargs)
            
        def patch(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.patch(*args, **kwargs)
            
        def delete(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.delete(*args, **kwargs)
    
    return CSRFClientWrapper(client, csrf_token)

@pytest.fixture
def authenticated_csrf_client(authenticated_client, csrf_token):
    """
    Create an authenticated client with CSRF token in headers.
    
    This fixture creates a wrapper around the authenticated Flask test client
    that automatically includes the CSRF token in all requests.
    
    Args:
        authenticated_client: Authenticated Flask test client
        csrf_token: CSRF token from the csrf_token fixture
        
    Returns:
        object: Authenticated client wrapper with CSRF token
    """
    class CSRFClientWrapper:
        def __init__(self, client, token):
            self.client = client
            self.token = token
            
        def get(self, *args, **kwargs):
            return self.client.get(*args, **kwargs)
            
        def post(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.post(*args, **kwargs)
            
        def put(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.put(*args, **kwargs)
            
        def patch(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.patch(*args, **kwargs)
            
        def delete(self, *args, **kwargs):
            headers = kwargs.get('headers', {})
            headers['X-CSRF-Token'] = self.token
            kwargs['headers'] = headers
            return self.client.delete(*args, **kwargs)
    
    return CSRFClientWrapper(authenticated_client, csrf_token)

@pytest.fixture
def disable_csrf(app):
    """
    Disable CSRF protection for testing.
    
    This fixture disables CSRF protection for tests that need to bypass it.
    
    Args:
        app: Flask application
        
    Returns:
        None
    """
    app.config['CSRF_DISABLED'] = True
    yield
    app.config['CSRF_DISABLED'] = False
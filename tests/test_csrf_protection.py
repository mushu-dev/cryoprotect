"""
Tests for CSRF protection in the CryoProtect API.

These tests verify that CSRF protection is properly implemented for all
state-changing API endpoints (POST, PUT, PATCH, DELETE).
"""

import pytest
import json
from flask import url_for, session

def test_csrf_token_endpoint(client):
    """Test that the CSRF token endpoint returns a valid token."""
    response = client.get('/api/v1/csrf-token')
    assert response.status_code == 200
    data = json.loads(response.data)
    assert 'csrf_token' in data
    assert len(data['csrf_token']) > 0
    
    # Check that the token is also set in a cookie
    assert 'csrf_token' in response.cookies
    assert response.cookies['csrf_token'] == data['csrf_token']

def test_csrf_token_in_session(client):
    """Test that the CSRF token is stored in the session."""
    client.get('/api/v1/csrf-token')
    with client.session_transaction() as sess:
        assert 'csrf_token' in sess
        assert len(sess['csrf_token']) > 0

def test_post_without_csrf_token(client):
    """Test that POST requests without CSRF token are rejected."""
    # Try to create a new molecule without CSRF token
    response = client.post('/api/v1/molecules', json={
        'name': 'Test Molecule',
        'smiles': 'CCO',
        'description': 'Test molecule for CSRF protection'
    })
    assert response.status_code == 403
    data = json.loads(response.data)
    assert 'error' in data
    assert 'CSRF token validation failed' in data['error']

def test_post_with_csrf_token_in_header(client):
    """Test that POST requests with CSRF token in header are accepted."""
    # Get CSRF token
    token_response = client.get('/api/v1/csrf-token')
    token_data = json.loads(token_response.data)
    csrf_token = token_data['csrf_token']
    
    # Create a new molecule with CSRF token in header
    response = client.post('/api/v1/molecules', json={
        'name': 'Test Molecule',
        'smiles': 'CCO',
        'description': 'Test molecule for CSRF protection'
    }, headers={'X-CSRF-Token': csrf_token})
    
    # The request might fail for other reasons (e.g., validation),
    # but it should not fail due to CSRF protection
    assert response.status_code != 403
    data = json.loads(response.data)
    assert 'CSRF token validation failed' not in str(data)

def test_post_with_csrf_token_in_json(client):
    """Test that POST requests with CSRF token in JSON body are accepted."""
    # Get CSRF token
    token_response = client.get('/api/v1/csrf-token')
    token_data = json.loads(token_response.data)
    csrf_token = token_data['csrf_token']
    
    # Create a new molecule with CSRF token in JSON body
    response = client.post('/api/v1/molecules', json={
        'name': 'Test Molecule',
        'smiles': 'CCO',
        'description': 'Test molecule for CSRF protection',
        'csrf_token': csrf_token
    })
    
    # The request might fail for other reasons (e.g., validation),
    # but it should not fail due to CSRF protection
    assert response.status_code != 403
    data = json.loads(response.data)
    assert 'CSRF token validation failed' not in str(data)

def test_put_without_csrf_token(client):
    """Test that PUT requests without CSRF token are rejected."""
    # Try to update a molecule without CSRF token
    response = client.put('/api/v1/molecules/test-id', json={
        'name': 'Updated Molecule',
        'description': 'Updated description'
    })
    assert response.status_code == 403
    data = json.loads(response.data)
    assert 'error' in data
    assert 'CSRF token validation failed' in data['error']

def test_delete_without_csrf_token(client):
    """Test that DELETE requests without CSRF token are rejected."""
    # Try to delete a molecule without CSRF token
    response = client.delete('/api/v1/molecules/test-id')
    assert response.status_code == 403
    data = json.loads(response.data)
    assert 'error' in data
    assert 'CSRF token validation failed' in data['error']

def test_get_without_csrf_token(client):
    """Test that GET requests without CSRF token are accepted."""
    # GET requests should not require CSRF token
    response = client.get('/api/v1/molecules')
    assert response.status_code != 403
    data = json.loads(response.data)
    assert 'CSRF token validation failed' not in str(data)

def test_csrf_exempt_endpoint(client, app):
    """Test that CSRF-exempt endpoints don't require CSRF token."""
    # Create a test endpoint that is exempt from CSRF protection
    @app.route('/api/v1/test/csrf-exempt', methods=['POST'])
    @csrf_exempt()
    def csrf_exempt_endpoint():
        return {'message': 'CSRF exempt endpoint'}, 200
    
    # POST to the exempt endpoint without CSRF token
    response = client.post('/api/v1/test/csrf-exempt', json={
        'test': 'data'
    })
    assert response.status_code == 200
    data = json.loads(response.data)
    assert data['message'] == 'CSRF exempt endpoint'

def test_csrf_token_rotation(client):
    """Test that CSRF tokens are rotated after expiry."""
    # Get initial CSRF token
    token_response = client.get('/api/v1/csrf-token')
    token_data = json.loads(token_response.data)
    initial_token = token_data['csrf_token']
    
    # Simulate token expiry by manipulating the session
    with client.session_transaction() as sess:
        sess['csrf_token_expiry'] = 0  # Set to expired
    
    # Get new CSRF token
    token_response = client.get('/api/v1/csrf-token')
    token_data = json.loads(token_response.data)
    new_token = token_data['csrf_token']
    
    # Tokens should be different
    assert initial_token != new_token

def test_csrf_token_validation_with_expired_token(client):
    """Test that expired CSRF tokens are rejected."""
    # Get CSRF token
    token_response = client.get('/api/v1/csrf-token')
    token_data = json.loads(token_response.data)
    csrf_token = token_data['csrf_token']
    
    # Simulate token expiry by manipulating the session
    with client.session_transaction() as sess:
        sess['csrf_token_expiry'] = 0  # Set to expired
    
    # Try to use the expired token
    response = client.post('/api/v1/molecules', json={
        'name': 'Test Molecule',
        'smiles': 'CCO',
        'description': 'Test molecule for CSRF protection'
    }, headers={'X-CSRF-Token': csrf_token})
    
    assert response.status_code == 403
    data = json.loads(response.data)
    assert 'error' in data
    assert 'CSRF token validation failed' in data['error']
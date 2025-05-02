# API Fixtures for CryoProtect v2

This directory contains fixtures for testing API endpoints in the CryoProtect v2 application. These fixtures provide standardized tools for testing API functionality with different authentication scenarios and user roles.

## Overview

The API fixtures are designed to:

1. Simplify testing of API endpoints with different authentication scenarios
2. Support testing with different user roles (admin, regular user, scientist)
3. Provide utilities for handling JSON requests and responses
4. Enable easy mocking of API responses for unit testing

## Available Fixtures

### Client Fixtures

- `api_client`: Unauthenticated API client for testing public endpoints
- `authenticated_client`: Authenticated API client with a regular user token
- `admin_client`: Authenticated API client with an admin user token
- `scientist_client`: Authenticated API client with a scientist user token
- `mock_api_client`: Mock API client for unit testing without a Flask app

### Authentication Fixtures

- `auth_token`: Regular user authentication token
- `admin_token`: Admin user authentication token
- `user_token`: Alias for auth_token
- `scientist_token`: Scientist user authentication token
- `expired_token`: Expired authentication token for testing token expiration
- `invalid_token`: Invalid authentication token for testing error handling

### Mock User Fixtures

- `mock_current_user`: Function to set the current user for testing
- `mock_admin_user`: Sets the current user as an admin
- `mock_regular_user`: Sets the current user as a regular user
- `mock_scientist_user`: Sets the current user as a scientist
- `mock_unauthenticated_user`: Sets no current user (unauthenticated)

## Usage Examples

### Testing Public Endpoints

```python
def test_public_endpoint(api_client):
    """Test accessing a public endpoint."""
    response = api_client.get('/api/public')
    assert response.status_code == 200
    data = api_client.parse_json(response)
    assert 'data' in data
```

### Testing Protected Endpoints

```python
def test_protected_endpoint(authenticated_client):
    """Test accessing a protected endpoint with authentication."""
    response = authenticated_client.get('/api/protected')
    assert response.status_code == 200
    data = authenticated_client.parse_json(response)
    assert 'data' in data
```

### Testing Admin-Only Endpoints

```python
def test_admin_endpoint(admin_client):
    """Test accessing an admin-only endpoint."""
    response = admin_client.get('/api/admin')
    assert response.status_code == 200
    data = admin_client.parse_json(response)
    assert 'data' in data
```

### Testing Scientist-Only Endpoints

```python
def test_scientist_endpoint(scientist_client):
    """Test accessing a scientist-only endpoint."""
    response = scientist_client.get('/api/scientist')
    assert response.status_code == 200
    data = scientist_client.parse_json(response)
    assert 'data' in data
```

### Testing JSON Request/Response

```python
def test_json_request(authenticated_client):
    """Test sending and receiving JSON data."""
    data = {
        'name': 'Test User',
        'email': 'test@example.com'
    }
    response = authenticated_client.post('/api/data', json_data=data)
    assert response.status_code == 200
    result = authenticated_client.parse_json(response)
    assert result['received'] == data
```

### Using Mock API Client

```python
def test_with_mock_client(mock_api_client):
    """Test using the mock API client."""
    # Configure mock responses
    mock_api_client.responses = {
        '/api/test': {'message': 'Test response'},
        '/api/data': lambda data, **kwargs: {'echo': data}
    }
    
    # Make request
    response = mock_api_client.get('/api/test')
    assert response.status_code == 200
    data = mock_api_client.parse_json(response)
    assert data['message'] == 'Test response'
    
    # Test with dynamic response
    test_data = {'test': 'data'}
    response = mock_api_client.post('/api/data', json_data=test_data)
    assert response.status_code == 200
    data = mock_api_client.parse_json(response)
    assert data['echo'] == test_data
    
    # Verify requests were recorded
    assert len(mock_api_client.requests) == 2
    assert mock_api_client.requests[0]['method'] == 'GET'
    assert mock_api_client.requests[1]['method'] == 'POST'
```

### Using Mock Authentication

```python
def test_with_mock_auth(mock_admin_user):
    """Test with mock authentication."""
    # The mock_admin_user fixture sets the current user as an admin
    # Now you can test functions that check the current user's role
    
    from api.utils import get_user_id
    user_id = get_user_id()
    assert user_id is not None
```

## Integration with Existing Tests

The API fixtures are designed to work with the existing testing framework, including the Database Fixtures and Mock Objects components.

### Combined Example

```python
def test_api_with_database(authenticated_client, mock_db):
    """Test API endpoint that interacts with the database."""
    # Add test data to the mock database
    mock_db.add_test_data('molecules', [
        {'id': 'mol-1', 'name': 'Test Molecule'}
    ])
    
    # Make API request
    response = authenticated_client.get('/api/molecules')
    assert response.status_code == 200
    data = authenticated_client.parse_json(response)
    assert len(data) == 1
    assert data[0]['name'] == 'Test Molecule'
```

## Extending the Fixtures

### Adding New Client Methods

To add new methods to the API client, extend the `APIClient` class in `client.py`:

```python
class ExtendedAPIClient(APIClient):
    def patch(self, endpoint, data=None, json_data=None, **kwargs):
        """Make a PATCH request to the API."""
        headers = kwargs.pop('headers', {})
        headers.update(self.headers)
        
        if json_data is not None:
            return self.client.patch(
                endpoint,
                json=json_data,
                headers=headers,
                **kwargs
            )
        else:
            return self.client.patch(
                endpoint,
                data=data,
                headers=headers,
                **kwargs
            )
```

### Adding New Authentication Roles

To add new authentication roles, extend the auth.py file:

```python
# Add new role constant
ROLE_MANAGER = "manager"

# Add new user ID
MANAGER_USER_ID = str(uuid.uuid4())

@pytest.fixture
def manager_token():
    """Generate a manager authentication token."""
    return generate_token(MANAGER_USER_ID, ROLE_MANAGER)

@pytest.fixture
def manager_client(client, manager_token):
    """Create an authenticated API client with a manager token."""
    from tests.fixtures.api.client import APIClient
    return APIClient(client, manager_token)
```

## Best Practices

1. **Use the appropriate client** for your test scenario:
   - `api_client` for public endpoints
   - `authenticated_client` for protected endpoints
   - `admin_client` for admin-only endpoints
   - `scientist_client` for scientist-only endpoints
   - `mock_api_client` for unit tests without a Flask app

2. **Combine with database fixtures** when testing endpoints that interact with the database.

3. **Use mock authentication** when testing functions that check the current user's role.

4. **Verify both success and error cases** to ensure proper error handling.

5. **Test with different user roles** to verify access control.

6. **Use the parse_json helper** to simplify working with JSON responses.
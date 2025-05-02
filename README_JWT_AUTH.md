# CryoProtect v2 - JWT Authentication System

This document provides information on the JWT-based authentication system implemented in CryoProtect v2, replacing the previous service role approach.

## Overview

The JWT authentication system provides secure, user-specific authentication with the following features:

- JWT token generation and validation
- Role-based access control (RBAC)
- Session management (timeout, refresh, revocation)
- Secure token storage using HTTP-only cookies

## Authentication Flow

1. **Login**: User submits credentials to the `/auth/login` endpoint
2. **Token Issuance**: Upon successful authentication, access and refresh tokens are issued
3. **Token Storage**: Tokens are stored securely (HTTP-only cookies or client-side storage)
4. **Authenticated Requests**: All API requests include the access token
5. **Token Validation**: Backend validates the token on every request
6. **Token Refresh**: When the access token expires, the refresh token is used to obtain a new one
7. **Logout**: Tokens are invalidated and removed from storage

## API Endpoints

### Authentication

- **POST /auth/login**: Authenticate with email and password
  ```json
  {
    "email": "user@example.com",
    "password": "your_password"
  }
  ```

- **POST /auth/logout**: Invalidate the current session (requires authentication)

- **POST /auth/refresh**: Refresh an expired access token
  ```json
  {
    "refresh_token": "your_refresh_token"
  }
  ```
  Note: If using HTTP-only cookies, the refresh token will be automatically included

- **GET /auth/validate**: Validate the current token and get user information (requires authentication)

### Multi-Factor Authentication (MFA)

- **POST /auth/mfa/initiate**: Initiate MFA challenge
  ```json
  {
    "email": "user@example.com"
  }
  ```

- **POST /auth/mfa/verify**: Verify MFA code
  ```json
  {
    "challenge_id": "challenge_id_from_initiate",
    "code": "verification_code"
  }
  ```

### User Management

- **POST /auth/register**: Register a new user
  ```json
  {
    "email": "user@example.com",
    "password": "your_password"
  }
  ```

- **POST /auth/reset-password**: Request password reset
  ```json
  {
    "email": "user@example.com"
  }
  ```

- **POST /auth/update-password**: Update password (requires authentication)
  ```json
  {
    "password": "new_password"
  }
  ```

- **POST /auth/update-profile**: Update user profile (requires authentication)
  ```json
  {
    "user_data": {
      "display_name": "John Doe",
      "avatar_url": "https://example.com/avatar.jpg"
    }
  }
  ```

## Using JWT Authentication in API Requests

### Including the Token in Requests

For authenticated requests, include the JWT token in the Authorization header:

```
Authorization: Bearer your_jwt_token
```

If using HTTP-only cookies (recommended), the token will be automatically included in requests to the same domain.

### Example API Request with Authentication

```python
import requests

# Login to get token
login_response = requests.post(
    'https://your-api.com/auth/login',
    json={
        'email': 'user@example.com',
        'password': 'your_password'
    }
)

# Extract token from response
token = login_response.json().get('token', {}).get('access_token')

# Make authenticated request
headers = {
    'Authorization': f'Bearer {token}'
}

response = requests.get('https://your-api.com/api/protected-resource', headers=headers)
```

## Role-Based Access Control (RBAC)

The system supports role-based access control with the following default roles:

- **user**: Regular user with access to their own data
- **admin**: Administrator with access to all data and administrative functions
- **curator**: User with special privileges for curating data
- **viewer**: User with read-only access to data

### Using Role-Based Decorators

In the backend code, you can use the `role_required` decorator to restrict access to specific roles:

```python
from api.jwt_auth import role_required

@app.route('/admin/dashboard')
@role_required('admin')
def admin_dashboard():
    # Only accessible to users with the 'admin' role
    return jsonify({'message': 'Admin dashboard'})

@app.route('/curator/edit')
@role_required(['curator', 'admin'])
def curator_edit():
    # Accessible to users with either 'curator' or 'admin' role
    return jsonify({'message': 'Curator edit page'})
```

## Configuration

Authentication settings can be configured in `auth_config.py`:

```python
# JWT Configuration
JWT_EXPIRY = 3600  # 1 hour in seconds
JWT_REFRESH_EXPIRY = 2592000  # 30 days in seconds

# Role Configuration
DEFAULT_ROLE = "user"
AVAILABLE_ROLES = ["user", "admin", "curator", "viewer"]

# Session Configuration
SESSION_TIMEOUT = 3600  # 1 hour in seconds
REFRESH_TOKEN_ROTATION = True  # Whether to rotate refresh tokens on use

# Security Configuration
SECURE_COOKIES = True  # Use secure cookies (HTTPS only)
HTTP_ONLY_COOKIES = True  # Use HTTP-only cookies for tokens
SAME_SITE_COOKIES = "Lax"  # SameSite cookie policy (Strict, Lax, None)
```

## Security Considerations

- **Token Storage**: Store tokens securely, preferably in HTTP-only cookies
- **HTTPS**: Always use HTTPS in production to protect tokens in transit
- **Token Expiry**: Keep access token expiry short (e.g., 1 hour) and refresh token expiry reasonable (e.g., 30 days)
- **Token Revocation**: Implement token revocation on logout and password change
- **CSRF Protection**: Use CSRF tokens for cookie-based authentication
- **XSS Protection**: Implement Content Security Policy (CSP) to prevent XSS attacks

## Migrating from Service Role Authentication

The service role approach has been completely removed from the codebase. All authentication now goes through the JWT-based system. If you were previously using the service role approach, you will need to:

1. Update your client code to use the new authentication endpoints
2. Include the JWT token in all authenticated requests
3. Update any code that was relying on the service role to use the new RBAC system

## Troubleshooting

### Common Issues

- **Token Expired**: The access token has expired. Use the refresh token to get a new one.
- **Invalid Token**: The token is malformed or has been tampered with.
- **Insufficient Permissions**: The user does not have the required role for the requested resource.
- **Token Not Found**: The token was not included in the request or could not be found in cookies.

### Debugging

For debugging authentication issues, check the application logs. Authentication-related errors are logged with the prefix "Authentication error" or "Token validation error".

## References

- [JWT.io](https://jwt.io/) - For more information on JWT tokens
- [OWASP Authentication Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Authentication_Cheat_Sheet.html) - For best practices on authentication
- [Supabase Auth Documentation](https://supabase.com/docs/guides/auth) - For Supabase-specific authentication information
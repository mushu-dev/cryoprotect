# CryoProtect v2 - Secure Session Management

This document provides information on the secure session management implementation for JWT authentication in CryoProtect v2.

## Overview

The session management system enhances the JWT authentication with the following features:

- **Token Revocation**: Invalidates tokens on logout and password change
- **Refresh Token Rotation**: Enhances security by rotating refresh tokens on use
- **Secure Session Storage**: Stores session information in the database
- **Audit Logging**: Tracks session activities for security monitoring
- **Session Verification**: Ensures tokens are associated with active sessions

## Architecture

### Components

1. **Session Manager**: Core component that handles session operations
   - Located in `api/session_management.py`
   - Provides methods for creating, validating, and revoking sessions

2. **Session Utilities**: Helper functions for session management
   - Located in `api/session_utils.py`
   - Includes background cleanup tasks and verification functions

3. **Enhanced JWT Authentication**: Extends JWT authentication with session verification
   - Located in `api/enhanced_jwt_auth.py`
   - Provides decorators that combine JWT validation with session verification

4. **Session API Routes**: Endpoints for managing sessions
   - Located in `api/session_routes.py`
   - Allows users to view and manage their active sessions

### Database Schema

The session management system uses two tables:

1. **sessions**: Stores session information
   ```sql
   CREATE TABLE sessions (
       id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
       user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
       refresh_token TEXT NOT NULL,
       refresh_token_hash TEXT NOT NULL,
       previous_refresh_token_hash TEXT,
       status TEXT NOT NULL,
       expires_at TIMESTAMP WITH TIME ZONE NOT NULL,
       created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
       updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
       ip_address TEXT,
       user_agent TEXT,
       device_info JSONB
   );
   ```

2. **session_audit_logs**: Tracks session activities
   ```sql
   CREATE TABLE session_audit_logs (
       id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
       user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
       session_id UUID REFERENCES sessions(id) ON DELETE SET NULL,
       action TEXT NOT NULL,
       status TEXT NOT NULL,
       ip_address TEXT,
       user_agent TEXT,
       details JSONB,
       created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
   );
   ```

## Session Lifecycle

### 1. Session Creation

When a user logs in:
1. The user is authenticated with Supabase Auth
2. Supabase Auth issues access and refresh tokens
3. The `SessionManager` creates a new session record in the database
4. The session is associated with the user's ID and refresh token
5. The activity is logged in the audit logs

### 2. Token Refresh

When an access token expires:
1. The client uses the refresh token to obtain a new access token
2. The `SessionManager` validates the refresh token against active sessions
3. If valid, the old refresh token is marked as rotated
4. A new session is created with the new refresh token
5. The activity is logged in the audit logs

### 3. Session Validation

On authenticated requests:
1. The JWT token is validated using standard JWT validation
2. The session is verified to ensure it's active and not expired
3. If the session is invalid, the request is rejected

### 4. Session Revocation

When a user logs out or changes password:
1. The `SessionManager` marks the session as revoked in the database
2. For password changes, all user sessions are revoked
3. The activity is logged in the audit logs

### 5. Session Cleanup

A background thread periodically:
1. Identifies expired sessions that are still marked as active
2. Updates their status to expired
3. Logs the cleanup activity

## API Endpoints

### Session Management

- **GET /api/sessions**: Get all active sessions for the current user
- **DELETE /api/sessions/{session_id}**: Revoke a specific session
- **DELETE /api/sessions/all**: Revoke all sessions except the current one
- **GET /api/sessions/audit**: Get audit logs for the current user's sessions

### Authentication

- **POST /auth/login**: Authenticate and create a new session
- **POST /auth/logout**: Revoke the current session
- **POST /auth/refresh**: Refresh tokens and rotate the session
- **POST /auth/update-password**: Update password and revoke all sessions

## Security Considerations

- **Token Storage**: Tokens are stored securely using HTTP-only cookies
- **Token Hashing**: Refresh tokens are hashed before storage in the database
- **Session Timeout**: Sessions have a configurable expiry time
- **Revocation**: Sessions can be revoked immediately when needed
- **Audit Trail**: All session activities are logged for security monitoring
- **Row Level Security**: Database tables use RLS to restrict access to user's own data

## Configuration

Session management settings can be configured in `auth_config.py`:

```python
# JWT Configuration
JWT_EXPIRY = 3600  # 1 hour in seconds
JWT_REFRESH_EXPIRY = 2592000  # 30 days in seconds

# Session Configuration
SESSION_TIMEOUT = 3600  # 1 hour in seconds
REFRESH_TOKEN_ROTATION = True  # Whether to rotate refresh tokens on use

# Security Configuration
SECURE_COOKIES = True  # Use secure cookies (HTTPS only)
HTTP_ONLY_COOKIES = True  # Use HTTP-only cookies for tokens
SAME_SITE_COOKIES = "Lax"  # SameSite cookie policy (Strict, Lax, None)
```

## Usage Examples

### Protecting Routes with Session Verification

```python
from api.enhanced_jwt_auth import session_verified_jwt_required

@app.route('/api/protected-resource')
@session_verified_jwt_required
def protected_resource():
    # This route requires both a valid JWT and an active session
    return jsonify({'message': 'Access granted'})
```

### Revoking a User's Sessions

```python
from api.session_management import get_session_manager

# Revoke all sessions for a user
session_manager = get_session_manager()
session_manager.revoke_all_user_sessions(user_id, reason='security_concern')
```

### Getting Active Sessions

```python
from api.session_utils import get_user_sessions

# Get all active sessions for a user
sessions = get_user_sessions(user_id)
```

## Troubleshooting

### Common Issues

1. **"Invalid session" error**: The session may have been revoked or expired. The user should log in again.
2. **"Failed to refresh session" error**: The refresh token may be invalid or the session may have been revoked. The user should log in again.
3. **"Session not found" error**: The session may have been deleted or never existed. The user should log in again.

### Debugging

For debugging session issues, check the application logs. Session-related errors are logged with the prefix "Session error" or "Session validation error".

## References

- [JWT.io](https://jwt.io/) - For more information on JWT tokens
- [OWASP Session Management Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Session_Management_Cheat_Sheet.html) - For best practices on session management
- [Supabase Auth Documentation](https://supabase.com/docs/guides/auth) - For Supabase-specific authentication information
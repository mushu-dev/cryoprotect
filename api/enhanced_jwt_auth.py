"""
CryoProtect Analyzer - Enhanced JWT Authentication

This module extends the JWT authentication with session verification,
ensuring that tokens are not only valid but also associated with active sessions.

Authentication Flow:
------------------
1. User Authentication:
   - User submits credentials to the /auth/login endpoint
   - Upon successful authentication, access and refresh tokens are issued
   - Access token contains user identity and permissions (short-lived, typically 1 hour)
   - Refresh token enables obtaining new access tokens (long-lived, typically 30 days)

2. Token Storage:
   - Access token is stored in memory or localStorage for API requests
   - Refresh token is stored in HTTP-only cookies for security
   - Both tokens are associated with a session record in the database

3. Request Authentication:
   - Client includes access token in Authorization header: "Bearer <token>"
   - Server validates token signature, expiration, and claims
   - Server verifies the token is associated with an active session
   - If valid, the request proceeds; otherwise, a 401 Unauthorized response is returned

4. Token Refresh:
   - When access token expires, client uses refresh token to obtain a new one
   - Old refresh token is invalidated and a new one is issued (token rotation)
   - This creates a new session record and marks the old one as rotated
   - If refresh token is invalid or expired, user must re-authenticate

5. Session Management:
   - Sessions track user login state across multiple devices
   - Sessions can be revoked individually or all at once (e.g., on password change)
   - Session status can be: active, revoked, expired, or rotated
   - Session verification adds an extra layer of security beyond JWT validation

6. Security Features:
   - Token signature verification using RSA public/private key pairs
   - Token expiration to limit the window of opportunity for token misuse
   - Session verification to detect revoked tokens
   - Refresh token rotation to prevent token reuse
   - Audit logging of all session-related activities

For detailed API usage, see README_JWT_AUTH.md and README_SESSION_MANAGEMENT.md
"""

import logging
from functools import wraps
from flask import request, g, jsonify

from .jwt_auth import jwt_required, get_token_from_request, extract_user_from_token
from .session_management import get_session_manager
from .session_utils import verify_session_status

# Setup logging
logger = logging.getLogger(__name__)

def session_verified_jwt_required(f):
    """
    Decorator to require a valid JWT token and verify the associated session.
    This extends the jwt_required decorator with session verification.
    
    This decorator provides an enhanced security layer by:
    1. Validating the JWT token signature, expiration, and claims
    2. Verifying that the token is associated with an active session
    3. Ensuring the session has not been revoked or expired
    
    The session verification step prevents the use of valid but revoked tokens,
    which is particularly important after password changes, account lockouts,
    or when a user explicitly logs out from all devices.
    
    Usage:
        @app.route('/api/protected')
        @session_verified_jwt_required
        def protected_endpoint():
            # This function will only execute if:
            # 1. A valid JWT token is provided
            # 2. The token is associated with an active session
            return {"message": "This is protected data"}
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function that returns a 401 Unauthorized response
        if the token is invalid or the session is not active
    """
    @wraps(f)
    @jwt_required
    def decorated(*args, **kwargs):
        try:
            # JWT validation is already done by jwt_required
            # Now verify the session status
            is_valid, error_message = verify_session_status()
            
            if not is_valid:
                return jsonify({'message': error_message or 'Invalid session'}), 401
            
            return f(*args, **kwargs)
        except Exception as e:
            logger.error(f"Session verification error: {str(e)}")
            return jsonify({'message': f'Session verification error: {str(e)}'}), 401
    
    return decorated

def validate_token_and_session(token):
    """
    Validate a JWT token and its associated session.
    
    This function performs a two-step validation process:
    1. Validates the JWT token's cryptographic signature, expiration time, and claims
    2. Verifies that the token is associated with an active session in the database
    
    The session validation step adds an important security layer beyond standard
    JWT validation, allowing for immediate token revocation through session management.
    This addresses one of the key limitations of stateless JWT authentication.
    
    If a refresh token is present in the request cookies, the function will also
    validate that the session associated with that refresh token is still active.
    This ensures that tokens cannot be used after a user has logged out or had
    their session revoked for security reasons.
    
    Args:
        token: The JWT token to validate (access token)
        
    Returns:
        Tuple of (is_valid, user_data, error_message) where:
        - is_valid: Boolean indicating if the token and session are valid
        - user_data: Dictionary containing user information extracted from the token
        - error_message: String containing error message if validation failed, None otherwise
    
    Security Notes:
    - This function is critical for maintaining secure authentication
    - It enables immediate session revocation despite JWT's stateless nature
    - It supports the security principle of defense in depth by adding session validation
    """
    try:
        # Validate the JWT token
        user_data, user_id = extract_user_from_token(token)
        
        # Get refresh token from cookies or request
        refresh_token = request.cookies.get('refresh_token')
        
        if not refresh_token:
            # No refresh token, can't verify session
            # This is still valid for API access with just an access token
            return True, user_data, None
        
        # Validate the session
        session_manager = get_session_manager()
        is_valid, session = session_manager.validate_session(refresh_token)
        
        if not is_valid:
            return False, None, "Session is invalid or expired"
        
        return True, user_data, None
    except Exception as e:
        logger.error(f"Token and session validation error: {str(e)}")
        return False, None, str(e)

def get_current_user_with_session_verification():
    """
    Get the current authenticated user with session verification.
    
    This function retrieves the current user's information from the JWT token
    while also verifying that the associated session is valid and active.
    It combines token extraction, validation, and session verification in one step.
    
    The function follows this process:
    1. Extract the JWT token from the request (Authorization header or cookies)
    2. If no token is found, return None (unauthenticated)
    3. Validate the token and verify the associated session
    4. If validation fails or session is invalid, return None
    5. If successful, return the user data extracted from the token
    
    This function is useful for routes that need user information but don't
    necessarily require the @session_verified_jwt_required decorator, such as
    routes that have optional authentication or need to display user-specific
    content without restricting access.
    
    Returns:
        User data dictionary or None if not authenticated or session is invalid
        
    Example Usage:
        @app.route('/api/profile')
        def get_profile():
            user = get_current_user_with_session_verification()
            if user:
                return jsonify({"profile": user})
            return jsonify({"message": "Not authenticated"}), 401
    """
    token = get_token_from_request()
    if not token:
        return None
    
    is_valid, user_data, _ = validate_token_and_session(token)
    if not is_valid:
        return None
    
    return user_data
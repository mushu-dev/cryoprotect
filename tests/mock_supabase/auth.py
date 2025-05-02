"""
CryoProtect Analyzer - Mock Supabase Auth

This module provides mock implementations of Supabase auth functionality.
"""

import uuid
from datetime import datetime, timedelta
from .data import _mock_data, DEFAULT_USER_ID
from .query import MockResponse

class MockAuth:
    """Mock implementation of Supabase auth."""
    
    def __init__(self):
        self.current_user = None
        self.current_session = None
    
    def sign_up(self, credentials):
        """
        Sign up a new user.
        
        Args:
            credentials: Dict with email and password
        
        Returns:
            MockResponse with user data
        """
        email = credentials.get('email')
        password = credentials.get('password')
        
        if not email or not password:
            return MockResponse(None, {"message": "Email and password are required"})
        
        # Check if user already exists
        existing_user = next((u for u in _mock_data['users'] if u.get('email') == email), None)
        if existing_user:
            return MockResponse(None, {"message": "User already exists"})
        
        # Create new user
        user_id = str(uuid.uuid4())
        now = datetime.now().isoformat()
        
        user = {
            'id': user_id,
            'email': email,
            'password': password,  # In a real system, this would be hashed
            'created_at': now,
            'updated_at': now
        }
        
        _mock_data['users'].append(user)
        
        # Create user profile
        profile = {
            'id': str(uuid.uuid4()),
            'user_id': user_id,
            'email': email,
            'name': None,
            'created_at': now,
            'updated_at': now
        }
        
        _mock_data['user_profiles'].append(profile)
        
        # Create session
        session = self._create_session(user)
        
        # Set current user and session
        self.current_user = user
        self.current_session = session
        
        # Return response with user
        response = MockResponse([user])
        response.user = user
        response.session = session
        
        return response
    
    def sign_in_with_password(self, credentials):
        """
        Sign in with email and password.
        
        Args:
            credentials: Dict with email and password
        
        Returns:
            MockResponse with user data
        """
        email = credentials.get('email')
        password = credentials.get('password')
        
        if not email or not password:
            return MockResponse(None, {"message": "Email and password are required"})
        
        # Find user
        user = next((u for u in _mock_data['users'] if u.get('email') == email), None)
        if not user or user.get('password') != password:
            return MockResponse(None, {"message": "Invalid credentials"})
        
        # Create session
        session = self._create_session(user)
        
        # Set current user and session
        self.current_user = user
        self.current_session = session
        
        # Return response with user
        response = MockResponse([user])
        response.user = user
        response.session = session
        
        return response
    
    def sign_out(self):
        """
        Sign out the current user.
        
        Returns:
            MockResponse
        """
        self.current_user = None
        self.current_session = None
        
        return MockResponse(None)
    
    def get_user(self, token):
        """
        Get user by token.
        
        Args:
            token: JWT token
        
        Returns:
            MockResponse with user data
        """
        # In a real implementation, this would validate the token
        # For testing, we'll just return the current user
        
        if not self.current_user:
            return MockResponse(None, {"message": "Invalid token"})
        
        response = MockResponse([self.current_user])
        response.user = self.current_user
        
        return response
    
    def update_user(self, updates):
        """
        Update user data.
        
        Args:
            updates: Dict with user data to update
        
        Returns:
            MockResponse with updated user data
        """
        if not self.current_user:
            return MockResponse(None, {"message": "No user is signed in"})
        
        # Update user
        for key, value in updates.items():
            if key == 'password':
                self.current_user['password'] = value
            elif key == 'email':
                self.current_user['email'] = value
            elif key == 'data':
                # User metadata
                if 'user_metadata' not in self.current_user:
                    self.current_user['user_metadata'] = {}
                
                for meta_key, meta_value in value.items():
                    self.current_user['user_metadata'][meta_key] = meta_value
        
        # Update timestamp
        self.current_user['updated_at'] = datetime.now().isoformat()
        
        # Update the user in the data store
        for i, user in enumerate(_mock_data['users']):
            if user['id'] == self.current_user['id']:
                _mock_data['users'][i] = self.current_user
                break
        
        # Return response with updated user
        response = MockResponse([self.current_user])
        response.user = self.current_user
        
        return response
    
    def reset_password_for_email(self, email, options=None):
        """
        Send password reset email.
        
        Args:
            email: User email
            options: Additional options
        
        Returns:
            MockResponse
        """
        # Find user
        user = next((u for u in _mock_data['users'] if u.get('email') == email), None)
        if not user:
            return MockResponse(None, {"message": "User not found"})
        
        # In a real implementation, this would send an email
        # For testing, we'll just return success
        
        return MockResponse(None)
    
    def refresh_session(self, refresh_token):
        """
        Refresh session with refresh token.
        
        Args:
            refresh_token: Refresh token
        
        Returns:
            MockResponse with new session
        """
        if not self.current_user:
            return MockResponse(None, {"message": "No user is signed in"})
        
        # Create new session
        session = self._create_session(self.current_user)
        
        # Update current session
        self.current_session = session
        
        # Return response with new session
        response = MockResponse(None)
        response.session = session
        
        return response
    
    def _create_session(self, user):
        """
        Create a new session for a user.
        
        Args:
            user: User data
        
        Returns:
            Session data
        """
        now = datetime.now()
        expires_at = now + timedelta(hours=1)
        
        session = {
            'access_token': f"mock_access_token_{user['id']}",
            'refresh_token': f"mock_refresh_token_{user['id']}",
            'expires_at': expires_at.isoformat(),
            'user': user
        }
        
        return session
    
    # MFA methods (simplified for testing)
    def mfa(self):
        """Get MFA handler."""
        return MockMFA()


class MockMFA:
    """Mock implementation of Supabase MFA."""
    
    def enroll(self, options):
        """
        Enroll in MFA.
        
        Args:
            options: MFA options
        
        Returns:
            MockResponse with MFA data
        """
        return MockResponse({
            'id': str(uuid.uuid4()),
            'factor_type': options.get('factor_type', 'totp'),
            'status': 'enrolled',
            'totp_secret': 'MOCK_TOTP_SECRET'
        })
    
    def challenge(self, options):
        """
        Challenge MFA.
        
        Args:
            options: Challenge options
        
        Returns:
            MockResponse with challenge result
        """
        factor_id = options.get('factor_id')
        code = options.get('code')
        
        if not factor_id or not code:
            return MockResponse(None, {"message": "Factor ID and code are required"})
        
        # For testing, accept any 6-digit code
        if len(str(code)) != 6:
            return MockResponse(None, {"message": "Invalid code"})
        
        return MockResponse({
            'id': str(uuid.uuid4()),
            'factor_id': factor_id,
            'status': 'verified'
        })
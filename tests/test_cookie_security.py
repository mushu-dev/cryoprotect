"""
Tests for cookie security in the CryoProtect API.

These tests verify that cookies are properly secured with the Secure, HttpOnly,
and SameSite attributes as required by the security audit remediation plan.
"""

import unittest
import json
from flask import session
from datetime import datetime, timedelta

# Import the base test case
from tests.base_test_case import BaseTestCase

class TestCookieSecurity(BaseTestCase):
    """Test case for cookie security."""
    
    def test_access_token_cookie_security(self):
        """Test that access_token cookie has proper security attributes."""
        # Mock login to get cookies
        response = self.client.post('/auth/login', json={
            'email': 'test@example.com',
            'password': 'password123'
        })
        
        # Check if access_token cookie exists
        self.assertIn('access_token', response.cookies)
        
        # Check security attributes
        cookie = response.cookies['access_token']
        self.assertTrue(cookie.secure)
        self.assertTrue(cookie.httponly)
        self.assertIn(cookie.samesite, ['Lax', 'Strict'])

    def test_refresh_token_cookie_security(self):
        """Test that refresh_token cookie has proper security attributes."""
        # Mock login to get cookies
        response = self.client.post('/auth/login', json={
            'email': 'test@example.com',
            'password': 'password123'
        })
        
        # Check if refresh_token cookie exists
        self.assertIn('refresh_token', response.cookies)
        
        # Check security attributes
        cookie = response.cookies['refresh_token']
        self.assertTrue(cookie.secure)
        self.assertTrue(cookie.httponly)
        self.assertIn(cookie.samesite, ['Lax', 'Strict'])

    def test_csrf_token_cookie_security(self):
        """Test that csrf_token cookie has proper security attributes."""
        # Get CSRF token to set cookie
        response = self.client.get('/api/v1/csrf-token')
        
        # Check if csrf_token cookie exists
        self.assertIn('csrf_token', response.cookies)
        
        # Check security attributes
        cookie = response.cookies['csrf_token']
        self.assertTrue(cookie.secure)
        self.assertFalse(cookie.httponly)  # Must be accessible to JavaScript
        self.assertIn(cookie.samesite, ['Lax', 'Strict'])

    def test_session_cookie_security(self):
        """Test that session cookie has proper security attributes."""
        # Access a route that uses session
        response = self.client.get('/api/v1/csrf-token')
        
        # Check if session cookie exists
        self.assertIn('session', response.cookies)
        
        # Check security attributes
        cookie = response.cookies['session']
        self.assertTrue(cookie.secure)
        self.assertTrue(cookie.httponly)
        self.assertIn(cookie.samesite, ['Lax', 'Strict'])

    def test_session_rotation_on_login(self):
        """Test that session is rotated on login."""
        # Get initial session
        self.client.get('/api/v1/csrf-token')
        with self.client.session_transaction() as sess:
            initial_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
        
        # Login to rotate session
        self.client.post('/auth/login', json={
            'email': 'test@example.com',
            'password': 'password123'
        })
        
        # Check that session has been rotated
        with self.client.session_transaction() as sess:
            new_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
        
        self.assertNotEqual(initial_session_id, new_session_id)

    def test_session_rotation_on_privilege_change(self):
        """Test that session is rotated on privilege change."""
        # Login as regular user
        self.client.post('/auth/login', json={
            'email': 'user@example.com',
            'password': 'password123'
        })
        
        with self.client.session_transaction() as sess:
            initial_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
            # Set user role to regular
            sess['user_role'] = 'user'
        
        # Simulate privilege escalation
        with self.client.session_transaction() as sess:
            sess['user_role'] = 'admin'
        
        # Access a protected route that checks privileges
        self.client.get('/api/v1/admin/users')
        
        # Check that session has been rotated
        with self.client.session_transaction() as sess:
            new_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
        
        self.assertNotEqual(initial_session_id, new_session_id)

    def test_session_rotation_on_password_change(self):
        """Test that session is rotated on password change."""
        # Login
        self.client.post('/auth/login', json={
            'email': 'test@example.com',
            'password': 'password123'
        })
        
        with self.client.session_transaction() as sess:
            initial_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
        
        # Change password
        self.client.post('/auth/update-password', json={
            'current_password': 'password123',
            'new_password': 'newpassword456'
        })
        
        # Check that session has been rotated
        with self.client.session_transaction() as sess:
            new_session_id = sess.sid if hasattr(sess, 'sid') else id(sess)
        
        self.assertNotEqual(initial_session_id, new_session_id)

    def test_session_expiry(self):
        """Test that session expires after the configured timeout."""
        # Login
        self.client.post('/auth/login', json={
            'email': 'test@example.com',
            'password': 'password123'
        })
        
        # Check that session has an expiry
        with self.client.session_transaction() as sess:
            self.assertTrue(any(key in sess for key in ['expires', '_expires', 'expiry', 'created_at']))

if __name__ == '__main__':
    unittest.main()
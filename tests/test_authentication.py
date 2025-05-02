"""
CryoProtect Analyzer - Authentication Tests

This module contains tests for the authentication system.
It tests user registration, login, logout, and access control with different user roles.
"""

import os
import sys
import json
import unittest
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

class TestAuthentication(MockSupabaseBaseTestCase):
    """Test cases for authentication system."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Sample user data for testing
        self.test_user = {
            'email': 'test@example.com',
            'password': 'Password123!',
            'display_name': 'Test User',
            'affiliation': 'Test Organization'
        }
        
        self.admin_user = {
            'email': 'admin@example.com',
            'password': 'AdminPass123!',
            'display_name': 'Admin User',
            'affiliation': 'Admin Organization',
            'is_admin': True
        }

    @patch_supabase(load_data=True)
    def test_user_registration(self, mock_client):
        """Test user registration."""
        # Mock the Supabase auth.sign_up method
        mock_client.auth.sign_up.return_value = MagicMock(
            user=MagicMock(id='test-user-id'),
            session=MagicMock(access_token='test-token')
        )
        
        # Make request to register endpoint
        response = self.client.post(
            '/auth/register',
            json=self.test_user
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertIn('user_id', data)
        self.assertIn('token', data)
        
        # Verify that the user profile was created
        mock_client.table.assert_called_with('user_profiles')
        mock_client.table().insert.assert_called_once()
        insert_data = mock_client.table().insert.call_args[0][0]
        self.assertEqual(insert_data['display_name'], self.test_user['display_name'])
        self.assertEqual(insert_data['affiliation'], self.test_user['affiliation'])

    @patch_supabase(load_data=True)
    def test_user_login(self, mock_client):
        """Test user login."""
        # Mock the Supabase auth.sign_in_with_password method
        mock_client.auth.sign_in_with_password.return_value = MagicMock(
            user=MagicMock(id='test-user-id'),
            session=MagicMock(access_token='test-token')
        )
        
        # Make request to login endpoint
        response = self.client.post(
            '/auth/login',
            json={
                'email': self.test_user['email'],
                'password': self.test_user['password']
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('user_id', data)
        self.assertIn('token', data)

    @patch_supabase(load_data=True)
    def test_user_logout(self, mock_client):
        """Test user logout."""
        # Mock the Supabase auth.sign_out method
        mock_client.auth.sign_out.return_value = None
        
        # Set up auth headers
        headers = {
            'Authorization': 'Bearer test-token'
        }
        
        # Make request to logout endpoint
        response = self.client.post(
            '/auth/logout',
            headers=headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Successfully logged out')

    @patch_supabase(load_data=True)
    def test_password_reset_request(self, mock_client):
        """Test password reset request."""
        # Mock the Supabase auth.reset_password_email method
        mock_client.auth.reset_password_email.return_value = None
        
        # Make request to password reset request endpoint
        response = self.client.post(
            '/auth/reset-password',
            json={
                'email': self.test_user['email']
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Password reset email sent')

    @patch_supabase(load_data=True)
    def test_access_control_anonymous_user(self, mock_client):
        """Test access control for anonymous users."""
        # Test access to public endpoints
        public_endpoints = [
            '/api/v1/molecules',
            '/api/v1/mixtures'
        ]
        
        for endpoint in public_endpoints:
            response = self.client.get(endpoint)
            self.assertEqual(response.status_code, 200, f"Anonymous user should be able to access {endpoint}")
        
        # Test access to protected endpoints
        protected_endpoints = [
            '/api/v1/molecules/import',
            '/api/v1/mixtures/create',
            '/api/v1/profile'
        ]
        
        for endpoint in protected_endpoints:
            response = self.client.get(endpoint)
            self.assertIn(response.status_code, [401, 403], 
                          f"Anonymous user should not be able to access {endpoint}")

    @patch_supabase(load_data=True)
    def test_access_control_authenticated_user(self, mock_client):
        """Test access control for authenticated users."""
        # Mock the Supabase auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id='test-user-id'
        )
        
        # Set up auth headers
        headers = {
            'Authorization': 'Bearer test-token'
        }
        
        # Test access to public endpoints
        public_endpoints = [
            '/api/v1/molecules',
            '/api/v1/mixtures'
        ]
        
        for endpoint in public_endpoints:
            response = self.client.get(endpoint, headers=headers)
            self.assertEqual(response.status_code, 200, 
                             f"Authenticated user should be able to access {endpoint}")
        
        # Test access to protected endpoints that require authentication
        protected_endpoints = [
            '/api/v1/profile',
            '/api/v1/molecules/import'
        ]
        
        for endpoint in protected_endpoints:
            response = self.client.get(endpoint, headers=headers)
            self.assertEqual(response.status_code, 200, 
                             f"Authenticated user should be able to access {endpoint}")

    @patch_supabase(load_data=True)
    def test_access_control_admin_user(self, mock_client):
        """Test access control for admin users."""
        # Mock the Supabase auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id='admin-user-id'
        )
        
        # Mock the user profile query to return admin status
        mock_client.table().select().eq().single.return_value = MagicMock(
            data={'is_admin': True}
        )
        
        # Set up auth headers
        headers = {
            'Authorization': 'Bearer admin-token'
        }
        
        # Test access to admin endpoints
        admin_endpoints = [
            '/api/v1/admin/users',
            '/api/v1/admin/settings'
        ]
        
        for endpoint in admin_endpoints:
            response = self.client.get(endpoint, headers=headers)
            self.assertEqual(response.status_code, 200, 
                             f"Admin user should be able to access {endpoint}")

    @patch_supabase(load_data=True)
    def test_user_cannot_access_admin_endpoints(self, mock_client):
        """Test that regular users cannot access admin endpoints."""
        # Mock the Supabase auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id='test-user-id'
        )
        
        # Mock the user profile query to return non-admin status
        mock_client.table().select().eq().single.return_value = MagicMock(
            data={'is_admin': False}
        )
        
        # Set up auth headers
        headers = {
            'Authorization': 'Bearer test-token'
        }
        
        # Test access to admin endpoints
        admin_endpoints = [
            '/api/v1/admin/users',
            '/api/v1/admin/settings'
        ]
        
        for endpoint in admin_endpoints:
            response = self.client.get(endpoint, headers=headers)
            self.assertIn(response.status_code, [401, 403], 
                          f"Regular user should not be able to access {endpoint}")

    @patch_supabase(load_data=True)
    def test_user_can_only_modify_own_data(self, mock_client):
        """Test that users can only modify their own data."""
        # Mock the Supabase auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id='test-user-id'
        )
        
        # Set up auth headers
        headers = {
            'Authorization': 'Bearer test-token'
        }
        
        # Test updating own profile
        response = self.client.put(
            '/api/v1/profile',
            headers=headers,
            json={
                'display_name': 'Updated Name',
                'affiliation': 'Updated Organization'
            }
        )
        self.assertEqual(response.status_code, 200, "User should be able to update own profile")
        
        # Test updating another user's profile
        response = self.client.put(
            '/api/v1/users/another-user-id',
            headers=headers,
            json={
                'display_name': 'Hacked Name',
                'affiliation': 'Hacked Organization'
            }
        )
        self.assertIn(response.status_code, [401, 403, 404], 
                      "User should not be able to update another user's profile")

if __name__ == '__main__':
    unittest.main()
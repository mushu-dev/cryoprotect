"""
CryoProtect Analyzer - API Endpoints Tests with Different User Roles

This module contains tests for the API endpoints with different user roles.
It tests request validation, response formats, and error handling for anonymous users,
authenticated users, and admin users.
"""

import os
import sys
import json
import uuid
from unittest.mock import patch, MagicMock
from flask import Flask
from flask_restful import Api

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

# Import the app and API resources
from app import create_app
from api.models import (
    Molecule, Mixture, Prediction, Experiment, Comparison, UserProfile
)

class TestAPIEndpointsWithRoles(MockSupabaseBaseTestCase):
    """Test cases for API endpoints with different user roles."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Sample data for testing
        self.sample_molecule_id = str(uuid.uuid4())
        self.sample_mixture_id = str(uuid.uuid4())
        
        # Sample user profiles
        self.anonymous_user = None
        
        self.regular_user = {
            'id': str(uuid.uuid4()),
            'auth_user_id': str(uuid.uuid4()),
            'display_name': 'Regular User',
            'email': 'regular@example.com',
            'affiliation': 'Regular Organization',
            'is_admin': False
        }
        
        self.admin_user = {
            'id': str(uuid.uuid4()),
            'auth_user_id': str(uuid.uuid4()),
            'display_name': 'Admin User',
            'email': 'admin@example.com',
            'affiliation': 'Admin Organization',
            'is_admin': True
        }
        
        # Sample molecule
        self.sample_molecule = {
            'id': self.sample_molecule_id,
            'cid': 123456,
            'name': 'Glycerol',
            'molecular_formula': 'C3H8O3',
            'smiles': 'C(C(CO)O)O',
            'inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
            'inchikey': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
            'properties': [],
            'created_by': self.regular_user['id']
        }
        
        # Sample mixture
        self.sample_mixture = {
            'id': self.sample_mixture_id,
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'components': [
                {
                    'molecule_id': self.sample_molecule_id,
                    'concentration': 100,
                    'concentration_unit': '%'
                }
            ],
            'created_by': self.regular_user['id']
        }
        
        # Auth headers
        self.regular_user_headers = {
            'Authorization': f'Bearer {str(uuid.uuid4())}'
        }
        
        self.admin_user_headers = {
            'Authorization': f'Bearer {str(uuid.uuid4())}'
        }

    # Anonymous User Tests
    
    @patch_supabase(load_data=True)
    def test_anonymous_user_can_get_molecules(self, mock_client):
        """Test that anonymous users can get molecules."""
        # Set up the mock Supabase client
        mock_client.table('molecules').select('*').execute.return_value = MagicMock(
            data=[self.sample_molecule]
        )
        
        # Make request
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertGreater(len(data), 0)

    @patch_supabase(load_data=True)
    def test_anonymous_user_can_get_mixtures(self, mock_client):
        """Test that anonymous users can get mixtures."""
        # Set up the mock Supabase client
        mock_client.table('mixture_with_components').select('*').execute.return_value = MagicMock(
            data=[self.sample_mixture]
        )
        
        # Make request
        response = self.client.get('/api/v1/mixtures')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertGreater(len(data), 0)

    @patch_supabase(load_data=True)
    def test_anonymous_user_cannot_create_mixture(self, mock_client):
        """Test that anonymous users cannot create mixtures."""
        # Make request
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': self.sample_molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            }
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])

    # Regular User Tests
    
    @patch_supabase(load_data=True)
    def test_regular_user_can_create_mixture(self, mock_client):
        """Test that regular users can create mixtures."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.regular_user['auth_user_id']
        )
        
        # Mock the table insert method
        mock_client.table('mixtures').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}]
        )
        
        mock_client.table('mixture_components').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}]
        )
        
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[self.sample_mixture]
        )
        
        # Make request
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': self.sample_molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['name'], 'Test Mixture')

    @patch_supabase(load_data=True)
    def test_regular_user_can_update_own_mixture(self, mock_client):
        """Test that regular users can update their own mixtures."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.regular_user['auth_user_id']
        )
        
        # Mock the table select method to return a mixture owned by the user
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'created_by': self.regular_user['id']}]
        )
        
        # Mock the table update method
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'name': 'Updated Mixture'}]
        )
        
        # Mock the mixture_with_components view
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[{**self.sample_mixture, 'name': 'Updated Mixture'}]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.sample_mixture_id}',
            json={
                'name': 'Updated Mixture',
                'description': 'An updated test mixture',
                'components': [
                    {
                        'molecule_id': self.sample_molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['name'], 'Updated Mixture')

    @patch_supabase(load_data=True)
    def test_regular_user_cannot_update_others_mixture(self, mock_client):
        """Test that regular users cannot update mixtures they don't own."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.regular_user['auth_user_id']
        )
        
        # Mock the table select method to return a mixture owned by another user
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'created_by': self.admin_user['id']}]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.sample_mixture_id}',
            json={
                'name': 'Updated Mixture',
                'description': 'An updated test mixture',
                'components': [
                    {
                        'molecule_id': self.sample_molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])

    # Admin User Tests
    
    @patch_supabase(load_data=True)
    def test_admin_user_can_update_any_mixture(self, mock_client):
        """Test that admin users can update any mixture."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.admin_user['auth_user_id']
        )
        
        # Mock the user profile query to return admin status
        mock_client.table('user_profiles').select().eq().single.return_value = MagicMock(
            data={'is_admin': True}
        )
        
        # Mock the table select method to return a mixture owned by another user
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'created_by': self.regular_user['id']}]
        )
        
        # Mock the table update method
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'name': 'Admin Updated Mixture'}]
        )
        
        # Mock the mixture_with_components view
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[{**self.sample_mixture, 'name': 'Admin Updated Mixture'}]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.sample_mixture_id}',
            json={
                'name': 'Admin Updated Mixture',
                'description': 'An admin updated test mixture',
                'components': [
                    {
                        'molecule_id': self.sample_molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.admin_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['name'], 'Admin Updated Mixture')

    @patch_supabase(load_data=True)
    def test_admin_user_can_access_admin_endpoints(self, mock_client):
        """Test that admin users can access admin endpoints."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.admin_user['auth_user_id']
        )
        
        # Mock the user profile query to return admin status
        mock_client.table('user_profiles').select().eq().single.return_value = MagicMock(
            data={'is_admin': True}
        )
        
        # Mock the admin endpoint response
        mock_client.table('user_profiles').select().execute.return_value = MagicMock(
            data=[self.regular_user, self.admin_user]
        )
        
        # Make request
        response = self.client.get(
            '/api/v1/admin/users',
            headers=self.admin_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 2)

    # Error Handling Tests
    
    @patch_supabase(load_data=True)
    def test_not_found_error(self, mock_client):
        """Test 404 error handling."""
        # Mock the table select method to return no data
        mock_client.table('molecule_with_properties').select().eq().execute.return_value = MagicMock(
            data=[]
        )
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{str(uuid.uuid4())}')
        
        # Assertions
        self.assertEqual(response.status_code, 404)
        data = json.loads(response.data)
        self.assertIn('message', data)

    @patch_supabase(load_data=True)
    def test_validation_error(self, mock_client):
        """Test validation error handling."""
        # Mock the auth.get_user method
        mock_client.auth.get_user.return_value = MagicMock(
            id=self.regular_user['auth_user_id']
        )
        
        # Make request with invalid data
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': []  # Invalid: empty components
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('message', data)

    @patch_supabase(load_data=True)
    def test_server_error(self, mock_client):
        """Test 500 error handling."""
        # Mock the table select method to raise an exception
        mock_client.table('molecule_with_properties').select().eq().execute.side_effect = Exception("Database error")
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{self.sample_molecule_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 500)
        data = json.loads(response.data)
        self.assertIn('message', data)

    # Performance Tests
    
    @patch_supabase(load_data=True)
    def test_pagination(self, mock_client):
        """Test pagination for large result sets."""
        # Create a large dataset
        molecules = []
        for i in range(100):
            molecules.append({
                'id': str(uuid.uuid4()),
                'name': f'Molecule {i}',
                'properties': []
            })
        
        # Mock the table select method to return paginated results
        mock_client.table('molecule_with_properties').select().range().execute.return_value = MagicMock(
            data=molecules[:10]
        )
        
        # Make request with pagination parameters
        response = self.client.get('/api/v1/molecules?limit=10&offset=0')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 10)
        
        # Mock the next page
        mock_client.table('molecule_with_properties').select().range().execute.return_value = MagicMock(
            data=molecules[10:20]
        )
        
        # Make request for the next page
        response = self.client.get('/api/v1/molecules?limit=10&offset=10')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 10)
        self.assertNotEqual(data[0]['name'], molecules[0]['name'])

if __name__ == '__main__':
    unittest.main()
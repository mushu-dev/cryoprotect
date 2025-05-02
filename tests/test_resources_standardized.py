"""
CryoProtect Analyzer - Standardized API Resources Tests

This module contains tests for the standardized error handling, authentication,
and response formatting patterns in api/resources.py.
"""

import os
import sys
import json
import uuid
from unittest.mock import patch, MagicMock
from datetime import datetime

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

# Import the app and API resources
from app import create_app
from api.models import (
    Molecule, Mixture, Prediction, Experiment
)
from api.utils import handle_supabase_error, handle_error, token_required

class TestStandardizedResources(MockSupabaseBaseTestCase):
    """Test cases for standardized API resources."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Sample data for testing
        self.sample_molecule_id = str(uuid.uuid4())
        self.sample_mixture_id = str(uuid.uuid4())
        self.sample_user_id = str(uuid.uuid4())
        self.sample_auth_id = str(uuid.uuid4())
        
        # Create a test user profile
        self.sample_user = {
            'id': self.sample_user_id,
            'auth_user_id': self.sample_auth_id,
            'display_name': 'Test User',
            'email': 'test@example.com',
            'affiliation': 'Test Organization',
            'is_admin': False
        }
        
        self.sample_molecule = {
            'id': self.sample_molecule_id,
            'cid': 123456,
            'name': 'Glycerol',
            'molecular_formula': 'C3H8O3',
            'smiles': 'C(C(CO)O)O',
            'inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
            'inchikey': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
            'created_by': self.sample_user_id,
            'created_at': datetime.now().isoformat()
        }
        
        self.sample_mixture = {
            'id': self.sample_mixture_id,
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_by': self.sample_user_id,
            'created_at': datetime.now().isoformat(),
            'components': [
                {
                    'molecule_id': self.sample_molecule_id,
                    'concentration': 100,
                    'concentration_unit': '%'
                }
            ]
        }
        
        # Auth headers
        self.auth_token = str(uuid.uuid4())
        self.auth_headers = {
            'Authorization': f'Bearer {self.auth_token}'
        }

    # Helper methods
    def mock_auth_user(self, mock_client):
        """Mock the auth.get_user method to return a specific user."""
        # Mock user object
        mock_user = MagicMock()
        mock_user.id = self.sample_auth_id
        
        # Mock auth.get_user response
        mock_client.auth.get_user.return_value = MagicMock(
            user=mock_user,
            error=None
        )
        
        # Mock user profile query
        mock_client.table('user_profiles').select().eq().single.return_value = MagicMock(
            data=self.sample_user
        )

    # Standardized Error Handling Tests
    
    @patch_supabase(load_data=True)
    def test_handle_supabase_error_pattern(self, mock_client):
        """Test the standardized handle_supabase_error pattern."""
        # Mock a Supabase error response
        mock_error_response = MagicMock()
        mock_error_response.error = {
            'message': 'Database error',
            'code': 'PGRST500'
        }
        
        # Mock the Supabase client to return the error
        mock_client.table('molecules').select('*').execute.return_value = mock_error_response
        
        # Make request
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 500)
        data = json.loads(response.data)
        self.assertEqual(data['status'], 'error')
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Database error')
        self.assertIn('meta', data)
        self.assertEqual(data['meta']['error_type'], 'Error')
        self.assertEqual(data['meta']['context'], 'Fetching molecules list')

    @patch_supabase(load_data=True)
    def test_handle_error_pattern(self, mock_client):
        """Test the standardized handle_error pattern."""
        # Mock the Supabase client to raise an exception
        mock_client.table('molecules').select('*').execute.side_effect = ValueError("Invalid parameter")
        
        # Make request
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertEqual(data['status'], 'error')
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Invalid parameter')
        self.assertIn('meta', data)
        self.assertEqual(data['meta']['error_type'], 'ValueError')
        self.assertEqual(data['meta']['context'], 'Fetching molecules list')

    @patch_supabase(load_data=True)
    def test_not_found_error_handling(self, mock_client):
        """Test 404 error handling with standardized pattern."""
        # Mock the Supabase response
        mock_client.table('molecules_with_properties').select('*').eq().execute.return_value = MagicMock(
            data=[],
            error=None
        )
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{str(uuid.uuid4())}')
        
        # Assertions
        self.assertEqual(response.status_code, 404)
        data = json.loads(response.data)
        self.assertEqual(data['status'], 'error')
        self.assertIn('message', data)
        self.assertIn('meta', data)
        self.assertEqual(data['meta']['error_type'], 'Error')
        self.assertIn('Fetching molecule', data['meta']['context'])

    @patch_supabase(load_data=True)
    def test_validation_error_handling(self, mock_client):
        """Test validation error handling with standardized pattern."""
        # Mock the auth user
        self.mock_auth_user(mock_client)
        
        # Make request with invalid data
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                # Missing required 'components' field
            },
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertEqual(data['status'], 'error')
        self.assertIn('message', data)
        self.assertIn('meta', data)
        self.assertEqual(data['meta']['error_type'], 'ValidationError')
        self.assertEqual(data['meta']['context'], 'Validating mixture data')

    # Authentication Tests
    
    @patch_supabase(load_data=True)
    def test_get_endpoint_no_auth_required(self, mock_client):
        """Test that GET endpoints don't require authentication."""
        # Mock the Supabase response
        mock_client.table('molecules_with_properties').select('*').range().execute.return_value = MagicMock(
            data=[self.sample_molecule],
            error=None
        )
        
        # Make request without auth header
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]['id'], self.sample_molecule_id)

    @patch_supabase(load_data=True)
    def test_post_endpoint_auth_required(self, mock_client):
        """Test that POST endpoints require authentication."""
        # Make request without auth header
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
        self.assertEqual(response.status_code, 401)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Authentication token is missing')

    @patch_supabase(load_data=True)
    def test_put_endpoint_auth_required(self, mock_client):
        """Test that PUT endpoints require authentication."""
        # Make request without auth header
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
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 401)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Authentication token is missing')

    @patch_supabase(load_data=True)
    def test_delete_endpoint_auth_required(self, mock_client):
        """Test that DELETE endpoints require authentication."""
        # Make request without auth header
        response = self.client.delete(f'/api/v1/mixtures/{self.sample_mixture_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 401)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertEqual(data['message'], 'Authentication token is missing')

    # Response Formatting Tests
    
    @patch_supabase(load_data=True)
    def test_get_response_formatting(self, mock_client):
        """Test response formatting for GET endpoints."""
        # Mock the Supabase response
        mock_client.table('molecules_with_properties').select('*').eq().execute.return_value = MagicMock(
            data=[self.sample_molecule],
            error=None
        )
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{self.sample_molecule_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.sample_molecule_id)
        self.assertEqual(data['name'], 'Glycerol')
        self.assertEqual(data['molecular_formula'], 'C3H8O3')
        self.assertEqual(data['smiles'], 'C(C(CO)O)O')

    @patch_supabase(load_data=True)
    def test_post_response_formatting(self, mock_client):
        """Test response formatting for POST endpoints."""
        # Mock the auth user
        self.mock_auth_user(mock_client)
        
        # Mock the Supabase response
        mock_client.table('mixtures').insert().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id}],
            error=None
        )
        
        mock_client.table('mixture_components').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}],
            error=None
        )
        
        mock_client.table('mixtures_with_components').select('*').eq().execute.return_value = MagicMock(
            data=[self.sample_mixture],
            error=None
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
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.sample_mixture_id)
        self.assertEqual(data['name'], 'Test Mixture')
        self.assertEqual(data['description'], 'A test mixture')
        self.assertIn('components', data)
        self.assertEqual(len(data['components']), 1)
        self.assertEqual(data['components'][0]['molecule_id'], self.sample_molecule_id)

    @patch_supabase(load_data=True)
    def test_put_response_formatting(self, mock_client):
        """Test response formatting for PUT endpoints."""
        # Mock the auth user
        self.mock_auth_user(mock_client)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select('*').eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'created_by': self.sample_user_id}],
            error=None
        )
        
        updated_mixture = {
            **self.sample_mixture,
            'name': 'Updated Mixture',
            'description': 'An updated test mixture'
        }
        
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id}],
            error=None
        )
        
        mock_client.table('mixture_components').delete().eq().execute.return_value = MagicMock(
            data=[],
            error=None
        )
        
        mock_client.table('mixture_components').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}],
            error=None
        )
        
        mock_client.table('mixtures_with_components').select('*').eq().execute.return_value = MagicMock(
            data=[updated_mixture],
            error=None
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
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.sample_mixture_id)
        self.assertEqual(data['name'], 'Updated Mixture')
        self.assertEqual(data['description'], 'An updated test mixture')
        self.assertIn('components', data)
        self.assertEqual(len(data['components']), 1)
        self.assertEqual(data['components'][0]['molecule_id'], self.sample_molecule_id)

    @patch_supabase(load_data=True)
    def test_delete_response_formatting(self, mock_client):
        """Test response formatting for DELETE endpoints."""
        # Mock the auth user
        self.mock_auth_user(mock_client)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select('*').eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id, 'created_by': self.sample_user_id}],
            error=None
        )
        
        mock_client.table('mixtures').delete().eq().execute.return_value = MagicMock(
            data=[{'id': self.sample_mixture_id}],
            error=None
        )
        
        # Make request
        response = self.client.delete(
            f'/api/v1/mixtures/{self.sample_mixture_id}',
            headers=self.auth_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('message', data)
        self.assertIn(self.sample_mixture_id, data['message'])
        self.assertIn('deleted successfully', data['message'])


if __name__ == '__main__':
    import unittest
    unittest.main()
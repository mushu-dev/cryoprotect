"""
CryoProtect Analyzer - Supabase API Integration Tests

This module contains comprehensive tests for API endpoints that interact with Supabase,
focusing on authentication scenarios and RLS policy enforcement.
"""

import os
import sys
import json
import uuid
import unittest
from unittest.mock import patch, MagicMock
from datetime import datetime, date

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase, reset_mock_data, load_test_data

# Import the app and API resources
from app import create_app
from api.models import (
    Molecule, Mixture, Prediction, Experiment, Comparison, UserProfile
)

class TestSupabaseAPIIntegration(MockSupabaseBaseTestCase):
    """Test cases for API endpoints with Supabase integration."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Create test users
        self.anonymous_user = None
        
        self.regular_user_id = str(uuid.uuid4())
        self.regular_auth_id = str(uuid.uuid4())
        self.regular_user = {
            'id': self.regular_user_id,
            'auth_user_id': self.regular_auth_id,
            'display_name': 'Regular User',
            'email': 'regular@example.com',
            'affiliation': 'Regular Organization',
            'is_admin': False
        }
        
        self.other_user_id = str(uuid.uuid4())
        self.other_auth_id = str(uuid.uuid4())
        self.other_user = {
            'id': self.other_user_id,
            'auth_user_id': self.other_auth_id,
            'display_name': 'Other User',
            'email': 'other@example.com',
            'affiliation': 'Other Organization',
            'is_admin': False
        }
        
        self.admin_user_id = str(uuid.uuid4())
        self.admin_auth_id = str(uuid.uuid4())
        self.admin_user = {
            'id': self.admin_user_id,
            'auth_user_id': self.admin_auth_id,
            'display_name': 'Admin User',
            'email': 'admin@example.com',
            'affiliation': 'Admin Organization',
            'is_admin': True
        }

        # Create test data
        self.molecule_id = str(uuid.uuid4())
        self.molecule = {
            'id': self.molecule_id,
            'cid': 123456,
            'name': 'Glycerol',
            'molecular_formula': 'C3H8O3',
            'smiles': 'C(C(CO)O)O',
            'inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
            'inchikey': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
            'created_by': self.regular_user_id,
            'created_at': datetime.now().isoformat()
        }
        
        self.regular_user_mixture_id = str(uuid.uuid4())
        self.regular_user_mixture = {
            'id': self.regular_user_mixture_id,
            'name': 'Regular User Mixture',
            'description': 'A test mixture owned by regular user',
            'created_by': self.regular_user_id,
            'created_at': datetime.now().isoformat(),
            'components': [
                {
                    'molecule_id': self.molecule_id,
                    'concentration': 100,
                    'concentration_unit': '%'
                }
            ]
        }
        
        self.other_user_mixture_id = str(uuid.uuid4())
        self.other_user_mixture = {
            'id': self.other_user_mixture_id,
            'name': 'Other User Mixture',
            'description': 'A test mixture owned by other user',
            'created_by': self.other_user_id,
            'created_at': datetime.now().isoformat(),
            'components': [
                {
                    'molecule_id': self.molecule_id,
                    'concentration': 100,
                    'concentration_unit': '%'
                }
            ]
        }
        
        self.prediction_id = str(uuid.uuid4())
        self.prediction = {
            'id': self.prediction_id,
            'mixture_id': self.regular_user_mixture_id,
            'property_name': 'Freezing Point',
            'numeric_value': -15.3,
            'confidence': 0.9,
            'calculation_method': 'CryoProtect Scoring',
            'created_by': self.regular_user_id,
            'created_at': datetime.now().isoformat()
        }
        
        self.experiment_id = str(uuid.uuid4())
        self.experiment = {
            'id': self.experiment_id,
            'mixture_id': self.regular_user_mixture_id,
            'property_name': 'Freezing Point',
            'numeric_value': -14.8,
            'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
            'date_performed': date.today().isoformat(),
            'created_by': self.regular_user_id,
            'created_at': datetime.now().isoformat()
        }
        
        # Auth headers
        self.regular_user_token = str(uuid.uuid4())
        self.regular_user_headers = {
            'Authorization': f'Bearer {self.regular_user_token}'
        }
        
        self.other_user_token = str(uuid.uuid4())
        self.other_user_headers = {
            'Authorization': f'Bearer {self.other_user_token}'
        }
        
        self.admin_user_token = str(uuid.uuid4())
        self.admin_user_headers = {
            'Authorization': f'Bearer {self.admin_user_token}'
        }
        
        self.service_role_token = str(uuid.uuid4())
        self.service_role_headers = {
            'Authorization': f'Bearer {self.service_role_token}'
        }

    # Helper methods
    def mock_auth_user(self, mock_client, user_id, auth_id, is_admin=False):
        """Mock the auth.get_user method to return a specific user."""
        # Mock user object
        mock_user = MagicMock()
        mock_user.id = auth_id
        
        # Mock auth.get_user response
        mock_client.auth.get_user.return_value = MagicMock(
            user=mock_user,
            error=None
        )
        
        # Mock user profile query
        mock_client.table('user_profiles').select().eq().single.return_value = MagicMock(
            data={'id': user_id, 'auth_user_id': auth_id, 'is_admin': is_admin}
        )

    # Unauthenticated Access Tests
    
    @patch_supabase(load_data=True)
    def test_unauthenticated_get_molecules(self, mock_client):
        """Test that unauthenticated users can get molecules."""
        # Mock the Supabase response
        mock_client.table('molecules').select('*').execute.return_value = MagicMock(
            data=[self.molecule]
        )
        
        # Make request
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]['id'], self.molecule_id)
    
    @patch_supabase(load_data=True)
    def test_unauthenticated_get_molecule_detail(self, mock_client):
        """Test that unauthenticated users can get molecule details."""
        # Mock the Supabase response
        mock_client.table('molecule_with_properties').select('*').eq().execute.return_value = MagicMock(
            data=[self.molecule]
        )
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{self.molecule_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], self.molecule_id)
    
    @patch_supabase(load_data=True)
    def test_unauthenticated_get_mixtures(self, mock_client):
        """Test that unauthenticated users can get public mixtures."""
        # Mock the Supabase response
        mock_client.table('mixture_with_components').select('*').execute.return_value = MagicMock(
            data=[self.regular_user_mixture, self.other_user_mixture]
        )
        
        # Make request
        response = self.client.get('/api/v1/mixtures')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 2)
    
    @patch_supabase(load_data=True)
    def test_unauthenticated_cannot_create_mixture(self, mock_client):
        """Test that unauthenticated users cannot create mixtures."""
        # Make request
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            }
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])
    
    # Authenticated User Tests - Own Data
    
    @patch_supabase(load_data=True)
    def test_authenticated_user_create_mixture(self, mock_client):
        """Test that authenticated users can create mixtures."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}]
        )
        
        mock_client.table('mixture_components').insert().execute.return_value = MagicMock(
            data=[{'id': str(uuid.uuid4())}]
        )
        
        new_mixture = {
            'id': str(uuid.uuid4()),
            'name': 'New Mixture',
            'description': 'A new test mixture',
            'created_by': self.regular_user_id,
            'components': [
                {
                    'molecule_id': self.molecule_id,
                    'concentration': 100,
                    'concentration_unit': '%'
                }
            ]
        }
        
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[new_mixture]
        )
        
        # Make request
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'New Mixture',
                'description': 'A new test mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
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
        self.assertEqual(data['name'], 'New Mixture')

    @patch_supabase(load_data=True)
    def test_authenticated_user_update_own_mixture(self, mock_client):
        """Test that authenticated users can update their own mixtures."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.regular_user_mixture_id, 'created_by': self.regular_user_id}]
        )
        
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.regular_user_mixture_id, 'name': 'Updated Mixture'}]
        )
        
        updated_mixture = {
            **self.regular_user_mixture,
            'name': 'Updated Mixture'
        }
        
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[updated_mixture]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.regular_user_mixture_id}',
            json={
                'name': 'Updated Mixture',
                'description': 'An updated test mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
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
    def test_authenticated_user_add_prediction(self, mock_client):
        """Test that authenticated users can add predictions to their own mixtures."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.regular_user_mixture_id, 'created_by': self.regular_user_id}]
        )
        
        new_prediction_id = str(uuid.uuid4())
        new_prediction = {
            'id': new_prediction_id,
            'mixture_id': self.regular_user_mixture_id,
            'property_name': 'Viscosity',
            'numeric_value': 1.5,
            'confidence': 0.8,
            'calculation_method': 'CryoProtect Scoring',
            'created_by': self.regular_user_id
        }
        
        mock_client.table('predictions').insert().execute.return_value = MagicMock(
            data=[new_prediction]
        )
        
        # Make request
        response = self.client.post(
            f'/api/v1/mixtures/{self.regular_user_mixture_id}/predictions',
            json={
                'property_name': 'Viscosity',
                'value': 1.5,
                'confidence': 0.8,
                'calculation_method': 'CryoProtect Scoring'
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['property_name'], 'Viscosity')
        self.assertEqual(data['numeric_value'], 1.5)
    
    @patch_supabase(load_data=True)
    def test_authenticated_user_add_experiment(self, mock_client):
        """Test that authenticated users can add experiments to their own mixtures."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.regular_user_mixture_id, 'created_by': self.regular_user_id}]
        )
        
        new_experiment_id = str(uuid.uuid4())
        new_experiment = {
            'id': new_experiment_id,
            'mixture_id': self.regular_user_mixture_id,
            'property_name': 'Viscosity',
            'numeric_value': 1.7,
            'experimental_conditions': 'Room temperature',
            'date_performed': date.today().isoformat(),
            'created_by': self.regular_user_id
        }
        
        mock_client.table('experiments').insert().execute.return_value = MagicMock(
            data=[new_experiment]
        )
        
        # Make request
        response = self.client.post(
            f'/api/v1/mixtures/{self.regular_user_mixture_id}/experiments',
            json={
                'property_name': 'Viscosity',
                'value': 1.7,
                'experimental_conditions': 'Room temperature',
                'date_performed': date.today().isoformat()
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['property_name'], 'Viscosity')
        self.assertEqual(data['numeric_value'], 1.7)

# Authenticated User Tests - Other Users' Data
    
    @patch_supabase(load_data=True)
    def test_authenticated_user_cannot_update_others_mixture(self, mock_client):
        """Test that authenticated users cannot update mixtures they don't own."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.other_user_mixture_id}',
            json={
                'name': 'Attempted Update',
                'description': 'An attempted update of another user\'s mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])
    
    @patch_supabase(load_data=True)
    def test_authenticated_user_cannot_add_prediction_to_others_mixture(self, mock_client):
        """Test that authenticated users cannot add predictions to mixtures they don't own."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        # Make request
        response = self.client.post(
            f'/api/v1/mixtures/{self.other_user_mixture_id}/predictions',
            json={
                'property_name': 'Viscosity',
                'value': 1.5,
                'confidence': 0.8,
                'calculation_method': 'CryoProtect Scoring'
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])
    
    @patch_supabase(load_data=True)
    def test_authenticated_user_cannot_add_experiment_to_others_mixture(self, mock_client):
        """Test that authenticated users cannot add experiments to mixtures they don't own."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        # Make request
        response = self.client.post(
            f'/api/v1/mixtures/{self.other_user_mixture_id}/experiments',
            json={
                'property_name': 'Viscosity',
                'value': 1.7,
                'experimental_conditions': 'Room temperature',
                'date_performed': date.today().isoformat()
            },
            headers=self.regular_user_headers
        )
        
        # Assertions
        self.assertIn(response.status_code, [401, 403])

# Admin User Tests
    
    @patch_supabase(load_data=True)
    def test_admin_user_can_update_any_mixture(self, mock_client):
        """Test that admin users can update any mixture."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.admin_user_id, self.admin_auth_id, is_admin=True)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'name': 'Admin Updated Mixture'}]
        )
        
        updated_mixture = {
            **self.other_user_mixture,
            'name': 'Admin Updated Mixture'
        }
        
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[updated_mixture]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.other_user_mixture_id}',
            json={
                'name': 'Admin Updated Mixture',
                'description': 'An admin updated test mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
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
    def test_admin_user_can_add_prediction_to_any_mixture(self, mock_client):
        """Test that admin users can add predictions to any mixture."""
        # Mock the auth user
        self.mock_auth_user(mock_client, self.admin_user_id, self.admin_auth_id, is_admin=True)
        
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        new_prediction_id = str(uuid.uuid4())
        new_prediction = {
            'id': new_prediction_id,
            'mixture_id': self.other_user_mixture_id,
            'property_name': 'Viscosity',
            'numeric_value': 1.5,
            'confidence': 0.8,
            'calculation_method': 'CryoProtect Scoring',
            'created_by': self.admin_user_id
        }
        
        mock_client.table('predictions').insert().execute.return_value = MagicMock(
            data=[new_prediction]
        )
        
        # Make request
        response = self.client.post(
            f'/api/v1/mixtures/{self.other_user_mixture_id}/predictions',
            json={
                'property_name': 'Viscosity',
                'value': 1.5,
                'confidence': 0.8,
                'calculation_method': 'CryoProtect Scoring'
            },
            headers=self.admin_user_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 201)
        data = json.loads(response.data)
        self.assertEqual(data['property_name'], 'Viscosity')
        self.assertEqual(data['numeric_value'], 1.5)
    
    # Service Role Tests
    
    @patch('api.utils.USE_SERVICE_ROLE', True)
    @patch('api.utils.USER_ID', 'service-role-user-id')
    @patch_supabase(load_data=True)
    def test_service_role_can_access_any_data(self, mock_client):
        """Test that service role can access any data."""
        # Mock the Supabase response
        mock_client.table('mixtures').select().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'created_by': self.other_user_id}]
        )
        
        mock_client.table('mixtures').update().eq().execute.return_value = MagicMock(
            data=[{'id': self.other_user_mixture_id, 'name': 'Service Role Updated Mixture'}]
        )
        
        updated_mixture = {
            **self.other_user_mixture,
            'name': 'Service Role Updated Mixture'
        }
        
        mock_client.table('mixture_with_components').select().eq().execute.return_value = MagicMock(
            data=[updated_mixture]
        )
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{self.other_user_mixture_id}',
            json={
                'name': 'Service Role Updated Mixture',
                'description': 'A service role updated test mixture',
                'components': [
                    {
                        'molecule_id': self.molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=self.service_role_headers
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['name'], 'Service Role Updated Mixture')

# Property Calculation Tests
    
    @patch_supabase(load_data=True)
    def test_calculate_molecular_properties(self, mock_client):
        """Test the molecular property calculation endpoint."""
        # Mock the RDKit utility function
        with patch('api.rdkit_utils.calculate_all_properties') as mock_calculate:
            # Set up mock
            mock_calculate.return_value = {
                'hydrogen_bonding': {'donors': 3, 'acceptors': 3, 'total': 6},
                'logp': -1.76,
                'tpsa': 60.69,
                'molecular_properties': {
                    'molecular_weight': 92.09,
                    'heavy_atom_count': 6
                },
                'functional_groups': {'alcohol': 3},
                'permeability': {'rule_of_5_violations': 0}
            }
            
            # Make request
            response = self.client.post(
                '/api/v1/rdkit/properties',
                json={
                    'molecule_data': 'C(C(CO)O)O',
                    'input_format': 'smiles'
                }
            )
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIn('hydrogen_bonding', data)
            self.assertIn('logp', data)
            self.assertIn('tpsa', data)
    
    # Search Tests
    
    @patch_supabase(load_data=True)
    def test_search_molecules_by_name(self, mock_client):
        """Test the molecule search endpoint."""
        # Mock the Supabase response
        mock_client.table('molecules').select().ilike().limit().execute.return_value = MagicMock(
            data=[self.molecule]
        )
        
        # Make request
        response = self.client.get('/api/v1/molecules/search?name=Glyc')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]['name'], 'Glycerol')
    
    # Error Handling Tests
    
    @patch_supabase(load_data=True)
    def test_not_found_error(self, mock_client):
        """Test 404 error handling."""
        # Mock the Supabase response
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
        # Mock the auth user
        self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
        
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
        # Mock the Supabase response to raise an exception
        mock_client.table('molecule_with_properties').select().eq().execute.side_effect = Exception("Database error")
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{self.molecule_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 500)
        data = json.loads(response.data)
        self.assertIn('message', data)


if __name__ == '__main__':
    unittest.main()

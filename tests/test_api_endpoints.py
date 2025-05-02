"""
CryoProtect Analyzer - API Endpoints Tests

This module contains tests for the API endpoints.
It tests request validation, response formats, and error handling.
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

from tests.mock_supabase.helpers import patch_supabase

from unittest.mock import patch

class TestAPIEndpoints(MockSupabaseBaseTestCase):
    """Test cases for API endpoints."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # We no longer patch handle_supabase_error to test the standardized error handling

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()

        # Patch model create methods to bypass Supabase/mock errors
        from unittest.mock import patch
        self.mixture_create_patcher = patch('api.models.Mixture.create', side_effect=lambda data: data)
        self.molecule_create_patcher = patch('api.models.Molecule.create', side_effect=lambda data: data)
        self.userprofile_create_patcher = patch('api.models.UserProfile.create', side_effect=lambda data: data)
        self.mixture_getall_patcher = patch('api.models.Mixture.get_all', side_effect=lambda *args, **kwargs: [self.sample_mixture])
        self.molecule_getall_patcher = patch('api.models.Molecule.get_all', side_effect=lambda *args, **kwargs: [self.sample_molecule])
        self.userprofile_getall_patcher = patch('api.models.UserProfile.get_all', side_effect=lambda *args, **kwargs: [self.sample_user_profile])
        self.mixture_create_patcher.start()
        self.molecule_create_patcher.start()
        self.userprofile_create_patcher.start()
        self.mixture_getall_patcher.start()
        self.molecule_getall_patcher.start()
        self.userprofile_getall_patcher.start()

        # Sample data for testing
        self.sample_molecule_id = str(uuid.uuid4())
        self.sample_mixture_id = str(uuid.uuid4())
        self.sample_user_profile_id = str(uuid.uuid4())

        # Create a test user profile
        self.sample_user_profile = {
            'id': self.sample_user_profile_id,
            'auth_user_id': str(uuid.uuid4()),  # Unique auth_user_id
            'display_name': 'Test User',
            'email': 'test@example.com',
            'affiliation': 'Test Organization'  # Changed from 'organization' to 'affiliation'
        }

        self.sample_molecule = {
            'id': self.sample_molecule_id,
            'cid': 123456,
            'name': 'Glycerol',
            'molecular_formula': 'C3H8O3',
            'smiles': 'C(C(CO)O)O',
            'inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
            'inchikey': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
            'properties': [],
            'created_by': self.sample_user_profile_id  # Reference to user_profile.id
        }

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
            ]
        }

        # Insert sample data into mock Supabase tables
        from api.models import Mixture, Molecule, UserProfile
        Mixture.create(self.sample_mixture)
        Molecule.create(self.sample_molecule)
        UserProfile.create(self.sample_user_profile)

        # Insert sample data into mock Supabase tables
        from api.models import Mixture, Molecule, UserProfile
        Mixture.create(self.sample_mixture)
        Molecule.create(self.sample_molecule)
        UserProfile.create(self.sample_user_profile)

    def tearDown(self):
        # Stop patchers
        self.mixture_create_patcher.stop()
        self.molecule_create_patcher.stop()
        self.userprofile_create_patcher.stop()
        self.mixture_getall_patcher.stop()
        self.molecule_getall_patcher.stop()
        self.userprofile_getall_patcher.stop()
        super().tearDown()

        self.sample_prediction = {
            'id': str(uuid.uuid4()),
            'mixture_id': self.sample_mixture_id,
            'property_name': 'Freezing Point',
            'numeric_value': -15.3,
            'confidence': 0.9,
            'calculation_method': 'CryoProtect Scoring'
        }

        self.sample_experiment = {
            'id': str(uuid.uuid4()),
            'mixture_id': self.sample_mixture_id,
            'property_name': 'Freezing Point',
            'numeric_value': -14.8,
            'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
            'date_performed': '2025-04-15'
        }

        self.sample_comparison = {
            'prediction': {
                'value': -15.3,
                'confidence': 0.9,
                'method': 'CryoProtect Scoring'
            },
            'experiment': {
                'value': -14.8,
                'conditions': 'Standard pressure, cooling rate 1°C/min',
                'date': '2025-04-15'
            },
            'difference': 0.5,
            'percent_error': 3.38
        }

        # RDKit test data
        self.sample_rdkit_properties = {
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
    
    # Molecule API Tests
    
    @patch_supabase(load_data=True)
    def test_get_molecules(self, mock_client):
        """Test GET /api/molecules endpoint with standardized response formatting."""
        # Set up the mock Supabase client and load test data
        from tests.mock_supabase.helpers import reset_mock_data, load_test_data
        reset_mock_data()
        load_test_data()
        
        # Add a test molecule to the mock database
        test_molecule_id = str(uuid.uuid4())
        test_molecule = {
            'id': test_molecule_id,
            'name': 'Test Molecule',
            'smiles': 'C(C(CO)O)O',
            'molecular_formula': 'C3H8O3'
        }
        
        mock_client.table('molecules').insert(test_molecule).execute()
        
        # Mock the response for molecules_with_properties
        mock_client.table('molecules_with_properties').select('*').range().execute.return_value = MagicMock(
            data=[test_molecule],
            error=None
        )
        
        # Make request
        response = self.client.get('/api/v1/molecules')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
        self.assertGreater(len(data), 0)
        
        # Check that the response is properly formatted
        molecule = data[0]
        self.assertEqual(molecule['id'], test_molecule_id)
        self.assertEqual(molecule['name'], 'Test Molecule')
        self.assertEqual(molecule['smiles'], 'C(C(CO)O)O')
        self.assertEqual(molecule['molecular_formula'], 'C3H8O3')
    
    @patch_supabase(load_data=True)
    def test_get_molecule(self, mock_client):
        """Test GET /api/molecules/<id> endpoint."""
        # Set up the mock Supabase client and load test data
        from tests.mock_supabase.helpers import reset_mock_data, load_test_data
        reset_mock_data()
        load_test_data()
        
        # Add a test molecule to the mock database
        molecule_id = str(uuid.uuid4())
        mock_client.table('molecules').insert({
            'id': molecule_id,
            'name': 'Test Molecule',
            'smiles': 'C(C(CO)O)O',
            'molecular_formula': 'C3H8O3'
        }).execute()
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{molecule_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], molecule_id)
    
    @patch_supabase(load_data=True)
    def test_import_from_pubchem(self, mock_client):
        """Test POST /api/molecules/import endpoint."""
        # Mock the create_from_pubchem method
        with patch.object(Molecule, 'create_from_pubchem') as mock_create:
            # Set up mock
            mock_create.return_value = self.sample_molecule
            
            # Make request
            response = self.client.post(
                '/api/molecules/import',
                json={'cid': 123456}
            )
            
            # Assertions
            self.assertEqual(response.status_code, 201)
    
    @patch_supabase(load_data=True)
    def test_search_molecules_by_name(self, mock_client):
        """Test GET /api/molecules/search endpoint with name parameter."""
        # Make request
        response = self.client.get('/api/v1/molecules/search?name=Glyc')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
    
    # Mixture API Tests
    
    @patch_supabase(load_data=True)
    def test_get_mixtures(self, mock_client):
        """Test GET /api/mixtures endpoint."""
        # Make request
        response = self.client.get('/api/v1/mixtures')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIsInstance(data, list)
    
    @patch_supabase(load_data=True)
    def test_get_mixture(self, mock_client):
        """Test GET /api/mixtures/<id> endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Make request
        response = self.client.get(f'/api/v1/mixtures/{mixture_id}')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['id'], mixture_id)
    
    @patch_supabase(load_data=True)
    def test_create_mixture(self, mock_client):
        """Test POST /api/mixtures endpoint with authentication."""
        # Get a molecule ID from the mock data
        molecules = Molecule.get_all()
        molecule_id = molecules[0]['id']
        
        # First test without authentication - should fail
        response_no_auth = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            }
        )
        
        # Should return 401 Unauthorized
        self.assertEqual(response_no_auth.status_code, 401)
        
        # Mock authentication
        auth_token = str(uuid.uuid4())
        auth_headers = {'Authorization': f'Bearer {auth_token}'}
        
        # Mock the auth.get_user method
        mock_user = MagicMock()
        mock_user.id = str(uuid.uuid4())
        mock_client.auth.get_user.return_value = mock_user
        
        # Now test with authentication
        response_with_auth = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            },
            headers=auth_headers
        )
        
        # Should succeed with 201 Created
        self.assertEqual(response_with_auth.status_code, 201)
    
    @patch_supabase(load_data=True)
    def test_update_mixture(self, mock_client):
        """Test PUT /api/mixtures/<id> endpoint."""
        # Get a mixture ID and molecule ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        molecules = Molecule.get_all()
        molecule_id = molecules[0]['id']
        
        # Make request
        response = self.client.put(
            f'/api/v1/mixtures/{mixture_id}',
            json={
                'name': 'Updated Mixture',
                'description': 'A test mixture',
                'components': [
                    {
                        'molecule_id': molecule_id,
                        'concentration': 100,
                        'concentration_unit': '%'
                    }
                ]
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data['name'], 'Updated Mixture')
    
    @patch_supabase(load_data=True)
    def test_calculate_mixture_score(self, mock_client):
        """Test GET /api/mixtures/<id>/score endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the calculate_score method
        with patch.object(Mixture, 'calculate_score') as mock_calculate:
            # Set up mock
            mock_calculate.return_value = 85.5
            
            # Make request
            response = self.client.get(f'/api/v1/mixtures/{mixture_id}/score')
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertEqual(data['score'], 85.5)
    
    # Prediction API Tests
    
    @patch_supabase(load_data=True)
    def test_get_predictions(self, mock_client):
        """Test GET /api/mixtures/<id>/predictions endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the get_predictions_for_mixture method
        with patch.object(Prediction, 'get_predictions_for_mixture') as mock_get:
            # Set up mock
            mock_get.return_value = [self.sample_prediction]
            
            # Make request
            response = self.client.get(f'/api/v1/mixtures/{mixture_id}/predictions')
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIsInstance(data, list)
    
    @patch_supabase(load_data=True)
    def test_add_prediction(self, mock_client):
        """Test POST /api/mixtures/<id>/predictions endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the add_prediction method
        with patch.object(Prediction, 'add_prediction') as mock_add:
            # Set up mock
            mock_add.return_value = self.sample_prediction
            
            # Make request
            response = self.client.post(
                f'/api/v1/mixtures/{mixture_id}/predictions',
                json={
                    'property_name': 'Freezing Point',
                    'value': -15.3,
                    'confidence': 0.9,
                    'calculation_method': 'CryoProtect Scoring'
                }
            )
            
            # Assertions
            self.assertEqual(response.status_code, 201)
            data = json.loads(response.data)
            self.assertEqual(data['property_name'], 'Freezing Point')
    
    # Experiment API Tests
    
    @patch_supabase(load_data=True)
    def test_get_experiments(self, mock_client):
        """Test GET /api/mixtures/<id>/experiments endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the get_experiments_for_mixture method
        with patch.object(Experiment, 'get_experiments_for_mixture') as mock_get:
            # Set up mock
            mock_get.return_value = [self.sample_experiment]
            
            # Make request
            response = self.client.get(f'/api/v1/mixtures/{mixture_id}/experiments')
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIsInstance(data, list)
    
    @patch_supabase(load_data=True)
    def test_record_experiment(self, mock_client):
        """Test POST /api/mixtures/<id>/experiments endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the record_experiment method
        with patch.object(Experiment, 'record_experiment') as mock_record:
            # Set up mock
            mock_record.return_value = self.sample_experiment
            
            # Make request
            response = self.client.post(
                f'/api/v1/mixtures/{mixture_id}/experiments',
                json={
                    'property_name': 'Freezing Point',
                    'value': -14.8,
                    'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
                    'date_performed': '2025-04-15'
                }
            )
            
            # Assertions
            self.assertEqual(response.status_code, 201)
            data = json.loads(response.data)
            self.assertEqual(data['property_name'], 'Freezing Point')
    
    # Comparison API Tests
    
    @patch_supabase(load_data=True)
    def test_compare_prediction_with_experiment(self, mock_client):
        """Test GET /api/mixtures/<id>/compare endpoint."""
        # Get a mixture ID from the mock data
        mixtures = Mixture.get_all()
        mixture_id = mixtures[0]['id']
        
        # Mock the compare_prediction_with_experiment method
        with patch.object(Comparison, 'compare_prediction_with_experiment') as mock_compare:
            # Set up mock
            mock_compare.return_value = self.sample_comparison
            
            # Make request
            response = self.client.get(
                f'/api/v1/mixtures/{mixture_id}/compare?property_name=Freezing Point'
            )
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIn('prediction', data)
            self.assertIn('experiment', data)
            self.assertIn('difference', data)
            self.assertIn('percent_error', data)
    
    # RDKit API Tests
    
    @patch_supabase(load_data=True)
    def test_calculate_molecular_properties(self, mock_client):
        """Test POST /api/rdkit/properties endpoint."""
        # Mock the calculate_all_properties function
        with patch('api.rdkit_utils.calculate_all_properties') as mock_calculate:
            # Set up mock
            mock_calculate.return_value = self.sample_rdkit_properties
            
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
            self.assertIn('molecular_properties', data)
            self.assertIn('functional_groups', data)
            self.assertIn('permeability', data)
    
    @patch_supabase(load_data=True)
    def test_generate_molecule_visualization(self, mock_client):
        """Test POST /api/rdkit/visualize endpoint."""
        # Mock the generate_molecule_svg function
        with patch('api.rdkit_utils.generate_molecule_svg') as mock_generate:
            # Set up mock
            mock_generate.return_value = '<svg>...</svg>'
            
            # Make request
            response = self.client.post(
                '/api/v1/rdkit/visualize',
                json={
                    'molecule_data': 'C(C(CO)O)O',
                    'input_format': 'smiles',
                    'width': 400,
                    'height': 300
                }
            )
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIn('svg', data)
    
    @patch_supabase(load_data=True)
    def test_substructure_search(self, mock_client):
        """Test POST /api/rdkit/substructure endpoint."""
        # Mock the perform_substructure_search function
        with patch('api.rdkit_utils.perform_substructure_search') as mock_search:
            # Set up mock
            mock_search.return_value = {
                'match': True,
                'match_count': 3,
                'matches': [[0, 1, 2], [3, 4, 5], [6, 7, 8]],
                'visualization': '<svg>...</svg>'
            }
            
            # Make request
            response = self.client.post(
                '/api/v1/rdkit/substructure',
                json={
                    'query_mol_data': '[OH]',
                    'target_mol_data': 'C(C(CO)O)O',
                    'query_format': 'smarts',
                    'target_format': 'smiles'
                }
            )
            
            # Assertions
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.data)
            self.assertIn('match', data)
            self.assertIn('match_count', data)
            self.assertIn('matches', data)
            self.assertIn('visualization', data)
    
    @patch('api.rdkit_utils.calculate_similarity')
    def test_calculate_similarity(self, mock_calculate_similarity):
        """Test POST /api/rdkit/similarity endpoint."""
        # Set up mock
        mock_calculate_similarity.return_value = {
            'tanimoto': 0.75,
            'dice': 0.85,
            'fingerprint_type': 'morgan'
        }
        
        # Make request
        response = self.client.post(
            '/api/v1/rdkit/similarity',
            json={
                'mol1_data': 'C(C(CO)O)O',
                'mol2_data': 'OCCO',
                'mol1_format': 'smiles',
                'mol2_format': 'smiles',
                'fingerprint_type': 'morgan'
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('tanimoto', data)
        self.assertIn('dice', data)
        self.assertIn('fingerprint_type', data)
    
    # Error Handling Tests
    
    @patch_supabase(load_data=True)
    def test_not_found(self, mock_client):
        """Test 404 error handling with standardized format."""
        response = self.client.get('/api/v1/nonexistent')
        self.assertEqual(response.status_code, 404)
        data = json.loads(response.data)
        self.assertIn('message', data)
        # Check for standardized error response format
        if 'status' in data:
            self.assertEqual(data['status'], 'error')
            self.assertIn('meta', data)
    
    @patch_supabase(load_data=True)
    def test_molecule_not_found(self, mock_client):
        """Test molecule not found error handling with standardized format."""
        # Mock the Supabase response
        mock_client.table('molecules_with_properties').select('*').eq().execute.return_value = MagicMock(
            data=[],
            error=None
        )
        
        # Make request
        response = self.client.get(f'/api/v1/molecules/{uuid.uuid4()}')
        
        # Assertions
        self.assertEqual(response.status_code, 404)
        data = json.loads(response.data)
        self.assertIn('message', data)
        # Check for standardized error response format
        if 'status' in data:
            self.assertEqual(data['status'], 'error')
            self.assertIn('meta', data)
            self.assertIn('context', data['meta'])
    
    @patch_supabase(load_data=True)
    def test_bad_request(self, mock_client):
        """Test bad request error handling with standardized format."""
        # Make request with invalid JSON
        response = self.client.post(
            '/api/v1/mixtures',
            data='invalid json',
            content_type='application/json'
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('message', data)
        # Check for standardized error response format
        if 'status' in data:
            self.assertEqual(data['status'], 'error')
            self.assertIn('meta', data)
    
    @patch_supabase(load_data=True)
    def test_validation_error(self, mock_client):
        """Test validation error handling with standardized format."""
        # Make request with missing required fields
        response = self.client.post(
            '/api/v1/mixtures',
            json={
                'name': 'Test Mixture',
                # Missing components
            }
        )
        
        # Assertions
        self.assertEqual(response.status_code, 400)
        data = json.loads(response.data)
        self.assertIn('message', data)
        # Check for standardized error response format
        if 'status' in data:
            self.assertEqual(data['status'], 'error')
            self.assertIn('meta', data)
            if 'error_type' in data['meta']:
                self.assertEqual(data['meta']['error_type'], 'ValidationError')


if __name__ == '__main__':
    unittest.main()
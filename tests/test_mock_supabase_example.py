"""
CryoProtect Analyzer - Mock Supabase Example Tests

This module demonstrates how to use the mock Supabase client in tests.
"""

import sys
import os
import uuid
from datetime import datetime
from unittest.mock import MagicMock, patch
from flask import current_app

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

# Import the mock Supabase helpers
from tests.mock_supabase.helpers import (
    patch_supabase,
    mock_rpc_function,
    configure_test_scenario
)

# Import the app and API resources
from api.models import Molecule, Mixture
from app import create_app


class TestMockSupabaseWithDecorator(BaseTestCase):
    """Test using the mock Supabase client with the decorator approach."""
    
    @patch_supabase(load_data=True)
    def test_get_molecules(self, mock_client):
        """Test getting molecules."""
        # Configure the mock client to return expected data
        expected_molecules = [
            {
                'id': str(uuid.uuid4()),
                'name': 'Test Molecule 1',
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'molecular_formula': 'C9H8O4',
                'created_at': datetime.now().isoformat(),
                'updated_at': datetime.now().isoformat(),
                'created_by': '00000000-0000-0000-0000-000000000001'
            },
            {
                'id': str(uuid.uuid4()),
                'name': 'Test Molecule 2',
                'smiles': 'C1=CC=C(C=C1)C(=O)O',
                'molecular_formula': 'C7H6O2',
                'created_at': datetime.now().isoformat(),
                'updated_at': datetime.now().isoformat(),
                'created_by': '00000000-0000-0000-0000-000000000001'
            }
        ]
        
        # Set up the mock to return our expected data
        mock_response = MagicMock()
        mock_response.data = expected_molecules
        mock_response.error = None
        
        # Configure the mock query builder with proper MagicMock chaining
        mock_execute = MagicMock()
        mock_execute.execute.return_value = mock_response
        
        mock_range = MagicMock()
        mock_range.range.return_value = mock_execute
        
        # Set up the mock client
        mock_table = MagicMock()
        mock_table.select.return_value = mock_range
        
        # Patch the Molecule.get_table method to use our mock
        with patch.object(Molecule, 'get_table', return_value=mock_table) as mock_get_table:
            # Get all molecules
            molecules = Molecule.get_all()
            
            # Verify we got the expected data
            self.assertEqual(len(molecules), len(expected_molecules))
            self.assertEqual(molecules, expected_molecules)
            self.assertIn('name', molecules[0])
            self.assertIn('smiles', molecules[0])
            
            # Verify the mock was called correctly
            mock_get_table.assert_called_once()
            mock_table.select.assert_called_with('*')
    
    @patch_supabase(load_data=True)
    def test_create_molecule(self, mock_client):
        """Test creating a molecule."""
        # Create a new molecule data
        molecule_data = {
            'name': 'Test Molecule',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'molecular_formula': 'C9H8O4',
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Expected response after creation
        expected_molecule = molecule_data.copy()
        expected_molecule['id'] = str(uuid.uuid4())
        expected_molecule['created_at'] = datetime.now().isoformat()
        expected_molecule['updated_at'] = datetime.now().isoformat()
        
        # Set up the mock to return our expected data
        mock_response = MagicMock()
        mock_response.data = [expected_molecule]
        mock_response.error = None
        
        # Configure the mock query builder with proper MagicMock chaining
        mock_execute = MagicMock()
        mock_execute.execute.return_value = mock_response
        
        # Set up the mock client
        mock_table = MagicMock()
        mock_table.insert.return_value = mock_execute
        
        # Patch the Molecule.get_table method to use our mock
        with patch.object(Molecule, 'get_table', return_value=mock_table):
            # Create the molecule
            molecule = Molecule.create(molecule_data)
            
            # Verify the molecule was created with expected values
            self.assertIsNotNone(molecule)
            self.assertEqual(molecule['name'], 'Test Molecule')
            self.assertEqual(molecule['smiles'], 'CC(=O)OC1=CC=CC=C1C(=O)O')
            self.assertEqual(molecule['id'], expected_molecule['id'])
            
            # Verify the mock was called correctly
            mock_table.insert.assert_called_once_with(molecule_data)
    
    @patch_supabase(load_data=True)
    def test_update_molecule(self, mock_client):
        """Test updating a molecule."""
        # Create a molecule to update
        molecule_id = str(uuid.uuid4())
        original_molecule = {
            'id': molecule_id,
            'name': 'Original Name',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'molecular_formula': 'C9H8O4',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Update data
        updated_data = {
            'name': 'Updated Molecule Name'
        }
        
        # Expected updated molecule
        updated_molecule = original_molecule.copy()
        updated_molecule['name'] = 'Updated Molecule Name'
        updated_molecule['updated_at'] = datetime.now().isoformat()
        
        # Set up mock for get_all to return our original molecule
        mock_get_all_response = MagicMock()
        mock_get_all_response.data = [original_molecule]
        mock_get_all_response.error = None
        
        # Set up mock for update
        mock_update_response = MagicMock()
        mock_update_response.data = [updated_molecule]
        mock_update_response.error = None
        
        # Set up mock for get after update
        mock_get_updated_response = MagicMock()
        mock_get_updated_response.data = [updated_molecule]
        mock_get_updated_response.error = None
        
        # Configure the mock query builder for get_all
        mock_get_all_execute = MagicMock()
        mock_get_all_execute.execute.return_value = mock_get_all_response
        
        mock_get_all_range = MagicMock()
        mock_get_all_range.range.return_value = mock_get_all_execute
        
        # Configure the mock query builder for update
        mock_update_execute = MagicMock()
        mock_update_execute.execute.return_value = mock_update_response
        
        mock_update_eq = MagicMock()
        mock_update_eq.eq.return_value = mock_update_execute
        
        # Configure the mock query builder for get after update
        mock_get_updated_execute = MagicMock()
        mock_get_updated_execute.execute.return_value = mock_get_updated_response
        
        mock_get_updated_eq = MagicMock()
        mock_get_updated_eq.eq.return_value = mock_get_updated_execute
        
        # Set up the mock client
        mock_table = MagicMock()
        mock_table.select.return_value = mock_get_all_range
        mock_table.update.return_value = mock_update_eq
        
        # Patch the Molecule.get_table method to use our mock
        with patch.object(Molecule, 'get_table', return_value=mock_table):
            # First, mock get_all to return our molecule
            molecules = Molecule.get_all()
            molecule = molecules[0]
            
            # Update the molecule
            updated = Molecule.update(molecule['id'], updated_data)
            
            # Verify the update
            self.assertEqual(updated['name'], 'Updated Molecule Name')
            
            # Verify the mock was called correctly
            mock_table.update.assert_called_once_with(updated_data)
            mock_update_eq.eq.assert_called_once_with('id', molecule_id)
            
            # For the get after update, we need to change the mock
            mock_table.select.return_value = mock_get_updated_eq
            
            # Verify we can retrieve it with the update
            retrieved = Molecule.get(molecule['id'])
            self.assertEqual(retrieved['name'], 'Updated Molecule Name')
    
    @patch_supabase(load_data=True)
    def test_delete_molecule(self, mock_client):
        """Test deleting a molecule."""
        # Create a molecule to delete
        molecule_id = str(uuid.uuid4())
        molecule = {
            'id': molecule_id,
            'name': 'Molecule to Delete',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'molecular_formula': 'C9H8O4',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Set up mock for get_all to return our molecule
        mock_get_all_response = MagicMock()
        mock_get_all_response.data = [molecule]
        mock_get_all_response.error = None
        
        # Set up mock for delete
        mock_delete_response = MagicMock()
        mock_delete_response.data = [molecule]  # Return the deleted molecule
        mock_delete_response.error = None
        
        # Set up mock for get after delete (should return None)
        mock_get_after_delete = MagicMock()
        mock_get_after_delete.data = []  # Empty list means not found
        mock_get_after_delete.error = None
        
        # Configure the mock query builder for get_all
        mock_get_all_execute = MagicMock()
        mock_get_all_execute.execute.return_value = mock_get_all_response
        
        mock_get_all_range = MagicMock()
        mock_get_all_range.range.return_value = mock_get_all_execute
        
        # Configure the mock query builder for delete
        mock_delete_execute = MagicMock()
        mock_delete_execute.execute.return_value = mock_delete_response
        
        mock_delete_eq = MagicMock()
        mock_delete_eq.eq.return_value = mock_delete_execute
        
        # Configure the mock query builder for get after delete
        mock_get_after_delete_execute = MagicMock()
        mock_get_after_delete_execute.execute.return_value = mock_get_after_delete
        
        mock_get_after_delete_eq = MagicMock()
        mock_get_after_delete_eq.eq.return_value = mock_get_after_delete_execute
        
        # Set up the mock client
        mock_table = MagicMock()
        mock_table.select.return_value = mock_get_all_range
        mock_table.delete.return_value = mock_delete_eq
        
        # Patch the Molecule.get_table method to use our mock
        with patch.object(Molecule, 'get_table', return_value=mock_table):
            # Get a molecule to delete
            molecules = Molecule.get_all()
            molecule = molecules[0]
            
            # Delete the molecule
            result = Molecule.delete(molecule['id'])
            
            # Verify the deletion
            self.assertTrue(result)
            
            # Verify the mock was called correctly
            mock_delete_eq.eq.assert_called_once_with('id', molecule_id)
            
            # For the get after delete, we need to change the mock
            mock_table.select.return_value = mock_get_after_delete_eq
            
            # Verify it's gone
            retrieved = Molecule.get(molecule['id'])
            self.assertIsNone(retrieved)
    
    @patch_supabase(load_data=True)
    def test_filter_molecules(self, mock_client):
        """Test filtering molecules."""
        # Create a molecule with specific attributes
        molecule_data = {
            'name': 'Filterable Molecule',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'molecular_formula': 'C9H8O4',
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Expected molecule after creation
        expected_molecule = molecule_data.copy()
        expected_molecule['id'] = str(uuid.uuid4())
        expected_molecule['created_at'] = datetime.now().isoformat()
        expected_molecule['updated_at'] = datetime.now().isoformat()
        
        # Set up mock for filter
        mock_filter_response = MagicMock()
        mock_filter_response.data = [expected_molecule]
        mock_filter_response.error = None
        
        # Configure the mock query builder for filter
        mock_filter_execute = MagicMock()
        mock_filter_execute.execute.return_value = mock_filter_response
        
        mock_filter_range = MagicMock()
        mock_filter_range.range.return_value = mock_filter_execute
        
        mock_filter_eq = MagicMock()
        mock_filter_eq.eq.return_value = mock_filter_range
        
        # Set up the mock client
        mock_table = MagicMock()
        mock_table.select.return_value = mock_filter_eq
        
        # Patch the Molecule.get_table method to use our mock
        with patch.object(Molecule, 'get_table', return_value=mock_table):
            # Filter molecules
            filtered = Molecule.filter({'name': 'Filterable Molecule'})
            
            # Verify the filter
            self.assertEqual(len(filtered), 1)
            self.assertEqual(filtered[0]['name'], 'Filterable Molecule')
            
            # Verify the mock was called correctly
            mock_table.select.assert_called_with('*')
            mock_filter_eq.eq.assert_called_with('name', 'Filterable Molecule')


class TestMockSupabaseWithBaseClass(MockSupabaseBaseTestCase):
    """Test using the mock Supabase client with the base class approach."""
    
    def setUp(self):
        """Set up the test case."""
        super().setUp()
        
        # Configure the mock client for this test class
        self.mixture_id = str(uuid.uuid4())
        self.test_mixture = {
            'id': self.mixture_id,
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Set up mock responses
        self.mock_create_response = MagicMock()
        self.mock_create_response.data = [self.test_mixture]
        self.mock_create_response.error = None
        
        self.mock_update_response = MagicMock()
        updated_mixture = self.test_mixture.copy()
        updated_mixture['name'] = 'Updated Mixture Name'
        self.mock_update_response.data = [updated_mixture]
        self.mock_update_response.error = None
        
        self.mock_delete_response = MagicMock()
        self.mock_delete_response.data = [self.test_mixture]
        self.mock_delete_response.error = None
        
        # Configure the mock client with proper MagicMock chaining
        # For insert
        mock_insert_execute = MagicMock()
        mock_insert_execute.execute.return_value = self.mock_create_response
        
        # For update
        mock_update_execute = MagicMock()
        mock_update_execute.execute.return_value = self.mock_update_response
        
        mock_update_eq = MagicMock()
        mock_update_eq.eq.return_value = mock_update_execute
        
        # For delete
        mock_delete_execute = MagicMock()
        mock_delete_execute.execute.return_value = self.mock_delete_response
        
        mock_delete_eq = MagicMock()
        mock_delete_eq.eq.return_value = mock_delete_execute
        
        # Set up the mock table
        self.mock_table = MagicMock()
        self.mock_table.insert.return_value = mock_insert_execute
        self.mock_table.update.return_value = mock_update_eq
        self.mock_table.delete.return_value = mock_delete_eq
        
        # Set up the mock client
        self.mock_client.table.return_value = self.mock_table
    
    def test_mixture_operations(self):
        """Test mixture operations."""
        # Create a new mixture
        mixture_data = {
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Patch the Mixture.get_table method to use our mock
        with patch.object(Mixture, 'get_table', return_value=self.mock_table):
            # Create the mixture
            mixture = Mixture.create(mixture_data)
            
            # Verify the mixture was created
            self.assertIsNotNone(mixture)
            self.assertEqual(mixture['name'], 'Test Mixture')
            self.assertEqual(mixture['id'], self.mixture_id)
            
            # Verify the mock was called correctly
            self.mock_table.insert.assert_called_once_with(mixture_data)
            
            # Update the mixture
            updated_data = {
                'name': 'Updated Mixture Name'
            }
            
            updated = Mixture.update(mixture['id'], updated_data)
            
            # Verify the update
            self.assertEqual(updated['name'], 'Updated Mixture Name')
            
            # Verify the mock was called correctly
            self.mock_table.update.assert_called_once_with(updated_data)
            self.mock_table.update.return_value.eq.assert_called_once_with('id', self.mixture_id)
            
            # Delete the mixture
            result = Mixture.delete(mixture['id'])
            
            # Verify the deletion
            self.assertTrue(result)
            
            # Verify the mock was called correctly
            self.mock_table.delete.assert_called_once()
            self.mock_table.delete.return_value.eq.assert_called_once_with('id', self.mixture_id)


class TestMockRPC(MockSupabaseBaseTestCase):
    """Test mocking RPC functions."""
    
    def setUp(self):
        """Set up the test case."""
        super().setUp()
        
        # Configure mock RPC responses
        self.test_function_response = MagicMock()
        self.test_function_response.data = {'success': True, 'data': 'test'}
        self.test_function_response.error = None
        
        self.mixture_score_response = MagicMock()
        self.mixture_score_response.data = 85.5  # A score between 0 and 100
        self.mixture_score_response.error = None
        
        # Configure the mock client for RPC
        mock_rpc_execute = MagicMock()
        mock_rpc_execute.execute.return_value = self.test_function_response
        
        # Set up the mock client
        self.mock_client.rpc.return_value = mock_rpc_execute
        
        # Set up a mixture for testing
        self.mixture_id = str(uuid.uuid4())
        self.mixture = {
            'id': self.mixture_id,
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
    
    def test_custom_rpc(self):
        """Test custom RPC function."""
        # Register a custom RPC function
        mock_rpc_function('test_function', result={'success': True, 'data': 'test'})
        
        # Call the RPC function
        response = self.mock_client.rpc('test_function', {'param': 'value'}).execute()
        
        # Verify the response
        self.assertIsNone(response.error)
        self.assertEqual(response.data, {'success': True, 'data': 'test'})
        
        # Verify the mock was called correctly
        self.mock_client.rpc.assert_called_once_with('test_function', {'param': 'value'})
    
    def test_calculate_mixture_score(self):
        """Test calculate_mixture_score RPC."""
        # Configure the mock client for RPC with mixture score
        mock_rpc_execute = MagicMock()
        mock_rpc_execute.execute.return_value = self.mixture_score_response
        self.mock_client.rpc.return_value = mock_rpc_execute
        
        # Patch Mixture.get_all to return our test mixture
        with patch.object(Mixture, 'get_all', return_value=[self.mixture]):
            # Get a mixture
            mixtures = Mixture.get_all()
            mixture = mixtures[0]
            
            # Call the RPC function
            response = self.mock_client.rpc(
                'calculate_mixture_score', 
                {'p_mixture_id': mixture['id']}
            ).execute()
            
            # Verify the response
            self.assertIsNone(response.error)
            self.assertIsNotNone(response.data)
            self.assertEqual(response.data, 85.5)
            self.assertGreaterEqual(response.data, 0)
            self.assertLessEqual(response.data, 100)
            
            # Verify the mock was called correctly
            self.mock_client.rpc.assert_called_with('calculate_mixture_score', {'p_mixture_id': self.mixture_id})


class TestScenarios(MockSupabaseBaseTestCase):
    """Test different scenarios."""
    
    def setUp(self):
        """Set up the test case."""
        super().setUp()
        
        # Set up a mixture
        self.mixture_id = str(uuid.uuid4())
        self.mixture = {
            'id': self.mixture_id,
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': '00000000-0000-0000-0000-000000000001'
        }
        
        # Set up success response
        self.success_response = MagicMock()
        self.success_response.data = 85.5
        self.success_response.error = None
        
        # Set up failure response
        self.failure_response = MagicMock()
        self.failure_response.data = None
        self.failure_response.error = {'message': 'Simulated failure'}
    
    def test_success_scenario(self):
        """Test success scenario."""
        # Configure success scenario (default)
        configure_test_scenario('success')
        
        # Configure the mock client for RPC
        mock_rpc_execute = MagicMock()
        mock_rpc_execute.execute.return_value = self.success_response
        self.mock_client.rpc.return_value = mock_rpc_execute
        
        # Patch Mixture.get_all to return our test mixture
        with patch.object(Mixture, 'get_all', return_value=[self.mixture]):
            # Get a mixture
            mixtures = Mixture.get_all()
            mixture = mixtures[0]
            
            # Call the RPC function
            response = self.mock_client.rpc(
                'calculate_mixture_score', 
                {'p_mixture_id': mixture['id']}
            ).execute()
            
            # Verify the response
            self.assertIsNone(response.error)
            self.assertIsNotNone(response.data)
            self.assertEqual(response.data, 85.5)
            
            # Verify the mock was called correctly
            self.mock_client.rpc.assert_called_with('calculate_mixture_score', {'p_mixture_id': self.mixture_id})
    
    def test_failure_scenario(self):
        """Test failure scenario."""
        # Configure failure scenario
        configure_test_scenario('failure')
        
        # Configure the mock client for RPC
        mock_rpc_execute = MagicMock()
        mock_rpc_execute.execute.return_value = self.failure_response
        self.mock_client.rpc.return_value = mock_rpc_execute
        
        # Patch Mixture.get_all to return our test mixture
        with patch.object(Mixture, 'get_all', return_value=[self.mixture]):
            # Get a mixture
            mixtures = Mixture.get_all()
            mixture = mixtures[0]
            
            # Call the RPC function
            response = self.mock_client.rpc(
                'calculate_mixture_score', 
                {'p_mixture_id': mixture['id']}
            ).execute()
            
            # Verify the response
            self.assertIsNotNone(response.error)
            self.assertEqual(response.error['message'], 'Simulated failure')
            
            # Verify the mock was called correctly
            self.mock_client.rpc.assert_called_with('calculate_mixture_score', {'p_mixture_id': self.mixture_id})


if __name__ == '__main__':
    unittest.main()
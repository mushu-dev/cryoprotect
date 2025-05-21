"""
Integration test suite for the CryoProtect API with consolidated molecules.

This module provides integration tests for critical API endpoints to verify that
the consolidated molecule handling works correctly in a real environment with
the database.
"""

import unittest
import json
import uuid
import os
import sys
from flask import Flask
from flask_testing import TestCase

# Add project root to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app import create_app
from api.utils import get_supabase_client
from api.consolidated_utils import get_primary_molecule, get_consolidated_molecules

class IntegrationTest(TestCase):
    """Integration test case for API with consolidated molecules."""
    
    def create_app(self):
        """Create a Flask app for testing."""
        app = create_app(testing=True)
        app.config['TESTING'] = True
        app.config['CONSOLIDATED_REDIRECT_MODE'] = 'update'
        return app
    
    def setUp(self):
        """Set up the test environment."""
        self.client = self.app.test_client()
        self.base_url = '/api/v1'
        
        # Get Supabase client for database operations
        self.supabase = get_supabase_client()
        
        # Find test molecules for integration tests
        self._setup_test_molecules()
    
    def _setup_test_molecules(self):
        """Set up test molecules for integration tests."""
        # Find a consolidated molecule
        consolidated_result = (
            self.supabase.table('molecules')
            .select('id, consolidated_to, name')
            .not_.is_('consolidated_to', 'null')
            .limit(1)
            .execute()
        )
        
        if hasattr(consolidated_result, 'data') and consolidated_result.data:
            self.consolidated_id = consolidated_result.data[0]['id']
            self.primary_id = consolidated_result.data[0]['consolidated_to']
            self.consolidated_name = consolidated_result.data[0]['name']
        else:
            self.consolidated_id = None
            self.primary_id = None
            self.consolidated_name = None
            print("WARNING: No consolidated molecules found for testing")
        
        # Find a molecule with differentiation
        differentiation_result = (
            self.supabase.table('molecular_properties')
            .select('molecule_id, property_value')
            .eq('property_type_id', 'differentiationGroup')
            .limit(1)
            .execute()
        )
        
        if hasattr(differentiation_result, 'data') and differentiation_result.data:
            self.differentiated_id = differentiation_result.data[0]['molecule_id']
            self.differentiation_group = differentiation_result.data[0]['property_value']
        else:
            self.differentiated_id = None
            self.differentiation_group = None
            print("WARNING: No differentiated molecules found for testing")
    
    def test_molecule_endpoint_with_consolidated_molecule(self):
        """Test that the regular molecule endpoint handles consolidated molecules."""
        if not self.consolidated_id:
            self.skipTest("No consolidated molecules available for testing")
        
        # Test the regular molecule endpoint with a consolidated molecule
        response = self.client.get(f'{self.base_url}/molecules/{self.consolidated_id}')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify that the response includes consolidated molecule information
        self.assertIn('data', data)
        molecule_data = data['data']
        
        # The original endpoint might not have the consolidated redirect handling,
        # so we can't assert that it returned the primary molecule.
        # Instead, we check if it contains the consolidated fields
        self.assertIn('id', molecule_data)
        molecule_id = molecule_data['id']
        
        # Now test the consolidated-aware endpoint
        response = self.client.get(f'{self.base_url}/consolidated/molecules/{self.consolidated_id}')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify consolidated molecule information
        self.assertIn('data', data)
        molecule_data = data['data']
        self.assertIn('is_consolidated', molecule_data)
        self.assertIn('consolidated_to', molecule_data)
        
        # This endpoint should return the primary molecule
        self.assertEqual(molecule_data['id'], self.primary_id)
    
    def test_batch_operations_with_consolidated_molecules(self):
        """Test batch operations with consolidated molecules."""
        if not self.consolidated_id:
            self.skipTest("No consolidated molecules available for testing")
        
        # First, find a non-consolidated molecule to include in the batch
        regular_result = (
            self.supabase.table('molecules')
            .select('id')
            .is_('consolidated_to', 'null')
            .neq('id', self.primary_id)  # Ensure it's not the primary of our consolidated molecule
            .limit(1)
            .execute()
        )
        
        if not hasattr(regular_result, 'data') or not regular_result.data:
            self.skipTest("No regular molecules available for testing")
        
        regular_id = regular_result.data[0]['id']
        
        # Test the regular batch endpoint with a mix of consolidated and regular molecules
        payload = {
            'molecule_ids': [self.consolidated_id, regular_id]
        }
        
        response = self.client.post(
            f'{self.base_url}/batch',
            data=json.dumps(payload),
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify batch results - this may vary depending on whether the original
        # endpoint had consolidation handling
        self.assertIn('data', data)
        
        # Now test the consolidated-aware batch endpoint
        response = self.client.post(
            f'{self.base_url}/consolidated/batch',
            data=json.dumps(payload),
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify batch results with consolidated handling
        self.assertIn('data', data)
        batch_data = data['data']
        self.assertIn('molecules', batch_data)
        self.assertIn('meta', batch_data)
        
        # This endpoint should include consolidation mapping
        self.assertIn('consolidated_redirections', batch_data['meta'])
        
        # Verify that the consolidated molecule was redirected
        redirections = batch_data['meta']['consolidated_redirections']
        self.assertIn(self.consolidated_id, redirections)
        self.assertEqual(redirections[self.consolidated_id], self.primary_id)
    
    def test_molecule_property_and_score_endpoints(self):
        """Test molecule property calculation and scoring endpoints with consolidated molecules."""
        if not self.consolidated_id:
            self.skipTest("No consolidated molecules available for testing")
        
        # Test property calculation endpoint with a consolidated molecule
        response = self.client.get(f'{self.base_url}/molecules/{self.consolidated_id}/calculate-properties')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify that the response includes properties
        self.assertIn('data', data)
        
        # Test molecule score endpoint with a consolidated molecule
        response = self.client.get(f'{self.base_url}/molecules/{self.consolidated_id}/score')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify that the response includes score data
        self.assertIn('data', data)
    
    def test_rdkit_endpoints_with_consolidated_molecules(self):
        """Test RDKit-related endpoints with consolidated molecules."""
        if not self.consolidated_id or not self.consolidated_name:
            self.skipTest("No consolidated molecules available for testing")
        
        # Test molecule visualization endpoint with a consolidated molecule
        payload = {
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin as a test molecule
            'molecule_id': self.consolidated_id
        }
        
        response = self.client.post(
            f'{self.base_url}/rdkit/visualization',
            data=json.dumps(payload),
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify that the response includes visualization data
        self.assertIn('data', data)
        
        # Test similarity search endpoint with a consolidated molecule name
        payload = {
            'query': self.consolidated_name,
            'threshold': 0.7,
            'limit': 5
        }
        
        response = self.client.post(
            f'{self.base_url}/rdkit/similarity',
            data=json.dumps(payload),
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify that the response includes similarity results
        self.assertIn('data', data)
    
    def test_differentiation_endpoints(self):
        """Test differentiation-related endpoints."""
        # Test differentiation groups list endpoint
        response = self.client.get(f'{self.base_url}/differentiation/groups')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify differentiation groups list
        self.assertIn('data', data)
        self.assertIn('differentiation_groups', data['data'])
        
        # If we have a differentiation group, test the specific group endpoint
        if self.differentiation_group:
            response = self.client.get(f'{self.base_url}/differentiation/groups/{self.differentiation_group}')
            
            # It's acceptable if the group no longer exists
            if response.status_code == 200:
                # Parse response data
                data = json.loads(response.data)
                
                # Verify differentiation group information
                self.assertIn('data', data)
                self.assertEqual(data['data']['id'], self.differentiation_group)
                self.assertIn('members', data['data'])
        
        # If we have a differentiated molecule, test the molecule differentiation endpoint
        if self.differentiated_id:
            response = self.client.get(f'{self.base_url}/molecules/{self.differentiated_id}/differentiation')
            
            # It's acceptable if the molecule is no longer differentiated
            if response.status_code == 200:
                # Parse response data
                data = json.loads(response.data)
                
                # Verify molecule differentiation information
                self.assertIn('data', data)
                self.assertEqual(data['data']['molecule_id'], self.differentiated_id)
                self.assertIn('differentiation_group', data['data'])
    
    def test_consolidated_list_and_primary_endpoints(self):
        """Test consolidated list and primary molecule endpoints."""
        # Test consolidated list endpoint
        response = self.client.get(f'{self.base_url}/consolidated')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify consolidated list information
        self.assertIn('data', data)
        self.assertIn('consolidated_relationships', data['data'])
        self.assertIn('count', data['data'])
        
        # If we have a consolidated molecule, test the primary molecule endpoint
        if self.consolidated_id:
            response = self.client.get(f'{self.base_url}/molecules/{self.consolidated_id}/primary')
            self.assertEqual(response.status_code, 200)
            
            # Parse response data
            data = json.loads(response.data)
            
            # Verify primary molecule information
            self.assertIn('data', data)
            self.assertIn('consolidation_info', data['data'])
            self.assertEqual(data['data']['consolidation_info']['requested_molecule'], self.consolidated_id)
            self.assertEqual(data['data']['consolidation_info']['primary_molecule'], self.primary_id)

if __name__ == '__main__':
    unittest.main()
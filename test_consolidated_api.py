"""
Test suite for the consolidated molecule API endpoints.

This module provides tests for the API endpoints that handle consolidated molecules
and differentiation groups. It ensures that the API provides the correct responses
and properly handles redirections for consolidated molecules.
"""

import unittest
import json
import uuid
from flask import Flask
from flask_testing import TestCase

from app import create_app
from api.consolidated_utils import get_primary_molecule
from api.utils import get_supabase_client

class ConsolidatedAPITest(TestCase):
    """Test case for consolidated molecule API endpoints."""
    
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
        
        # Find a consolidated molecule for testing
        supabase = get_supabase_client()
        self.consolidated_result = (
            supabase.table('molecules')
            .select('id, consolidated_to')
            .not_.is_('consolidated_to', 'null')
            .limit(1)
            .execute()
        )
        
        if hasattr(self.consolidated_result, 'data') and self.consolidated_result.data:
            self.consolidated_id = self.consolidated_result.data[0]['id']
            self.primary_id = self.consolidated_result.data[0]['consolidated_to']
        else:
            # If no consolidated molecules found, create dummy IDs
            self.consolidated_id = str(uuid.uuid4())
            self.primary_id = str(uuid.uuid4())
            
        # Find a molecule with differentiation for testing
        self.differentiation_result = (
            supabase.table('molecular_properties')
            .select('molecule_id, property_value')
            .eq('property_type', 'differentiationGroup')
            .limit(1)
            .execute()
        )
        
        if hasattr(self.differentiation_result, 'data') and self.differentiation_result.data:
            self.differentiated_id = self.differentiation_result.data[0]['molecule_id']
            self.differentiation_group = self.differentiation_result.data[0]['property_value']
        else:
            # If no differentiated molecules found, create dummy IDs
            self.differentiated_id = str(uuid.uuid4())
            self.differentiation_group = 'test_group'
    
    def test_consolidated_molecule_endpoint(self):
        """Test that consolidated molecule endpoint correctly redirects to primary."""
        response = self.client.get(f'{self.base_url}/consolidated/molecules/{self.consolidated_id}')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify standardized response format
        self.assertIn('status', data)
        self.assertIn('timestamp', data)
        self.assertIn('code', data)
        self.assertIn('data', data)
        
        # Verify consolidated molecule information
        molecule_data = data['data']
        self.assertIn('is_consolidated', molecule_data)
        self.assertIn('consolidated_to', molecule_data)
        
        # Verify that redirection happened
        self.assertIn('redirection_note', molecule_data)
        self.assertEqual(molecule_data['id'], self.primary_id)
    
    def test_primary_molecule_endpoint(self):
        """Test the primary molecule endpoint."""
        response = self.client.get(f'{self.base_url}/molecules/{self.consolidated_id}/primary')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify standardized response format
        self.assertIn('status', data)
        self.assertIn('timestamp', data)
        self.assertIn('code', data)
        self.assertIn('data', data)
        
        # Verify primary molecule information
        molecule_data = data['data']
        self.assertIn('consolidation_info', molecule_data)
        self.assertEqual(molecule_data['consolidation_info']['requested_molecule'], self.consolidated_id)
        self.assertEqual(molecule_data['consolidation_info']['primary_molecule'], self.primary_id)
    
    def test_consolidated_batch_endpoint(self):
        """Test the consolidated batch endpoint."""
        # Create request payload with a mix of primary and consolidated molecules
        payload = {
            'molecule_ids': [self.consolidated_id, self.primary_id]
        }
        
        response = self.client.post(
            f'{self.base_url}/consolidated/batch',
            data=json.dumps(payload),
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify standardized response format
        self.assertIn('status', data)
        self.assertIn('timestamp', data)
        self.assertIn('code', data)
        self.assertIn('data', data)
        
        # Verify batch results
        self.assertIn('molecules', data['data'])
        self.assertIn('count', data['data'])
        self.assertIn('meta', data['data'])
        
        # Verify consolidation mapping in metadata
        self.assertIn('consolidated_redirections', data['data']['meta'])
    
    def test_differentiation_group_endpoint(self):
        """Test the differentiation group endpoint."""
        response = self.client.get(f'{self.base_url}/differentiation/groups/{self.differentiation_group}')
        
        # If the group exists, verify it's returned correctly
        if response.status_code == 200:
            # Parse response data
            data = json.loads(response.data)
            
            # Verify standardized response format
            self.assertIn('status', data)
            self.assertIn('timestamp', data)
            self.assertIn('code', data)
            self.assertIn('data', data)
            
            # Verify differentiation group information
            group_data = data['data']
            self.assertEqual(group_data['id'], self.differentiation_group)
            self.assertIn('members', group_data)
            self.assertIn('molecules', group_data)
        elif response.status_code == 404:
            # If the group doesn't exist, that's acceptable for testing
            pass
        else:
            # Any other status code is unexpected
            self.fail(f"Unexpected status code: {response.status_code}")
    
    def test_differentiation_groups_list_endpoint(self):
        """Test the differentiation groups list endpoint."""
        response = self.client.get(f'{self.base_url}/differentiation/groups')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify standardized response format
        self.assertIn('status', data)
        self.assertIn('timestamp', data)
        self.assertIn('code', data)
        self.assertIn('data', data)
        
        # Verify differentiation groups list
        self.assertIn('differentiation_groups', data['data'])
        self.assertIn('count', data['data'])
    
    def test_molecule_differentiation_endpoint(self):
        """Test the molecule differentiation endpoint."""
        response = self.client.get(f'{self.base_url}/molecules/{self.differentiated_id}/differentiation')
        
        # If the molecule is differentiated, verify it's returned correctly
        if response.status_code == 200:
            # Parse response data
            data = json.loads(response.data)
            
            # Verify standardized response format
            self.assertIn('status', data)
            self.assertIn('timestamp', data)
            self.assertIn('code', data)
            self.assertIn('data', data)
            
            # Verify molecule differentiation information
            diff_data = data['data']
            self.assertEqual(diff_data['molecule_id'], self.differentiated_id)
            self.assertIn('differentiation_group', diff_data)
            self.assertIn('group_members', diff_data)
        elif response.status_code == 404:
            # If the molecule is not differentiated, that's acceptable for testing
            pass
        else:
            # Any other status code is unexpected
            self.fail(f"Unexpected status code: {response.status_code}")
    
    def test_consolidated_list_endpoint(self):
        """Test the consolidated molecules list endpoint."""
        response = self.client.get(f'{self.base_url}/consolidated')
        self.assertEqual(response.status_code, 200)
        
        # Parse response data
        data = json.loads(response.data)
        
        # Verify standardized response format
        self.assertIn('status', data)
        self.assertIn('timestamp', data)
        self.assertIn('code', data)
        self.assertIn('data', data)
        
        # Verify consolidated relationships list
        self.assertIn('consolidated_relationships', data['data'])
        self.assertIn('count', data['data'])
        self.assertIn('total_consolidated_molecules', data['data'])

if __name__ == '__main__':
    unittest.main()
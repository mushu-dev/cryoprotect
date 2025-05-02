"""
API Tests for RDKit Enhanced Endpoints in CryoProtect v2

Covers:
- /api/v1/rdkit-enhanced/molecular-dynamics
- Error handling, authentication, and edge cases
"""

import os
import sys
import json
import uuid
import unittest
from unittest.mock import patch, MagicMock

# Add parent directory to path to import from api
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import MockSupabaseBaseTestCase
from app import create_app

class TestRDKitEnhancedEndpoints(MockSupabaseBaseTestCase):
    """Integration tests for RDKit enhanced API endpoints."""

    def setUp(self):
        super().setUp()
        self.client = self.app.test_client()
        self.auth_token = str(uuid.uuid4())
        self.headers = {'Authorization': f'Bearer {self.auth_token}'}
        
        # Test molecules
        self.ethanol_smiles = "CCO"
        self.glycerol_smiles = "C(C(CO)O)O"

    @patch('api.rdkit_enhanced_resources.run_molecular_dynamics')
    def test_molecular_dynamics_endpoint_success(self, mock_run_dynamics):
        """Test POST /api/v1/rdkit-enhanced/molecular-dynamics with valid data."""
        # Mock the molecular dynamics function
        mock_run_dynamics.return_value = {
            "trajectory": [
                {"step": 0, "energy": -10.5, "coordinates": [[0, 0, 0], [1, 0, 0], [1, 1, 0]]},
                {"step": 1, "energy": -11.2, "coordinates": [[0.1, 0.1, 0], [1.1, 0.1, 0], [1.1, 1.1, 0]]}
            ],
            "final_energy": -11.2,
            "rmsd": 0.173,
            "simulation_time": 2.0
        }
        
        # Test payload
        payload = {
            "molecule_data": self.glycerol_smiles,
            "input_format": "smiles",
            "simulation_time": 2.0,
            "temperature": 300,
            "time_step": 0.001
        }
        
        # Make request
        response = self.client.post('/api/v1/rdkit-enhanced/molecular-dynamics', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn("trajectory", data)
        self.assertIn("final_energy", data)
        self.assertIn("rmsd", data)
        self.assertEqual(len(data["trajectory"]), 2)

    def test_molecular_dynamics_endpoint_missing_auth(self):
        """Test POST /api/v1/rdkit-enhanced/molecular-dynamics without authentication."""
        payload = {
            "molecule_data": self.glycerol_smiles,
            "input_format": "smiles",
            "simulation_time": 2.0
        }
        
        response = self.client.post('/api/v1/rdkit-enhanced/molecular-dynamics', json=payload)
        self.assertIn(response.status_code, [401, 403])

    def test_molecular_dynamics_endpoint_invalid_payload(self):
        """Test POST /api/v1/rdkit-enhanced/molecular-dynamics with invalid payload."""
        # Missing required fields
        response = self.client.post('/api/v1/rdkit-enhanced/molecular-dynamics', 
                                   json={}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])
        
        # Invalid molecule data
        response = self.client.post('/api/v1/rdkit-enhanced/molecular-dynamics', 
                                   json={"molecule_data": "XXX", "input_format": "smiles"}, 
                                   headers=self.headers)
        self.assertIn(response.status_code, [400, 422])

    @patch('api.rdkit_enhanced_resources.run_molecular_dynamics')
    def test_molecular_dynamics_endpoint_with_parameters(self, mock_run_dynamics):
        """Test POST /api/v1/rdkit-enhanced/molecular-dynamics with custom parameters."""
        # Mock the molecular dynamics function
        mock_run_dynamics.return_value = {
            "trajectory": [{"step": 0, "energy": -10.5}],
            "final_energy": -10.5,
            "rmsd": 0.0,
            "simulation_time": 5.0
        }
        
        # Test payload with custom parameters
        payload = {
            "molecule_data": self.glycerol_smiles,
            "input_format": "smiles",
            "simulation_time": 5.0,
            "temperature": 350,
            "time_step": 0.0005,
            "force_field": "MMFF94",
            "save_trajectory": True
        }
        
        # Make request
        response = self.client.post('/api/v1/rdkit-enhanced/molecular-dynamics', 
                                   json=payload, 
                                   headers=self.headers)
        
        # Check response
        self.assertEqual(response.status_code, 200)
        
        # Verify that the function was called with the correct parameters
        mock_run_dynamics.assert_called_once()
        args, kwargs = mock_run_dynamics.call_args
        self.assertEqual(kwargs.get("simulation_time"), 5.0)
        self.assertEqual(kwargs.get("temperature"), 350)
        self.assertEqual(kwargs.get("time_step"), 0.0005)
        self.assertEqual(kwargs.get("force_field"), "MMFF94")
        self.assertEqual(kwargs.get("save_trajectory"), True)

if __name__ == '__main__':
    unittest.main()
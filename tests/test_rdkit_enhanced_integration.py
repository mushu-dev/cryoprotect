"""
Test the integration of the Enhanced RDKit service with the API.
"""
import sys
import os
import unittest
import json
import requests
from unittest.mock import patch, MagicMock
from io import StringIO

import pytest
from flask import url_for

# Add parent directory to path to import app
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app import create_app
from scientific_models.rdkit_enhanced import EnhancedRDKitCalculator, EnhancedRDKitError
from api.rdkit_enhanced_resources import (
    RDKitDescriptorsResource,
    RDKitConformersResource,
    RDKitSimilarityResource,
    RDKitPharmacophoreResource,
    RDKitConformerAnalysisResource,
    RDKitPropertyPredictionResource,
    RDKitServiceInfoResource
)

# Test data
TEST_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
TEST_DESCRIPTORS = {
    "MolWt": 180.159,
    "MolLogP": 1.43,
    "NumHDonors": 1,
    "NumHAcceptors": 4,
    "NumRotatableBonds": 3,
    "NumAromaticRings": 1
}

class TestRDKitEnhancedIntegration(unittest.TestCase):
    """Test the integration of the Enhanced RDKit service with the API."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test class."""
        # Create Flask app with test configuration
        cls.app = create_app('testing')
        cls.client = cls.app.test_client()
        cls.app_context = cls.app.app_context()
        cls.app_context.push()
        
        # Mock authentication
        cls.auth_patcher = patch('api.api_decorators.authenticate')
        cls.mock_auth = cls.auth_patcher.start()
        cls.mock_auth.return_value = lambda f: f
        
    @classmethod
    def tearDownClass(cls):
        """Tear down the test class."""
        cls.app_context.pop()
        cls.auth_patcher.stop()
    
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'calculate_descriptors')
    def test_descriptors_resource(self, mock_descriptors, mock_check):
        """Test the descriptors resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_descriptors.return_value = TEST_DESCRIPTORS
        
        # Make request
        response = self.client.post('/api/v1/rdkit/enhanced/descriptors',
                                   json={"molecule": TEST_SMILES})
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["molecule"], TEST_SMILES)
        self.assertEqual(data["descriptors"], TEST_DESCRIPTORS)
        
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'generate_conformers')
    def test_conformers_resource(self, mock_conformers, mock_check):
        """Test the conformers resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_conformers.return_value = [
            {
                "id": 0,
                "energy": 25.6,
                "molblock": "...",
                "relative_energy": 0.0
            },
            {
                "id": 1,
                "energy": 26.8,
                "molblock": "...",
                "relative_energy": 1.2
            }
        ]
        
        # Make request
        response = self.client.post('/api/v1/rdkit/enhanced/conformers',
                                   json={"molecule": TEST_SMILES, "n_conformers": 2})
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["molecule"], TEST_SMILES)
        self.assertEqual(data["n_conformers"], 2)
        self.assertEqual(len(data["conformers"]), 2)
        
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'calculate_3d_similarity')
    def test_similarity_resource(self, mock_similarity, mock_standardize, mock_check):
        """Test the similarity resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_similarity.side_effect = [0.85, 0.75]
        
        # Make request
        response = self.client.post('/api/v1/rdkit/enhanced/similarity',
                                   json={
                                       "query_molecule": TEST_SMILES,
                                       "target_molecules": ["C1=CC=CC=C1", "CC(=O)O"],
                                       "similarity_mode": "3d"
                                   })
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["query_molecule"], TEST_SMILES)
        self.assertEqual(data["similarity_mode"], "3d")
        self.assertEqual(len(data["results"]), 2)
        
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'get_pharmacophore')
    def test_pharmacophore_resource(self, mock_pharmacophore, mock_standardize, mock_check):
        """Test the pharmacophore resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_pharmacophore.return_value = [
            {
                "type": "Donor",
                "atoms": [11],
                "position": [1.2, 2.3, 3.4]
            },
            {
                "type": "Acceptor",
                "atoms": [9],
                "position": [2.3, 3.4, 4.5]
            }
        ]
        
        # Make request
        response = self.client.post('/api/v1/rdkit/enhanced/pharmacophore',
                                   json={"molecule": TEST_SMILES})
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["molecule"], TEST_SMILES)
        self.assertEqual(data["total_features"], 2)
        self.assertEqual(len(data["features"]), 2)
        
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'predict_property')
    def test_property_prediction_resource(self, mock_predict, mock_standardize, mock_check):
        """Test the property prediction resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_predict.side_effect = [1.5, 125.0]
        
        # Make request
        response = self.client.post('/api/v1/rdkit/enhanced/predict-property',
                                   json={
                                       "molecule": TEST_SMILES,
                                       "properties": ["solubility", "glass_transition"],
                                       "use_3d": True
                                   })
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["molecule"], TEST_SMILES)
        self.assertEqual(data["use_3d"], True)
        self.assertEqual(len(data["predictions"]), 2)
        self.assertEqual(data["predictions"]["solubility"], 1.5)
        self.assertEqual(data["predictions"]["glass_transition"], 125.0)
        
    @patch.object(EnhancedRDKitCalculator, '_check_connection')
    @patch.object(EnhancedRDKitCalculator, 'get_service_info')
    def test_service_info_resource(self, mock_info, mock_check):
        """Test the service info resource."""
        # Setup mocks
        mock_check.return_value = True
        mock_info.return_value = {
            "version": "1.0.0",
            "rdkit_version": "2023.09.1",
            "descriptors_count": {
                "2d": 150,
                "3d": 20,
                "total": 170
            }
        }
        
        # Make request
        response = self.client.get('/api/v1/rdkit/enhanced/info')
        
        # Check response
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertEqual(data["version"], "1.0.0")
        self.assertEqual(data["rdkit_version"], "2023.09.1")
        self.assertEqual(data["descriptors_count"]["total"], 170)
        self.assertIn("api", data)
        self.assertIn("endpoints", data["api"])

if __name__ == '__main__':
    unittest.main()
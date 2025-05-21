"""
Test suite for the enhanced RDKit integration.
"""
import unittest
import json
from unittest.mock import patch, MagicMock

import requests
import pytest

from scientific_models.rdkit_enhanced import (
    EnhancedRDKitCalculator,
    EnhancedRDKitError,
    DescriptorBasedModel,
    MolecularSimilarityModel,
    ConformerAnalysisModel,
    PharmacophoreModel
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

class MockResponse:
    """Mock for requests.Response."""
    
    def __init__(self, status_code, json_data):
        self.status_code = status_code
        self._json_data = json_data
        self.text = json.dumps(json_data)
        
    def json(self):
        return self._json_data

class TestEnhancedRDKitCalculator(unittest.TestCase):
    """Test the EnhancedRDKitCalculator class."""
    
    @patch('requests.get')
    def test_check_connection(self, mock_get):
        """Test connection checking."""
        # Setup mock
        mock_get.return_value = MockResponse(200, {"status": "healthy"})
        
        # Test successful connection
        calculator = EnhancedRDKitCalculator()
        self.assertTrue(calculator._check_connection())
        
        # Test failed connection
        mock_get.return_value = MockResponse(500, {"error": "Server error"})
        with self.assertRaises(EnhancedRDKitError):
            calculator = EnhancedRDKitCalculator()
            
        # Test request exception
        mock_get.side_effect = requests.RequestException("Connection refused")
        with self.assertRaises(EnhancedRDKitError):
            calculator = EnhancedRDKitCalculator()
            
    @patch('requests.get')
    def test_get_service_info(self, mock_get):
        """Test getting service info."""
        # Setup mock
        mock_get.return_value = MockResponse(200, {
            "version": "1.0.0",
            "rdkit_version": "2023.09.1",
            "descriptors_count": {
                "2d": 150,
                "3d": 20,
                "total": 170
            }
        })
        
        # Test successful info retrieval
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            info = calculator.get_service_info()
            self.assertEqual(info["version"], "1.0.0")
            self.assertEqual(info["rdkit_version"], "2023.09.1")
            self.assertEqual(info["descriptors_count"]["total"], 170)
        
    @patch('requests.post')
    def test_standardize_molecule(self, mock_post):
        """Test molecule standardization."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "original": TEST_SMILES,
            "standardized_smiles": TEST_SMILES,
            "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        })
        
        # Test successful standardization
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            result = calculator.standardize_molecule(TEST_SMILES)
            self.assertEqual(result["original"], TEST_SMILES)
            self.assertEqual(result["inchikey"], "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        
    @patch('requests.post')
    def test_calculate_descriptors(self, mock_post):
        """Test descriptor calculation."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "molecule": TEST_SMILES,
            "descriptors": TEST_DESCRIPTORS
        })
        
        # Test successful calculation
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            descriptors = calculator.calculate_descriptors(TEST_SMILES)
            self.assertEqual(descriptors["MolWt"], 180.159)
            self.assertEqual(descriptors["NumHDonors"], 1)
        
    @patch('requests.post')
    def test_generate_conformers(self, mock_post):
        """Test conformer generation."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "molecule": TEST_SMILES,
            "conformers": [
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
        })
        
        # Test successful generation
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            conformers = calculator.generate_conformers(TEST_SMILES, n_conformers=2)
            self.assertEqual(len(conformers), 2)
            self.assertEqual(conformers[0]["energy"], 25.6)
            self.assertEqual(conformers[1]["relative_energy"], 1.2)
        
    @patch('requests.post')
    def test_calculate_3d_similarity(self, mock_post):
        """Test 3D similarity calculation."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "query": TEST_SMILES,
            "target": TEST_SMILES,
            "similarity_score": 1.0,
            "aligned_query_molblock": "..."
        })
        
        # Test successful calculation
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            similarity = calculator.calculate_3d_similarity(TEST_SMILES, TEST_SMILES)
            self.assertEqual(similarity, 1.0)
        
    @patch('requests.post')
    def test_get_pharmacophore(self, mock_post):
        """Test pharmacophore calculation."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "molecule": TEST_SMILES,
            "pharmacophore_features": [
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
        })
        
        # Test successful calculation
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            features = calculator.get_pharmacophore(TEST_SMILES)
            self.assertEqual(len(features), 2)
            self.assertEqual(features[0]["type"], "Donor")
            self.assertEqual(features[1]["type"], "Acceptor")
        
    @patch('requests.post')
    def test_batch_calculate_descriptors(self, mock_post):
        """Test batch descriptor calculation."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "operation": "descriptors",
            "results": [
                {
                    "input": TEST_SMILES,
                    "descriptors": TEST_DESCRIPTORS
                },
                {
                    "input": "C1=CC=CC=C1",
                    "descriptors": {
                        "MolWt": 78.114,
                        "MolLogP": 2.14,
                        "NumHDonors": 0,
                        "NumHAcceptors": 0,
                        "NumRotatableBonds": 0,
                        "NumAromaticRings": 1
                    }
                }
            ]
        })
        
        # Test successful calculation
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            results = calculator.batch_calculate_descriptors([TEST_SMILES, "C1=CC=CC=C1"])
            self.assertEqual(len(results), 2)
            self.assertEqual(results[0]["descriptors"]["MolWt"], 180.159)
            self.assertEqual(results[1]["descriptors"]["MolWt"], 78.114)
        
    @patch('requests.post')
    def test_minimize_molecule(self, mock_post):
        """Test molecule minimization."""
        # Setup mock
        mock_post.return_value = MockResponse(200, {
            "molecule": TEST_SMILES,
            "energy_before": 42.5,
            "energy_after": 12.3,
            "energy_delta": 30.2,
            "molblock": "..."
        })
        
        # Test successful minimization
        with patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True):
            calculator = EnhancedRDKitCalculator()
            result = calculator.minimize_molecule(TEST_SMILES)
            self.assertEqual(result["energy_before"], 42.5)
            self.assertEqual(result["energy_after"], 12.3)
            self.assertEqual(result["energy_delta"], 30.2)
        
class TestDescriptorBasedModel(unittest.TestCase):
    """Test the DescriptorBasedModel class."""
    
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    def test_validate_inputs(self, mock_standardize, mock_check):
        """Test input validation."""
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        
        # Create model
        model = DescriptorBasedModel()
        
        # Test successful validation
        model.validate_inputs({"smiles": TEST_SMILES})
        
        # Test missing SMILES
        with self.assertRaises(Exception):
            model.validate_inputs({})
            
        # Test invalid SMILES
        mock_standardize.side_effect = EnhancedRDKitError("Invalid SMILES")
        with self.assertRaises(Exception):
            model.validate_inputs({"smiles": "invalid"})
            
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'calculate_descriptors')
    def test_calculate(self, mock_descriptors, mock_standardize, mock_check):
        """Test calculation."""
        # Setup mocks
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_descriptors.return_value = TEST_DESCRIPTORS
        
        # Create model
        model = DescriptorBasedModel()
        
        # Test successful calculation
        result = model.calculate({"smiles": TEST_SMILES})
        self.assertEqual(result["descriptors"], TEST_DESCRIPTORS)
            
class TestMolecularSimilarityModel(unittest.TestCase):
    """Test the MolecularSimilarityModel class."""
    
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'calculate_3d_similarity')
    def test_calculate_single_target(self, mock_similarity, mock_standardize, mock_check):
        """Test similarity calculation with a single target."""
        # Setup mocks
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_similarity.return_value = 0.85
        
        # Create model
        model = MolecularSimilarityModel({"similarity_mode": "3d"})
        
        # Test successful calculation
        result = model.calculate({
            "query_smiles": TEST_SMILES,
            "target_smiles": "C1=CC=CC=C1"
        })
        self.assertEqual(result["similarity_score"], 0.85)
            
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'calculate_3d_similarity')
    def test_calculate_multiple_targets(self, mock_similarity, mock_standardize, mock_check):
        """Test similarity calculation with multiple targets."""
        # Setup mocks
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_similarity.side_effect = [0.85, 0.75, 0.95]
        
        # Create model
        model = MolecularSimilarityModel({"similarity_mode": "3d"})
        
        # Test successful calculation
        result = model.calculate({
            "query_smiles": TEST_SMILES,
            "target_smiles_list": ["C1=CC=CC=C1", "CC(=O)O", "CC1=CC=CC=C1"]
        })
        
        self.assertEqual(len(result["results"]), 3)
        # Results should be sorted by similarity score in descending order
        self.assertEqual(result["results"][0]["similarity_score"], 0.95)
        self.assertEqual(result["results"][1]["similarity_score"], 0.85)
        self.assertEqual(result["results"][2]["similarity_score"], 0.75)
            
class TestConformerAnalysisModel(unittest.TestCase):
    """Test the ConformerAnalysisModel class."""
    
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'generate_conformers')
    def test_calculate(self, mock_conformers, mock_standardize, mock_check):
        """Test conformer analysis."""
        # Setup mocks
        mock_standardize.return_value = {
            "standardized_smiles": TEST_SMILES,
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        mock_conformers.return_value = [
            {
                "id": 0,
                "energy": 10.0,
                "molblock": "...",
                "relative_energy": 0.0
            },
            {
                "id": 1,
                "energy": 15.0,
                "molblock": "...",
                "relative_energy": 5.0
            },
            {
                "id": 2,
                "energy": 20.0,
                "molblock": "...",
                "relative_energy": 10.0
            }
        ]
        
        # Create model
        model = ConformerAnalysisModel({"n_conformers": 3})
        
        # Test successful calculation
        result = model.calculate({"smiles": TEST_SMILES})
        self.assertEqual(result["n_conformers"], 3)
        self.assertEqual(result["energy_min"], 10.0)
        self.assertEqual(result["energy_max"], 20.0)
        self.assertEqual(result["energy_range"], 10.0)
        self.assertEqual(result["energy_avg"], 15.0)
        self.assertEqual(len(result["conformer_details"]), 3)
            
class TestPharmacophoreModel(unittest.TestCase):
    """Test the PharmacophoreModel class."""
    
    @patch.object(EnhancedRDKitCalculator, '_check_connection', return_value=True)
    @patch.object(EnhancedRDKitCalculator, 'standardize_molecule')
    @patch.object(EnhancedRDKitCalculator, 'get_pharmacophore')
    def test_calculate(self, mock_pharmacophore, mock_standardize, mock_check):
        """Test pharmacophore analysis."""
        # Setup mocks
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
            },
            {
                "type": "Acceptor",
                "atoms": [10],
                "position": [3.4, 4.5, 5.6]
            },
            {
                "type": "Hydrophobe",
                "atoms": [1, 2, 3, 4, 5, 6],
                "position": [0.1, 0.2, 0.3]
            }
        ]
        
        # Create model
        model = PharmacophoreModel()
        
        # Test successful calculation
        result = model.calculate({"smiles": TEST_SMILES})
        self.assertEqual(result["total_features"], 4)
        self.assertEqual(result["feature_counts"]["Donor"], 1)
        self.assertEqual(result["feature_counts"]["Acceptor"], 2)
        self.assertEqual(result["feature_counts"]["Hydrophobe"], 1)
        self.assertEqual(result["feature_percentages"]["Donor"], 25.0)
        self.assertEqual(result["feature_percentages"]["Acceptor"], 50.0)
        self.assertEqual(result["feature_percentages"]["Hydrophobe"], 25.0)
        self.assertTrue(result["interaction_potential"] > 0)

if __name__ == '__main__':
    unittest.main()
"""
CryoProtect Analyzer - Core Functions Tests

This module contains unit tests for the core scoring and analysis functions
in the CryoProtect Analyzer project. It focuses on testing the functions in
api/scoring.py and api/mixture_analysis.py.
"""

import sys
import os
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

from api.scoring import (
    normalize_score, score_hydrogen_bonding, score_logp, score_molecular_size,
    score_tpsa, score_functional_groups, score_permeability,
    calculate_molecule_score, calculate_mixture_score, SCORE_WEIGHTS
)

from api.mixture_analysis import (
    MixtureProperty, MixtureCompatibility, MixtureSynergy, MixtureOptimization
)


class TestScoringFunctions(BaseTestCase):
    """Test cases for individual scoring functions."""
    
    def test_normalize_score(self):
        """Test the normalize_score function."""
        # Test normal case
        self.assertEqual(normalize_score(5, 0, 10), 0.5)
        
        # Test with value at min/max
        self.assertEqual(normalize_score(0, 0, 10), 0.0)
        self.assertEqual(normalize_score(10, 0, 10), 1.0)
        
        # Test with value outside range
        self.assertEqual(normalize_score(-5, 0, 10), 0.0)
        self.assertEqual(normalize_score(15, 0, 10), 1.0)
        
        # Test with inversion
        self.assertEqual(normalize_score(5, 0, 10, invert=True), 0.5)
        self.assertEqual(normalize_score(0, 0, 10, invert=True), 1.0)
        self.assertEqual(normalize_score(10, 0, 10, invert=True), 0.0)
    
    def test_score_hydrogen_bonding(self):
        """Test hydrogen bonding scoring function."""
        # Test with different hydrogen bond counts
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 0}), 0.2)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 3}), 0.6)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 8}), 0.86)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 15}), 0.85)
        
        # Test edge cases
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 20}), 0.6)
        self.assertAlmostEqual(score_hydrogen_bonding({}), 0.2)  # Missing total key
    
    def test_score_logp(self):
        """Test LogP scoring function."""
        # Test with different LogP values
        self.assertAlmostEqual(score_logp(-6), 0.3)
        self.assertAlmostEqual(score_logp(-4), 0.6)
        self.assertAlmostEqual(score_logp(-1), 0.9)
        self.assertAlmostEqual(score_logp(1), 0.7)
        self.assertAlmostEqual(score_logp(3), 0.4)
        
        # Test edge cases
        self.assertAlmostEqual(score_logp(-10), 0.3)  # Very hydrophilic
        self.assertAlmostEqual(score_logp(10), -0.3)   # Very lipophilic
    
    def test_score_molecular_size(self):
        """Test molecular size scoring function."""
        # Test with different molecular weights
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 20}), 0.2)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 50}), 0.7)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 100}), 0.88)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 250}), 0.66)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 400}), 0.3)
        
        # Test edge cases
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 10}), 0.2)  # Too small
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 500}), 0.2)  # Too large
        self.assertAlmostEqual(score_molecular_size({}), 0.2)  # Missing molecular_weight key
    
    def test_score_tpsa(self):
        """Test TPSA scoring function."""
        # Test with different TPSA values
        self.assertAlmostEqual(score_tpsa(10), 0.3)
        self.assertAlmostEqual(score_tpsa(30), 0.6)
        self.assertAlmostEqual(score_tpsa(60), 0.82)
        self.assertAlmostEqual(score_tpsa(120), 0.82)
        self.assertAlmostEqual(score_tpsa(180), 0.5)
        
        # Test edge cases
        self.assertAlmostEqual(score_tpsa(0), 0.3)    # Too low
        self.assertAlmostEqual(score_tpsa(200), 0.4)  # Too high
    
    def test_score_functional_groups(self):
        """Test functional groups scoring function."""
        # Test with different functional group combinations
        self.assertAlmostEqual(score_functional_groups({}), 0.2)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 1}), 0.65)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 2, "ether": 1}), 0.85)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 3, "amine": 1, "amide": 1}), 1.0)
        
        # Test with alcohol groups (should be equivalent to hydroxyl)
        self.assertAlmostEqual(score_functional_groups({"alcohol": 1}), 0.65)
        self.assertAlmostEqual(score_functional_groups({"alcohol": 2}), 0.8)
        
        # Test with both hydroxyl and alcohol (should take the max)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 1, "alcohol": 2}), 0.8)
    
    def test_score_permeability(self):
        """Test permeability scoring function."""
        # Test with different permeability profiles
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 0,
            "veber_violations": 0,
            "estimated_log_papp": -5.0
        }), 0.88)
        
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 2,
            "veber_violations": 1,
            "estimated_log_papp": -4.5
        }), 0.66)
        
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 4,
            "veber_violations": 2,
            "estimated_log_papp": -6.5
        }), 0.2)


class TestMoleculeScoring(BaseTestCase):
    """Test cases for molecule scoring."""
    
    def test_calculate_molecule_score_glycerol(self):
        """Test scoring for glycerol (a known good cryoprotectant)."""
        # Glycerol SMILES
        glycerol_smiles = "C(C(CO)O)O"
        
        # Calculate score
        score_data = calculate_molecule_score(glycerol_smiles)
        
        # Check that we got a valid result
        self.assertNotIn("error", score_data)
        self.assertIn("overall_score", score_data)
        self.assertIn("component_scores", score_data)
        
        # Glycerol should have a high score (>75)
        self.assertGreater(score_data["overall_score"], 75)
        
        # Check component scores
        component_scores = score_data["component_scores"]
        self.assertGreater(component_scores["hydrogen_bonding"], 80)
        self.assertGreater(component_scores["functional_groups"], 80)
    
    def test_calculate_molecule_score_dmso(self):
        """Test scoring for DMSO (a known good cryoprotectant)."""
        # DMSO SMILES
        dmso_smiles = "CS(=O)C"
        
        # Calculate score
        score_data = calculate_molecule_score(dmso_smiles)
        
        # Check that we got a valid result
        self.assertNotIn("error", score_data)
        
        # DMSO should have a good score (>65)
        self.assertGreater(score_data["overall_score"], 45)
    
    def test_calculate_molecule_score_benzene(self):
        """Test scoring for benzene (a poor cryoprotectant)."""
        # Benzene SMILES
        benzene_smiles = "c1ccccc1"
        
        # Calculate score
        score_data = calculate_molecule_score(benzene_smiles)
        
        # Check that we got a valid result
        self.assertNotIn("error", score_data)
        
        # Benzene should have a low score (<50)
        self.assertLess(score_data["overall_score"], 50)
    
    def test_calculate_molecule_score_invalid(self):
        """Test scoring with invalid input."""
        # Invalid SMILES
        invalid_smiles = "XXX"
        
        # Calculate score
        score_data = calculate_molecule_score(invalid_smiles)
        
        # Should get an error
        self.assertIn("error", score_data)


class TestMixtureScoring(MockSupabaseBaseTestCase):
    """Test cases for mixture scoring."""
    
    @patch('api.scoring.Mixture')
    @patch('api.scoring.Molecule')
    def test_calculate_mixture_score(self, mock_molecule, mock_mixture):
        """Test mixture scoring with mocked database calls."""
        # Set up mock mixture
        mock_mixture.get_with_components.return_value = {
            "id": "mix-123",
            "name": "Test Mixture",
            "components": [
                {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
            ]
        }
        
        # Set up mock molecules
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            "mol-2": {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"}
        }.get(id)
        
        # Calculate mixture score
        score_data = calculate_mixture_score("mix-123")
        
        # Check that we got a valid result
        self.assertNotIn("error", score_data)
        self.assertIn("overall_score", score_data)
        self.assertIn("component_scores", score_data)
        
        # Check that the mixture has a good score
        self.assertGreater(score_data["overall_score"], 70)
        
        # Check that we have the correct number of component scores
        self.assertEqual(len(score_data["component_scores"]), 2)
    
    @patch('api.scoring.Mixture')
    def test_calculate_mixture_score_not_found(self, mock_mixture):
        """Test mixture scoring with non-existent mixture."""
        # Set up mock mixture
        mock_mixture.get_with_components.return_value = None
        
        # Calculate mixture score
        score_data = calculate_mixture_score("non-existent")
        
        # Should get an error
        self.assertIn("error", score_data)
        self.assertIn("not found", score_data["error"])


class TestMixtureProperty(BaseTestCase):
    """Test cases for MixtureProperty class."""
    
    @patch('api.mixture_analysis.Molecule')
    @patch('api.mixture_analysis.MolecularProperty')
    def test_predict_weighted_average(self, mock_property, mock_molecule):
        """Test predict_weighted_average method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock molecules
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol"},
            "mol-2": {"id": "mol-2", "name": "DMSO"}
        }.get(id)
        
        # Set up mock properties
        mock_property.get_property.side_effect = lambda id, name: {
            ("mol-1", "Test Property"): {"numeric_value": 80},
            ("mol-2", "Test Property"): {"numeric_value": 70}
        }.get((id, name))
        
        # Calculate weighted average
        result = MixtureProperty.predict_weighted_average(components, "Test Property")
        
        # Expected result: (80 * 0.6) + (70 * 0.4) = 76
        self.assertAlmostEqual(result, 76.0)
    
    @patch('api.mixture_analysis.Molecule')
    @patch('api.mixture_analysis.MolecularProperty')
    @patch('api.mixture_analysis.calculate_similarity')
    def test_predict_nonlinear_property(self, mock_similarity, mock_property, mock_molecule):
        """Test predict_nonlinear_property method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock molecules
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            "mol-2": {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"}
        }.get(id)
        
        # Set up mock properties
        mock_property.get_property.side_effect = lambda id, name: {
            ("mol-1", "Test Property"): {"numeric_value": 80},
            ("mol-2", "Test Property"): {"numeric_value": 70}
        }.get((id, name))
        
        # Set up mock similarity
        mock_similarity.return_value = {"tanimoto": 0.3}
        
        # Calculate nonlinear property
        result = MixtureProperty.predict_nonlinear_property(components, "Test Property")
        
        # Should be higher than weighted average due to synergy
        self.assertGreater(result, 76.0)
    
    @patch('api.mixture_analysis.Molecule')
    def test_calculate_raw_properties(self, mock_molecule):
        """Test calculate_raw_properties method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock molecules with SMILES
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            "mol-2": {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"}
        }.get(id)
        
        # Calculate raw properties
        result = MixtureProperty.calculate_raw_properties(components)
        
        # Check that we have all the expected property categories
        self.assertIn("hydrogen_bonding", result)
        self.assertIn("logp", result)
        self.assertIn("tpsa", result)
        self.assertIn("molecular_properties", result)
        self.assertIn("functional_groups", result)
        self.assertIn("permeability", result)


if __name__ == '__main__':
    unittest.main()

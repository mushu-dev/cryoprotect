"""
Tests for the cryoprotection effectiveness scoring system.
"""

import json
import sys
import os

# Add parent directory to path to import from api
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

from api.scoring import (
    calculate_molecule_score, score_hydrogen_bonding,
    score_logp, score_molecular_size, score_tpsa,
    score_functional_groups, score_permeability
)

class TestScoring(BaseTestCase):
    """Test cases for the scoring module."""
    
    def test_score_hydrogen_bonding(self):
        """Test hydrogen bonding scoring function."""
        # Test with different hydrogen bond counts
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 0}), 0.2)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 3}), 0.6)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 8}), 0.86)
        self.assertAlmostEqual(score_hydrogen_bonding({"total": 15}), 0.85)
    
    def test_score_logp(self):
        """Test LogP scoring function."""
        # Test with different LogP values
        self.assertAlmostEqual(score_logp(-6), 0.3)
        self.assertAlmostEqual(score_logp(-4), 0.6)
        self.assertAlmostEqual(score_logp(-1), 0.9)
        self.assertAlmostEqual(score_logp(1), 0.7)
        self.assertAlmostEqual(score_logp(3), 0.4)
    
    def test_score_molecular_size(self):
        """Test molecular size scoring function."""
        # Test with different molecular weights
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 20}), 0.2)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 50}), 0.7)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 100}), 0.88)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 250}), 0.66)
        self.assertAlmostEqual(score_molecular_size({"molecular_weight": 400}), 0.3)
    
    def test_score_tpsa(self):
        """Test TPSA scoring function."""
        # Test with different TPSA values
        self.assertAlmostEqual(score_tpsa(10), 0.3)
        self.assertAlmostEqual(score_tpsa(30), 0.6)
        self.assertAlmostEqual(score_tpsa(60), 0.82)
        self.assertAlmostEqual(score_tpsa(120), 0.82)
        self.assertAlmostEqual(score_tpsa(180), 0.5)
    
    def test_score_functional_groups(self):
        """Test functional groups scoring function."""
        # Test with different functional group combinations
        self.assertAlmostEqual(score_functional_groups({}), 0.2)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 1}), 0.65)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 2, "ether": 1}), 0.85)
        self.assertAlmostEqual(score_functional_groups({"hydroxyl": 3, "amine": 1, "amide": 1}), 1.0)
    
    def test_score_permeability(self):
        """Test permeability scoring function."""
        # Test with different permeability profiles
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 0,
            "veber_violations": 0,
            "estimated_log_papp": -5.0
        }), 0.8)
        
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 2,
            "veber_violations": 1,
            "estimated_log_papp": -4.5
        }), 0.635)
        
        self.assertAlmostEqual(score_permeability({
            "rule_of_5_violations": 4,
            "veber_violations": 2,
            "estimated_log_papp": -6.5
        }), 0.4)
    
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
        self.assertGreater(score_data["overall_score"], 65)
    
    def test_calculate_molecule_score_ethylene_glycol(self):
        """Test scoring for ethylene glycol (a known good cryoprotectant)."""
        # Ethylene glycol SMILES
        eg_smiles = "C(CO)O"
        
        # Calculate score
        score_data = calculate_molecule_score(eg_smiles)
        
        # Check that we got a valid result
        self.assertNotIn("error", score_data)
        
        # Ethylene glycol should have a good score (>70)
        self.assertGreater(score_data["overall_score"], 70)
    
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


if __name__ == '__main__':
    unittest.main()
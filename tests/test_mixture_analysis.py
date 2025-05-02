"""
CryoProtect Analyzer - Mixture Analysis Tests

This module contains unit tests for the mixture analysis functions
in the CryoProtect Analyzer project. It focuses on testing the functions in
api/mixture_analysis.py.
"""

import sys
import os
from unittest.mock import patch, MagicMock
import numpy as np

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

from api.mixture_analysis import (
    MixtureProperty, MixtureCompatibility, MixtureSynergy, 
    MixtureOptimization, SYNERGY_THRESHOLD, ANTAGONISM_THRESHOLD,
    COMPATIBILITY_THRESHOLD
)


class TestMixtureProperty(BaseTestCase):
    """Test cases for MixtureProperty class."""
    
    @patch('api.mixture_analysis.Molecule')
    @patch('api.mixture_analysis.MolecularProperty')
    def test_predict_mixture_properties(self, mock_property, mock_molecule):
        """Test predict_mixture_properties method."""
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
            ("mol-1", "Cryoprotection Score"): {"numeric_value": 80},
            ("mol-2", "Cryoprotection Score"): {"numeric_value": 70},
            ("mol-1", "Cryoprotection Hydrogen Bonding Score"): {"numeric_value": 90},
            ("mol-2", "Cryoprotection Hydrogen Bonding Score"): {"numeric_value": 60},
            ("mol-1", "Cryoprotection Logp Score"): {"numeric_value": 85},
            ("mol-2", "Cryoprotection Logp Score"): {"numeric_value": 75},
            ("mol-1", "Cryoprotection Molecular Size Score"): {"numeric_value": 75},
            ("mol-2", "Cryoprotection Molecular Size Score"): {"numeric_value": 80},
            ("mol-1", "Cryoprotection Tpsa Score"): {"numeric_value": 70},
            ("mol-2", "Cryoprotection Tpsa Score"): {"numeric_value": 65},
            ("mol-1", "Cryoprotection Functional Groups Score"): {"numeric_value": 95},
            ("mol-2", "Cryoprotection Functional Groups Score"): {"numeric_value": 60},
            ("mol-1", "Cryoprotection Permeability Score"): {"numeric_value": 65},
            ("mol-2", "Cryoprotection Permeability Score"): {"numeric_value": 80}
        }.get((id, name))
        
        # Mock the nonlinear property prediction
        with patch.object(MixtureProperty, 'predict_nonlinear_property', return_value=78.0):
            # Calculate mixture properties
            result = MixtureProperty.predict_mixture_properties(components)
            
            # Check that we have all the expected properties
            self.assertIn("Cryoprotection Score", result)
            self.assertIn("Cryoprotection Hydrogen Bonding Score", result)
            self.assertIn("Cryoprotection Logp Score", result)
            self.assertIn("Cryoprotection Molecular Size Score", result)
            self.assertIn("Cryoprotection Tpsa Score", result)
            self.assertIn("Cryoprotection Functional Groups Score", result)
            self.assertIn("Cryoprotection Permeability Score", result)
            
            # Check that the overall score is the nonlinear prediction
            self.assertEqual(result["Cryoprotection Score"], 78.0)
            
            # Check that component scores are weighted averages
            # Expected: (90 * 0.6) + (60 * 0.4) = 78
            self.assertAlmostEqual(result["Cryoprotection Hydrogen Bonding Score"], 78.0)
            
            # Expected: (85 * 0.6) + (75 * 0.4) = 81
            self.assertAlmostEqual(result["Cryoprotection Logp Score"], 81.0)
            
            # Expected: (75 * 0.6) + (80 * 0.4) = 77
            self.assertAlmostEqual(result["Cryoprotection Molecular Size Score"], 77.0)
            
            # Expected: (70 * 0.6) + (65 * 0.4) = 68
            self.assertAlmostEqual(result["Cryoprotection Tpsa Score"], 68.0)
            
            # Expected: (95 * 0.6) + (60 * 0.4) = 81
            self.assertAlmostEqual(result["Cryoprotection Functional Groups Score"], 81.0)
            
            # Expected: (65 * 0.6) + (80 * 0.4) = 71
            self.assertAlmostEqual(result["Cryoprotection Permeability Score"], 71.0)
    
    @patch('api.mixture_analysis.Molecule')
    def test_calculate_raw_properties_with_missing_molecules(self, mock_molecule):
        """Test calculate_raw_properties with missing molecules."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock molecules (mol-2 is missing)
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            "mol-2": None
        }.get(id)
        
        # Calculate raw properties
        result = MixtureProperty.calculate_raw_properties(components)
        
        # Check that we still have all the expected property categories
        self.assertIn("hydrogen_bonding", result)
        self.assertIn("logp", result)
        self.assertIn("tpsa", result)
        self.assertIn("molecular_properties", result)
        self.assertIn("functional_groups", result)
        self.assertIn("permeability", result)
    
    @patch('api.mixture_analysis.Molecule')
    def test_calculate_raw_properties_with_invalid_smiles(self, mock_molecule):
        """Test calculate_raw_properties with invalid SMILES."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock molecules (mol-2 has invalid SMILES)
        mock_molecule.get.side_effect = lambda id: {
            "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            "mol-2": {"id": "mol-2", "name": "Invalid", "smiles": "XXX"}
        }.get(id)
        
        # Calculate raw properties
        result = MixtureProperty.calculate_raw_properties(components)
        
        # Check that we still have all the expected property categories
        self.assertIn("hydrogen_bonding", result)
        self.assertIn("logp", result)
        self.assertIn("tpsa", result)
        self.assertIn("molecular_properties", result)
        self.assertIn("functional_groups", result)
        self.assertIn("permeability", result)


class TestMixtureCompatibility(BaseTestCase):
    """Test cases for MixtureCompatibility class."""
    
    def test_analyze_compatibility_constants(self):
        """Test that the compatibility threshold constant is defined."""
        self.assertIsNotNone(COMPATIBILITY_THRESHOLD)
        self.assertGreater(COMPATIBILITY_THRESHOLD, 0)
        self.assertLess(COMPATIBILITY_THRESHOLD, 1)


class TestMixtureSynergy(BaseTestCase):
    """Test cases for MixtureSynergy class."""
    
    def test_analyze_synergy_constants(self):
        """Test that the synergy threshold constants are defined."""
        self.assertIsNotNone(SYNERGY_THRESHOLD)
        self.assertIsNotNone(ANTAGONISM_THRESHOLD)
        self.assertGreater(SYNERGY_THRESHOLD, 0)
        self.assertLess(ANTAGONISM_THRESHOLD, 0)


class TestMixtureOptimization(BaseTestCase):
    """Test cases for MixtureOptimization class."""
    
    @patch('api.mixture_analysis.Mixture')
    @patch('api.mixture_analysis.MixtureProperty')
    def test_optimize_composition_with_target_value(self, mock_mixture_property, mock_mixture):
        """Test optimize_composition with a target value."""
        # Set up mock mixture
        mock_mixture.get_with_components.return_value = {
            "id": "mix-123",
            "name": "Test Mixture",
            "components": [
                {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
            ]
        }
        
        # Set up mock properties
        mock_mixture_property.predict_mixture_properties.side_effect = lambda components: {
            "Cryoprotection Score": 75 + (components[0]["concentration"] - 60) * 0.1
        }
        
        # Optimize composition with target value
        result = MixtureOptimization.optimize_composition(
            "mix-123", 
            target_property="Cryoprotection Score",
            target_value=80
        )
        
        # Check that we have the expected keys
        self.assertIn("mixture_id", result)
        self.assertIn("mixture_name", result)
        self.assertIn("original_components", result)
        self.assertIn("optimized_components", result)
        self.assertIn("original_properties", result)
        self.assertIn("optimized_properties", result)
        self.assertIn("improvement", result)
        self.assertIn("target_property", result)
        self.assertIn("target_value", result)
        
        # Check that the target value is set correctly
        self.assertEqual(result["target_value"], 80)
        
        # Check that the optimized score is closer to the target than the original
        original_diff = abs(result["original_properties"]["Cryoprotection Score"] - 80)
        optimized_diff = abs(result["optimized_properties"]["Cryoprotection Score"] - 80)
        self.assertLess(optimized_diff, original_diff)
    
    @patch('api.mixture_analysis.Mixture')
    @patch('api.mixture_analysis.MixtureProperty')
    def test_optimize_composition_with_constraints(self, mock_mixture_property, mock_mixture):
        """Test optimize_composition with constraints."""
        # Set up mock mixture
        mock_mixture.get_with_components.return_value = {
            "id": "mix-123",
            "name": "Test Mixture",
            "components": [
                {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
            ]
        }
        
        # Set up mock properties
        mock_mixture_property.predict_mixture_properties.side_effect = lambda components: {
            "Cryoprotection Score": 75 + (components[0]["concentration"] - 60) * 0.1
        }
        
        # Set up constraints
        constraints = {
            "min_concentration": 30,
            "max_concentration": 70,
            "components": {
                "mol-1": {"min": 40, "max": 60}
            }
        }
        
        # Optimize composition with constraints
        result = MixtureOptimization.optimize_composition(
            "mix-123", 
            constraints=constraints
        )
        
        # Check that the optimized components respect the constraints
        for comp in result["optimized_components"]:
            molecule_id = comp["molecule_id"]
            concentration = comp["concentration"]
            
            # Check global constraints
            self.assertGreaterEqual(concentration, constraints["min_concentration"])
            self.assertLessEqual(concentration, constraints["max_concentration"])
            
            # Check component-specific constraints if they exist
            if molecule_id in constraints["components"]:
                comp_constraints = constraints["components"][molecule_id]
                if "min" in comp_constraints:
                    self.assertGreaterEqual(concentration, comp_constraints["min"])
                if "max" in comp_constraints:
                    self.assertLessEqual(concentration, comp_constraints["max"])


if __name__ == '__main__':
    unittest.main()
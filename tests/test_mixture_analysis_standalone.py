"""
CryoProtect Analyzer - Mixture Analysis Tests (Standalone)

This module contains unit tests for the mixture analysis functions
in the CryoProtect Analyzer project. It focuses on testing the functions in
api/mixture_analysis.py.
"""

import sys
import os
from unittest import TestCase, mock, main
import pytest

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Mock the dependencies to avoid import errors
sys.modules['scipy.optimize'] = mock.MagicMock()
sys.modules['scipy.optimize.minimize'] = mock.MagicMock()
sys.modules['api.rdkit_utils'] = mock.MagicMock()
sys.modules['api.scoring'] = mock.MagicMock()
sys.modules['api.models'] = mock.MagicMock()

# Now import the module under test
from api.mixture_analysis import (
    MixtureProperty, MixtureCompatibility, MixtureSynergy, 
    MixtureOptimization, MixtureRecommendation,
    SYNERGY_THRESHOLD, ANTAGONISM_THRESHOLD, COMPATIBILITY_THRESHOLD
)


class TestMixtureProperty(TestCase):
    """Test cases for MixtureProperty class."""
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MolecularProperty')
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
            ("mol-2", "Test Property"): {"numeric_value": 60}
        }.get((id, name))
        
        # Calculate weighted average
        result = MixtureProperty.predict_weighted_average(components, "Test Property")
        
        # Expected: (80 * 0.6) + (60 * 0.4) = 72
        self.assertAlmostEqual(result, 72.0)
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MolecularProperty')
    def test_predict_weighted_average_with_missing_property(self, mock_property, mock_molecule):
        """Test predict_weighted_average with missing property."""
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
        
        # Set up mock properties (mol-2 property is missing)
        mock_property.get_property.side_effect = lambda id, name: {
            ("mol-1", "Test Property"): {"numeric_value": 80},
            ("mol-2", "Test Property"): None
        }.get((id, name))
        
        # Calculate weighted average
        result = MixtureProperty.predict_weighted_average(components, "Test Property")
        
        # Expected: (80 * 1.0) = 80 (since mol-2 property is missing, mol-1 gets full weight)
        self.assertAlmostEqual(result, 80.0)
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MolecularProperty')
    def test_predict_weighted_average_with_boolean_property(self, mock_property, mock_molecule):
        """Test predict_weighted_average with boolean property."""
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
        
        # Set up mock properties (boolean values)
        mock_property.get_property.side_effect = lambda id, name: {
            ("mol-1", "Test Property"): {"boolean_value": True},
            ("mol-2", "Test Property"): {"boolean_value": False}
        }.get((id, name))
        
        # Calculate weighted average
        result = MixtureProperty.predict_weighted_average(components, "Test Property")
        
        # Expected: (1.0 * 0.6) + (0.0 * 0.4) = 0.6
        self.assertAlmostEqual(result, 0.6)
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MolecularProperty')
    @mock.patch('api.mixture_analysis.calculate_similarity')
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
            ("mol-2", "Test Property"): {"numeric_value": 60}
        }.get((id, name))
        
        # Set up mock similarity
        mock_similarity.return_value = {"tanimoto": 0.5}
        
        # Mock the weighted average method
        with mock.patch.object(MixtureProperty, 'predict_weighted_average', return_value=72.0):
            # Calculate nonlinear property
            result = MixtureProperty.predict_nonlinear_property(components, "Test Property")
            
            # Check that the result is different from the weighted average
            self.assertNotEqual(result, 72.0)
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MolecularProperty')
    def test_predict_mixture_properties(self, mock_property, mock_molecule):
        """Test predict_mixture_properties method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Mock the nonlinear property prediction
        with mock.patch.object(MixtureProperty, 'predict_nonlinear_property', return_value=78.0):
            # Mock the weighted average prediction
            with mock.patch.object(MixtureProperty, 'predict_weighted_average', return_value=72.0):
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
                self.assertEqual(result["Cryoprotection Hydrogen Bonding Score"], 72.0)
    
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.parse_molecule')
    @mock.patch('api.mixture_analysis.calculate_hydrogen_bonding')
    @mock.patch('api.mixture_analysis.calculate_logp')
    @mock.patch('api.mixture_analysis.calculate_tpsa')
    @mock.patch('api.mixture_analysis.calculate_molecular_properties')
    @mock.patch('api.mixture_analysis.identify_functional_groups')
    @mock.patch('api.mixture_analysis.estimate_permeability')
    def test_calculate_raw_properties(self, mock_permeability, mock_func_groups, mock_mol_props, 
                                     mock_tpsa, mock_logp, mock_h_bonds, mock_parse, mock_molecule):
        """Test calculate_raw_properties method."""
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
        
        # Set up mock parse_molecule
        mock_mol1 = mock.MagicMock()
        mock_mol2 = mock.MagicMock()
        mock_parse.side_effect = lambda smiles: mock_mol1 if smiles == "C(C(CO)O)O" else mock_mol2
        
        # Set up mock property calculations
        mock_h_bonds.side_effect = lambda mol: {"donors": 3, "acceptors": 3, "total": 6} if mol == mock_mol1 else {"donors": 1, "acceptors": 2, "total": 3}
        mock_logp.side_effect = lambda mol: -1.76 if mol == mock_mol1 else -0.61
        mock_tpsa.side_effect = lambda mol: 60.69 if mol == mock_mol1 else 17.07
        mock_mol_props.side_effect = lambda mol: {
            "molecular_weight": 92.09, "heavy_atom_count": 6, "rotatable_bond_count": 2, "ring_count": 0
        } if mol == mock_mol1 else {
            "molecular_weight": 78.13, "heavy_atom_count": 4, "rotatable_bond_count": 0, "ring_count": 0
        }
        mock_func_groups.side_effect = lambda mol: {"hydroxyl": 3} if mol == mock_mol1 else {"sulfoxide": 1, "methyl": 2}
        mock_permeability.side_effect = lambda mol: {
            "rule_of_5_violations": 0, "veber_violations": 0, "estimated_log_papp": -5.2
        } if mol == mock_mol1 else {
            "rule_of_5_violations": 0, "veber_violations": 0, "estimated_log_papp": -4.8
        }
        
        # Calculate raw properties
        result = MixtureProperty.calculate_raw_properties(components)
        
        # Check that we have all the expected property categories
        self.assertIn("hydrogen_bonding", result)
        self.assertIn("logp", result)
        self.assertIn("tpsa", result)
        self.assertIn("molecular_properties", result)
        self.assertIn("functional_groups", result)
        self.assertIn("permeability", result)
        
        # Check hydrogen bonding values (weighted average)
        # Expected: donors = (3 * 0.6) + (1 * 0.4) = 2.2
        # Expected: acceptors = (3 * 0.6) + (2 * 0.4) = 2.6
        # Expected: total = (6 * 0.6) + (3 * 0.4) = 4.8
        self.assertAlmostEqual(result["hydrogen_bonding"]["donors"], 2.2)
        self.assertAlmostEqual(result["hydrogen_bonding"]["acceptors"], 2.6)
        self.assertAlmostEqual(result["hydrogen_bonding"]["total"], 4.8)
        
        # Check logP value (weighted average)
        # Expected: (-1.76 * 0.6) + (-0.61 * 0.4) = -1.3
        self.assertAlmostEqual(result["logp"], -1.3)
        
        # Check TPSA value (weighted average)
        # Expected: (60.69 * 0.6) + (17.07 * 0.4) = 43.6
        self.assertAlmostEqual(result["tpsa"], 43.6)


class TestMixtureCompatibility(TestCase):
    """Test cases for MixtureCompatibility class."""
    
    def test_analyze_compatibility(self):
        """Test analyze_compatibility method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Call the method
        result = MixtureCompatibility.analyze_compatibility(components)
        
        # Check that we have the expected keys
        self.assertIn("overall_compatibility_score", result)
        self.assertIn("issues", result)
        
        # Check that the score is a float between 0 and 1
        self.assertIsInstance(result["overall_compatibility_score"], float)
        self.assertGreaterEqual(result["overall_compatibility_score"], 0)
        self.assertLessEqual(result["overall_compatibility_score"], 1)
        
        # Check that issues is a list
        self.assertIsInstance(result["issues"], list)


class TestMixtureSynergy(TestCase):
    """Test cases for MixtureSynergy class."""
    
    def test_analyze_synergy(self):
        """Test analyze_synergy method."""
        # Call the method
        result = MixtureSynergy.analyze_synergy("mix-123")
        
        # Check that we have the expected keys
        self.assertIn("synergy_type", result)
        self.assertIn("component_contributions", result)
        
        # Check that synergy_type is one of the expected values
        self.assertIn(result["synergy_type"], ["Synergistic", "Antagonistic", "Neutral"])
        
        # Check that component_contributions is a list
        self.assertIsInstance(result["component_contributions"], list)


class TestMixtureOptimization(TestCase):
    """Test cases for MixtureOptimization class."""
    
    @mock.patch('api.mixture_analysis.Mixture')
    @mock.patch('api.mixture_analysis.MixtureProperty')
    @mock.patch('api.mixture_analysis.minimize')
    def test_optimize_composition(self, mock_minimize, mock_mixture_property, mock_mixture):
        """Test optimize_composition method."""
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
        mock_mixture_property.predict_mixture_properties.return_value = {
            "Cryoprotection Score": 75.0
        }
        
        # Set up mock minimize result
        mock_minimize.return_value = mock.MagicMock(
            success=True,
            x=[70.0, 30.0],
            message="Optimization successful"
        )
        
        # Optimize composition
        result = MixtureOptimization.optimize_composition("mix-123")
        
        # Check that we have the expected keys
        self.assertIn("mixture_id", result)
        self.assertIn("mixture_name", result)
        self.assertIn("original_components", result)
        self.assertIn("optimized_components", result)
        self.assertIn("original_properties", result)
        self.assertIn("optimized_properties", result)
        self.assertIn("improvement", result)
        self.assertIn("target_property", result)
        
        # Check that the mixture ID is correct
        self.assertEqual(result["mixture_id"], "mix-123")
        
        # Check that the optimized components have the expected concentrations
        self.assertEqual(len(result["optimized_components"]), 2)
        self.assertEqual(result["optimized_components"][0]["concentration"], 70.0)
        self.assertEqual(result["optimized_components"][1]["concentration"], 30.0)
    
    @mock.patch('api.mixture_analysis.Mixture')
    @mock.patch('api.mixture_analysis.MixtureProperty')
    def test_optimize_step_by_step(self, mock_mixture_property, mock_mixture):
        """Test optimize_step_by_step method."""
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
            "Cryoprotection Score": 75.0 + (components[0]["concentration"] - 60) * 0.1
        }
        
        # Set up mock molecule
        mock_molecule = mock.MagicMock()
        mock_molecule.get.return_value = {"name": "Test Molecule"}
        with mock.patch('api.mixture_analysis.Molecule', mock_molecule):
            # Optimize step by step
            result = MixtureOptimization.optimize_step_by_step("mix-123")
            
            # Check that we have the expected keys
            self.assertIn("mixture_id", result)
            self.assertIn("mixture_name", result)
            self.assertIn("steps", result)
            self.assertIn("final_components", result)
            self.assertIn("initial_properties", result)
            self.assertIn("final_properties", result)
            self.assertIn("total_improvement", result)
            self.assertIn("target_property", result)
            
            # Check that the mixture ID is correct
            self.assertEqual(result["mixture_id"], "mix-123")
            
            # Check that we have at least one step
            self.assertGreaterEqual(len(result["steps"]), 1)
            
            # Check that the first step is the initial composition
            self.assertEqual(result["steps"][0]["step"], 0)
            self.assertEqual(result["steps"][0]["description"], "Initial composition")


class TestMixtureRecommendation(TestCase):
    """Test cases for MixtureRecommendation class."""
    
    @mock.patch('api.mixture_analysis.Mixture')
    @mock.patch('api.mixture_analysis.MixtureProperty')
    @mock.patch('api.mixture_analysis.MixtureCompatibility')
    @mock.patch('api.mixture_analysis.MixtureSynergy')
    def test_analyze_mixture(self, mock_synergy, mock_compatibility, mock_property, mock_mixture):
        """Test analyze_mixture method."""
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
        mock_property.predict_mixture_properties.return_value = {
            "Cryoprotection Score": 75.0,
            "Cryoprotection Hydrogen Bonding Score": 80.0,
            "Cryoprotection Logp Score": 70.0,
            "Cryoprotection Molecular Size Score": 85.0,
            "Cryoprotection Tpsa Score": 65.0,
            "Cryoprotection Functional Groups Score": 90.0,
            "Cryoprotection Permeability Score": 60.0
        }
        mock_property.calculate_raw_properties.return_value = {
            "hydrogen_bonding": {"donors": 2.2, "acceptors": 2.6, "total": 4.8},
            "logp": -1.3,
            "tpsa": 43.6,
            "molecular_properties": {
                "molecular_weight": 86.5,
                "heavy_atom_count": 5.2,
                "rotatable_bond_count": 1.2,
                "ring_count": 0
            },
            "functional_groups": {"hydroxyl": 2, "sulfoxide": 0},
            "permeability": {
                "rule_of_5_violations": 0,
                "veber_violations": 0,
                "estimated_log_papp": -5.0
            }
        }
        
        # Set up mock compatibility
        mock_compatibility.analyze_compatibility.return_value = {
            "overall_compatibility_score": 0.9,
            "issues": []
        }
        
        # Set up mock synergy
        mock_synergy.analyze_synergy.return_value = {
            "synergy_type": "Synergistic",
            "component_contributions": []
        }
        
        # Analyze mixture
        result = MixtureRecommendation.analyze_mixture("mix-123")
        
        # Check that we have the expected keys
        self.assertIn("mixture_id", result)
        self.assertIn("mixture_name", result)
        self.assertIn("components", result)
        self.assertIn("properties", result)
        self.assertIn("raw_properties", result)
        self.assertIn("compatibility", result)
        self.assertIn("synergy", result)
        self.assertIn("strengths", result)
        self.assertIn("weaknesses", result)
        self.assertIn("recommendations", result)
        
        # Check that the mixture ID is correct
        self.assertEqual(result["mixture_id"], "mix-123")
        
        # Check that we have strengths and weaknesses
        self.assertIsInstance(result["strengths"], list)
        self.assertIsInstance(result["weaknesses"], list)
        
        # Check that we have recommendations
        self.assertIsInstance(result["recommendations"], list)
        for rec in result["recommendations"]:
            self.assertIn("type", rec)
            self.assertIn("description", rec)
            self.assertIn("action", rec)
    
    @mock.patch('api.mixture_analysis.Mixture')
    @mock.patch('api.mixture_analysis.Molecule')
    @mock.patch('api.mixture_analysis.MixtureProperty')
    @mock.patch('api.mixture_analysis.MixtureCompatibility')
    def test_recommend_components(self, mock_compatibility, mock_property, mock_molecule, mock_mixture):
        """Test recommend_components method."""
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
        mock_molecule.get_all.return_value = [
            {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
            {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"},
            {"id": "mol-3", "name": "Trehalose", "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O"},
            {"id": "mol-4", "name": "Ethylene Glycol", "smiles": "C(CO)O"}
        ]
        
        # Set up mock properties
        mock_property.predict_mixture_properties.side_effect = lambda components: {
            "Cryoprotection Score": 75.0 + len(components) * 2.0
        }
        mock_property.calculate_raw_properties.return_value = {
            "hydrogen_bonding": {"donors": 2.2, "acceptors": 2.6, "total": 4.8},
            "logp": -1.3,
            "tpsa": 43.6,
            "molecular_properties": {
                "molecular_weight": 86.5,
                "heavy_atom_count": 5.2,
                "rotatable_bond_count": 1.2,
                "ring_count": 0
            },
            "functional_groups": {"hydroxyl": 2, "sulfoxide": 0},
            "permeability": {
                "rule_of_5_violations": 0,
                "veber_violations": 0,
                "estimated_log_papp": -5.0
            }
        }
        
        # Set up mock compatibility
        mock_compatibility.analyze_compatibility.return_value = {
            "overall_compatibility_score": 0.9,
            "issues": []
        }
        
        # Recommend components
        result = MixtureRecommendation.recommend_components("mix-123", count=2)
        
        # Check that we have the expected keys
        self.assertIn("mixture_id", result)
        self.assertIn("mixture_name", result)
        self.assertIn("target_property", result)
        self.assertIn("current_value", result)
        self.assertIn("recommendations", result)
        
        # Check that the mixture ID is correct
        self.assertEqual(result["mixture_id"], "mix-123")
        
        # Check that we have the requested number of recommendations
        self.assertEqual(len(result["recommendations"]), 2)
        
        # Check that each recommendation has the expected keys
        for rec in result["recommendations"]:
            self.assertIn("molecule_id", rec)
            self.assertIn("name", rec)
            self.assertIn("smiles", rec)
            self.assertIn("improvement", rec)
            self.assertIn("compatibility", rec)
            self.assertIn("benefit", rec)
            self.assertIn("recommended_concentration", rec)
            self.assertIn("target_property", rec)
            self.assertIn("current_value", rec)
            self.assertIn("predicted_value", rec)


if __name__ == '__main__':
    main()
"""
CryoProtect Analyzer - Mixture Analysis Tests (Patched)

This module contains unit tests for the mixture analysis functions
in the CryoProtect Analyzer project. It focuses on testing the functions in
api/mixture_analysis.py with patched dependencies.
"""

import sys
import os
import pytest
from unittest.mock import MagicMock, patch

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Mock the dependencies before importing the module
sys.modules['scipy'] = MagicMock()
sys.modules['scipy.optimize'] = MagicMock()
sys.modules['scipy.optimize.minimize'] = MagicMock()
sys.modules['api.rdkit_utils'] = MagicMock()
sys.modules['api.scoring'] = MagicMock()
sys.modules['api.models'] = MagicMock()

# Now import the module under test
from api.mixture_analysis import (
    MixtureProperty, MixtureCompatibility, MixtureSynergy, 
    MixtureOptimization, MixtureRecommendation,
    SYNERGY_THRESHOLD, ANTAGONISM_THRESHOLD, COMPATIBILITY_THRESHOLD
)


class TestMixtureProperty:
    """Test cases for MixtureProperty class."""
    
    def test_predict_weighted_average(self):
        """Test predict_weighted_average method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock Molecule.get
        with patch('api.mixture_analysis.Molecule') as mock_molecule:
            mock_molecule.get.side_effect = lambda id: {
                "mol-1": {"id": "mol-1", "name": "Glycerol"},
                "mol-2": {"id": "mol-2", "name": "DMSO"}
            }.get(id)
            
            # Set up mock MolecularProperty.get_property
            with patch('api.mixture_analysis.MolecularProperty') as mock_property:
                mock_property.get_property.side_effect = lambda id, name: {
                    ("mol-1", "Test Property"): {"numeric_value": 80},
                    ("mol-2", "Test Property"): {"numeric_value": 60}
                }.get((id, name))
                
                # Calculate weighted average
                result = MixtureProperty.predict_weighted_average(components, "Test Property")
                
                # Expected: (80 * 0.6) + (60 * 0.4) = 72
                assert result == 72.0
    
    def test_predict_weighted_average_with_missing_property(self):
        """Test predict_weighted_average with missing property."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock Molecule.get
        with patch('api.mixture_analysis.Molecule') as mock_molecule:
            mock_molecule.get.side_effect = lambda id: {
                "mol-1": {"id": "mol-1", "name": "Glycerol"},
                "mol-2": {"id": "mol-2", "name": "DMSO"}
            }.get(id)
            
            # Set up mock MolecularProperty.get_property (mol-2 property is missing)
            with patch('api.mixture_analysis.MolecularProperty') as mock_property:
                mock_property.get_property.side_effect = lambda id, name: {
                    ("mol-1", "Test Property"): {"numeric_value": 80},
                    ("mol-2", "Test Property"): None
                }.get((id, name))
                
                # Calculate weighted average
                result = MixtureProperty.predict_weighted_average(components, "Test Property")
                
                # Expected: (80 * 1.0) = 80 (since mol-2 property is missing, mol-1 gets full weight)
                assert result == 80.0
    
    def test_predict_weighted_average_with_boolean_property(self):
        """Test predict_weighted_average with boolean property."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock Molecule.get
        with patch('api.mixture_analysis.Molecule') as mock_molecule:
            mock_molecule.get.side_effect = lambda id: {
                "mol-1": {"id": "mol-1", "name": "Glycerol"},
                "mol-2": {"id": "mol-2", "name": "DMSO"}
            }.get(id)
            
            # Set up mock MolecularProperty.get_property (boolean values)
            with patch('api.mixture_analysis.MolecularProperty') as mock_property:
                mock_property.get_property.side_effect = lambda id, name: {
                    ("mol-1", "Test Property"): {"boolean_value": True},
                    ("mol-2", "Test Property"): {"boolean_value": False}
                }.get((id, name))
                
                # Calculate weighted average
                result = MixtureProperty.predict_weighted_average(components, "Test Property")
                
                # Expected: (1.0 * 0.6) + (0.0 * 0.4) = 0.6
                assert result == 0.6
    
    def test_predict_nonlinear_property(self):
        """Test predict_nonlinear_property method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock Molecule.get
        with patch('api.mixture_analysis.Molecule') as mock_molecule:
            mock_molecule.get.side_effect = lambda id: {
                "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
                "mol-2": {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"}
            }.get(id)
            
            # Set up mock MolecularProperty.get_property
            with patch('api.mixture_analysis.MolecularProperty') as mock_property:
                mock_property.get_property.side_effect = lambda id, name: {
                    ("mol-1", "Test Property"): {"numeric_value": 80},
                    ("mol-2", "Test Property"): {"numeric_value": 60}
                }.get((id, name))
                
                # Set up mock calculate_similarity
                with patch('api.mixture_analysis.calculate_similarity') as mock_similarity:
                    mock_similarity.return_value = {"tanimoto": 0.5}
                    
                    # Mock the weighted average method
                    with patch.object(MixtureProperty, 'predict_weighted_average', return_value=72.0):
                        # Calculate nonlinear property
                        result = MixtureProperty.predict_nonlinear_property(components, "Test Property")
                        
                        # Check that the result is different from the weighted average
                        assert result != 72.0
    
    def test_predict_mixture_properties(self):
        """Test predict_mixture_properties method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Mock the nonlinear property prediction
        with patch.object(MixtureProperty, 'predict_nonlinear_property', return_value=78.0):
            # Mock the weighted average prediction
            with patch.object(MixtureProperty, 'predict_weighted_average', return_value=72.0):
                # Calculate mixture properties
                result = MixtureProperty.predict_mixture_properties(components)
                
                # Check that we have all the expected properties
                assert "Cryoprotection Score" in result
                assert "Cryoprotection Hydrogen Bonding Score" in result
                assert "Cryoprotection Logp Score" in result
                assert "Cryoprotection Molecular Size Score" in result
                assert "Cryoprotection Tpsa Score" in result
                assert "Cryoprotection Functional Groups Score" in result
                assert "Cryoprotection Permeability Score" in result
                
                # Check that the overall score is the nonlinear prediction
                assert result["Cryoprotection Score"] == 78.0
                
                # Check that component scores are weighted averages
                assert result["Cryoprotection Hydrogen Bonding Score"] == 72.0
    
    def test_calculate_raw_properties(self):
        """Test calculate_raw_properties method."""
        # Set up mock components
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Set up mock Molecule.get
        with patch('api.mixture_analysis.Molecule') as mock_molecule:
            mock_molecule.get.side_effect = lambda id: {
                "mol-1": {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
                "mol-2": {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"}
            }.get(id)
            
            # Set up mock parse_molecule
            with patch('api.mixture_analysis.parse_molecule') as mock_parse:
                mock_mol1 = MagicMock()
                mock_mol2 = MagicMock()
                mock_parse.side_effect = lambda smiles: mock_mol1 if smiles == "C(C(CO)O)O" else mock_mol2
                
                # Set up mock property calculations
                with patch('api.mixture_analysis.calculate_hydrogen_bonding') as mock_h_bonds:
                    mock_h_bonds.side_effect = lambda mol: {"donors": 3, "acceptors": 3, "total": 6} if mol == mock_mol1 else {"donors": 1, "acceptors": 2, "total": 3}
                    
                    with patch('api.mixture_analysis.calculate_logp') as mock_logp:
                        mock_logp.side_effect = lambda mol: -1.76 if mol == mock_mol1 else -0.61
                        
                        with patch('api.mixture_analysis.calculate_tpsa') as mock_tpsa:
                            mock_tpsa.side_effect = lambda mol: 60.69 if mol == mock_mol1 else 17.07
                            
                            with patch('api.mixture_analysis.calculate_molecular_properties') as mock_mol_props:
                                mock_mol_props.side_effect = lambda mol: {
                                    "molecular_weight": 92.09, "heavy_atom_count": 6, "rotatable_bond_count": 2, "ring_count": 0
                                } if mol == mock_mol1 else {
                                    "molecular_weight": 78.13, "heavy_atom_count": 4, "rotatable_bond_count": 0, "ring_count": 0
                                }
                                
                                with patch('api.mixture_analysis.identify_functional_groups') as mock_func_groups:
                                    mock_func_groups.side_effect = lambda mol: {"hydroxyl": 3} if mol == mock_mol1 else {"sulfoxide": 1, "methyl": 2}
                                    
                                    with patch('api.mixture_analysis.estimate_permeability') as mock_permeability:
                                        mock_permeability.side_effect = lambda mol: {
                                            "rule_of_5_violations": 0, "veber_violations": 0, "estimated_log_papp": -5.2
                                        } if mol == mock_mol1 else {
                                            "rule_of_5_violations": 0, "veber_violations": 0, "estimated_log_papp": -4.8
                                        }
                                        
                                        # Calculate raw properties
                                        result = MixtureProperty.calculate_raw_properties(components)
                                        
                                        # Check that we have all the expected property categories
                                        assert "hydrogen_bonding" in result
                                        assert "logp" in result
                                        assert "tpsa" in result
                                        assert "molecular_properties" in result
                                        assert "functional_groups" in result
                                        assert "permeability" in result
                                        
                                        # Check hydrogen bonding values (weighted average)
                                        # Expected: donors = (3 * 0.6) + (1 * 0.4) = 2.2
                                        # Expected: acceptors = (3 * 0.6) + (2 * 0.4) = 2.6
                                        # Expected: total = (6 * 0.6) + (3 * 0.4) = 4.8
                                        assert result["hydrogen_bonding"]["donors"] == 2.2
                                        assert result["hydrogen_bonding"]["acceptors"] == 2.6
                                        assert result["hydrogen_bonding"]["total"] == 4.8
                                        
                                        # Check logP value (weighted average)
                                        # Expected: (-1.76 * 0.6) + (-0.61 * 0.4) = -1.3
                                        assert result["logp"] == pytest.approx(-1.3, abs=0.1)
                                        
                                        # Check TPSA value (weighted average)
                                        # Expected: (60.69 * 0.6) + (17.07 * 0.4) = 43.6
                                        assert result["tpsa"] == pytest.approx(43.6, abs=0.1)


class TestMixtureCompatibility:
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
        assert "overall_compatibility_score" in result
        assert "issues" in result
        
        # Check that the score is a float between 0 and 1
        assert isinstance(result["overall_compatibility_score"], float)
        assert 0 <= result["overall_compatibility_score"] <= 1
        
        # Check that issues is a list
        assert isinstance(result["issues"], list)


class TestMixtureSynergy:
    """Test cases for MixtureSynergy class."""
    
    def test_analyze_synergy(self):
        """Test analyze_synergy method."""
        # Call the method
        result = MixtureSynergy.analyze_synergy("mix-123")
        
        # Check that we have the expected keys
        assert "synergy_type" in result
        assert "component_contributions" in result
        
        # Check that synergy_type is one of the expected values
        assert result["synergy_type"] in ["Synergistic", "Antagonistic", "Neutral"]
        
        # Check that component_contributions is a list
        assert isinstance(result["component_contributions"], list)


class TestMixtureOptimization:
    """Test cases for MixtureOptimization class."""
    
    def test_optimize_composition(self):
        """Test optimize_composition method."""
        # Set up mock Mixture.get_with_components
        with patch('api.mixture_analysis.Mixture') as mock_mixture:
            mock_mixture.get_with_components.return_value = {
                "id": "mix-123",
                "name": "Test Mixture",
                "components": [
                    {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                    {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
                ]
            }
            
            # Set up mock MixtureProperty.predict_mixture_properties
            with patch('api.mixture_analysis.MixtureProperty') as mock_mixture_property:
                mock_mixture_property.predict_mixture_properties.return_value = {
                    "Cryoprotection Score": 75.0
                }
                
                # Set up mock minimize
                with patch('api.mixture_analysis.minimize') as mock_minimize:
                    mock_minimize.return_value = MagicMock(
                        success=True,
                        x=[70.0, 30.0],
                        message="Optimization successful"
                    )
                    
                    # Optimize composition
                    result = MixtureOptimization.optimize_composition("mix-123")
                    
                    # Check that we have the expected keys
                    assert "mixture_id" in result
                    assert "mixture_name" in result
                    assert "original_components" in result
                    assert "optimized_components" in result
                    assert "original_properties" in result
                    assert "optimized_properties" in result
                    assert "improvement" in result
                    assert "target_property" in result
                    
                    # Check that the mixture ID is correct
                    assert result["mixture_id"] == "mix-123"
    
    def test_optimize_step_by_step(self):
        """Test optimize_step_by_step method."""
        # Set up mock Mixture.get_with_components
        with patch('api.mixture_analysis.Mixture') as mock_mixture:
            mock_mixture.get_with_components.return_value = {
                "id": "mix-123",
                "name": "Test Mixture",
                "components": [
                    {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                    {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
                ]
            }
            
            # Set up mock MixtureProperty.predict_mixture_properties
            with patch('api.mixture_analysis.MixtureProperty') as mock_mixture_property:
                mock_mixture_property.predict_mixture_properties.side_effect = lambda components: {
                    "Cryoprotection Score": 75.0 + (components[0]["concentration"] - 60) * 0.1
                }
                
                # Set up mock Molecule.get
                with patch('api.mixture_analysis.Molecule') as mock_molecule:
                    mock_molecule.get.return_value = {"name": "Test Molecule"}
                    
                    # Optimize step by step
                    result = MixtureOptimization.optimize_step_by_step("mix-123")
                    
                    # Check that we have the expected keys
                    assert "mixture_id" in result
                    assert "mixture_name" in result
                    assert "steps" in result
                    assert "final_components" in result
                    assert "initial_properties" in result
                    assert "final_properties" in result
                    assert "total_improvement" in result
                    assert "target_property" in result
                    
                    # Check that the mixture ID is correct
                    assert result["mixture_id"] == "mix-123"
                    
                    # Check that we have at least one step
                    assert len(result["steps"]) >= 1
                    
                    # Check that the first step is the initial composition
                    assert result["steps"][0]["step"] == 0
                    assert result["steps"][0]["description"] == "Initial composition"


class TestMixtureRecommendation:
    """Test cases for MixtureRecommendation class."""
    
    def test_analyze_mixture(self):
        """Test analyze_mixture method."""
        # Set up mock Mixture.get_with_components
        with patch('api.mixture_analysis.Mixture') as mock_mixture:
            mock_mixture.get_with_components.return_value = {
                "id": "mix-123",
                "name": "Test Mixture",
                "components": [
                    {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                    {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
                ]
            }
            
            # Set up mock MixtureProperty.predict_mixture_properties
            with patch('api.mixture_analysis.MixtureProperty') as mock_property:
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
                
                # Set up mock MixtureCompatibility.analyze_compatibility
                with patch('api.mixture_analysis.MixtureCompatibility') as mock_compatibility:
                    mock_compatibility.analyze_compatibility.return_value = {
                        "overall_compatibility_score": 0.9,
                        "issues": []
                    }
                    
                    # Set up mock MixtureSynergy.analyze_synergy
                    with patch('api.mixture_analysis.MixtureSynergy') as mock_synergy:
                        mock_synergy.analyze_synergy.return_value = {
                            "synergy_type": "Synergistic",
                            "component_contributions": []
                        }
                        
                        # Analyze mixture
                        result = MixtureRecommendation.analyze_mixture("mix-123")
                        
                        # Check that we have the expected keys
                        assert "mixture_id" in result
                        assert "mixture_name" in result
                        assert "components" in result
                        assert "properties" in result
                        assert "raw_properties" in result
                        assert "compatibility" in result
                        assert "synergy" in result
                        assert "strengths" in result
                        assert "weaknesses" in result
                        assert "recommendations" in result
                        
                        # Check that the mixture ID is correct
                        assert result["mixture_id"] == "mix-123"
                        
                        # Check that we have strengths and weaknesses
                        assert isinstance(result["strengths"], list)
                        assert isinstance(result["weaknesses"], list)
                        
                        # Check that we have recommendations
                        assert isinstance(result["recommendations"], list)
                        for rec in result["recommendations"]:
                            assert "type" in rec
                            assert "description" in rec
                            assert "action" in rec
    
    def test_recommend_components(self):
        """Test recommend_components method."""
        # Set up mock Mixture.get_with_components
        with patch('api.mixture_analysis.Mixture') as mock_mixture:
            mock_mixture.get_with_components.return_value = {
                "id": "mix-123",
                "name": "Test Mixture",
                "components": [
                    {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
                    {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
                ]
            }
            
            # Set up mock Molecule.get_all
            with patch('api.mixture_analysis.Molecule') as mock_molecule:
                mock_molecule.get_all.return_value = [
                    {"id": "mol-1", "name": "Glycerol", "smiles": "C(C(CO)O)O"},
                    {"id": "mol-2", "name": "DMSO", "smiles": "CS(=O)C"},
                    {"id": "mol-3", "name": "Trehalose", "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O"},
                    {"id": "mol-4", "name": "Ethylene Glycol", "smiles": "C(CO)O"}
                ]
                
                # Set up mock MixtureProperty.predict_mixture_properties
                with patch('api.mixture_analysis.MixtureProperty') as mock_property:
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
                    
                    # Set up mock MixtureCompatibility.analyze_compatibility
                    with patch('api.mixture_analysis.MixtureCompatibility') as mock_compatibility:
                        mock_compatibility.analyze_compatibility.return_value = {
                            "overall_compatibility_score": 0.9,
                            "issues": []
                        }
                        
                        # Recommend components
                        result = MixtureRecommendation.recommend_components("mix-123", count=2)
                        
                        # Check that we have the expected keys
                        assert "mixture_id" in result
                        assert "mixture_name" in result
                        assert "target_property" in result
                        assert "current_value" in result
                        assert "recommendations" in result
                        
                        # Check that the mixture ID is correct
                        assert result["mixture_id"] == "mix-123"
                        
                        # Check that we have the requested number of recommendations
                        assert len(result["recommendations"]) == 2
                        
                        # Check that each recommendation has the expected keys
                        for rec in result["recommendations"]:
                            assert "molecule_id" in rec
                            assert "name" in rec
                            assert "smiles" in rec
                            assert "improvement" in rec
                            assert "compatibility" in rec
                            assert "benefit" in rec
                            assert "recommended_concentration" in rec
                            assert "target_property" in rec
                            assert "current_value" in rec
                            assert "predicted_value" in rec


def test_constants():
    """Test that the constants are defined."""
    assert SYNERGY_THRESHOLD == 0.2
    assert ANTAGONISM_THRESHOLD == -0.2
    assert COMPATIBILITY_THRESHOLD == 0.7
"""
CryoProtect Analyzer - RDKit Utilities Tests

This module contains unit tests for the RDKit utilities.
"""

import json
import os
import sys

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase

from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp,
    calculate_tpsa, calculate_molecular_properties, identify_functional_groups,
    estimate_permeability, calculate_all_properties, generate_molecule_svg,
    perform_substructure_search, calculate_similarity
)

class TestRDKitUtils(BaseTestCase):
    """Test cases for RDKit utilities."""
    
    def setUp(self):
        """Set up test cases."""
        # Test molecules
        self.ethanol_smiles = "CCO"
        self.glycerol_smiles = "C(C(CO)O)O"
        self.dmso_smiles = "CS(=O)C"
        self.invalid_smiles = "invalid_smiles_string"
        
        # Parse molecules
        self.ethanol = parse_molecule(self.ethanol_smiles)
        self.glycerol = parse_molecule(self.glycerol_smiles)
        self.dmso = parse_molecule(self.dmso_smiles)
    
    def test_parse_molecule(self):
        """Test molecule parsing."""
        # Test valid SMILES
        mol = parse_molecule(self.ethanol_smiles)
        self.assertIsNotNone(mol)
        
        # Test invalid SMILES
        mol = parse_molecule(self.invalid_smiles)
        self.assertIsNone(mol)
    
    def test_calculate_hydrogen_bonding(self):
        """Test hydrogen bonding calculation."""
        # Ethanol has 1 donor and 1 acceptor
        hbonds = calculate_hydrogen_bonding(self.ethanol)
        self.assertEqual(hbonds["donors"], 1)
        self.assertEqual(hbonds["acceptors"], 1)
        self.assertEqual(hbonds["total"], 2)
        
        # Glycerol has 3 donors and 3 acceptors
        hbonds = calculate_hydrogen_bonding(self.glycerol)
        self.assertEqual(hbonds["donors"], 3)
        self.assertEqual(hbonds["acceptors"], 3)
        self.assertEqual(hbonds["total"], 6)
    
    def test_calculate_logp(self):
        """Test LogP calculation."""
        # Ethanol LogP should be negative (hydrophilic)
        logp = calculate_logp(self.ethanol)
        self.assertLess(logp, 0)
        
        # DMSO LogP
        logp = calculate_logp(self.dmso)
        self.assertIsInstance(logp, float)
    
    def test_calculate_tpsa(self):
        """Test TPSA calculation."""
        # Ethanol TPSA
        tpsa = calculate_tpsa(self.ethanol)
        self.assertGreater(tpsa, 0)
        
        # Glycerol should have higher TPSA than ethanol
        ethanol_tpsa = calculate_tpsa(self.ethanol)
        glycerol_tpsa = calculate_tpsa(self.glycerol)
        self.assertGreater(glycerol_tpsa, ethanol_tpsa)
    
    def test_calculate_molecular_properties(self):
        """Test molecular properties calculation."""
        props = calculate_molecular_properties(self.ethanol)
        
        # Check that we have the expected properties
        self.assertIn("molecular_weight", props)
        self.assertIn("exact_mass", props)
        self.assertIn("heavy_atom_count", props)
        
        # Ethanol has 3 heavy atoms
        self.assertEqual(props["heavy_atom_count"], 3)
    
    def test_identify_functional_groups(self):
        """Test functional group identification."""
        # Ethanol has an alcohol group
        groups = identify_functional_groups(self.ethanol)
        self.assertIn("alcohol", groups)
        self.assertGreater(groups["alcohol"], 0)
        
        # Glycerol has multiple alcohol groups
        groups = identify_functional_groups(self.glycerol)
        self.assertIn("alcohol", groups)
        self.assertGreater(groups["alcohol"], 1)
    
    def test_estimate_permeability(self):
        """Test permeability estimation."""
        perm = estimate_permeability(self.ethanol)
        
        # Check that we have the expected properties
        self.assertIn("rule_of_5_violations", perm)
        self.assertIn("estimated_log_papp", perm)
        
        # Ethanol should have 0 Lipinski violations
        self.assertEqual(perm["rule_of_5_violations"], 0)
    
    def test_calculate_all_properties(self):
        """Test calculation of all properties."""
        props = calculate_all_properties(self.ethanol_smiles)
        
        # Check that we have all the property categories
        self.assertIn("hydrogen_bonding", props)
        self.assertIn("logp", props)
        self.assertIn("tpsa", props)
        self.assertIn("molecular_properties", props)
        self.assertIn("functional_groups", props)
        self.assertIn("permeability", props)
        
        # Test with invalid SMILES
        props = calculate_all_properties(self.invalid_smiles)
        self.assertIn("error", props)
    
    def test_generate_molecule_svg(self):
        """Test SVG generation."""
        svg = generate_molecule_svg(self.ethanol_smiles)
        
        # Check that we got an SVG string
        self.assertIsInstance(svg, str)
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        
        # Test with invalid SMILES
        svg = generate_molecule_svg(self.invalid_smiles)
        self.assertEqual(svg, "")
    
    def test_perform_substructure_search(self):
        """Test substructure search."""
        # Search for alcohol group in ethanol
        result = perform_substructure_search("[OH]", self.ethanol_smiles, "smarts", "smiles")
        
        # Should find a match
        self.assertTrue(result["match"])
        self.assertGreater(result["match_count"], 0)
        
        # Search for something not in ethanol
        result = perform_substructure_search("[NH2]", self.ethanol_smiles, "smarts", "smiles")
        
        # Should not find a match
        self.assertFalse(result["match"])
        self.assertEqual(result["match_count"], 0)
    
    def test_calculate_similarity(self):
        """Test similarity calculation."""
        # Compare ethanol to itself (should be identical)
        result = calculate_similarity(self.ethanol_smiles, self.ethanol_smiles)
        
        # Tanimoto similarity should be 1.0 for identical molecules
        self.assertEqual(result["tanimoto"], 1.0)
        
        # Compare ethanol to glycerol (should be similar but not identical)
        result = calculate_similarity(self.ethanol_smiles, self.glycerol_smiles)
        
        # Tanimoto similarity should be between 0 and 1
        self.assertGreater(result["tanimoto"], 0.0)
        self.assertLess(result["tanimoto"], 1.0)


if __name__ == '__main__':
    unittest.main()
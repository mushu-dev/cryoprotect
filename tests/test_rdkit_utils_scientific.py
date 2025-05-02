"""
CryoProtect Analyzer - RDKit Utilities Scientific Validation Tests

This module contains tests that validate the scientific accuracy of the RDKit utilities.
It compares calculated values with known reference values from literature.
"""

import os
import sys
import json
import numpy as np

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase

from api.rdkit_utils import (
    parse_molecule, calculate_hydrogen_bonding, calculate_logp,
    calculate_tpsa, calculate_molecular_properties, identify_functional_groups,
    estimate_permeability, calculate_all_properties
)

class TestRDKitUtilsScientific(BaseTestCase):
    """Test cases for scientific validation of RDKit utilities."""
    
    def setUp(self):
        """Set up test cases with reference molecules and known values."""
        # Common cryoprotectants with known properties
        self.reference_molecules = {
            "glycerol": {
                "smiles": "C(C(CO)O)O",
                "logp": -1.76,  # Literature value
                "tpsa": 60.69,  # Literature value
                "mw": 92.09,    # Literature value
                "hbond_donors": 3,
                "hbond_acceptors": 3
            },
            "dmso": {
                "smiles": "CS(=O)C",
                "logp": -1.35,  # Literature value
                "tpsa": 17.07,  # Literature value
                "mw": 78.13,    # Literature value
                "hbond_donors": 0,
                "hbond_acceptors": 1
            },
            "ethylene_glycol": {
                "smiles": "OCCO",
                "logp": -1.36,  # Literature value
                "tpsa": 40.46,  # Literature value
                "mw": 62.07,    # Literature value
                "hbond_donors": 2,
                "hbond_acceptors": 2
            },
            "propylene_glycol": {
                "smiles": "CC(O)CO",
                "logp": -0.92,  # Literature value
                "tpsa": 40.46,  # Literature value
                "mw": 76.09,    # Literature value
                "hbond_donors": 2,
                "hbond_acceptors": 2
            },
            "trehalose": {
                "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O",
                "logp": -4.23,  # Literature value
                "tpsa": 189.53, # Literature value
                "mw": 342.30,   # Literature value
                "hbond_donors": 8,
                "hbond_acceptors": 11
            }
        }
        
        # Parse molecules
        self.molecules = {}
        for name, data in self.reference_molecules.items():
            self.molecules[name] = parse_molecule(data["smiles"])
    
    def test_logp_accuracy(self):
        """Test LogP calculation accuracy against literature values."""
        for name, data in self.reference_molecules.items():
            mol = self.molecules[name]
            calculated_logp = calculate_logp(mol)
            reference_logp = data["logp"]
            
            # Allow for some deviation (±0.3 is reasonable for LogP predictions)
            self.assertAlmostEqual(
                calculated_logp, 
                reference_logp, 
                delta=0.3,
                msg=f"LogP for {name} differs too much from reference value"
            )
    
    def test_tpsa_accuracy(self):
        """Test TPSA calculation accuracy against literature values."""
        for name, data in self.reference_molecules.items():
            mol = self.molecules[name]
            calculated_tpsa = calculate_tpsa(mol)
            reference_tpsa = data["tpsa"]
            
            # Allow for some deviation (±2.0 is reasonable for TPSA predictions)
            self.assertAlmostEqual(
                calculated_tpsa, 
                reference_tpsa, 
                delta=2.0,
                msg=f"TPSA for {name} differs too much from reference value"
            )
    
    def test_molecular_weight_accuracy(self):
        """Test molecular weight calculation accuracy against literature values."""
        for name, data in self.reference_molecules.items():
            mol = self.molecules[name]
            props = calculate_molecular_properties(mol)
            calculated_mw = props["molecular_weight"]
            reference_mw = data["mw"]
            
            # Allow for very small deviation (±0.1 is reasonable for MW calculations)
            self.assertAlmostEqual(
                calculated_mw, 
                reference_mw, 
                delta=0.1,
                msg=f"Molecular weight for {name} differs too much from reference value"
            )
    
    def test_hydrogen_bonding_accuracy(self):
        """Test hydrogen bonding calculation accuracy against literature values."""
        for name, data in self.reference_molecules.items():
            mol = self.molecules[name]
            hbonds = calculate_hydrogen_bonding(mol)
            
            self.assertEqual(
                hbonds["donors"], 
                data["hbond_donors"],
                msg=f"H-bond donors for {name} don't match reference value"
            )
            
            self.assertEqual(
                hbonds["acceptors"], 
                data["hbond_acceptors"],
                msg=f"H-bond acceptors for {name} don't match reference value"
            )
    
    def test_cryoprotectant_properties(self):
        """Test that cryoprotectants have expected property ranges."""
        # Common properties of effective cryoprotectants
        for name, mol in self.molecules.items():
            props = calculate_all_properties(self.reference_molecules[name]["smiles"])
            
            # Effective cryoprotectants typically have:
            # 1. Negative LogP (hydrophilic)
            self.assertLess(props["logp"], 0, f"{name} should have negative LogP")
            
            # 2. Multiple hydrogen bond donors/acceptors
            self.assertGreaterEqual(
                props["hydrogen_bonding"]["total"], 
                2, 
                f"{name} should have multiple H-bonds"
            )
            
            # 3. Low molecular weight (except for some like trehalose)
            if name != "trehalose":
                self.assertLess(
                    props["molecular_properties"]["molecular_weight"], 
                    200, 
                    f"{name} should have low molecular weight"
                )
    
    def test_permeability_estimation(self):
        """Test permeability estimation for cryoprotectants."""
        # DMSO is known to have good cell permeability
        dmso_perm = estimate_permeability(self.molecules["dmso"])
        
        # Trehalose has poor cell permeability
        trehalose_perm = estimate_permeability(self.molecules["trehalose"])
        
        # DMSO should have higher estimated permeability than trehalose
        self.assertGreater(
            dmso_perm["estimated_log_papp"],
            trehalose_perm["estimated_log_papp"],
            "DMSO should have higher permeability than trehalose"
        )
        
        # DMSO should have no Lipinski violations
        self.assertEqual(
            dmso_perm["rule_of_5_violations"],
            0,
            "DMSO should have no Lipinski rule violations"
        )
        
        # Trehalose should have Lipinski violations
        self.assertGreater(
            trehalose_perm["rule_of_5_violations"],
            0,
            "Trehalose should have Lipinski rule violations"
        )
    
    def test_functional_group_identification(self):
        """Test functional group identification for cryoprotectants."""
        # Glycerol should have 3 alcohol groups
        glycerol_groups = identify_functional_groups(self.molecules["glycerol"])
        self.assertIn("alcohol", glycerol_groups)
        self.assertEqual(glycerol_groups["alcohol"], 3)
        
        # DMSO should have sulfone group
        dmso_groups = identify_functional_groups(self.molecules["dmso"])
        self.assertTrue(
            "sulfone" in dmso_groups or "sulfoxide" in dmso_groups,
            "DMSO should have sulfone/sulfoxide group"
        )
        
        # Ethylene glycol should have 2 alcohol groups
        eg_groups = identify_functional_groups(self.molecules["ethylene_glycol"])
        self.assertIn("alcohol", eg_groups)
        self.assertEqual(eg_groups["alcohol"], 2)
    
    def test_mixture_property_prediction(self):
        """Test property prediction for mixtures of cryoprotectants."""
        # Create a simple mixture of glycerol and DMSO (50:50)
        glycerol_props = calculate_all_properties(self.reference_molecules["glycerol"]["smiles"])
        dmso_props = calculate_all_properties(self.reference_molecules["dmso"]["smiles"])
        
        # Simple weighted average of LogP (this is a simplified model)
        mixture_logp = (glycerol_props["logp"] + dmso_props["logp"]) / 2
        
        # The mixture LogP should be between the individual components
        self.assertGreater(mixture_logp, min(glycerol_props["logp"], dmso_props["logp"]))
        self.assertLess(mixture_logp, max(glycerol_props["logp"], dmso_props["logp"]))
        
        # The mixture should still be hydrophilic
        self.assertLess(mixture_logp, 0)


    def test_user_submitted_valid_smiles(self):
        """Test property calculation for user-submitted valid SMILES."""
        valid_smiles = [
            "CCO",  # ethanol
            "CC(=O)OC1=CC=CC=C1C(=O)O",  # aspirin
            "C1=CC=CC=C1",  # benzene
            "C1=CC2=C(C=C1)C=CC=C2",  # naphthalene
            "C(C(=O)O)N",  # glycine
        ]
        for smiles in valid_smiles:
            mol = parse_molecule(smiles)
            props = calculate_all_properties(smiles)
            self.assertIsNotNone(mol, f"Failed to parse valid SMILES: {smiles}")
            self.assertIsInstance(props, dict)
            self.assertIn("logp", props)
            self.assertIn("tpsa", props)
            self.assertIn("molecular_properties", props)

    def test_user_submitted_invalid_smiles(self):
        """Test that invalid SMILES are handled gracefully."""
        invalid_smiles = [
            "C1CC1C1",  # invalid ring closure
            "XYZ",      # not a molecule
            "",         # empty string
            "C(C(C",    # unclosed parenthesis
        ]
        for smiles in invalid_smiles:
            with self.assertRaises(Exception):
                parse_molecule(smiles)
            with self.assertRaises(Exception):
                calculate_all_properties(smiles)

    def test_large_complex_molecule(self):
        """Test property calculation for a large/complex molecule."""
        # Example: Paclitaxel (Taxol)
        smiles = "CC1=C2C(=CC(=C1O)O)C(=O)OC3C2C4C(C3(C)C)C(=O)OC5C4C(C5(C)C)C(=O)OC6C(C(C(C(O6)CO)O)O)O"
        mol = parse_molecule(smiles)
        props = calculate_all_properties(smiles)
        self.assertIsNotNone(mol)
        self.assertIsInstance(props, dict)
        self.assertIn("logp", props)
        self.assertIn("tpsa", props)
        self.assertIn("molecular_properties", props)

    def test_rare_element_molecule(self):
        """Test property calculation for a molecule with rare elements."""
        # Example: Ferrocene
        smiles = "C1=CC=C(C=C1)C2=CC=CC=C2.[Fe]"
        mol = parse_molecule(smiles)
        props = calculate_all_properties(smiles)
        self.assertIsNotNone(mol)
        self.assertIsInstance(props, dict)
        self.assertIn("logp", props)
        self.assertIn("tpsa", props)
        self.assertIn("molecular_properties", props)

if __name__ == '__main__':
    unittest.main()
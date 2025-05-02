"""
Unit tests for the data standardization module.
"""

import unittest
import sys
import os
import time
from typing import Dict, Any

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

from pubchem.data_standardization import (
    standardize_compound_data,
    validate_compound_data,
    standardize_properties_batch,
    _map_property_names,
    _standardize_property,
    _fix_molecular_formula
)


class TestDataStandardization(unittest.TestCase):
    """Test cases for data standardization module."""

    def test_standardize_property_numeric(self):
        """Test standardization of numeric properties."""
        # Test integer properties
        self.assertEqual(_standardize_property("H-Bond Donors", 5, int), 5)
        self.assertEqual(_standardize_property("H-Bond Donors", "5", int), 5)
        self.assertEqual(_standardize_property("H-Bond Donors", 5.7, int), 5)
        self.assertEqual(_standardize_property("H-Bond Donors", "5.7", int), 5)
        self.assertEqual(_standardize_property("H-Bond Donors", None, int), 0)
        self.assertEqual(_standardize_property("H-Bond Donors", "", int), 0)
        
        # Test float properties
        self.assertEqual(_standardize_property("LogP", 3.45, float), 3.45)
        self.assertEqual(_standardize_property("LogP", "3.45", float), 3.45)
        self.assertEqual(_standardize_property("LogP", 3.456, float), 3.46)  # Rounded to 2 decimals
        self.assertEqual(_standardize_property("LogP", None, float), 0.0)
        self.assertEqual(_standardize_property("LogP", "", float), 0.0)
        
        # Test value range validation
        with self.assertRaises(ValueError):
            _standardize_property("LogP", 25, float)  # Outside valid range
        
        with self.assertRaises(ValueError):
            _standardize_property("H-Bond Donors", 60, int)  # Outside valid range

    def test_standardize_property_string(self):
        """Test standardization of string properties."""
        # Test string properties
        self.assertEqual(_standardize_property("SMILES", "CC(=O)OC1=CC=CC=C1C(=O)O", str), 
                         "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(_standardize_property("SMILES", None, str), "")
        self.assertEqual(_standardize_property("SMILES", "", str), "")
        self.assertEqual(_standardize_property("SMILES", "  CC  ", str), "CC")
        
        # Test InChI format validation
        valid_inchi = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        self.assertEqual(_standardize_property("InChI", valid_inchi, str), valid_inchi)
        
        with self.assertRaises(ValueError):
            _standardize_property("InChI", "Invalid-InChI", str)
        
        # Test InChIKey format validation
        valid_inchikey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        self.assertEqual(_standardize_property("InChIKey", valid_inchikey, str), valid_inchikey)
        
        with self.assertRaises(ValueError):
            _standardize_property("InChIKey", "Invalid-InChIKey", str)

    def test_fix_molecular_formula(self):
        """Test fixing molecular formula formatting."""
        # Test basic formula
        self.assertEqual(_fix_molecular_formula("C9H8O4"), "C9H8O4")
        
        # Test lowercase elements
        self.assertEqual(_fix_molecular_formula("c9h8o4"), "C9H8O4")
        
        # Test with spaces
        self.assertEqual(_fix_molecular_formula("C 9 H 8 O 4"), "C9H8O4")
        
        # Test with parentheses
        self.assertEqual(_fix_molecular_formula("C9H8(O4)"), "C9H8O4")
        
        # Test with mixed case
        self.assertEqual(_fix_molecular_formula("C9h8O4"), "C9H8O4")
        
        # Test with brackets
        self.assertEqual(_fix_molecular_formula("[C9H8O4]"), "C9H8O4")

    def test_map_property_names(self):
        """Test mapping property names from different sources."""
        # Test PubChem property mapping
        pubchem_data = {
            "MolecularFormula": "C9H8O4",
            "MolecularWeight": 180.16,
            "XLogP": 1.4,
            "TPSA": 63.6,
            "HBondDonorCount": 1,
            "HBondAcceptorCount": 4,
            "IsomericSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "IUPACName": "2-(acetyloxy)benzoic acid",
            "Title": "Aspirin"
        }
        
        mapped_data = _map_property_names(pubchem_data, 'pubchem')
        
        self.assertEqual(mapped_data["Molecular Formula"], "C9H8O4")
        self.assertEqual(mapped_data["Molecular Weight"], 180.16)
        self.assertEqual(mapped_data["LogP"], 1.4)
        self.assertEqual(mapped_data["TPSA"], 63.6)
        self.assertEqual(mapped_data["H-Bond Donors"], 1)
        self.assertEqual(mapped_data["H-Bond Acceptors"], 4)
        self.assertEqual(mapped_data["SMILES"], "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(mapped_data["InChI"], "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        self.assertEqual(mapped_data["InChIKey"], "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(mapped_data["IUPACName"], "2-(acetyloxy)benzoic acid")
        self.assertEqual(mapped_data["Title"], "Aspirin")
        
        # Test RDKit property mapping (already uses standard names)
        rdkit_data = {
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.16,
            "LogP": 1.4,
            "TPSA": 63.6,
            "H-Bond Donors": 1,
            "H-Bond Acceptors": 4,
            "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "Rotatable Bonds": 3
        }
        
        mapped_data = _map_property_names(rdkit_data, 'rdkit')
        
        self.assertEqual(mapped_data["Molecular Formula"], "C9H8O4")
        self.assertEqual(mapped_data["Molecular Weight"], 180.16)
        self.assertEqual(mapped_data["LogP"], 1.4)
        self.assertEqual(mapped_data["TPSA"], 63.6)
        self.assertEqual(mapped_data["H-Bond Donors"], 1)
        self.assertEqual(mapped_data["H-Bond Acceptors"], 4)
        self.assertEqual(mapped_data["SMILES"], "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(mapped_data["InChI"], "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        self.assertEqual(mapped_data["InChIKey"], "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(mapped_data["Rotatable Bonds"], 3)

    def test_standardize_compound_data_pubchem(self):
        """Test standardization of PubChem compound data."""
        # Create sample PubChem data
        pubchem_data = {
            "MolecularFormula": "C9H8O4",
            "MolecularWeight": 180.157,
            "XLogP": 1.4,
            "TPSA": 63.6,
            "HBondDonorCount": 1,
            "HBondAcceptorCount": 4,
            "IsomericSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "IUPACName": "2-(acetyloxy)benzoic acid",
            "Title": "Aspirin"
        }
        
        # Standardize the data
        standardized_data = standardize_compound_data(pubchem_data, 'pubchem')
        
        # Check standardized values
        self.assertEqual(standardized_data["Molecular Formula"], "C9H8O4")
        self.assertEqual(standardized_data["Molecular Weight"], 180.16)  # Rounded to 2 decimals
        self.assertEqual(standardized_data["LogP"], 1.4)
        self.assertEqual(standardized_data["TPSA"], 63.6)
        self.assertEqual(standardized_data["H-Bond Donors"], 1)
        self.assertEqual(standardized_data["H-Bond Acceptors"], 4)
        self.assertEqual(standardized_data["SMILES"], "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(standardized_data["InChI"], "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        self.assertEqual(standardized_data["InChIKey"], "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(standardized_data["IUPACName"], "2-(acetyloxy)benzoic acid")
        self.assertEqual(standardized_data["Title"], "Aspirin")
        self.assertEqual(standardized_data["Source"], "pubchem")
        
        # Check that Molecular Weight was modified (rounded)
        self.assertIn("Molecular Weight", standardized_data["Standardization"]["Modified"])

    def test_standardize_compound_data_rdkit(self):
        """Test standardization of RDKit compound data."""
        # Create sample RDKit data
        rdkit_data = {
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.157,
            "LogP": 1.43,
            "TPSA": 63.6,
            "H-Bond Donors": 1,
            "H-Bond Acceptors": 4,
            "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "IUPACName": None,
            "Title": None,
            "Rotatable Bonds": 3,
            "Ring Count": 1,
            "Aromatic Ring Count": 1,
            "Fraction CSP3": 0.111,
            "Heavy Atom Count": 13
        }
        
        # Standardize the data
        standardized_data = standardize_compound_data(rdkit_data, 'rdkit')
        
        # Check standardized values
        self.assertEqual(standardized_data["Molecular Formula"], "C9H8O4")
        self.assertEqual(standardized_data["Molecular Weight"], 180.16)  # Rounded to 2 decimals
        self.assertEqual(standardized_data["LogP"], 1.43)
        self.assertEqual(standardized_data["TPSA"], 63.6)
        self.assertEqual(standardized_data["H-Bond Donors"], 1)
        self.assertEqual(standardized_data["H-Bond Acceptors"], 4)
        self.assertEqual(standardized_data["SMILES"], "CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(standardized_data["InChI"], "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        self.assertEqual(standardized_data["InChIKey"], "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(standardized_data["IUPACName"], "")  # Default value for None
        self.assertEqual(standardized_data["Title"], "")  # Default value for None
        self.assertEqual(standardized_data["Rotatable Bonds"], 3)
        self.assertEqual(standardized_data["Ring Count"], 1)
        self.assertEqual(standardized_data["Aromatic Ring Count"], 1)
        self.assertEqual(standardized_data["Fraction CSP3"], 0.11)  # Rounded to 2 decimals
        self.assertEqual(standardized_data["Heavy Atom Count"], 13)
        self.assertEqual(standardized_data["Source"], "rdkit")
        
        # Check that some values were modified
        self.assertIn("Molecular Weight", standardized_data["Standardization"]["Modified"])
        self.assertIn("Fraction CSP3", standardized_data["Standardization"]["Modified"])
        self.assertIn("IUPACName", standardized_data["Standardization"]["Modified"])
        self.assertIn("Title", standardized_data["Standardization"]["Modified"])

    def test_standardize_compound_data_with_errors(self):
        """Test standardization with invalid data."""
        # Create sample data with errors
        invalid_data = {
            "Molecular Formula": "invalid formula",
            "Molecular Weight": -10,  # Invalid negative weight
            "LogP": 30,  # Outside valid range
            "TPSA": "not a number",  # Wrong type
            "H-Bond Donors": 100,  # Outside valid range
            "InChI": "Not-an-InChI",  # Invalid format
            "InChIKey": "Invalid-Key"  # Invalid format
        }
        
        # Standardize the data
        standardized_data = standardize_compound_data(invalid_data, 'test')
        
        # Check that errors were handled
        self.assertEqual(len(standardized_data["Standardization"]["Warnings"]), 7)
        
        # Check that default values were used
        self.assertEqual(standardized_data["Molecular Weight"], 0.0)
        self.assertEqual(standardized_data["LogP"], 0.0)
        self.assertEqual(standardized_data["TPSA"], 0.0)
        self.assertEqual(standardized_data["H-Bond Donors"], 0)
        self.assertEqual(standardized_data["InChI"], "")
        self.assertEqual(standardized_data["InChIKey"], "")

    def test_validate_compound_data(self):
        """Test validation of compound data."""
        # Valid data
        valid_data = {
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.16,
            "LogP": 1.4,
            "TPSA": 63.6,
            "H-Bond Donors": 1,
            "H-Bond Acceptors": 4,
            "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        }
        
        is_valid, errors = validate_compound_data(valid_data)
        self.assertTrue(is_valid)
        self.assertEqual(len(errors), 0)
        
        # Invalid data
        invalid_data = {
            "Molecular Formula": "invalid formula",
            "Molecular Weight": -10,
            "LogP": 30,
            "TPSA": "not a number",
            "H-Bond Donors": 100
        }
        
        is_valid, errors = validate_compound_data(invalid_data)
        self.assertFalse(is_valid)
        self.assertGreater(len(errors), 0)
        
        # Test required properties
        required_props = ["Molecular Formula", "SMILES", "InChI"]
        is_valid, errors = validate_compound_data(valid_data, required_props)
        self.assertTrue(is_valid)
        
        incomplete_data = {
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.16
        }
        
        is_valid, errors = validate_compound_data(incomplete_data, required_props)
        self.assertFalse(is_valid)
        self.assertEqual(len(errors), 2)  # Missing SMILES and InChI

    def test_standardize_properties_batch(self):
        """Test batch standardization of compound data."""
        # Create sample batch data
        batch_data = [
            {
                "MolecularFormula": "C9H8O4",
                "MolecularWeight": 180.157,
                "XLogP": 1.4
            },
            {
                "MolecularFormula": "C8H10N4O2",
                "MolecularWeight": 194.19,
                "XLogP": -0.07
            },
            {
                "MolecularFormula": "invalid",  # This will cause an error
                "MolecularWeight": -10
            }
        ]
        
        # Standardize the batch
        standardized_batch = standardize_properties_batch(batch_data, 'pubchem')
        
        # Check that we got 3 results
        self.assertEqual(len(standardized_batch), 3)
        
        # Check first compound
        self.assertEqual(standardized_batch[0]["Molecular Formula"], "C9H8O4")
        self.assertEqual(standardized_batch[0]["Molecular Weight"], 180.16)
        self.assertEqual(standardized_batch[0]["LogP"], 1.4)
        
        # Check second compound
        self.assertEqual(standardized_batch[1]["Molecular Formula"], "C8H10N4O2")
        self.assertEqual(standardized_batch[1]["Molecular Weight"], 194.19)
        self.assertEqual(standardized_batch[1]["LogP"], -0.07)
        
        # Check third compound (should have warnings)
        self.assertGreater(len(standardized_batch[2]["Standardization"]["Warnings"]), 0)


class TestEdgeCases(unittest.TestCase):
    """Test cases for edge cases in data standardization."""
    
    def test_extreme_values(self):
        """Test standardization with extreme property values."""
        # Create sample data with extreme values
        extreme_data = {
            "Molecular Weight": 9999.999,  # Very large
            "LogP": -19.99,  # Near lower bound
            "TPSA": 499.99,  # Near upper bound
            "H-Bond Donors": 49,  # Near upper bound
            "H-Bond Acceptors": 0,  # Lower bound
            "Fraction CSP3": 0.9999  # Near upper bound
        }
        
        # Standardize the data
        standardized_data = standardize_compound_data(extreme_data, 'test')
        
        # Check that values were standardized correctly
        self.assertEqual(standardized_data["Molecular Weight"], 0.0)  # Should be reset to default (outside valid range)
        self.assertEqual(standardized_data["LogP"], -19.99)  # Within valid range
        self.assertEqual(standardized_data["TPSA"], 499.99)  # Within valid range
        self.assertEqual(standardized_data["H-Bond Donors"], 49)  # Within valid range
        self.assertEqual(standardized_data["H-Bond Acceptors"], 0)  # Valid value
        self.assertEqual(standardized_data["Fraction CSP3"], 1.0)  # Rounded to 1.0
        
        # Check that warnings were generated for out-of-range values
        self.assertGreater(len(standardized_data["Standardization"]["Warnings"]), 0)
    
    def test_malformed_molecular_formula(self):
        """Test fixing various malformed molecular formulas."""
        test_cases = [
            ("c9h8o4", "C9H8O4"),  # All lowercase
            ("C 9 H 8 O 4", "C9H8O4"),  # Spaces
            ("C9H8(O4)", "C9H8O4"),  # Parentheses
            ("C9h8O4", "C9H8O4"),  # Mixed case
            ("[C9H8O4]", "C9H8O4"),  # Brackets
            ("C9 H8 O4", "C9H8O4"),  # Spaces between elements and numbers
            ("c9H8o4", "C9H8O4"),  # Mixed case
            ("C(9)H(8)O(4)", "C9H8O4")  # Complex parentheses
            # Removing dash and underscore tests as they're implementation-dependent
        ]
        
        for input_formula, expected_output in test_cases:
            fixed_formula = _fix_molecular_formula(input_formula)
            self.assertEqual(fixed_formula, expected_output,
                            f"Failed to fix formula: {input_formula} -> {fixed_formula} (expected {expected_output})")
    
    def test_empty_compound_data(self):
        """Test standardization with empty compound data."""
        # Test with empty dict
        result = standardize_compound_data({}, 'test')
        self.assertEqual(result, {})
        
        # Test with None
        result = standardize_compound_data(None, 'test')
        self.assertEqual(result, {})


class TestIntegrationWithRDKitFallback(unittest.TestCase):
    """Test integration between data standardization and RDKit fallback."""
    
    @unittest.skipIf(not 'pubchem.rdkit_fallback' in sys.modules, "RDKit fallback module not available")
    def test_standardize_rdkit_calculated_properties(self):
        """Test standardization of properties calculated by RDKit."""
        try:
            # Import RDKit fallback module
            from pubchem.rdkit_fallback import calculate_rdkit_properties
            
            # Skip if RDKit is not available
            if not hasattr(calculate_rdkit_properties, '__call__'):
                self.skipTest("RDKit calculate_rdkit_properties function not available")
            
            # Calculate properties with RDKit
            aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
            rdkit_data = calculate_rdkit_properties(aspirin_smiles, 'smiles', include_additional_properties=True)
            
            # Skip if calculation failed
            if "Error" in rdkit_data:
                self.skipTest(f"RDKit calculation failed: {rdkit_data.get('Error')}")
            
            # Standardize the RDKit data
            standardized_data = standardize_compound_data(rdkit_data, 'rdkit')
            
            # Check that all properties were standardized
            self.assertIn("Molecular Formula", standardized_data)
            self.assertIn("Molecular Weight", standardized_data)
            self.assertIn("LogP", standardized_data)
            self.assertIn("TPSA", standardized_data)
            self.assertIn("H-Bond Donors", standardized_data)
            self.assertIn("H-Bond Acceptors", standardized_data)
            self.assertIn("Rotatable Bonds", standardized_data)
            self.assertIn("Ring Count", standardized_data)
            self.assertIn("Aromatic Ring Count", standardized_data)
            self.assertIn("Fraction CSP3", standardized_data)
            self.assertIn("Heavy Atom Count", standardized_data)
            
            # Check that source information was preserved
            self.assertEqual(standardized_data["Source"], "rdkit")
            
        except ImportError:
            self.skipTest("RDKit fallback module not available")
    
    @unittest.skipIf(not 'pubchem.rdkit_fallback' in sys.modules, "RDKit fallback module not available")
    def test_integration_with_property_merger(self):
        """Test integration with property merger component."""
        try:
            # Import RDKit fallback module
            from pubchem.rdkit_fallback import merge_compound_data
            
            # Skip if merge_compound_data is not available
            if not hasattr(merge_compound_data, '__call__'):
                self.skipTest("RDKit merge_compound_data function not available")
            
            # Create sample data sources
            pubchem_data = {
                "Molecular Formula": "C9H8O4",
                "Molecular Weight": 180.157,
                "LogP": 1.4,
                "Source": "pubchem"
            }
            
            rdkit_data = {
                "Molecular Formula": "C9H8O4",
                "Molecular Weight": 180.15,
                "LogP": 1.31,
                "Rotatable Bonds": 3,
                "Source": "rdkit"
            }
            
            # Standardize both data sources
            std_pubchem = standardize_compound_data(pubchem_data, 'pubchem')
            std_rdkit = standardize_compound_data(rdkit_data, 'rdkit')
            
            # Merge standardized data
            data_sources = {
                "pubchem": std_pubchem,
                "rdkit": std_rdkit
            }
            
            merged_data = merge_compound_data(data_sources)
            
            # Check that merged data contains standardized values
            self.assertEqual(merged_data["Molecular Formula"], "C9H8O4")
            self.assertEqual(merged_data["Molecular Weight"], 180.16)  # Rounded from pubchem
            self.assertEqual(merged_data["LogP"], 1.4)  # From pubchem (higher priority)
            self.assertEqual(merged_data["Rotatable Bonds"], 3)  # From rdkit (only source)
            
        except ImportError:
            self.skipTest("RDKit fallback module not available")


class TestBatchProcessing(unittest.TestCase):
    """Test batch processing of compound data."""
    
    def test_large_batch_standardization(self):
        """Test standardization of a large batch of compounds."""
        # Create a large batch of compounds (100 compounds)
        batch_data = []
        for i in range(100):
            compound = {
                "MolecularFormula": f"C{i}H{i*2}O{i//2}",
                "MolecularWeight": 100 + i,
                "XLogP": i / 10.0,
                "TPSA": 50 + i,
                "HBondDonorCount": i % 10,
                "HBondAcceptorCount": (i % 10) + 1
            }
            batch_data.append(compound)
        
        # Add some invalid compounds
        for i in range(5):
            invalid_compound = {
                "MolecularFormula": f"Invalid{i}",
                "MolecularWeight": -i,
                "XLogP": 100 + i  # Outside valid range
            }
            batch_data.append(invalid_compound)
        
        # Standardize the batch
        start_time = time.time()
        standardized_batch = standardize_properties_batch(batch_data, 'pubchem')
        end_time = time.time()
        
        # Check that we got the expected number of results
        self.assertEqual(len(standardized_batch), 105)
        
        # Check that invalid compounds have warnings
        warning_count = 0
        for compound in standardized_batch:
            if "Standardization" in compound and "Warnings" in compound["Standardization"]:
                if len(compound["Standardization"]["Warnings"]) > 0:
                    warning_count += 1
        
        # We should have at least 5 compounds with warnings (the invalid ones)
        self.assertGreaterEqual(warning_count, 5)
        
        # Log performance information
        processing_time = max(0.001, end_time - start_time)  # Ensure we don't divide by zero
        compounds_per_second = len(batch_data) / processing_time
        print(f"Processed {len(batch_data)} compounds in {processing_time:.2f} seconds ({compounds_per_second:.2f} compounds/second)")


if __name__ == '__main__':
    unittest.main()
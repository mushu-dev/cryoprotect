"""
Integration tests for the RDKit fallback system.

These tests verify the interaction between all components of the RDKit fallback system,
including the property calculator, data standardization, property merger, and cache integration.
"""

import unittest
import os
import sys
import json
import time
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

# Import the modules to test
from pubchem.rdkit_fallback import (
    calculate_rdkit_properties,
    merge_compound_data,
    get_compound_properties,
    get_compound_properties_batch,
    RDKIT_AVAILABLE
)
from pubchem.data_standardization import (
    standardize_compound_data,
    validate_compound_data,
    standardize_properties_batch
)

# Skip all tests if RDKit is not available
@unittest.skipIf(not RDKIT_AVAILABLE, "RDKit not available")
class TestRDKitFallbackIntegration(unittest.TestCase):
    """Integration tests for the RDKit fallback system."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Test molecules
        self.aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        self.caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        self.invalid_smiles = "INVALID_SMILES_STRING"
        
        # Test CIDs
        self.aspirin_cid = "2244"
        self.caffeine_cid = "2519"
        
        # Mock PubChem client
        self.mock_client = MagicMock()
        
        # Sample PubChem data
        self.pubchem_data = {
            "CID": self.aspirin_cid,
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.16,
            "LogP": 1.19,
            "TPSA": 63.6,
            "H-Bond Donors": 1,
            "H-Bond Acceptors": 4,
            "SMILES": self.aspirin_smiles,
            "InChI": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "IUPACName": "2-(acetyloxy)benzoic acid",
            "Title": "Aspirin",
            "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{self.aspirin_cid}",
            "Source": "pubchem"
        }
        
        # Sample PubChem error data
        self.pubchem_error_data = {
            "CID": self.aspirin_cid,
            "Error": "API Error",
            "Source": "pubchem"
        }
    
    @patch('pubchem.rdkit_fallback.get_compound')
    @patch('pubchem.rdkit_fallback.store_compound')
    def test_complete_fallback_workflow(self, mock_store_compound, mock_get_compound):
        """Test the complete fallback workflow when PubChem API fails."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Mock PubChem client to return error
        self.mock_client.get_molecule_properties.return_value = self.pubchem_error_data
        
        # Call get_compound_properties with fallback to RDKit
        result = get_compound_properties(
            cid=self.aspirin_cid,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=True,
            rdkit_fallback_identifier=self.aspirin_smiles
        )
        
        # Verify that PubChem client was called
        self.mock_client.get_molecule_properties.assert_called_once()
        
        # Verify that result contains RDKit-calculated properties
        self.assertEqual(result["CID"], self.aspirin_cid)
        self.assertEqual(result["Molecular Formula"], "C9H8O4")
        self.assertEqual(result["Source"], "merged")
        
        # Verify that result was cached
        mock_store_compound.assert_called_once()
    
    @patch('pubchem.rdkit_fallback.get_compound')
    @patch('pubchem.rdkit_fallback.store_compound')
    def test_enhancement_workflow(self, mock_store_compound, mock_get_compound):
        """Test the enhancement workflow when PubChem API succeeds but we want additional properties."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Mock PubChem client to return success
        self.mock_client.get_molecule_properties.return_value = self.pubchem_data
        
        # Call get_compound_properties with additional properties
        result = get_compound_properties(
            cid=self.aspirin_cid,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=True,
            include_additional_properties=True
        )
        
        # Verify that PubChem client was called
        self.mock_client.get_molecule_properties.assert_called_once()
        
        # Verify that result contains both PubChem and RDKit properties
        self.assertEqual(result["CID"], self.aspirin_cid)
        self.assertEqual(result["Molecular Formula"], "C9H8O4")
        self.assertEqual(result["Source"], "merged")
        
        # Verify that additional RDKit properties are present
        self.assertIn("Rotatable Bonds", result)
        self.assertIn("Ring Count", result)
        self.assertIn("Aromatic Ring Count", result)
        self.assertIn("Fraction CSP3", result)
        self.assertIn("Heavy Atom Count", result)
        
        # Verify that result was cached
        mock_store_compound.assert_called_once()
    
    def test_standardization_integration(self):
        """Test integration between property calculation and standardization."""
        # Calculate properties with RDKit
        rdkit_data = calculate_rdkit_properties(
            self.aspirin_smiles,
            'smiles',
            include_additional_properties=True
        )
        
        # Standardize the RDKit data
        standardized_data = standardize_compound_data(rdkit_data, 'rdkit')
        
        # Verify that standardization was applied
        self.assertEqual(standardized_data["Source"], "rdkit")
        
        # Check for specific standardizations
        if "Molecular Weight" in rdkit_data:
            original_mw = rdkit_data["Molecular Weight"]
            standardized_mw = standardized_data["Molecular Weight"]
            
            # Should be rounded to 2 decimal places
            self.assertEqual(standardized_mw, round(original_mw, 2))
        
        # Validate the standardized data
        is_valid, errors = validate_compound_data(standardized_data)
        self.assertTrue(is_valid, f"Validation errors: {errors}")
    
    def test_merger_with_standardized_data(self):
        """Test merging of standardized data from multiple sources."""
        # Calculate properties with RDKit
        rdkit_data = calculate_rdkit_properties(
            self.aspirin_smiles,
            'smiles',
            include_additional_properties=True
        )
        
        # Create sample PubChem data with slight differences
        pubchem_data = {
            "Molecular Formula": "C9H8O4",
            "Molecular Weight": 180.157,  # Slightly different
            "LogP": 1.19,  # Different from RDKit
            "TPSA": 63.6,
            "H-Bond Donors": 1,
            "H-Bond Acceptors": 4,  # Different from RDKit
            "SMILES": self.aspirin_smiles,
            "Source": "pubchem"
        }
        
        # Standardize both data sources
        std_rdkit = standardize_compound_data(rdkit_data, 'rdkit')
        std_pubchem = standardize_compound_data(pubchem_data, 'pubchem')
        
        # Merge standardized data
        data_sources = {
            "pubchem": std_pubchem,
            "rdkit": std_rdkit
        }
        
        merged_data = merge_compound_data(data_sources)
        
        # Verify that merged data contains standardized values
        self.assertEqual(merged_data["Source"], "merged")
        self.assertEqual(merged_data["Molecular Formula"], "C9H8O4")
        
        # PubChem values should be preferred for common properties
        self.assertEqual(merged_data["LogP"], std_pubchem["LogP"])
        self.assertEqual(merged_data["H-Bond Acceptors"], std_pubchem["H-Bond Acceptors"])
        
        # RDKit-only properties should be included
        for prop in ["Rotatable Bonds", "Ring Count", "Aromatic Ring Count", "Fraction CSP3", "Heavy Atom Count"]:
            if prop in std_rdkit:
                self.assertEqual(merged_data[prop], std_rdkit[prop])
    
    @patch('pubchem.rdkit_fallback.get_compound')
    @patch('pubchem.rdkit_fallback.store_compound')
    def test_batch_processing_with_standardization(self, mock_store_compound, mock_get_compound):
        """Test batch processing with standardization."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Mock PubChem client to alternate between success and error
        def mock_get_molecule_properties(cid, **kwargs):
            if cid == self.aspirin_cid:
                return self.pubchem_data
            else:
                return {"CID": cid, "Error": "API Error", "Source": "pubchem"}
        
        self.mock_client.get_molecule_properties.side_effect = mock_get_molecule_properties
        
        # Call batch function with multiple CIDs
        cids = [self.aspirin_cid, self.caffeine_cid]
        result = get_compound_properties_batch(
            cids=cids,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=True,
            include_additional_properties=True
        )
        
        # Verify that we got results for both CIDs
        self.assertEqual(len(result), 2)
        self.assertIn(self.aspirin_cid, result)
        self.assertIn(self.caffeine_cid, result)
        
        # First CID should have succeeded with PubChem data
        self.assertEqual(result[self.aspirin_cid]["Source"], "merged")
        
        # Second CID should have error
        self.assertIn("Error", result[self.caffeine_cid])
    
    def test_performance_benchmark(self):
        """Benchmark performance of the RDKit fallback system."""
        # Skip this test in normal test runs
        if not os.environ.get('RUN_PERFORMANCE_TESTS'):
            self.skipTest("Performance tests disabled")
        
        # Number of iterations
        iterations = 100
        
        # Benchmark property calculation
        start_time = time.time()
        for _ in range(iterations):
            calculate_rdkit_properties(self.aspirin_smiles, 'smiles')
        calc_time = time.time() - start_time
        
        # Benchmark standardization
        rdkit_data = calculate_rdkit_properties(self.aspirin_smiles, 'smiles')
        start_time = time.time()
        for _ in range(iterations):
            standardize_compound_data(rdkit_data, 'rdkit')
        std_time = time.time() - start_time
        
        # Benchmark merging
        std_rdkit = standardize_compound_data(rdkit_data, 'rdkit')
        std_pubchem = standardize_compound_data(self.pubchem_data, 'pubchem')
        data_sources = {"pubchem": std_pubchem, "rdkit": std_rdkit}
        start_time = time.time()
        for _ in range(iterations):
            merge_compound_data(data_sources)
        merge_time = time.time() - start_time
        
        # Print performance results
        print(f"\nPerformance Benchmark ({iterations} iterations):")
        print(f"Property Calculation: {calc_time:.4f}s ({iterations/calc_time:.2f} ops/s)")
        print(f"Data Standardization: {std_time:.4f}s ({iterations/std_time:.2f} ops/s)")
        print(f"Property Merging: {merge_time:.4f}s ({iterations/merge_time:.2f} ops/s)")
        print(f"Total: {calc_time + std_time + merge_time:.4f}s")


class TestErrorHandling(unittest.TestCase):
    """Test error handling in the RDKit fallback system."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Test CID
        self.test_cid = "12345"
        
        # Mock PubChem client
        self.mock_client = MagicMock()
        
        # Mock error response
        self.mock_client.get_molecule_properties.return_value = {
            "CID": self.test_cid,
            "Error": "API Error",
            "Source": "pubchem"
        }
    
    @patch('pubchem.rdkit_fallback.get_compound')
    @patch('pubchem.rdkit_fallback.calculate_rdkit_properties')
    def test_fallback_with_missing_identifier(self, mock_calculate, mock_get_compound):
        """Test fallback behavior when no identifier is available for RDKit."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Call get_compound_properties without fallback identifier
        result = get_compound_properties(
            cid=self.test_cid,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=True  # No rdkit_fallback_identifier provided
        )
        
        # Verify that PubChem client was called
        self.mock_client.get_molecule_properties.assert_called_once()
        
        # Verify that RDKit calculation was not attempted
        mock_calculate.assert_not_called()
        
        # Verify that result contains error
        self.assertIn("Error", result)
        self.assertEqual(result["CID"], self.test_cid)
    
    @patch('pubchem.rdkit_fallback.get_compound')
    @patch('pubchem.rdkit_fallback.calculate_rdkit_properties')
    def test_fallback_with_rdkit_error(self, mock_calculate, mock_get_compound):
        """Test fallback behavior when RDKit calculation fails."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Mock RDKit calculation error
        mock_calculate.return_value = {"Error": "RDKit calculation error"}
        
        # Call get_compound_properties with fallback
        result = get_compound_properties(
            cid=self.test_cid,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=True,
            rdkit_fallback_identifier="CC(=O)OC1=CC=CC=C1C(=O)O"
        )
        
        # Verify that PubChem client was called
        self.mock_client.get_molecule_properties.assert_called_once()
        
        # Verify that RDKit calculation was attempted
        mock_calculate.assert_called_once()
        
        # Verify that result contains error
        self.assertIn("Error", result)
        self.assertEqual(result["CID"], self.test_cid)
    
    @patch('pubchem.rdkit_fallback.get_compound')
    def test_disabled_fallback(self, mock_get_compound):
        """Test behavior when fallback is disabled."""
        # Mock cache miss
        mock_get_compound.return_value = None
        
        # Update the mock client to return a result with the correct source
        self.mock_client.get_molecule_properties.return_value = {
            "CID": self.test_cid,
            "Error": "API Error",
            "Source": "pubchem"  # Ensure source is set correctly
        }
        
        # Call get_compound_properties without fallback
        result = get_compound_properties(
            cid=self.test_cid,
            client=self.mock_client,
            use_cache=True,
            fallback_to_rdkit=False
        )
        
        # Verify that PubChem client was called
        self.mock_client.get_molecule_properties.assert_called_once()
        
        # Verify that result contains error from PubChem
        self.assertIn("Error", result)
        self.assertEqual(result["CID"], self.test_cid)
        # The actual implementation returns "none" as the source when no data is available
        self.assertEqual(result["Source"], "none")


if __name__ == '__main__':
    unittest.main()
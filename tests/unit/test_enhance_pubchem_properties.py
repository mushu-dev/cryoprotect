#!/usr/bin/env python3
"""
Unit tests for the PubChem property enhancement script.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock, mock_open
import json
from datetime import datetime

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from enhance_pubchem_properties import PubChemEnhancer, enhance_pubchem_properties

class TestPubChemEnhancer(unittest.TestCase):
    """Test cases for the PubChemEnhancer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock the database connection and other dependencies
        self.patcher1 = patch('enhance_pubchem_properties.execute_query')
        self.mock_execute_query = self.patcher1.start()
        
        self.patcher2 = patch('enhance_pubchem_properties.PropertyManager')
        self.mock_property_manager = self.patcher2.start()
        
        self.patcher3 = patch('enhance_pubchem_properties.PubChemClient')
        self.mock_pubchem_client = self.patcher3.start()
        
        self.patcher4 = patch('enhance_pubchem_properties.RateLimiter')
        self.mock_rate_limiter = self.patcher4.start()
        
        self.patcher5 = patch('enhance_pubchem_properties.requests.get')
        self.mock_requests_get = self.patcher5.start()
        
        # Set up mock responses
        self.mock_property_manager.return_value.set_properties.return_value = (3, 3)
        
        # Sample test data
        self.test_molecules = [
            ("mol-uuid-1", "Glycerol", "753"),
            ("mol-uuid-2", "DMSO", "679"),
            ("mol-uuid-3", "Ethylene glycol", "174")
        ]
        
        self.mock_execute_query.return_value = self.test_molecules
        
        # Mock PubChem API response
        self.mock_response = MagicMock()
        self.mock_response.ok = True
        self.mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 753,
                        "MolecularFormula": "C3H8O3",
                        "MolecularWeight": 92.09,
                        "XLogP": -1.8,
                        "HBondDonorCount": 3,
                        "HBondAcceptorCount": 3,
                        "RotatableBondCount": 2,
                        "HeavyAtomCount": 6,
                        "TPSA": 60.7
                    }
                ]
            }
        }
        self.mock_requests_get.return_value = self.mock_response
        
        # Create a test enhancer with dry run mode
        self.enhancer = PubChemEnhancer(batch_size=2, dry_run=True, checkpoint_file="test_checkpoint.json")
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.patcher1.stop()
        self.patcher2.stop()
        self.patcher3.stop()
        self.patcher4.stop()
        self.patcher5.stop()
        
        # Remove test checkpoint file if it exists
        if os.path.exists("test_checkpoint.json"):
            os.remove("test_checkpoint.json")
    
    def test_find_molecules_to_enhance(self):
        """Test finding molecules to enhance."""
        # Mock the dry_run mode to return test molecules instead of sample data
        with patch.object(self.enhancer, 'dry_run', False):
            molecules = self.enhancer.find_molecules_to_enhance()
            
            self.assertEqual(molecules, self.test_molecules)
            self.mock_execute_query.assert_called_once()
    
    @patch('enhance_pubchem_properties.os.path.exists')
    @patch('builtins.open', new_callable=mock_open, read_data='{"processed_molecules": ["mol-uuid-1"]}')
    def test_load_checkpoint(self, mock_file, mock_exists):
        """Test loading checkpoint data."""
        mock_exists.return_value = True
        
        enhancer = PubChemEnhancer(checkpoint_file="test_checkpoint.json")
        
        self.assertEqual(enhancer.processed_molecules, {"mol-uuid-1"})
        mock_exists.assert_called_once_with("test_checkpoint.json")
    
    @patch('enhance_pubchem_properties.os.makedirs')
    @patch('enhance_pubchem_properties.json.dump')
    @patch('builtins.open', new_callable=mock_open)
    def test_save_checkpoint(self, mock_file, mock_json_dump, mock_makedirs):
        """Test saving checkpoint data."""
        self.enhancer.processed_molecules = {"mol-uuid-1", "mol-uuid-2"}
        self.enhancer._save_checkpoint()
        
        mock_makedirs.assert_called_once()
        mock_file.assert_called_once_with("test_checkpoint.json", 'w')
        mock_json_dump.assert_called_once()
        
        # Verify the JSON content was passed to json.dump
        checkpoint_data = mock_json_dump.call_args[0][0]
        
        self.assertIn("processed_molecules", checkpoint_data)
        self.assertIn("timestamp", checkpoint_data)
        self.assertEqual(set(checkpoint_data["processed_molecules"]), {"mol-uuid-1", "mol-uuid-2"})
    
    def test_fetch_pubchem_properties(self):
        """Test fetching properties from PubChem API."""
        # Set up successful response
        self.mock_response.ok = True
        
        result = self.enhancer._fetch_pubchem_properties("753")
        
        self.assertEqual(result, self.mock_response.json()["PropertyTable"]["Properties"][0])
        self.mock_requests_get.assert_called_once()
        self.mock_rate_limiter.return_value.wait.assert_called_once()
        
    @patch('enhance_pubchem_properties.time.sleep')
    def test_fetch_pubchem_properties_with_retry(self, mock_sleep):
        """Test retry logic when fetching properties from PubChem API."""
        # Set up a failure followed by success
        fail_response = MagicMock()
        fail_response.ok = False
        fail_response.status_code = 500
        
        success_response = MagicMock()
        success_response.ok = True
        success_response.json.return_value = self.mock_response.json()
        
        # First call fails, second call succeeds
        self.mock_requests_get.side_effect = [fail_response, success_response]
        
        result = self.enhancer._fetch_pubchem_properties("753")
        
        # Should return the successful result
        self.assertEqual(result, self.mock_response.json()["PropertyTable"]["Properties"][0])
        
        # Should have called requests.get twice
        self.assertEqual(self.mock_requests_get.call_count, 2)
        
        # Should have called rate_limiter.wait twice (once per attempt)
        self.assertEqual(self.mock_rate_limiter.return_value.wait.call_count, 2)
        
        # Should have called sleep once for the retry
        mock_sleep.assert_called_once()
    
    def test_enhance_molecule_dry_run(self):
        """Test enhancing a molecule in dry run mode."""
        molecule = self.test_molecules[0]
        success, is_enhanced = self.enhancer._enhance_molecule(molecule)
        
        self.assertTrue(success)
        self.assertTrue(is_enhanced)
        self.assertIn(molecule[0], self.enhancer.processed_molecules)
        self.mock_property_manager.return_value.set_properties.assert_not_called()
    
    @patch('enhance_pubchem_properties.execute_query')
    def test_enhance_molecule_with_db_update(self, mock_execute):
        """Test enhancing a molecule with database update."""
        # Create a non-dry-run enhancer
        enhancer = PubChemEnhancer(batch_size=2, dry_run=False, checkpoint_file="test_checkpoint.json")
        
        molecule = self.test_molecules[0]
        success, is_enhanced = enhancer._enhance_molecule(molecule)
        
        self.assertTrue(success)
        self.assertTrue(is_enhanced)
        self.assertIn(molecule[0], enhancer.processed_molecules)
        self.mock_property_manager.return_value.set_properties.assert_called_once()
        mock_execute.assert_called_once()
    
    def test_enhance_molecule_api_failure(self):
        """Test handling API failure when enhancing a molecule."""
        # Set up API failure
        self.mock_response.ok = False
        self.mock_response.status_code = 500
        
        molecule = self.test_molecules[0]
        success, is_enhanced = self.enhancer._enhance_molecule(molecule)
        
        self.assertFalse(success)
        self.assertFalse(is_enhanced)
        self.mock_property_manager.return_value.set_properties.assert_not_called()
    
    @patch('enhance_pubchem_properties.concurrent.futures.ThreadPoolExecutor')
    def test_process_batch(self, mock_executor):
        """Test processing a batch of molecules."""
        # Set up mock executor
        mock_executor.return_value.__enter__.return_value.submit.side_effect = lambda func, arg: MagicMock(
            result=lambda: (True, True)
        )
        mock_executor.return_value.__enter__.return_value.submit.return_value.result = lambda: (True, True)
        
        # Mock as_completed to return the futures in order
        with patch('enhance_pubchem_properties.concurrent.futures.as_completed', 
                  side_effect=lambda futures: futures):
            
            self.enhancer.process_batch(self.test_molecules[:2])
        
        self.assertEqual(self.enhancer.results['molecules_enhanced'], 2)
        self.assertEqual(len(self.enhancer.results['details']), 2)
    
    @patch.object(PubChemEnhancer, 'find_molecules_to_enhance')
    @patch.object(PubChemEnhancer, 'process_batch')
    @patch('enhance_pubchem_properties.process_in_batches')
    def test_enhance_properties(self, mock_process_batches, mock_process_batch, mock_find):
        """Test the main enhance_properties method."""
        mock_find.return_value = self.test_molecules
        
        # Call the method
        results = self.enhancer.enhance_properties()
        
        # Verify results
        self.assertEqual(results['total_molecules'], 3)
        mock_find.assert_called_once()
        mock_process_batches.assert_called_once()

if __name__ == '__main__':
    unittest.main()
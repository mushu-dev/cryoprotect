#!/usr/bin/env python3
"""
Test suite for operational safeguards in ChEMBL data population.

This module tests the rate limiting, batch processing, and memory safeguards
implemented in the ChEMBL_CryoProtectants_Supabase.py script.
"""

import os
import sys
import unittest
import tempfile
import json
import time
from unittest.mock import patch, MagicMock, call

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the modules to test
from chembl.rate_limiter import AdaptiveRateLimiter
from ChEMBL_CryoProtectants_Supabase import process_batches, verify_configuration

class TestRateLimiter(unittest.TestCase):
    """Test the AdaptiveRateLimiter class."""

    def test_rate_limiter_initialization(self):
        """Test that the rate limiter initializes with correct parameters."""
        limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=5.0,
            weekend_requests_per_second=5.0,
            min_delay=0.1,
            max_delay=5.0,
            memory_threshold_percent=80.0
        )
        
        # Check default values
        self.assertEqual(limiter.weekday_delay, 0.2)  # 1.0 / 5.0
        self.assertEqual(limiter.weekend_delay, 0.2)  # 1.0 / 5.0
        self.assertEqual(limiter.min_delay, 0.1)
        self.assertEqual(limiter.max_delay, 5.0)
        self.assertEqual(limiter.memory_threshold_percent, 80.0)
        self.assertEqual(limiter.memory_slowdowns, 0)

    @patch('time.sleep')
    @patch('time.time')
    def test_rate_limiter_wait(self, mock_time, mock_sleep):
        """Test that the rate limiter enforces delays between requests."""
        # Mock time.time() to return consistent values
        mock_time.side_effect = [0, 0.1, 0.1]  # Initial, elapsed time check, after wait
        
        limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=5.0,
            weekend_requests_per_second=5.0
        )
        
        # First call should set last_request_time but not sleep
        limiter.wait()
        mock_sleep.assert_not_called()
        
        # Second call should sleep for the remaining delay
        limiter.wait()
        # Should sleep for 0.2 - 0.1 = 0.1 seconds
        mock_sleep.assert_called_once_with(0.1)

    @patch('psutil.Process')
    @patch('time.sleep')
    def test_memory_threshold_slows_requests(self, mock_sleep, mock_process):
        """Test that high memory usage increases delay."""
        # Mock process.memory_percent to return high memory usage
        process_mock = MagicMock()
        process_mock.memory_percent.return_value = 85.0  # Above threshold
        mock_process.return_value = process_mock
        
        limiter = AdaptiveRateLimiter(
            weekday_requests_per_second=5.0,
            memory_threshold_percent=80.0
        )
        
        # Wait should increase delay due to high memory usage
        limiter.wait()
        
        # Check that memory_slowdowns was incremented
        self.assertEqual(limiter.memory_slowdowns, 1)
        
        # Delay should be increased by a factor based on memory usage
        # 85% is 5% above threshold, so factor should be 1.0 + (5/20) = 1.25
        # Original delay is 0.2, so new delay should be around 0.25
        self.assertGreater(limiter.current_delay, 0.2)


class TestBatchProcessing(unittest.TestCase):
    """Test the batch processing functionality."""

    @patch('ChEMBL_CryoProtectants_Supabase.get_molecule_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.filter_molecule')
    @patch('ChEMBL_CryoProtectants_Supabase.check_molecule_exists')
    @patch('ChEMBL_CryoProtectants_Supabase.get_additional_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.score_molecule')
    @patch('ChEMBL_CryoProtectants_Supabase.bulk_insert_molecules')
    @patch('ChEMBL_CryoProtectants_Supabase.bulk_insert_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.write_checkpoint')
    @patch('ChEMBL_CryoProtectants_Supabase.fetch_property_types')
    @patch('ChEMBL_CryoProtectants_Supabase.supabase')
    @patch('ChEMBL_CryoProtectants_Supabase.CHEMBL_API_DELAY', 0)  # Skip delays for testing
    def test_batch_processing(self, mock_supabase, mock_fetch_property_types, 
                             mock_write_checkpoint, mock_bulk_insert_properties,
                             mock_bulk_insert_molecules, mock_score_molecule,
                             mock_get_additional_properties, mock_check_molecule_exists,
                             mock_filter_molecule, mock_get_molecule_properties):
        """Test that batch processing works with the correct batch size."""
        # Setup mocks
        mock_supabase.auth.current_user = MagicMock()
        mock_supabase.auth.current_user.id = "test-user-id"
        
        mock_fetch_property_types.return_value = [
            {"id": 1, "name": "LogP", "data_type": "numeric"},
            {"id": 2, "name": "TPSA", "data_type": "numeric"},
            {"id": 3, "name": "H-Bond Donors", "data_type": "numeric"},
            {"id": 4, "name": "H-Bond Acceptors", "data_type": "numeric"},
            {"id": 5, "name": "Toxicity", "data_type": "text"},
            {"id": 6, "name": "Stability", "data_type": "text"},
            {"id": 7, "name": "Environmental Safety", "data_type": "text"},
            {"id": 8, "name": "Total Score", "data_type": "numeric"},
            {"id": 9, "name": "ChEMBL ID", "data_type": "text"}
        ]
        
        mock_get_molecule_properties.return_value = {
            "ChEMBL ID": "CHEMBL123",
            "Name": "Test Molecule",
            "Molecular Formula": "C10H20O2",
            "Molecular Weight": "172.26",
            "LogP": "2.5",
            "TPSA": "40.5",
            "H-Bond Donors": "1",
            "H-Bond Acceptors": "2",
            "SMILES": "CCCCCCCCCC(=O)O",
            "InChI": "InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)",
            "InChIKey": "KJTXGVHSILVNPO-UHFFFAOYSA-N"
        }
        
        mock_filter_molecule.return_value = True
        mock_check_molecule_exists.return_value = None
        mock_get_additional_properties.return_value = {
            "Toxicity": None,
            "Stability": None,
            "Environmental Safety": None
        }
        mock_score_molecule.return_value = 150
        mock_bulk_insert_molecules.return_value = ["mol-id-1", "mol-id-2", "mol-id-3"]
        mock_bulk_insert_properties.return_value = True
        mock_write_checkpoint.return_value = True
        
        # Create a temporary directory for checkpoints
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test with a small batch size for testing
            batch_size = 3
            chembl_ids = ["CHEMBL123", "CHEMBL456", "CHEMBL789", "CHEMBL101", "CHEMBL202"]
            
            # Call the function
            process_batches(
                chembl_ids=chembl_ids,
                batch_size=batch_size,
                checkpoint_path=temp_dir,
                resume=False,
                reset=True,
                skipped_log_path=os.path.join(temp_dir, "skipped.log"),
                checkpoint_frequency=1,
                memory_check_frequency=1  # Check memory for every molecule
            )
            
            # Verify batch processing
            # Should process in 2 batches: [CHEMBL123, CHEMBL456, CHEMBL789] and [CHEMBL101, CHEMBL202]
            self.assertEqual(mock_get_molecule_properties.call_count, 5)
            self.assertEqual(mock_bulk_insert_molecules.call_count, 2)
            self.assertEqual(mock_write_checkpoint.call_count, 2)  # One for each batch

    @patch('psutil.Process')
    @patch('gc.collect')
    @patch('ChEMBL_CryoProtectants_Supabase.get_molecule_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.filter_molecule')
    @patch('ChEMBL_CryoProtectants_Supabase.check_molecule_exists')
    @patch('ChEMBL_CryoProtectants_Supabase.get_additional_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.score_molecule')
    @patch('ChEMBL_CryoProtectants_Supabase.bulk_insert_molecules')
    @patch('ChEMBL_CryoProtectants_Supabase.bulk_insert_properties')
    @patch('ChEMBL_CryoProtectants_Supabase.write_checkpoint')
    @patch('ChEMBL_CryoProtectants_Supabase.fetch_property_types')
    @patch('ChEMBL_CryoProtectants_Supabase.supabase')
    @patch('ChEMBL_CryoProtectants_Supabase.CHEMBL_API_DELAY', 0)  # Skip delays for testing
    def test_memory_monitoring(self, mock_supabase, mock_fetch_property_types, 
                              mock_write_checkpoint, mock_bulk_insert_properties,
                              mock_bulk_insert_molecules, mock_score_molecule,
                              mock_get_additional_properties, mock_check_molecule_exists,
                              mock_filter_molecule, mock_get_molecule_properties,
                              mock_gc_collect, mock_process):
        """Test that memory monitoring triggers garbage collection when memory usage is high."""
        # Setup mocks
        mock_supabase.auth.current_user = MagicMock()
        mock_supabase.auth.current_user.id = "test-user-id"
        
        # Mock process.memory_percent to return high memory usage
        process_mock = MagicMock()
        process_mock.memory_percent.return_value = 90.0  # Very high memory usage
        process_mock.memory_info.return_value = MagicMock()
        process_mock.memory_info.return_value.rss = 1000000000  # 1GB
        mock_process.return_value = process_mock
        
        mock_fetch_property_types.return_value = [
            {"id": 1, "name": "LogP", "data_type": "numeric"}
        ]
        
        mock_get_molecule_properties.return_value = {
            "ChEMBL ID": "CHEMBL123",
            "Name": "Test Molecule",
            "Molecular Formula": "C10H20O2",
            "Molecular Weight": "172.26",
            "LogP": "2.5",
            "SMILES": "CCCCCCCCCC(=O)O",
            "InChI": "InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)",
            "InChIKey": "KJTXGVHSILVNPO-UHFFFAOYSA-N"
        }
        
        mock_filter_molecule.return_value = True
        mock_check_molecule_exists.return_value = None
        mock_get_additional_properties.return_value = {}
        mock_score_molecule.return_value = 150
        mock_bulk_insert_molecules.return_value = ["mol-id-1", "mol-id-2"]
        mock_bulk_insert_properties.return_value = True
        mock_write_checkpoint.return_value = True
        
        # Create a temporary directory for checkpoints
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test with a small batch size for testing
            batch_size = 2
            chembl_ids = ["CHEMBL123", "CHEMBL456"]
            
            # Call the function
            process_batches(
                chembl_ids=chembl_ids,
                batch_size=batch_size,
                checkpoint_path=temp_dir,
                resume=False,
                reset=True,
                skipped_log_path=os.path.join(temp_dir, "skipped.log"),
                checkpoint_frequency=1,
                memory_check_frequency=1  # Check memory for every molecule
            )
            
            # Verify that garbage collection was called due to high memory usage
            self.assertTrue(mock_gc_collect.called)


if __name__ == '__main__':
    unittest.main()
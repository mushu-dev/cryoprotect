#!/usr/bin/env python3
"""
Unit tests for progress tracking and resumability in ChEMBL_Integrated_Import.py
"""

import unittest
from unittest.mock import patch, MagicMock, call
import json
import os
import tempfile
from datetime import datetime
import shutil

# Import the module to test
import ChEMBL_Integrated_Import as chembl_import
from supabase_direct import SupabaseDirectConnection


class TestChEMBLProgressTracking(unittest.TestCase):
    """Test cases for progress tracking and resumability in ChEMBL_Integrated_Import.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a temporary directory for checkpoints
        self.temp_dir = tempfile.mkdtemp()
        self.checkpoint_file = os.path.join(self.temp_dir, "test_checkpoint.json")
        self.progress_file = os.path.join(self.temp_dir, "test_progress.json")
        
        # Create a mock for SupabaseDirectConnection
        self.mock_db = MagicMock(spec=SupabaseDirectConnection)
        
        # Mock the get_db_connection function to return our mock
        patcher = patch.object(chembl_import, 'get_db_connection', return_value=self.mock_db)
        self.addCleanup(patcher.stop)
        self.mock_get_db = patcher.start()
        
        # Mock user profile ID
        self.user_profile_id = 'test-user-profile-id'
        
        # Mock get_user_id and ensure_user_profile
        patcher = patch.object(chembl_import, 'get_user_id', return_value='test-user-id')
        self.addCleanup(patcher.stop)
        self.mock_get_user_id = patcher.start()
        
        patcher = patch.object(chembl_import, 'ensure_user_profile', return_value=self.user_profile_id)
        self.addCleanup(patcher.stop)
        self.mock_ensure_user_profile = patcher.start()
        
        # Create mock compounds
        self.mock_compounds = []
        for i in range(5):
            self.mock_compounds.append({
                'molecule_chembl_id': f'CHEMBL{1000+i}',
                'pref_name': f'Test Compound {i}',
                'molecule_structures': {
                    'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                    'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                    'standard_inchi_key': f'INCHIKEY{1000+i}'
                },
                'molecule_properties': {
                    'full_molformula': 'C9H8O4',
                    'full_mwt': 180.16,
                    'alogp': 1.31,
                    'hba': 4,
                    'hbd': 1,
                    'psa': 63.6,
                    'rtb': 3
                }
            })

    def tearDown(self):
        """Clean up test fixtures"""
        # Remove temporary directory
        shutil.rmtree(self.temp_dir)

    def test_progress_tracker_initialization(self):
        """Test that the progress tracker initializes correctly"""
        # Create a progress tracker
        progress_tracker = chembl_import.ChEMBLProgressTracker(
            total_compounds=100,
            batch_size=10,
            checkpoint_file=self.progress_file
        )
        
        # Verify initial state
        self.assertEqual(progress_tracker.total_compounds, 100)
        self.assertEqual(progress_tracker.batch_size, 10)
        self.assertEqual(progress_tracker.total_batches, 10)
        self.assertEqual(progress_tracker.total_processed, 0)
        self.assertEqual(progress_tracker.total_imported, 0)
        self.assertEqual(progress_tracker.total_skipped, 0)
        self.assertEqual(progress_tracker.total_errors, 0)
        self.assertEqual(progress_tracker.status, "Running")
        
        # Verify progress data
        progress_data = progress_tracker.get_progress_data()
        self.assertEqual(progress_data["total_compounds"], 100)
        self.assertEqual(progress_data["progress_percentage"], 0)
        self.assertEqual(progress_data["status"], "Running")

    def test_progress_tracker_batch_processing(self):
        """Test that the progress tracker correctly tracks batch processing"""
        # Create a progress tracker
        progress_tracker = chembl_import.ChEMBLProgressTracker(
            total_compounds=100,
            batch_size=10,
            checkpoint_file=self.progress_file
        )
        
        # Process a batch
        progress_tracker.start_batch(0, 10)
        progress_tracker.end_batch(10, 8, 1, 1)
        
        # Verify state after batch
        self.assertEqual(progress_tracker.total_processed, 10)
        self.assertEqual(progress_tracker.total_imported, 8)
        self.assertEqual(progress_tracker.total_skipped, 1)
        self.assertEqual(progress_tracker.total_errors, 1)
        self.assertEqual(progress_tracker.current_batch, 1)
        
        # Verify progress data
        progress_data = progress_tracker.get_progress_data()
        self.assertEqual(progress_data["total_processed"], 10)
        self.assertEqual(progress_data["progress_percentage"], 10)
        
        # Process another batch
        progress_tracker.start_batch(1, 10)
        progress_tracker.end_batch(10, 9, 0, 1)
        
        # Verify state after second batch
        self.assertEqual(progress_tracker.total_processed, 20)
        self.assertEqual(progress_tracker.total_imported, 17)
        self.assertEqual(progress_tracker.total_skipped, 1)
        self.assertEqual(progress_tracker.total_errors, 2)
        self.assertEqual(progress_tracker.current_batch, 2)
        
        # Verify progress data
        progress_data = progress_tracker.get_progress_data()
        self.assertEqual(progress_data["total_processed"], 20)
        self.assertEqual(progress_data["progress_percentage"], 20)

    def test_checkpoint_manager_save_load(self):
        """Test that the checkpoint manager can save and load checkpoints"""
        # Create a checkpoint manager
        checkpoint_manager = chembl_import.ChEMBLCheckpointManager(self.checkpoint_file)
        
        # Save a checkpoint
        compounds = self.mock_compounds[:2]
        next_index = 2
        status = "Running"
        
        result = checkpoint_manager.save(compounds, next_index, status)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(self.checkpoint_file))
        
        # Load the checkpoint
        checkpoint = checkpoint_manager.load()
        self.assertIsNotNone(checkpoint)
        self.assertEqual(len(checkpoint["compounds"]), 2)
        self.assertEqual(checkpoint["next_index"], 2)
        self.assertEqual(checkpoint["status"], "Running")
        
        # Save another checkpoint
        compounds = self.mock_compounds[:3]
        next_index = 3
        status = "Running"
        
        result = checkpoint_manager.save(compounds, next_index, status)
        self.assertTrue(result)
        
        # Load the updated checkpoint
        checkpoint = checkpoint_manager.load()
        self.assertIsNotNone(checkpoint)
        self.assertEqual(len(checkpoint["compounds"]), 3)
        self.assertEqual(checkpoint["next_index"], 3)

    def test_checkpoint_manager_backup_recovery(self):
        """Test that the checkpoint manager can recover from backup"""
        # Create a checkpoint manager
        checkpoint_manager = chembl_import.ChEMBLCheckpointManager(self.checkpoint_file)
        
        # Save a checkpoint
        compounds = self.mock_compounds[:2]
        next_index = 2
        status = "Running"
        
        result = checkpoint_manager.save(compounds, next_index, status)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(self.checkpoint_file))
        self.assertTrue(os.path.exists(f"{self.checkpoint_file}.bak"))
        
        # Corrupt the main checkpoint file
        with open(self.checkpoint_file, 'w') as f:
            f.write("corrupted data")
        
        # Load the checkpoint (should recover from backup)
        checkpoint = checkpoint_manager.load()
        self.assertIsNotNone(checkpoint)
        self.assertEqual(len(checkpoint["compounds"]), 2)
        self.assertEqual(checkpoint["next_index"], 2)
        self.assertEqual(checkpoint["status"], "Running")

    @patch('ChEMBL_Integrated_Import.transform_chembl_to_molecule')
    @patch('ChEMBL_Integrated_Import.transform_chembl_to_properties')
    def test_import_with_progress_tracking(self, mock_transform_properties, mock_transform_molecule):
        """Test that import_compounds_to_database uses progress tracking"""
        # Mock transform_chembl_to_molecule to return a valid molecule
        mock_transform_molecule.return_value = {
            'name': 'Test Molecule',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
            'inchikey': 'INCHIKEY1000',
            'formula': 'C9H8O4',
            'molecular_weight': 180.16,
            'created_by': self.user_profile_id,
            'data_source': 'ChEMBL Test',
            'version': 1
        }
        
        # Mock transform_chembl_to_properties to return a list of properties
        mock_transform_properties.return_value = [
            {
                'molecule_id': 'test-molecule-id',
                'property_type_id': 'test-property-type-id',
                'numeric_value': 1.31,
                'text_value': None,
                'boolean_value': None,
                'created_by': self.user_profile_id,
                'data_source': 'ChEMBL Test',
                'version': 1
            }
        ]
        
        # Mock database operations
        self.mock_db.execute_query.side_effect = [
            # get_property_types
            [{'id': 'prop-type-1', 'name': 'logp'}],
            # check_molecule_exists (first compound)
            [],
            # insert_molecule (first compound)
            [{'id': 'new-molecule-id-1'}],
            # check_molecule_exists (second compound)
            [],
            # insert_molecule (second compound)
            [{'id': 'new-molecule-id-2'}]
        ]
        
        # Mock execute_batch to return True (success)
        self.mock_db.execute_batch.return_value = True
        
        # Call the function with progress tracking
        result = chembl_import.import_compounds_to_database(
            self.mock_compounds[:2],
            batch_size=1,
            dry_run=False,
            checkpoint_file=self.progress_file
        )
        
        # Verify the result
        self.assertEqual(result["total_compounds"], 2)
        self.assertEqual(result["processed"], 2)
        self.assertTrue("progress_data" in result)
        
        # Verify progress data
        progress_data = result["progress_data"]
        self.assertEqual(progress_data["total_compounds"], 2)
        self.assertEqual(progress_data["total_processed"], 2)
        self.assertEqual(progress_data["status"], "Completed")


if __name__ == '__main__':
    unittest.main()
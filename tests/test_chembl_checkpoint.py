import unittest
import os
import shutil
import tempfile
import json
import time
from datetime import datetime, timedelta
from chembl.checkpoint import CheckpointManager


class TestCheckpointManager(unittest.TestCase):
    """Test cases for ChEMBL checkpoint system."""

    def setUp(self):
        """Set up a temporary directory for checkpoint files"""
        self.temp_dir = tempfile.mkdtemp()
        self.checkpoint_manager = CheckpointManager(self.temp_dir)
        
    def tearDown(self):
        """Remove temporary directory after tests"""
        shutil.rmtree(self.temp_dir)
        
    def test_initialization(self):
        """Test checkpoint manager initialization"""
        self.assertEqual(self.checkpoint_manager.checkpoint_dir, self.temp_dir)
        self.assertEqual(self.checkpoint_manager.prefix, "chembl_import")
        
        # Check that the state is initialized correctly
        self.assertIn("processed_compounds", self.checkpoint_manager.state)
        self.assertIn("successful_compounds", self.checkpoint_manager.state)
        self.assertIn("failed_compounds", self.checkpoint_manager.state)
        self.assertIn("current_batch", self.checkpoint_manager.state)
        self.assertIn("total_processed", self.checkpoint_manager.state)
        self.assertIn("success_count", self.checkpoint_manager.state)
        self.assertIn("error_count", self.checkpoint_manager.state)
        self.assertIn("start_time", self.checkpoint_manager.state)
        self.assertIn("last_updated", self.checkpoint_manager.state)
        self.assertIn("config", self.checkpoint_manager.state)
        
        # Check initial values
        self.assertEqual(self.checkpoint_manager.state["processed_compounds"], [])
        self.assertEqual(self.checkpoint_manager.state["successful_compounds"], [])
        self.assertEqual(self.checkpoint_manager.state["failed_compounds"], {})
        self.assertEqual(self.checkpoint_manager.state["current_batch"], 0)
        self.assertEqual(self.checkpoint_manager.state["total_processed"], 0)
        self.assertEqual(self.checkpoint_manager.state["success_count"], 0)
        self.assertEqual(self.checkpoint_manager.state["error_count"], 0)
        
    def test_custom_prefix(self):
        """Test initialization with custom prefix"""
        custom_prefix = "custom_prefix"
        manager = CheckpointManager(self.temp_dir, prefix=custom_prefix)
        self.assertEqual(manager.prefix, custom_prefix)
        
    def test_get_checkpoint_path(self):
        """Test getting checkpoint path"""
        path = self.checkpoint_manager.get_checkpoint_path()
        self.assertTrue(path.startswith(os.path.join(self.temp_dir, "chembl_import_")))
        self.assertTrue(path.endswith(".json"))
        
    def test_save_load_checkpoint(self):
        """Test saving and loading checkpoints"""
        # Update state
        self.checkpoint_manager.state["total_processed"] = 10
        self.checkpoint_manager.state["success_count"] = 8
        self.checkpoint_manager.state["error_count"] = 2
        
        # Save checkpoint
        checkpoint_path = self.checkpoint_manager.save_checkpoint()
        self.assertTrue(os.path.exists(checkpoint_path))
        
        # Create new manager and load checkpoint
        new_manager = CheckpointManager(self.temp_dir)
        success = new_manager.load_checkpoint()
        
        # Verify state was loaded
        self.assertTrue(success)
        self.assertEqual(new_manager.state["total_processed"], 10)
        self.assertEqual(new_manager.state["success_count"], 8)
        self.assertEqual(new_manager.state["error_count"], 2)
        
    def test_load_nonexistent_checkpoint(self):
        """Test loading when no checkpoint exists"""
        # Create a new manager in an empty directory
        empty_dir = tempfile.mkdtemp()
        try:
            manager = CheckpointManager(empty_dir)
            success = manager.load_checkpoint()
            self.assertFalse(success)
        finally:
            shutil.rmtree(empty_dir)
            
    def test_get_latest_checkpoint(self):
        """Test getting the latest checkpoint"""
        # No checkpoints initially
        self.assertIsNone(self.checkpoint_manager.get_latest_checkpoint())
        
        # Create multiple checkpoints
        self.checkpoint_manager.state["total_processed"] = 5
        first_path = self.checkpoint_manager.save_checkpoint()
        
        time.sleep(0.1)  # Ensure different modification times
        
        self.checkpoint_manager.state["total_processed"] = 10
        second_path = self.checkpoint_manager.save_checkpoint()
        
        # Latest should be the second one
        latest = self.checkpoint_manager.get_latest_checkpoint()
        self.assertEqual(latest, second_path)
        
    def test_cleanup_old_checkpoints(self):
        """Test cleanup of old checkpoints"""
        # Create 10 checkpoints
        for i in range(10):
            self.checkpoint_manager.state["total_processed"] = i
            self.checkpoint_manager.save_checkpoint()
            time.sleep(0.1)  # Ensure different modification times
            
        # Should only keep 5 most recent
        checkpoint_files = [
            f for f in os.listdir(self.temp_dir)
            if f.startswith("chembl_import_") and f.endswith(".json")
        ]
        self.assertEqual(len(checkpoint_files), 5)
        
    def test_update_progress(self):
        """Test updating progress for compounds"""
        # Update progress for a successful compound
        self.checkpoint_manager.update_progress("CHEMBL123", success=True)
        
        # Check state updates
        self.assertIn("CHEMBL123", self.checkpoint_manager.state["processed_compounds"])
        self.assertIn("CHEMBL123", self.checkpoint_manager.state["successful_compounds"])
        self.assertEqual(self.checkpoint_manager.state["total_processed"], 1)
        self.assertEqual(self.checkpoint_manager.state["success_count"], 1)
        self.assertEqual(self.checkpoint_manager.state["error_count"], 0)
        
        # Update progress for a failed compound
        error_msg = "API error"
        self.checkpoint_manager.update_progress("CHEMBL456", success=False, error=error_msg)
        
        # Check state updates
        self.assertIn("CHEMBL456", self.checkpoint_manager.state["processed_compounds"])
        self.assertNotIn("CHEMBL456", self.checkpoint_manager.state["successful_compounds"])
        self.assertEqual(self.checkpoint_manager.state["failed_compounds"]["CHEMBL456"], error_msg)
        self.assertEqual(self.checkpoint_manager.state["total_processed"], 2)
        self.assertEqual(self.checkpoint_manager.state["success_count"], 1)
        self.assertEqual(self.checkpoint_manager.state["error_count"], 1)
        
        # Update an already processed compound
        self.checkpoint_manager.update_progress("CHEMBL123", success=True)
        
        # Counts should not change for already processed compounds
        self.assertEqual(self.checkpoint_manager.state["total_processed"], 2)
        self.assertEqual(self.checkpoint_manager.state["success_count"], 1)
        
    def test_batch_processing(self):
        """Test batch processing methods"""
        # Start a batch
        batch_num = 1
        batch_size = 5
        self.checkpoint_manager.start_batch(batch_num, batch_size)
        
        # Check batch initialization
        self.assertEqual(self.checkpoint_manager.state["current_batch"], batch_num)
        self.assertEqual(self.checkpoint_manager.state["current_batch_size"], batch_size)
        self.assertIn("current_batch_start", self.checkpoint_manager.state)
        self.assertEqual(self.checkpoint_manager.state["current_batch_processed"], 0)
        
        # Process some compounds
        self.checkpoint_manager.update_progress("CHEMBL1", success=True)
        self.checkpoint_manager.update_progress("CHEMBL2", success=True)
        self.checkpoint_manager.update_progress("CHEMBL3", success=False, error="Error")
        
        # End the batch
        self.checkpoint_manager.end_batch(processed=3, success=2, errors=1)
        
        # Check batch statistics
        self.assertIn("batches", self.checkpoint_manager.state)
        self.assertIn(str(batch_num), self.checkpoint_manager.state["batches"])
        
        batch_stats = self.checkpoint_manager.state["batches"][str(batch_num)]
        self.assertEqual(batch_stats["processed"], 3)
        self.assertEqual(batch_stats["success"], 2)
        self.assertEqual(batch_stats["errors"], 1)
        self.assertIn("start_time", batch_stats)
        self.assertIn("end_time", batch_stats)
        self.assertIn("duration_seconds", batch_stats)
        self.assertIn("items_per_second", batch_stats)
        
        # Check performance metrics
        self.assertIn("performance", self.checkpoint_manager.state)
        self.assertEqual(self.checkpoint_manager.state["performance"]["batch_count"], 1)
        self.assertGreaterEqual(self.checkpoint_manager.state["performance"]["total_duration_seconds"], 0)
        
    def test_elapsed_time_calculation(self):
        """Test elapsed time calculation during checkpoint save"""
        # Set a start time in the past
        past_time = datetime.now() - timedelta(minutes=5)
        self.checkpoint_manager.state["start_time"] = past_time.isoformat()
        
        # Save checkpoint
        self.checkpoint_manager.save_checkpoint()
        
        # Check elapsed time calculation
        self.assertIn("elapsed_seconds", self.checkpoint_manager.state)
        self.assertGreaterEqual(self.checkpoint_manager.state["elapsed_seconds"], 300)  # At least 5 minutes
        
    def test_invalid_json_handling(self):
        """Test handling of invalid JSON during load"""
        # Create an invalid JSON file
        invalid_path = os.path.join(self.temp_dir, "chembl_import_invalid.json")
        with open(invalid_path, 'w') as f:
            f.write("This is not valid JSON")
        
        # Attempt to load it
        original_get_latest = self.checkpoint_manager.get_latest_checkpoint
        try:
            # Mock get_latest_checkpoint to return our invalid file
            self.checkpoint_manager.get_latest_checkpoint = lambda: invalid_path
            success = self.checkpoint_manager.load_checkpoint()
            self.assertFalse(success)
        finally:
            # Restore original method
            self.checkpoint_manager.get_latest_checkpoint = original_get_latest


if __name__ == '__main__':
    unittest.main()
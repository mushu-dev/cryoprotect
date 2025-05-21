"""
Test cases for the checkpoint system.

This module tests the functionality of the checkpoint system that enables
resumable imports of molecular data from various sources.
"""

import os
import json
import tempfile
import shutil
import pytest
import time
import logging
from typing import Dict, Any, Set

from ..core.checkpoint import CheckpointManager


class TestCheckpointManager:
    """Test the checkpoint manager functionality."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for checkpoint files."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def checkpoint_file(self, temp_dir):
        """Create a temporary checkpoint file path."""
        return os.path.join(temp_dir, "test_checkpoint.json")
    
    @pytest.fixture
    def logger(self):
        """Set up a logger for testing."""
        logger = logging.getLogger("test_checkpoint")
        logger.setLevel(logging.DEBUG)
        return logger
    
    @pytest.fixture
    def checkpoint_manager(self, checkpoint_file, logger):
        """Create a checkpoint manager for testing."""
        return CheckpointManager(checkpoint_file, backup_interval=2, logger=logger)
    
    def test_initialization(self, checkpoint_manager):
        """Test checkpoint manager initialization."""
        # Check initial state
        assert checkpoint_manager.data['progress']['total_items'] == 0
        assert checkpoint_manager.data['progress']['processed_items'] == 0
        assert checkpoint_manager.data['progress']['successful_items'] == 0
        assert checkpoint_manager.data['progress']['failed_items'] == 0
        assert checkpoint_manager.data['progress']['skipped_items'] == 0
        
        # Check ID tracking sets
        assert isinstance(checkpoint_manager.data['processed_ids'], set)
        assert isinstance(checkpoint_manager.data['failed_ids'], set)
        assert isinstance(checkpoint_manager.data['skipped_ids'], set)
        assert len(checkpoint_manager.data['processed_ids']) == 0
        assert len(checkpoint_manager.data['failed_ids']) == 0
        assert len(checkpoint_manager.data['skipped_ids']) == 0
        
        # Check metadata
        assert 'created_at' in checkpoint_manager.data['metadata']
        assert 'last_updated' in checkpoint_manager.data['metadata']
        assert checkpoint_manager.data['metadata']['save_count'] == 0
        assert checkpoint_manager.data['metadata']['version'] == '1.0'
    
    def test_save_and_load(self, checkpoint_manager, checkpoint_file):
        """Test saving and loading checkpoint data."""
        # Mark some items
        checkpoint_manager.mark_processed("item1")
        checkpoint_manager.mark_failed("item2")
        checkpoint_manager.mark_skipped("item3")
        
        # Set total items and position
        checkpoint_manager.set_total_items(100)
        checkpoint_manager.set_position({"page": 5, "offset": 50})
        
        # Save custom state
        checkpoint_manager.set_custom_state("source", "ChEMBL")
        checkpoint_manager.set_custom_state("query", "aspirin")
        
        # Save the checkpoint
        checkpoint_manager.save()
        
        # Check that the file exists
        assert os.path.exists(checkpoint_file)
        
        # Create a new checkpoint manager to load the data
        new_manager = CheckpointManager(checkpoint_file)
        
        # Check loaded data
        assert new_manager.is_processed("item1")
        assert new_manager.is_failed("item2")
        assert new_manager.is_skipped("item3")
        assert new_manager.data['progress']['total_items'] == 100
        assert new_manager.get_position() == {"page": 5, "offset": 50}
        assert new_manager.get_custom_state("source") == "ChEMBL"
        assert new_manager.get_custom_state("query") == "aspirin"
    
    def test_backup_creation(self, checkpoint_manager, checkpoint_file):
        """Test backup file creation."""
        backup_file = f"{checkpoint_file}.bak"
        
        # Save multiple times to trigger backup
        for i in range(3):  # Backup interval is 2
            checkpoint_manager.mark_processed(f"item{i}")
            checkpoint_manager.save()
        
        # Check that the backup file exists
        assert os.path.exists(backup_file)
        
        # Verify backup content
        with open(backup_file, 'r') as f:
            backup_data = json.load(f)
            assert len(backup_data['processed_ids']) > 0
    
    def test_backup_recovery(self, checkpoint_manager, checkpoint_file, temp_dir):
        """Test recovery from backup when main file is corrupted."""
        # Mark some items and save
        checkpoint_manager.mark_processed("item1")
        checkpoint_manager.mark_processed("item2")
        checkpoint_manager.save()
        checkpoint_manager.save()  # Second save to create backup
        
        # Corrupt the main file
        with open(checkpoint_file, 'w') as f:
            f.write("invalid json")
        
        # Create a new checkpoint manager, should recover from backup
        recovered_manager = CheckpointManager(checkpoint_file)
        
        # Check recovery
        assert recovered_manager.is_processed("item1")
        assert recovered_manager.is_processed("item2")
    
    def test_marking_items(self, checkpoint_manager):
        """Test marking items as processed, failed, or skipped."""
        # Mark an item as processed
        checkpoint_manager.mark_processed("item1")
        assert checkpoint_manager.is_processed("item1")
        assert not checkpoint_manager.is_failed("item1")
        assert not checkpoint_manager.is_skipped("item1")
        assert checkpoint_manager.data['progress']['processed_items'] == 1
        assert checkpoint_manager.data['progress']['successful_items'] == 1
        
        # Mark an item as failed
        checkpoint_manager.mark_failed("item2")
        assert not checkpoint_manager.is_processed("item2")
        assert checkpoint_manager.is_failed("item2")
        assert not checkpoint_manager.is_skipped("item2")
        assert checkpoint_manager.data['progress']['processed_items'] == 2
        assert checkpoint_manager.data['progress']['failed_items'] == 1
        
        # Mark an item as skipped
        checkpoint_manager.mark_skipped("item3")
        assert not checkpoint_manager.is_processed("item3")
        assert not checkpoint_manager.is_failed("item3")
        assert checkpoint_manager.is_skipped("item3")
        assert checkpoint_manager.data['progress']['skipped_items'] == 1
        
        # Test changing item state
        checkpoint_manager.mark_failed("item1")  # Was processed, now failed
        assert not checkpoint_manager.is_processed("item1")
        assert checkpoint_manager.is_failed("item1")
        assert checkpoint_manager.data['progress']['successful_items'] == 0
        assert checkpoint_manager.data['progress']['failed_items'] == 2
        
        checkpoint_manager.mark_processed("item2")  # Was failed, now processed
        assert checkpoint_manager.is_processed("item2")
        assert not checkpoint_manager.is_failed("item2")
        assert checkpoint_manager.data['progress']['successful_items'] == 1
        assert checkpoint_manager.data['progress']['failed_items'] == 1
    
    def test_batch_management(self, checkpoint_manager):
        """Test batch management functionality."""
        # Set a batch
        batch = ["item1", "item2", "item3", "item4"]
        checkpoint_manager.update_batch(batch)
        
        # Mark some items as processed
        checkpoint_manager.mark_processed("item1")
        checkpoint_manager.mark_processed("item3")
        
        # Get unprocessed items from batch
        unprocessed = checkpoint_manager.get_unprocessed_from_batch()
        assert len(unprocessed) == 2
        assert "item2" in unprocessed
        assert "item4" in unprocessed
    
    def test_custom_state(self, checkpoint_manager):
        """Test custom state storage and retrieval."""
        # Store custom state values
        checkpoint_manager.set_custom_state("source", "PubChem")
        checkpoint_manager.set_custom_state("query", "glucose")
        checkpoint_manager.set_custom_state("options", {"limit": 100, "batch_size": 10})
        
        # Retrieve custom state
        assert checkpoint_manager.get_custom_state("source") == "PubChem"
        assert checkpoint_manager.get_custom_state("query") == "glucose"
        assert checkpoint_manager.get_custom_state("options") == {"limit": 100, "batch_size": 10}
        
        # Test default value
        assert checkpoint_manager.get_custom_state("missing_key", "default") == "default"
    
    def test_progress_tracking(self, checkpoint_manager):
        """Test progress tracking functionality."""
        # Set total items
        checkpoint_manager.set_total_items(100)
        
        # Mark some items
        for i in range(30):
            checkpoint_manager.mark_processed(f"success{i}")
        
        for i in range(10):
            checkpoint_manager.mark_failed(f"fail{i}")
        
        for i in range(5):
            checkpoint_manager.mark_skipped(f"skip{i}")
        
        # Check counts
        assert checkpoint_manager.get_processed_count() == 30
        assert checkpoint_manager.get_failed_count() == 10
        assert checkpoint_manager.get_skipped_count() == 5
        
        # Check progress data
        progress = checkpoint_manager.get_progress()
        assert progress['processed_items'] == 40  # 30 success + 10 failed
        assert progress['successful_items'] == 30
        assert progress['failed_items'] == 10
        assert progress['skipped_items'] == 5
        assert progress['total_items'] == 100
    
    def test_summary_generation(self, checkpoint_manager):
        """Test summary generation functionality."""
        # Set up some data
        checkpoint_manager.set_total_items(200)
        
        for i in range(50):
            checkpoint_manager.mark_processed(f"item{i}")
        
        for i in range(10):
            checkpoint_manager.mark_failed(f"fail{i}")
        
        checkpoint_manager.set_position({"page": 3})
        
        # Get summary
        summary = checkpoint_manager.get_summary()
        
        # Check summary content
        assert summary['total_items'] == 200
        assert summary['processed_items'] == 60  # 50 success + 10 failed
        assert summary['successful_items'] == 50
        assert summary['failed_items'] == 10
        assert summary['completion_percentage'] == 30.0  # (60/200)*100
        assert summary['current_position'] == {"page": 3}
    
    def test_clear(self, checkpoint_manager):
        """Test clearing checkpoint data."""
        # Add some data
        checkpoint_manager.mark_processed("item1")
        checkpoint_manager.set_total_items(100)
        checkpoint_manager.set_custom_state("source", "ChEMBL")
        
        # Clear the data
        checkpoint_manager.clear()
        
        # Check that data is reset
        assert not checkpoint_manager.is_processed("item1")
        assert checkpoint_manager.data['progress']['total_items'] == 0
        assert checkpoint_manager.data['progress']['processed_items'] == 0
        assert checkpoint_manager.get_custom_state("source") is None
    
    def test_integer_ids(self, checkpoint_manager):
        """Test handling of integer item IDs."""
        # Mark items with integer IDs
        checkpoint_manager.mark_processed(1)
        checkpoint_manager.mark_failed(2)
        checkpoint_manager.mark_skipped(3)
        
        # Check that the IDs are properly converted to strings and tracked
        assert checkpoint_manager.is_processed("1")
        assert checkpoint_manager.is_failed("2")
        assert checkpoint_manager.is_skipped("3")
        
        # Check directly with integer arguments
        assert checkpoint_manager.is_processed(1)
        assert checkpoint_manager.is_failed(2)
        assert checkpoint_manager.is_skipped(3)


class TestRealWorldUsage:
    """Test real-world usage scenarios for the checkpoint system."""
    
    @pytest.fixture
    def temp_checkpoint_file(self):
        """Create a temporary checkpoint file for testing."""
        # Create a temporary file
        fd, path = tempfile.mkstemp()
        os.close(fd)
        
        yield path
        
        # Clean up temporary file
        if os.path.exists(path):
            os.unlink(path)
    
    def test_resumable_import_simulation(self, temp_checkpoint_file):
        """Simulate a resumable import process."""
        # Create initial checkpoint manager
        manager = CheckpointManager(temp_checkpoint_file)
        
        # Set up import parameters
        manager.set_total_items(1000)
        manager.set_custom_state("source", "ChEMBL")
        manager.set_custom_state("query", "cryoprotectant")
        manager.set_position({"offset": 0})
        
        # Process first batch
        batch1 = [f"item{i}" for i in range(1, 101)]
        manager.update_batch(batch1)
        
        # Mark some as processed, some as failed
        for i in range(1, 81):
            manager.mark_processed(f"item{i}")
        
        for i in range(81, 91):
            manager.mark_failed(f"item{i}")
        
        for i in range(91, 101):
            manager.mark_skipped(f"item{i}")
        
        # Update position and save
        manager.set_position({"offset": 100})
        manager.save()
        
        # Simulate process restart by creating new manager
        resumed_manager = CheckpointManager(temp_checkpoint_file)
        
        # Check that we can continue where we left off
        assert resumed_manager.get_position() == {"offset": 100}
        assert resumed_manager.get_custom_state("source") == "ChEMBL"
        assert resumed_manager.get_custom_state("query") == "cryoprotectant"
        
        # Items from first batch should be marked correctly
        assert resumed_manager.is_processed("item50")
        assert resumed_manager.is_failed("item85")
        assert resumed_manager.is_skipped("item95")
        
        # Process next batch
        batch2 = [f"item{i}" for i in range(101, 201)]
        resumed_manager.update_batch(batch2)
        
        # Check for already processed items (there shouldn't be any in the new batch)
        unprocessed = resumed_manager.get_unprocessed_from_batch()
        assert len(unprocessed) == 100
        
        # Process some more items
        for i in range(101, 151):
            resumed_manager.mark_processed(f"item{i}")
        
        # Update position and save
        resumed_manager.set_position({"offset": 200})
        resumed_manager.save()
        
        # Verify progress
        summary = resumed_manager.get_summary()
        assert summary['successful_items'] == 130  # 80 from first batch + 50 from second
        assert summary['failed_items'] == 10  # All from first batch
        assert summary['skipped_items'] == 10  # All from first batch
        assert summary['completion_percentage'] == 15.0  # (150/1000)*100
    
    def test_error_recovery(self, temp_checkpoint_file):
        """Test recovery from processing errors."""
        # Create checkpoint manager
        manager = CheckpointManager(temp_checkpoint_file)
        
        # Set up import
        manager.set_total_items(100)
        
        # Process items with some failures
        for i in range(1, 31):
            manager.mark_processed(f"item{i}")
        
        for i in range(31, 41):
            manager.mark_failed(f"item{i}")
        
        manager.save()
        
        # Simulate restart and retry failed items
        retry_manager = CheckpointManager(temp_checkpoint_file)
        
        # Get failed items
        failed_ids = retry_manager.get_failed_ids()
        
        # Retry failed items
        for item_id in failed_ids:
            # Simulate successful retry
            retry_manager.mark_processed(item_id)
        
        # Check that all items are now processed
        for i in range(1, 41):
            assert retry_manager.is_processed(f"item{i}")
        
        # Check progress
        assert retry_manager.get_failed_count() == 0
        assert retry_manager.get_processed_count() == 40
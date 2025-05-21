"""
Integration tests for the checkpoint and progress tracking systems.

This module tests the combined functionality of the checkpoint and progress
tracking systems in a realistic import simulation.
"""

import os
import tempfile
import shutil
import asyncio
import pytest
import time
import logging
from typing import Dict, Any, List, Set, Tuple

from ..core.checkpoint import CheckpointManager
from ..core.progress import ProgressTracker


class MockDataSource:
    """
    Mock data source for testing import processes.
    
    Simulates importing data with controlled success, failure, and delay rates.
    """
    
    def __init__(
        self,
        data_items: List[str],
        checkpoint_manager: CheckpointManager,
        progress_tracker: ProgressTracker,
        success_rate: float = 0.9,
        failure_rate: float = 0.05,
        skip_rate: float = 0.05,
        delay_per_item: float = 0.01
    ):
        """
        Initialize mock data source.
        
        Args:
            data_items: List of item IDs to process
            checkpoint_manager: Checkpoint manager instance
            progress_tracker: Progress tracker instance
            success_rate: Percentage of items to succeed (0.0-1.0)
            failure_rate: Percentage of items to fail (0.0-1.0)
            skip_rate: Percentage of items to skip (0.0-1.0)
            delay_per_item: Processing delay per item in seconds
        """
        self.data_items = data_items
        self.checkpoint_manager = checkpoint_manager
        self.progress_tracker = progress_tracker
        self.success_rate = success_rate
        self.failure_rate = failure_rate
        self.skip_rate = skip_rate
        self.delay_per_item = delay_per_item
        
        # Set total items in tracker and checkpoint manager
        self.progress_tracker.total_items = len(data_items)
        self.checkpoint_manager.set_total_items(len(data_items))
        
        # Counters for verification
        self.successful_count = 0
        self.failed_count = 0
        self.skipped_count = 0
    
    async def process_item(self, item_id: str) -> Tuple[bool, str, Dict[str, Any]]:
        """
        Process a single data item.
        
        Args:
            item_id: ID of the item to process
            
        Returns:
            Tuple of (success, error_message, result_data)
        """
        # Add some delay to simulate processing
        await asyncio.sleep(self.delay_per_item)
        
        # Determine outcome based on configured rates
        outcome = self._determine_outcome()
        
        if outcome == "success":
            self.successful_count += 1
            self.checkpoint_manager.mark_processed(item_id)
            self.progress_tracker.update(processed=1, successful=1)
            return True, None, {"id": item_id, "status": "success"}
            
        elif outcome == "fail":
            self.failed_count += 1
            self.checkpoint_manager.mark_failed(item_id)
            self.progress_tracker.update(processed=1, failed=1, errors={"processing_error": 1})
            return False, "Simulated processing error", None
            
        else:  # skip
            self.skipped_count += 1
            self.checkpoint_manager.mark_skipped(item_id)
            self.progress_tracker.update(skipped=1)
            return True, None, {"id": item_id, "status": "skipped"}
    
    def _determine_outcome(self) -> str:
        """
        Determine the outcome of processing based on configured rates.
        
        Returns:
            Outcome string: "success", "fail", or "skip"
        """
        import random
        
        # Generate a random value between 0 and 1
        r = random.random()
        
        if r < self.success_rate:
            return "success"
        elif r < self.success_rate + self.failure_rate:
            return "fail"
        else:
            return "skip"
    
    async def process_batch(self, batch: List[str]) -> Tuple[int, int, List[Tuple[str, str]]]:
        """
        Process a batch of data items.
        
        Args:
            batch: List of item IDs to process
            
        Returns:
            Tuple of (success_count, failure_count, failures_with_reasons)
        """
        # Update checkpoint with current batch
        self.checkpoint_manager.update_batch(batch)
        
        # Create tasks for each item
        tasks = []
        for item_id in batch:
            # Skip already processed items
            if self.checkpoint_manager.is_processed(item_id) or self.checkpoint_manager.is_failed(item_id) or self.checkpoint_manager.is_skipped(item_id):
                continue
                
            tasks.append(self.process_item(item_id))
        
        # Process items concurrently
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # Process results
        success_count = 0
        failure_count = 0
        failures = []
        
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                failure_count += 1
                if i < len(batch):
                    failures.append((batch[i], str(result)))
            else:
                success, error, _ = result
                if success:
                    success_count += 1
                else:
                    failure_count += 1
                    if i < len(batch) and error:
                        failures.append((batch[i], error))
        
        # Save checkpoint after batch
        self.checkpoint_manager.save()
        
        # Log progress
        self.progress_tracker.log_progress()
        
        return success_count, failure_count, failures
    
    async def process_all(self, batch_size: int = 10) -> Dict[str, Any]:
        """
        Process all data items in batches.
        
        Args:
            batch_size: Number of items to process in each batch
            
        Returns:
            Results summary
        """
        # Process in batches
        for i in range(0, len(self.data_items), batch_size):
            batch = self.data_items[i:i + batch_size]
            await self.process_batch(batch)
            
            # Update position after each batch
            self.checkpoint_manager.set_position(i + batch_size)
        
        # Save final checkpoint and report
        self.checkpoint_manager.save()
        report = self.progress_tracker.generate_report()
        
        return {
            "total_items": len(self.data_items),
            "successful": self.successful_count,
            "failed": self.failed_count,
            "skipped": self.skipped_count,
            "processed": self.successful_count + self.failed_count,
            "tracking_report": report
        }


@pytest.mark.asyncio
class TestCheckpointProgressIntegration:
    """Integration tests for checkpoint and progress systems."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def checkpoint_file(self, temp_dir):
        """Create a checkpoint file path."""
        return os.path.join(temp_dir, "test_checkpoint.json")
    
    @pytest.fixture
    def report_file(self, temp_dir):
        """Create a report file path."""
        return os.path.join(temp_dir, "progress_report.json")
    
    @pytest.fixture
    def logger(self):
        """Set up a logger for testing."""
        logger = logging.getLogger("test_integration")
        logger.setLevel(logging.DEBUG)
        return logger
    
    @pytest.fixture
    def checkpoint_manager(self, checkpoint_file, logger):
        """Create a checkpoint manager."""
        return CheckpointManager(checkpoint_file, backup_interval=5, logger=logger)
    
    @pytest.fixture
    def progress_tracker(self, logger):
        """Create a progress tracker."""
        return ProgressTracker(window_size=20, logger=logger)
    
    @pytest.fixture
    def test_data(self):
        """Generate test data items."""
        return [f"item{i}" for i in range(100)]
    
    @pytest.fixture
    def mock_source(self, test_data, checkpoint_manager, progress_tracker):
        """Create a mock data source."""
        return MockDataSource(
            test_data,
            checkpoint_manager,
            progress_tracker,
            success_rate=0.8,
            failure_rate=0.1,
            skip_rate=0.1,
            delay_per_item=0.01
        )
    
    async def test_full_import_process(self, mock_source, checkpoint_manager, progress_tracker, report_file):
        """Test a full import process with checkpointing and progress tracking."""
        # Process all items
        results = await mock_source.process_all(batch_size=20)
        
        # Save final report
        progress_tracker.save_report(report_file)
        
        # Verify results
        assert results["total_items"] == 100
        assert results["processed"] + results["skipped"] == 100
        
        # Check that stats from progress tracker match
        assert progress_tracker.successful_items == results["successful"]
        assert progress_tracker.failed_items == results["failed"]
        assert progress_tracker.skipped_items == results["skipped"]
        
        # Check that stats from checkpoint manager match
        assert checkpoint_manager.get_processed_count() == results["successful"]
        assert checkpoint_manager.get_failed_count() == results["failed"]
        assert checkpoint_manager.get_skipped_count() == results["skipped"]
        
        # Verify the report file exists
        assert os.path.exists(report_file)
    
    async def test_interrupted_import_resume(self, test_data, checkpoint_manager, progress_tracker, report_file):
        """Test resuming an interrupted import process."""
        # Create first data source
        first_source = MockDataSource(
            test_data,
            checkpoint_manager,
            progress_tracker,
            success_rate=0.8,
            failure_rate=0.1,
            skip_rate=0.1,
            delay_per_item=0.01
        )
        
        # Process first half of the data
        for i in range(0, 50, 10):
            batch = test_data[i:i + 10]
            await first_source.process_batch(batch)
            checkpoint_manager.set_position(i + 10)
        
        # Save checkpoint
        checkpoint_manager.save()
        
        # Check progress
        first_half_stats = {
            "successful": first_source.successful_count,
            "failed": first_source.failed_count,
            "skipped": first_source.skipped_count
        }
        
        # Create a new data source to resume from checkpoint
        second_source = MockDataSource(
            test_data,
            checkpoint_manager,
            progress_tracker,
            success_rate=0.8,
            failure_rate=0.1,
            skip_rate=0.1,
            delay_per_item=0.01
        )
        
        # Resume from the current position (50)
        current_position = checkpoint_manager.get_position()
        
        # Process second half of the data
        for i in range(current_position, 100, 10):
            batch = test_data[i:i + 10]
            await second_source.process_batch(batch)
            checkpoint_manager.set_position(i + 10)
        
        # Save final report
        progress_tracker.save_report(report_file)
        
        # Verify combined results
        combined_successful = first_half_stats["successful"] + second_source.successful_count
        combined_failed = first_half_stats["failed"] + second_source.failed_count
        combined_skipped = first_half_stats["skipped"] + second_source.skipped_count
        
        # Check that progress tracker has the combined results
        assert progress_tracker.successful_items == combined_successful
        assert progress_tracker.failed_items == combined_failed
        assert progress_tracker.skipped_items == combined_skipped
        assert progress_tracker.processed_items == combined_successful + combined_failed
    
    async def test_retry_failed_items(self, test_data, checkpoint_manager, progress_tracker):
        """Test retrying failed items."""
        # Create first data source with high failure rate
        first_source = MockDataSource(
            test_data,
            checkpoint_manager,
            progress_tracker,
            success_rate=0.6,
            failure_rate=0.3,  # High failure rate
            skip_rate=0.1,
            delay_per_item=0.01
        )
        
        # Process all items
        await first_source.process_all(batch_size=20)
        
        # Get the ids of failed items
        failed_ids = list(checkpoint_manager.get_failed_ids())
        assert len(failed_ids) > 0  # Should have some failures
        
        # Save failure count
        first_failed_count = len(failed_ids)
        
        # Create a recovery data source with higher success rate
        recovery_source = MockDataSource(
            failed_ids,  # Only process the failed items
            checkpoint_manager,
            progress_tracker,
            success_rate=0.9,  # Higher success rate for retry
            failure_rate=0.1,
            skip_rate=0.0,
            delay_per_item=0.01
        )
        
        # Reset progress tracker counts for retry operation
        original_processed = progress_tracker.processed_items
        original_successful = progress_tracker.successful_items
        original_failed = progress_tracker.failed_items
        original_skipped = progress_tracker.skipped_items
        
        # Adjust progress tracker to only count the retries
        recovery_tracker = ProgressTracker(total_items=len(failed_ids))
        
        # Process all failed items
        await recovery_source.process_all(batch_size=10)
        
        # Get the new count of failed items
        new_failed_ids = checkpoint_manager.get_failed_ids()
        
        # Should have fewer failures after retry
        assert len(new_failed_ids) < first_failed_count
        
        # Success count should have increased
        assert recovery_source.successful_count > 0
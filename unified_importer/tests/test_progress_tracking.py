"""
Unit tests for the progress tracking system.
"""

import unittest
import os
import time
import tempfile
import json
import logging
from unittest.mock import Mock, patch

from unified_importer.core.progress_tracking import (
    ProgressTracker, ProgressStats, ReportFormat,
    format_time_duration, format_timestamp
)

class TestProgressTracking(unittest.TestCase):
    """Tests for the progress tracking system."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger("test")
        self.logger.setLevel(logging.DEBUG)
    
    def test_progress_stats_initialization(self):
        """Test initialization of ProgressStats."""
        stats = ProgressStats(total_items=100)
        
        self.assertEqual(stats.total_items, 100)
        self.assertEqual(stats.processed_items, 0)
        self.assertEqual(stats.success_count, 0)
        self.assertEqual(stats.error_count, 0)
        self.assertEqual(stats.warning_count, 0)
        self.assertEqual(stats.skipped_count, 0)
        self.assertEqual(stats.error_categories, {})
        self.assertIsNotNone(stats.start_time)
        self.assertIsNone(stats.end_time)
        self.assertIsNotNone(stats.last_update_time)
        self.assertEqual(stats.elapsed_time, 0.0)
        self.assertEqual(stats.items_per_second, 0.0)
        self.assertEqual(stats.average_process_time, 0.0)
        self.assertEqual(stats.recent_process_times, [])
        self.assertIsNone(stats.estimated_completion_time)
        self.assertEqual(stats.estimated_time_remaining, 0.0)
        self.assertIsNone(stats.estimated_completion_datetime)
        self.assertIsNotNone(stats.eta_last_calculated)
        self.assertEqual(stats.eta_smoothing_factor, 0.2)
        self.assertEqual(stats.eta_samples, [])
        self.assertEqual(stats.eta_max_samples, 10)
        self.assertEqual(stats.eta_variance, 0.0)
        self.assertEqual(stats.checkpoint_times, [])
    
    def test_progress_stats_update_elapsed_time(self):
        """Test updating elapsed time in ProgressStats."""
        stats = ProgressStats(total_items=100)
        
        # Set a specific start time for testing
        stats.start_time = time.time() - 10.0  # 10 seconds ago
        
        # Update elapsed time
        stats.update_elapsed_time()
        
        # Elapsed time should be approximately 10 seconds
        self.assertAlmostEqual(stats.elapsed_time, 10.0, delta=0.1)
        
        # Set end time and check again
        stats.end_time = time.time()
        stats.update_elapsed_time()
        
        # Elapsed time should now be based on end time
        self.assertAlmostEqual(stats.elapsed_time, 10.0, delta=0.2)
    
    def test_progress_stats_update_processing_rate(self):
        """Test updating processing rate in ProgressStats."""
        stats = ProgressStats(total_items=100)
        
        # Set start time and processed items
        stats.start_time = time.time() - 10.0  # 10 seconds ago
        stats.processed_items = 50
        
        # Update elapsed time and processing rate
        stats.update_elapsed_time()
        stats.update_processing_rate()
        
        # Items per second should be approximately 5
        self.assertAlmostEqual(stats.items_per_second, 5.0, delta=0.5)
        
        # Average process time should be approximately 0.2 seconds
        self.assertAlmostEqual(stats.average_process_time, 0.2, delta=0.02)
    
    def test_progress_stats_calculate_eta(self):
        """Test ETA calculation in ProgressStats."""
        stats = ProgressStats(total_items=100)
        
        # Set start time and processed items
        stats.start_time = time.time() - 10.0  # 10 seconds ago
        stats.processed_items = 50
        
        # Calculate ETA
        stats.calculate_eta()
        
        # ETA should be approximately another 10 seconds
        self.assertIsNotNone(stats.estimated_completion_time)
        self.assertIsNotNone(stats.estimated_completion_datetime)
        self.assertAlmostEqual(stats.estimated_time_remaining, 10.0, delta=2.0)
        
        # Completion time should be approximately 10 seconds from now
        current_time = time.time()
        self.assertAlmostEqual(stats.estimated_completion_time, current_time + 10.0, delta=2.0)
    
    def test_progress_stats_add_checkpoint(self):
        """Test adding checkpoints in ProgressStats."""
        stats = ProgressStats(total_items=100)
        
        # Add a checkpoint
        stats.add_checkpoint("test_checkpoint", {"test": "data"})
        
        # Check that checkpoint was added
        self.assertEqual(len(stats.checkpoint_times), 1)
        self.assertEqual(stats.checkpoint_times[0]["name"], "test_checkpoint")
        self.assertEqual(stats.checkpoint_times[0]["data"], {"test": "data"})
        self.assertIn("time", stats.checkpoint_times[0])
        self.assertIn("elapsed", stats.checkpoint_times[0])
        self.assertEqual(stats.checkpoint_times[0]["processed_items"], 0)
        
        # Add another checkpoint with updated stats
        stats.processed_items = 50
        stats.success_count = 40
        stats.error_count = 10
        stats.add_checkpoint("second_checkpoint")
        
        # Check that second checkpoint was added with updated stats
        self.assertEqual(len(stats.checkpoint_times), 2)
        self.assertEqual(stats.checkpoint_times[1]["name"], "second_checkpoint")
        self.assertEqual(stats.checkpoint_times[1]["processed_items"], 50)
        self.assertEqual(stats.checkpoint_times[1]["success_count"], 40)
        self.assertEqual(stats.checkpoint_times[1]["error_count"], 10)
    
    def test_progress_stats_to_dict(self):
        """Test converting ProgressStats to dictionary."""
        stats = ProgressStats(total_items=100)
        
        # Set some values
        stats.processed_items = 50
        stats.success_count = 40
        stats.error_count = 10
        stats.start_time = time.time() - 10.0  # 10 seconds ago
        
        # Update time-based metrics
        stats.update_elapsed_time()
        stats.update_processing_rate()
        stats.calculate_eta()
        
        # Add a checkpoint
        stats.add_checkpoint("test_checkpoint")
        
        # Convert to dictionary
        data = stats.to_dict()
        
        # Check values
        self.assertEqual(data["total_items"], 100)
        self.assertEqual(data["processed_items"], 50)
        self.assertAlmostEqual(data["percent_complete"], 50.0)
        self.assertEqual(data["success_count"], 40)
        self.assertEqual(data["error_count"], 10)
        self.assertAlmostEqual(data["elapsed_time"], 10.0, delta=0.2)
        self.assertIn("elapsed_time_formatted", data)
        self.assertAlmostEqual(data["items_per_second"], 5.0, delta=0.5)
        self.assertAlmostEqual(data["average_process_time"], 0.2, delta=0.05)
        self.assertIn("start_time", data)
        self.assertIn("start_time_formatted", data)
        
        # Check ETA-related values
        self.assertIn("estimated_time_remaining", data)
        self.assertIn("estimated_time_remaining_formatted", data)
        self.assertIn("estimated_completion_time", data)
        self.assertIn("estimated_completion_time_formatted", data)
        
        # Check checkpoint data
        self.assertIn("checkpoints", data)
        self.assertEqual(len(data["checkpoints"]), 1)
        self.assertEqual(data["checkpoints"][0]["name"], "test_checkpoint")
    
    def test_progress_stats_mark_complete(self):
        """Test marking ProgressStats as complete."""
        stats = ProgressStats(total_items=100)
        stats.processed_items = 100
        
        # Mark as complete
        stats.mark_complete()
        
        # Check values
        self.assertIsNotNone(stats.end_time)
        self.assertEqual(stats.estimated_time_remaining, 0)
        self.assertEqual(stats.estimated_completion_time, stats.end_time)
        self.assertIsNotNone(stats.estimated_completion_datetime)
    
    def test_progress_stats_to_json(self):
        """Test converting ProgressStats to JSON."""
        stats = ProgressStats(total_items=100)
        stats.processed_items = 50
        
        # Convert to JSON
        json_str = stats.to_json()
        
        # Parse and verify
        data = json.loads(json_str)
        self.assertEqual(data["total_items"], 100)
        self.assertEqual(data["processed_items"], 50)
        self.assertAlmostEqual(data["percent_complete"], 50.0)
    
    def test_progress_stats_to_csv(self):
        """Test converting ProgressStats to CSV."""
        stats = ProgressStats(total_items=100)
        stats.processed_items = 50
        
        # Convert to CSV
        csv_str = stats.to_csv()
        
        # Verify format
        self.assertIn("total_items,processed_items,percent_complete", csv_str)
        self.assertIn("100,50,50.0", csv_str)
    
    def test_progress_stats_to_console(self):
        """Test converting ProgressStats to console output."""
        stats = ProgressStats(total_items=100)
        stats.processed_items = 50
        stats.success_count = 40
        stats.error_count = 10
        stats.start_time = time.time() - 10.0  # 10 seconds ago
        
        # Update time-based metrics
        stats.update_elapsed_time()
        stats.update_processing_rate()
        stats.calculate_eta()
        
        # Get minimal console output
        minimal = stats.to_console(detailed=False)
        
        # Verify minimal format
        self.assertIn("50/100", minimal)
        self.assertIn("40 ok", minimal)
        self.assertIn("10 err", minimal)
        
        # Get detailed console output
        detailed = stats.to_console(detailed=True)
        
        # Verify detailed format
        self.assertIn("Progress: 50/100", detailed)
        self.assertIn("Status: 40 successful", detailed)
        self.assertIn("Time:", detailed)
    
    def test_progress_tracker_initialization(self):
        """Test initialization of ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger,
            update_interval=2.0,
            report_format=ReportFormat.CONSOLE
        )
        
        self.assertEqual(tracker.description, "Test Tracker")
        self.assertEqual(tracker.update_interval, 2.0)
        self.assertEqual(tracker.report_format, ReportFormat.CONSOLE)
        self.assertEqual(tracker.stats.total_items, 100)
    
    def test_progress_tracker_start(self):
        """Test starting ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger
        )
        
        # Start tracking
        tracker.start()
        
        # Check initial values
        self.assertEqual(tracker.stats.processed_items, 0)
        self.assertEqual(tracker.stats.success_count, 0)
        self.assertEqual(tracker.stats.error_count, 0)
        self.assertEqual(tracker.stats.warning_count, 0)
        self.assertEqual(tracker.stats.skipped_count, 0)
        self.assertEqual(tracker.stats.error_categories, {})
        
        # Check that start checkpoint was added
        self.assertEqual(len(tracker.stats.checkpoint_times), 1)
        self.assertEqual(tracker.stats.checkpoint_times[0]["name"], "start")
        
        # Test start with updated total items
        tracker.start(total_items=200)
        self.assertEqual(tracker.stats.total_items, 200)
    
    def test_progress_tracker_update(self):
        """Test updating ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger
        )
        
        # Start tracking
        tracker.start()
        
        # Update with success
        tracker.update(increment=10, success=10)
        self.assertEqual(tracker.stats.processed_items, 10)
        self.assertEqual(tracker.stats.success_count, 10)
        
        # Update with errors
        tracker.update(increment=5, error=5, error_category="NETWORK")
        self.assertEqual(tracker.stats.processed_items, 15)
        self.assertEqual(tracker.stats.error_count, 5)
        self.assertEqual(tracker.stats.error_categories["NETWORK"], 5)
        
        # Update with warnings and skipped
        tracker.update(increment=8, warning=3, skipped=5)
        self.assertEqual(tracker.stats.processed_items, 23)
        self.assertEqual(tracker.stats.warning_count, 3)
        self.assertEqual(tracker.stats.skipped_count, 5)
        
        # Force a progress report
        with patch.object(tracker, 'report_progress') as mock_report:
            tracker.update(increment=1, force_report=True)
            mock_report.assert_called_once()
    
    def test_progress_tracker_add_checkpoint(self):
        """Test adding checkpoints in ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger
        )
        
        # Start tracking
        tracker.start()
        
        # Update progress
        tracker.update(increment=25, success=25)
        
        # Add checkpoint
        tracker.add_checkpoint("quarter_done", {"percent": 25})
        
        # Check that checkpoint was added
        self.assertEqual(len(tracker.stats.checkpoint_times), 2)  # start + new checkpoint
        self.assertEqual(tracker.stats.checkpoint_times[1]["name"], "quarter_done")
        self.assertEqual(tracker.stats.checkpoint_times[1]["data"], {"percent": 25})
        self.assertEqual(tracker.stats.checkpoint_times[1]["processed_items"], 25)
    
    def test_progress_tracker_checkpoints_file(self):
        """Test saving and loading checkpoints from file."""
        with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as temp_file:
            checkpoint_path = temp_file.name
        
        try:
            # Create tracker with checkpoint file
            tracker = ProgressTracker(
                total_items=100,
                description="Test Tracker",
                logger=self.logger,
                checkpoint_file=checkpoint_path
            )
            
            # Start tracking
            tracker.start()
            
            # Update progress
            tracker.update(increment=25, success=25)
            
            # Add checkpoint
            tracker.add_checkpoint("quarter_done", {"percent": 25})
            
            # Check that file exists
            self.assertTrue(os.path.exists(checkpoint_path))
            
            # Create new tracker and load checkpoint
            tracker2 = ProgressTracker(
                logger=self.logger,
                checkpoint_file=checkpoint_path
            )
            
            # Load checkpoint
            success = tracker2.load_checkpoint()
            
            # Check that checkpoint was loaded
            self.assertTrue(success)
            self.assertEqual(tracker2.description, "Test Tracker")
            self.assertEqual(tracker2.stats.total_items, 100)
            self.assertEqual(tracker2.stats.processed_items, 25)
            self.assertEqual(tracker2.stats.success_count, 25)
            
            # Check for resume checkpoint
            self.assertEqual(tracker2.stats.checkpoint_times[-1]["name"], "resume")
        finally:
            # Clean up
            if os.path.exists(checkpoint_path):
                os.remove(checkpoint_path)
    
    def test_progress_tracker_report_progress(self):
        """Test progress reporting in ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger
        )
        
        # Start tracking
        tracker.start()
        
        # Update progress
        tracker.update(increment=25, success=25)
        
        # Test different report formats
        console_report = tracker.report_progress(format=ReportFormat.CONSOLE)
        json_report = tracker.report_progress(format=ReportFormat.JSON)
        csv_report = tracker.report_progress(format=ReportFormat.CSV)
        detailed_report = tracker.report_progress(format=ReportFormat.DETAILED)
        
        # Verify formats
        self.assertIn("25/100", console_report)
        self.assertIn("\"total_items\": 100", json_report)
        self.assertIn("total_items,processed_items", csv_report)
        self.assertIn("Progress: 25/100", detailed_report)
    
    def test_progress_tracker_complete(self):
        """Test completing ProgressTracker."""
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger
        )
        
        # Start tracking
        tracker.start()
        
        # Update progress
        tracker.update(increment=100, success=100)
        
        # Complete tracking
        tracker.complete()
        
        # Check that stats were updated
        self.assertIsNotNone(tracker.stats.end_time)
        self.assertEqual(tracker.stats.estimated_time_remaining, 0)
        
        # Check that complete checkpoint was added
        self.assertEqual(tracker.stats.checkpoint_times[-1]["name"], "complete")
    
    def test_format_time_duration(self):
        """Test formatting time durations."""
        # Test microseconds
        self.assertEqual(format_time_duration(0.0000005), "0.50μs")
        
        # Test milliseconds
        self.assertEqual(format_time_duration(0.005), "5.00ms")
        
        # Test seconds
        self.assertEqual(format_time_duration(5), "5.00s")
        
        # Test minutes
        self.assertEqual(format_time_duration(65), "1m 5.0s")
        
        # Test hours
        self.assertEqual(format_time_duration(3665), "1h 1m 5.0s")
        
        # Test days
        self.assertEqual(format_time_duration(90000), "1d 1h 0m")
    
    def test_format_timestamp(self):
        """Test formatting timestamps."""
        # Use a fixed timestamp for testing
        timestamp = 1609459200.0  # 2021-01-01 00:00:00 UTC
        
        formatted = format_timestamp(timestamp)
        
        # Format should be YYYY-MM-DD HH:MM:SS
        self.assertRegex(formatted, r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}')
    
    def test_on_update_callback(self):
        """Test the on_update callback function."""
        # Create a mock callback
        mock_callback = Mock()
        
        tracker = ProgressTracker(
            total_items=100,
            description="Test Tracker",
            logger=self.logger,
            on_update=mock_callback
        )
        
        # Start tracking
        tracker.start()
        
        # Callback should be called on start
        mock_callback.assert_called_once()
        
        # Reset mock
        mock_callback.reset_mock()
        
        # Update progress
        tracker.update(increment=25, success=25)
        
        # Callback should be called on update
        mock_callback.assert_called_once()
        
        # Reset mock
        mock_callback.reset_mock()
        
        # Add checkpoint
        tracker.add_checkpoint("quarter_done")
        
        # Callback should be called on checkpoint
        mock_callback.assert_called_once()
        
        # Reset mock
        mock_callback.reset_mock()
        
        # Complete tracking
        tracker.complete()
        
        # Callback should be called on complete
        mock_callback.assert_called_once()


if __name__ == '__main__':
    unittest.main()
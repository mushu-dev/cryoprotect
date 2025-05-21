"""
Test cases for the progress tracking system.

This module tests the functionality of the progress tracker that monitors
import operations and provides metrics and reporting.
"""

import os
import json
import time
import tempfile
import logging
import pytest
from datetime import datetime, timedelta
from typing import Dict, Any

from ..core.progress import ProgressTracker, ConsoleProgressReporter


class TestProgressTracker:
    """Test the progress tracker functionality."""
    
    @pytest.fixture
    def logger(self):
        """Set up a logger for testing."""
        logger = logging.getLogger("test_progress")
        logger.setLevel(logging.DEBUG)
        return logger
    
    @pytest.fixture
    def tracker(self, logger):
        """Create a progress tracker for testing."""
        return ProgressTracker(total_items=1000, window_size=10, logger=logger)
    
    def test_initialization(self, tracker):
        """Test progress tracker initialization."""
        assert tracker.total_items == 1000
        assert tracker.processed_items == 0
        assert tracker.successful_items == 0
        assert tracker.failed_items == 0
        assert tracker.skipped_items == 0
        assert tracker.window_size == 10
        assert len(tracker.processing_times) == 0
        assert len(tracker.batch_sizes) == 0
        assert tracker.error_counts == {}
    
    def test_update(self, tracker):
        """Test updating progress metrics."""
        # Single update
        tracker.update(processed=10, successful=8, failed=2)
        
        # Check counters
        assert tracker.processed_items == 10
        assert tracker.successful_items == 8
        assert tracker.failed_items == 2
        assert tracker.skipped_items == 0
        
        # Add skipped items
        tracker.update(skipped=5)
        assert tracker.skipped_items == 5
        
        # Add errors
        tracker.update(errors={"API Error": 2, "Timeout": 1})
        assert tracker.error_counts == {"API Error": 2, "Timeout": 1}
        
        # Cumulative updates
        tracker.update(processed=20, successful=18, failed=2)
        assert tracker.processed_items == 30  # 10 + 20
        assert tracker.successful_items == 26  # 8 + 18
        assert tracker.failed_items == 4  # 2 + 2
        
        # Additional errors
        tracker.update(errors={"API Error": 1, "Database Error": 3})
        assert tracker.error_counts == {"API Error": 3, "Timeout": 1, "Database Error": 3}
    
    def test_processing_rate(self, tracker):
        """Test processing rate calculation."""
        # Add some processing data
        for _ in range(5):
            time.sleep(0.01)  # Small delay
            tracker.update(processed=10, successful=10)
        
        # Check that rate is calculated
        rate = tracker.get_processing_rate()
        assert rate > 0  # Should be non-zero
        
        # Check items per minute
        items_per_minute = tracker.get_processing_rate() * 60
        assert items_per_minute > 0
    
    def test_progress_percentage(self, tracker):
        """Test progress percentage calculation."""
        # No progress yet
        assert tracker.get_progress_percentage() == 0
        
        # Some progress
        tracker.update(processed=250)
        assert tracker.get_progress_percentage() == 25.0
        
        # Complete
        tracker.update(processed=750)
        assert tracker.get_progress_percentage() == 100.0
        
        # Test with zero total
        zero_tracker = ProgressTracker(total_items=0)
        assert zero_tracker.get_progress_percentage() == 0
    
    def test_estimated_completion_time(self, tracker):
        """Test completion time estimation."""
        # No progress yet, should return None
        assert tracker.get_estimated_completion_time() is None
        
        # Some progress - need to simulate actual processing time
        start_time = time.time()
        
        # Process 100 items in about 0.1 seconds (very fast for testing)
        time.sleep(0.1)
        tracker.update(processed=100)
        
        # Get estimated time
        eta = tracker.get_estimated_completion_time()
        assert eta is not None
        
        # Should be in the future
        assert eta > time.time()
        
        # Calculate expected time based on processing rate
        elapsed = time.time() - start_time
        rate = 100 / elapsed
        expected_remaining = (1000 - 100) / rate
        expected_eta = time.time() + expected_remaining
        
        # Check that our estimate is close to expected
        assert abs(eta - expected_eta) < 1.0  # Within 1 second
    
    def test_elapsed_time(self, tracker):
        """Test elapsed time calculation."""
        # Get initial elapsed time
        initial_elapsed = tracker.get_elapsed_time()
        assert initial_elapsed >= 0
        
        # Wait a bit
        time.sleep(0.1)
        
        # Get new elapsed time
        new_elapsed = tracker.get_elapsed_time()
        
        # Should be greater than initial
        assert new_elapsed > initial_elapsed
        assert new_elapsed - initial_elapsed >= 0.1
    
    def test_progress_data(self, tracker):
        """Test progress data collection."""
        # Add some progress
        tracker.update(processed=300, successful=280, failed=20)
        
        # Get progress data
        data = tracker.get_progress_data()
        
        # Check data fields
        assert 'timestamp' in data
        assert 'elapsed_seconds' in data
        assert 'elapsed_formatted' in data
        assert 'total_items' in data
        assert 'processed_items' in data
        assert 'successful_items' in data
        assert 'failed_items' in data
        assert 'skipped_items' in data
        assert 'progress_percentage' in data
        assert 'processing_rate' in data
        assert 'items_per_minute' in data
        
        # Check values
        assert data['total_items'] == 1000
        assert data['processed_items'] == 300
        assert data['successful_items'] == 280
        assert data['failed_items'] == 20
        assert data['progress_percentage'] == 30.0
    
    def test_generate_report(self, tracker):
        """Test report generation."""
        # Add some progress with mixed results
        tracker.update(processed=500, successful=400, failed=50, skipped=50)
        tracker.update(errors={"API Error": 30, "Database Error": 20})
        
        # Generate report
        report = tracker.generate_report()
        
        # Check additional report fields
        assert 'success_rate' in report
        assert 'error_rate' in report
        assert 'skipped_rate' in report
        assert 'detailed_error_breakdown' in report
        
        # Check calculations
        assert report['success_rate'] == 80.0  # (400/500)*100
        assert report['error_rate'] == 10.0  # (50/500)*100
        assert report['skipped_rate'] == 10.0  # (50/500)*100
        assert report['detailed_error_breakdown'] == {"API Error": 30, "Database Error": 20}
    
    def test_save_report(self, tracker):
        """Test saving report to file."""
        # Create a temporary file
        fd, path = tempfile.mkstemp()
        os.close(fd)
        
        try:
            # Add some progress
            tracker.update(processed=100, successful=90, failed=10)
            
            # Save report
            tracker.save_report(path)
            
            # Check that file exists
            assert os.path.exists(path)
            
            # Load the file and check content
            with open(path, 'r') as f:
                report = json.load(f)
                
            assert report['total_items'] == 1000
            assert report['processed_items'] == 100
            assert report['successful_items'] == 90
            assert report['failed_items'] == 10
        finally:
            # Clean up
            if os.path.exists(path):
                os.unlink(path)
    
    def test_reset(self, tracker):
        """Test resetting progress metrics."""
        # Add some progress
        tracker.update(processed=200, successful=180, failed=20)
        tracker.update(errors={"API Error": 10, "Timeout": 10})
        
        # Reset tracker
        tracker.reset()
        
        # Check that metrics are reset
        assert tracker.processed_items == 0
        assert tracker.successful_items == 0
        assert tracker.failed_items == 0
        assert tracker.skipped_items == 0
        assert len(tracker.processing_times) == 0
        assert len(tracker.batch_sizes) == 0
        assert tracker.error_counts == {}
        
        # Start time should be updated
        assert tracker.start_time >= time.time() - 1  # Within the last second
    
    def test_callbacks(self, tracker):
        """Test progress callback system."""
        # Create a callback function
        callback_data = []
        
        def test_callback(data):
            callback_data.append(data)
        
        # Register callback
        tracker.register_callback(test_callback)
        
        # Update progress
        tracker.update(processed=10)
        
        # Check that callback was called
        assert len(callback_data) == 1
        assert callback_data[0]['processed_items'] == 10
        
        # Another update
        tracker.update(processed=20)
        assert len(callback_data) == 2
        assert callback_data[1]['processed_items'] == 30
        
        # Unregister callback
        tracker.unregister_callback(test_callback)
        
        # Update again
        tracker.update(processed=10)
        
        # Should still have only 2 callback calls
        assert len(callback_data) == 2


class TestConsoleProgressReporter:
    """Test the console progress reporter functionality."""
    
    @pytest.fixture
    def logger(self):
        """Set up a logger for testing."""
        logger = logging.getLogger("test_reporter")
        logger.setLevel(logging.DEBUG)
        return logger
    
    @pytest.fixture
    def tracker(self, logger):
        """Create a progress tracker for testing."""
        return ProgressTracker(total_items=1000, logger=logger)
    
    def test_reporter_initialization(self, tracker, logger):
        """Test reporter initialization."""
        # Create reporter
        reporter = ConsoleProgressReporter(tracker, update_interval=0.1, logger=logger)
        
        # Check that it's registered with tracker
        assert reporter._progress_callback in tracker.progress_callbacks
        
        # Clean up
        reporter.stop()
        
        # Check that it's unregistered
        assert reporter._progress_callback not in tracker.progress_callbacks
    
    def test_reporter_update_interval(self, tracker, logger):
        """Test reporter update interval."""
        # Create reporter with 1 second interval
        reporter = ConsoleProgressReporter(tracker, update_interval=1.0, logger=logger)
        
        try:
            # Define a mock progress callback counter
            display_count = [0]
            original_display = reporter._display_progress
            
            def mock_display(data):
                display_count[0] += 1
                original_display(data)
            
            reporter._display_progress = mock_display
            
            # Send multiple updates quickly
            for _ in range(5):
                tracker.update(processed=10)
                time.sleep(0.1)
            
            # Should only have displayed once due to interval
            assert display_count[0] <= 2
            
            # Wait longer than interval
            time.sleep(1.0)
            
            # Update again
            tracker.update(processed=10)
            
            # Should have displayed again
            assert display_count[0] <= 3
        finally:
            # Clean up
            reporter.stop()
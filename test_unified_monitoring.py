#!/usr/bin/env python3
"""
Unit tests for the Unified Monitoring Module

This script tests the key functionality of the unified_monitoring module:
1. MonitoringService initialization and singleton pattern
2. Performance metrics recording and retrieval
3. Progress tracking and status reporting
4. Operation tracking with context manager
5. Error recording and retrieval
6. Alert generation and handling
"""

import unittest
import time
import tempfile
import os
import json
import threading
from unittest.mock import patch, MagicMock

# Import the module to test
from unified_monitoring import (
    MonitoringService,
    PerformanceMetrics,
    ProgressTracker,
    track_operation,
    track_database_operation,
    start_monitoring
)

class TestMonitoringService(unittest.TestCase):
    """Tests for the MonitoringService class."""
    
    def setUp(self):
        """Set up test case."""
        # Reset singleton instance for each test
        MonitoringService._instance = None
        
        # Create a monitoring service
        self.monitor = MonitoringService()
    
    def test_singleton_pattern(self):
        """Test that MonitoringService is a singleton."""
        # Get the instance
        instance1 = MonitoringService.get_instance()
        instance2 = MonitoringService.get_instance()
        
        # Verify that they are the same instance
        self.assertIs(instance1, instance2)
        
        # Verify that attempting to create a new instance raises an error
        with self.assertRaises(RuntimeError):
            MonitoringService()
    
    def test_start_stop_monitoring(self):
        """Test starting and stopping monitoring."""
        # Start monitoring
        self.monitor.start()
        
        # Verify that monitoring is active
        self.assertTrue(self.monitor.monitoring_active)
        
        # Verify that monitoring threads are running
        self.assertTrue(len(self.monitor.monitoring_threads) > 0)
        
        # Stop monitoring
        self.monitor.stop()
        
        # Verify that monitoring is inactive
        self.assertFalse(self.monitor.monitoring_active)
        
        # Verify that monitoring threads are stopped
        self.assertEqual(len(self.monitor.monitoring_threads), 0)
    
    def test_record_error(self):
        """Test error recording."""
        # Record an error
        error = self.monitor.record_error(
            "test_error",
            "Test error message",
            "Test stack trace",
            {"test_context": "context_value"}
        )
        
        # Verify that the error was recorded
        self.assertEqual(len(self.monitor.errors), 1)
        
        # Verify error details
        recorded_error = self.monitor.errors[0]
        self.assertEqual(recorded_error["error_type"], "test_error")
        self.assertEqual(recorded_error["message"], "Test error message")
        self.assertEqual(recorded_error["stack_trace"], "Test stack trace")
        self.assertEqual(recorded_error["context"]["test_context"], "context_value")
        
        # Verify error count
        self.assertEqual(self.monitor.error_counts["test_error"], 1)
    
    def test_trigger_alert(self):
        """Test alert triggering."""
        # Create a mock alert handler
        mock_handler = MagicMock()
        
        # Register the mock handler
        self.monitor.add_alert_handler(mock_handler)
        
        # Trigger an alert
        alert = self.monitor._trigger_alert(
            "test_alert",
            "Test alert message",
            "test_source",
            "warning",
            {"test_context": "context_value"}
        )
        
        # Verify that the alert was recorded
        self.assertEqual(len(self.monitor.alerts), 1)
        
        # Verify alert details
        recorded_alert = self.monitor.alerts[0]
        self.assertEqual(recorded_alert["alert_type"], "test_alert")
        self.assertEqual(recorded_alert["message"], "Test alert message")
        self.assertEqual(recorded_alert["source"], "test_source")
        self.assertEqual(recorded_alert["severity"], "warning")
        
        # Verify that the handler was called
        mock_handler.assert_called_once()
        args, _ = mock_handler.call_args
        called_alert = args[0]
        self.assertEqual(called_alert["alert_type"], "test_alert")
    
    def test_track_operation(self):
        """Test operation tracking with context manager."""
        # Use the context manager
        with self.monitor.track_operation("test_operation", items_count=5):
            # Simulate some work
            time.sleep(0.1)
        
        # Verify that the operation was recorded in the appropriate metrics
        api_metrics = self.monitor.performance_metrics["api"].get_metrics()
        self.assertEqual(api_metrics["operations"]["total"], 1)
        self.assertEqual(api_metrics["throughput"]["items_processed"], 5)
    
    def test_silence_alert(self):
        """Test alert silencing."""
        # Silence an alert
        self.monitor.silence_alert("test_alert", 60)
        
        # Verify that the alert is silenced
        self.assertIn("test_alert", self.monitor.alert_silence_periods)
        
        # Trigger the silenced alert
        self.monitor._trigger_alert(
            "test_alert",
            "This should be silenced",
            "test_source"
        )
        
        # Verify that no alert was recorded
        self.assertEqual(len(self.monitor.alerts), 0)
        
        # Unsilence the alert
        self.monitor.unsilence_alert("test_alert")
        
        # Verify that the alert is no longer silenced
        self.assertNotIn("test_alert", self.monitor.alert_silence_periods)
        
        # Trigger the unsilenced alert
        self.monitor._trigger_alert(
            "test_alert",
            "This should not be silenced",
            "test_source"
        )
        
        # Verify that the alert was recorded
        self.assertEqual(len(self.monitor.alerts), 1)


class TestPerformanceMetrics(unittest.TestCase):
    """Tests for the PerformanceMetrics class."""
    
    def setUp(self):
        """Set up test case."""
        # Create a temporary file for metrics
        self.temp_dir = tempfile.mkdtemp()
        self.metrics_file = os.path.join(self.temp_dir, "test_metrics.json")
        
        # Create a performance metrics collector
        self.metrics = PerformanceMetrics("test_metrics", self.metrics_file)
    
    def tearDown(self):
        """Clean up after test case."""
        # Remove temporary files
        if os.path.exists(self.metrics_file):
            os.remove(self.metrics_file)
        os.rmdir(self.temp_dir)
    
    def test_record_operation(self):
        """Test recording operations."""
        # Record a successful operation
        self.metrics.record_operation(
            operation_type="test_query",
            success=True,
            execution_time=0.1,
            items_processed=5
        )
        
        # Verify operation counts
        metrics = self.metrics.get_metrics()
        self.assertEqual(metrics["operations"]["total"], 1)
        self.assertEqual(metrics["operations"]["successful"], 1)
        self.assertEqual(metrics["operations"]["failed"], 0)
        
        # Verify timing metrics
        self.assertEqual(metrics["timing"]["total_execution_time"], 0.1)
        self.assertEqual(metrics["timing"]["min_execution_time"], 0.1)
        self.assertEqual(metrics["timing"]["max_execution_time"], 0.1)
        self.assertEqual(metrics["timing"]["avg_execution_time"], 0.1)
        
        # Verify throughput metrics
        self.assertEqual(metrics["throughput"]["items_processed"], 5)
        
        # Record a failed operation
        self.metrics.record_operation(
            operation_type="test_query",
            success=False,
            execution_time=0.2,
            items_processed=3
        )
        
        # Verify updated metrics
        metrics = self.metrics.get_metrics()
        self.assertEqual(metrics["operations"]["total"], 2)
        self.assertEqual(metrics["operations"]["successful"], 1)
        self.assertEqual(metrics["operations"]["failed"], 1)
        self.assertEqual(metrics["timing"]["total_execution_time"], 0.3)
        self.assertEqual(metrics["timing"]["min_execution_time"], 0.1)
        self.assertEqual(metrics["timing"]["max_execution_time"], 0.2)
        self.assertEqual(metrics["timing"]["avg_execution_time"], 0.15)
        self.assertEqual(metrics["throughput"]["items_processed"], 8)
    
    def test_record_custom_metric(self):
        """Test recording custom metrics."""
        # Record a custom metric
        self.metrics.record_custom_metric("test_metric", 42)
        
        # Verify custom metric
        metrics = self.metrics.get_metrics()
        self.assertEqual(metrics["custom_metrics"]["test_metric"], 42)
    
    def test_metrics_persistence(self):
        """Test metrics persistence."""
        # Record an operation
        self.metrics.record_operation(
            operation_type="test_query",
            success=True,
            execution_time=0.1,
            items_processed=5
        )
        
        # Verify that the metrics file was created
        self.assertTrue(os.path.exists(self.metrics_file))
        
        # Load the metrics file
        with open(self.metrics_file, 'r') as f:
            saved_metrics = json.load(f)
        
        # Verify saved metrics
        self.assertEqual(saved_metrics["name"], "test_metrics")
        self.assertEqual(saved_metrics["operations"]["total"], 1)
        self.assertEqual(saved_metrics["throughput"]["items_processed"], 5)
    
    def test_reset(self):
        """Test resetting metrics."""
        # Record an operation
        self.metrics.record_operation(
            operation_type="test_query",
            success=True,
            execution_time=0.1,
            items_processed=5
        )
        
        # Reset metrics
        self.metrics.reset()
        
        # Verify that metrics are reset
        metrics = self.metrics.get_metrics()
        self.assertEqual(metrics["operations"]["total"], 0)
        self.assertEqual(metrics["operations"]["successful"], 0)
        self.assertEqual(metrics["operations"]["failed"], 0)
        self.assertEqual(metrics["timing"]["total_execution_time"], 0)
        self.assertEqual(metrics["throughput"]["items_processed"], 0)


class TestProgressTracker(unittest.TestCase):
    """Tests for the ProgressTracker class."""
    
    def setUp(self):
        """Set up test case."""
        # Create a temporary file for checkpoint
        self.temp_dir = tempfile.mkdtemp()
        self.checkpoint_file = os.path.join(self.temp_dir, "test_checkpoint.json")
        
        # Create a progress tracker
        self.tracker = ProgressTracker("test_progress", 100, self.checkpoint_file)
    
    def tearDown(self):
        """Clean up after test case."""
        # Remove temporary files
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)
        os.rmdir(self.temp_dir)
    
    def test_update(self):
        """Test progress updates."""
        # Update progress
        status = self.tracker.update(10, successful=9, failed=1, batch_size=10)
        
        # Verify progress status
        self.assertEqual(status["processed_items"], 10)
        self.assertEqual(status["successful_items"], 9)
        self.assertEqual(status["failed_items"], 1)
        self.assertEqual(status["progress_percentage"], 10.0)
        
        # Update progress again
        status = self.tracker.update(20, successful=19, failed=1, batch_size=20)
        
        # Verify updated progress status
        self.assertEqual(status["processed_items"], 30)
        self.assertEqual(status["successful_items"], 28)
        self.assertEqual(status["failed_items"], 2)
        self.assertEqual(status["progress_percentage"], 30.0)
    
    def test_checkpoint(self):
        """Test checkpoint saving and loading."""
        # Update progress
        self.tracker.update(10, successful=9, failed=1, batch_size=10)
        
        # Verify that checkpoint file was created
        self.assertTrue(os.path.exists(self.checkpoint_file))
        
        # Load the checkpoint file
        with open(self.checkpoint_file, 'r') as f:
            checkpoint = json.load(f)
        
        # Verify checkpoint data
        self.assertEqual(checkpoint["processed"], 10)
        self.assertEqual(checkpoint["successful"], 9)
        self.assertEqual(checkpoint["failed"], 1)
        
        # Create a new progress tracker with the same checkpoint file
        new_tracker = ProgressTracker("test_progress", 100, self.checkpoint_file)
        
        # Verify that progress was loaded from checkpoint
        self.assertEqual(new_tracker.processed_items, 10)
        self.assertEqual(new_tracker.successful_items, 9)
        self.assertEqual(new_tracker.failed_items, 1)
    
    def test_reset(self):
        """Test resetting progress."""
        # Update progress
        self.tracker.update(10, successful=9, failed=1, batch_size=10)
        
        # Reset progress
        self.tracker.reset()
        
        # Verify that progress is reset
        self.assertEqual(self.tracker.processed_items, 0)
        self.assertEqual(self.tracker.successful_items, 0)
        self.assertEqual(self.tracker.failed_items, 0)
        
        # Verify that checkpoint file is updated
        with open(self.checkpoint_file, 'r') as f:
            checkpoint = json.load(f)
        
        self.assertEqual(checkpoint["processed"], 0)
        self.assertEqual(checkpoint["successful"], 0)
        self.assertEqual(checkpoint["failed"], 0)


class TestDecorators(unittest.TestCase):
    """Tests for the decorator functions."""
    
    def setUp(self):
        """Set up test case."""
        # Reset singleton instance
        MonitoringService._instance = None
        
        # Create a monitoring service
        self.monitor = MonitoringService()
    
    def test_track_operation_decorator(self):
        """Test the track_operation decorator."""
        # Define a function with the decorator
        @track_operation("test_operation", items_count=5)
        def test_function():
            time.sleep(0.1)
            return 42
        
        # Call the function
        result = test_function()
        
        # Verify that the operation was recorded
        metrics = self.monitor.performance_metrics["api"].get_metrics()
        self.assertEqual(metrics["operations"]["total"], 1)
        self.assertEqual(metrics["throughput"]["items_processed"], 5)
        
        # Verify that the function returned the correct result
        self.assertEqual(result, 42)
    
    def test_track_database_operation_decorator(self):
        """Test the track_database_operation decorator."""
        # Define a function with the decorator
        @track_database_operation(items_count=10)
        def test_db_function():
            time.sleep(0.1)
            return {"result": "success"}
        
        # Call the function
        result = test_db_function()
        
        # Verify that the operation was recorded
        metrics = self.monitor.performance_metrics["database"].get_metrics()
        self.assertEqual(metrics["operations"]["total"], 1)
        self.assertEqual(metrics["throughput"]["items_processed"], 10)
        
        # Verify that the function returned the correct result
        self.assertEqual(result, {"result": "success"})


class TestStartMonitoring(unittest.TestCase):
    """Tests for the start_monitoring function."""
    
    def setUp(self):
        """Set up test case."""
        # Reset singleton instance
        MonitoringService._instance = None
    
    def test_start_monitoring(self):
        """Test the start_monitoring function."""
        # Create a mock Flask app
        mock_app = MagicMock()
        
        # Start monitoring
        monitor = start_monitoring(mock_app, dashboard_port=5002)
        
        # Verify that monitoring is started
        self.assertTrue(monitor.monitoring_active)
        self.assertEqual(monitor.config["dashboard_port"], 5002)
        
        # Verify that the app was initialized
        self.assertEqual(monitor.app, mock_app)
        
        # Clean up
        monitor.stop()


if __name__ == '__main__':
    unittest.main()
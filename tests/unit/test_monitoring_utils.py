#!/usr/bin/env python3
"""
Unit tests for monitoring_utils.py
"""

import os
import time
import json
import unittest
import tempfile
from unittest.mock import patch, MagicMock, mock_open

# Import the module to test
import monitoring_utils
from monitoring_utils import (
    PerformanceMetrics,
    ConnectionHealthMonitor,
    ProgressTracker,
    create_monitoring_directory,
    monitor_connection_health,
    track_progress,
    collect_performance_metrics,
    get_connection_health_status,
    generate_monitoring_report,
    timed_operation
)


class TestPerformanceMetrics(unittest.TestCase):
    """Test cases for PerformanceMetrics class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.metrics_file = os.path.join(self.temp_dir.name, "test_metrics.json")
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_init(self):
        """Test initialization of PerformanceMetrics."""
        metrics = PerformanceMetrics("test_metrics")
        self.assertEqual(metrics.name, "test_metrics")
        self.assertIsNone(metrics.metrics_file)
        
        metrics = PerformanceMetrics("test_metrics", self.metrics_file)
        self.assertEqual(metrics.name, "test_metrics")
        self.assertEqual(metrics.metrics_file, self.metrics_file)
    
    @patch('monitoring_utils.time.time')
    def test_record_operation(self, mock_time):
        """Test recording an operation."""
        mock_time.return_value = 1000.0
        
        metrics = PerformanceMetrics("test_metrics")
        metrics.record_operation(
            operation_type="test",
            success=True,
            execution_time=0.5,
            items_processed=10,
            bytes_processed=1024
        )
        
        self.assertEqual(metrics.metrics["operations"]["total"], 1)
        self.assertEqual(metrics.metrics["operations"]["successful"], 1)
        self.assertEqual(metrics.metrics["operations"]["failed"], 0)
        self.assertEqual(metrics.metrics["timing"]["total_execution_time"], 0.5)
        self.assertEqual(metrics.metrics["timing"]["min_execution_time"], 0.5)
        self.assertEqual(metrics.metrics["timing"]["max_execution_time"], 0.5)
        self.assertEqual(metrics.metrics["timing"]["avg_execution_time"], 0.5)
        self.assertEqual(metrics.metrics["throughput"]["items_processed"], 10)
        self.assertEqual(metrics.metrics["throughput"]["bytes_processed"], 1024)
        
        # Test failed operation
        metrics.record_operation(
            operation_type="test",
            success=False,
            execution_time=1.0,
            items_processed=5,
            bytes_processed=512
        )
        
        self.assertEqual(metrics.metrics["operations"]["total"], 2)
        self.assertEqual(metrics.metrics["operations"]["successful"], 1)
        self.assertEqual(metrics.metrics["operations"]["failed"], 1)
        self.assertEqual(metrics.metrics["timing"]["total_execution_time"], 1.5)
        self.assertEqual(metrics.metrics["timing"]["min_execution_time"], 0.5)
        self.assertEqual(metrics.metrics["timing"]["max_execution_time"], 1.0)
        self.assertEqual(metrics.metrics["timing"]["avg_execution_time"], 0.75)
        self.assertEqual(metrics.metrics["throughput"]["items_processed"], 15)
        self.assertEqual(metrics.metrics["throughput"]["bytes_processed"], 1536)
    
    def test_record_custom_metric(self):
        """Test recording a custom metric."""
        metrics = PerformanceMetrics("test_metrics")
        metrics.record_custom_metric("test_metric", "test_value")
        
        self.assertEqual(metrics.metrics["custom_metrics"]["test_metric"], "test_value")
    
    @patch('monitoring_utils.time.time')
    def test_get_metrics(self, mock_time):
        """Test getting metrics."""
        mock_time.return_value = 1000.0
        
        metrics = PerformanceMetrics("test_metrics")
        metrics.record_operation(
            operation_type="test",
            success=True,
            execution_time=0.5,
            items_processed=10,
            bytes_processed=1024
        )
        
        result = metrics.get_metrics()
        
        self.assertEqual(result["name"], "test_metrics")
        self.assertEqual(result["operations"]["total"], 1)
        self.assertEqual(result["operations"]["successful"], 1)
        self.assertEqual(result["operations"]["failed"], 0)
        self.assertEqual(result["timing"]["total_execution_time"], 0.5)
        self.assertEqual(result["timing"]["min_execution_time"], 0.5)
        self.assertEqual(result["timing"]["max_execution_time"], 0.5)
        self.assertEqual(result["timing"]["avg_execution_time"], 0.5)
        self.assertEqual(result["throughput"]["items_processed"], 10)
        self.assertEqual(result["throughput"]["bytes_processed"], 1024)
    
    def test_generate_report(self):
        """Test generating a report."""
        metrics = PerformanceMetrics("test_metrics")
        metrics.record_operation(
            operation_type="test",
            success=True,
            execution_time=0.5,
            items_processed=10,
            bytes_processed=1024
        )
        
        report = metrics.generate_report()
        
        self.assertIn("Performance Report for test_metrics", report)
        self.assertIn("Operations:", report)
        self.assertIn("Total: 1", report)
        self.assertIn("Successful: 1", report)
        self.assertIn("Failed: 0", report)
        self.assertIn("Timing:", report)
        self.assertIn("Average Execution Time:", report)
        self.assertIn("Throughput:", report)
        self.assertIn("Items Processed: 10", report)
    
    def test_reset(self):
        """Test resetting metrics."""
        metrics = PerformanceMetrics("test_metrics")
        metrics.record_operation(
            operation_type="test",
            success=True,
            execution_time=0.5,
            items_processed=10,
            bytes_processed=1024
        )
        
        metrics.reset()
        
        self.assertEqual(metrics.metrics["operations"]["total"], 0)
        self.assertEqual(metrics.metrics["operations"]["successful"], 0)
        self.assertEqual(metrics.metrics["operations"]["failed"], 0)
        self.assertEqual(metrics.metrics["timing"]["total_execution_time"], 0)
        self.assertEqual(metrics.metrics["timing"]["min_execution_time"], float('inf'))
        self.assertEqual(metrics.metrics["timing"]["max_execution_time"], 0)
        self.assertEqual(metrics.metrics["timing"]["avg_execution_time"], 0)
        self.assertEqual(metrics.metrics["throughput"]["items_processed"], 0)
        self.assertEqual(metrics.metrics["throughput"]["bytes_processed"], 0)


class TestProgressTracker(unittest.TestCase):
    """Test cases for ProgressTracker class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.checkpoint_file = os.path.join(self.temp_dir.name, "test_checkpoint.json")
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_init(self):
        """Test initialization of ProgressTracker."""
        tracker = ProgressTracker("test_tracker", 100)
        self.assertEqual(tracker.name, "test_tracker")
        self.assertEqual(tracker.total_items, 100)
        self.assertIsNone(tracker.checkpoint_file)
        self.assertEqual(tracker.processed_items, 0)
        self.assertEqual(tracker.successful_items, 0)
        self.assertEqual(tracker.failed_items, 0)
        
        # Test with checkpoint file
        with open(self.checkpoint_file, 'w') as f:
            json.dump({"processed": 50}, f)
        
        tracker = ProgressTracker("test_tracker", 100, self.checkpoint_file)
        self.assertEqual(tracker.processed_items, 50)
    
    @patch('monitoring_utils.time.time')
    def test_update(self, mock_time):
        """Test updating progress."""
        mock_time.return_value = 1000.0
        
        tracker = ProgressTracker("test_tracker", 100)
        status = tracker.update(10)
        
        self.assertEqual(tracker.processed_items, 10)
        self.assertEqual(tracker.successful_items, 10)
        self.assertEqual(tracker.failed_items, 0)
        self.assertEqual(status["processed_items"], 10)
        self.assertEqual(status["progress_percentage"], 10.0)
        
        # Test with successful and failed counts
        status = tracker.update(20, 15, 5)
        
        self.assertEqual(tracker.processed_items, 30)
        self.assertEqual(tracker.successful_items, 25)
        self.assertEqual(tracker.failed_items, 5)
        self.assertEqual(status["processed_items"], 30)
        self.assertEqual(status["progress_percentage"], 30.0)
    
    @patch('monitoring_utils.time.time')
    def test_get_status(self, mock_time):
        """Test getting status."""
        mock_time.return_value = 1000.0
        
        tracker = ProgressTracker("test_tracker", 100)
        tracker.update(10)
        
        status = tracker.get_status()
        
        self.assertEqual(status["name"], "test_tracker")
        self.assertEqual(status["total_items"], 100)
        self.assertEqual(status["processed_items"], 10)
        self.assertEqual(status["successful_items"], 10)
        self.assertEqual(status["failed_items"], 0)
        self.assertEqual(status["progress_percentage"], 10.0)
    
    def test_generate_report(self):
        """Test generating a report."""
        tracker = ProgressTracker("test_tracker", 100)
        tracker.update(10)
        
        report = tracker.generate_report()
        
        self.assertIn("Progress Report for test_tracker", report)
        self.assertIn("Status: 10/100 items processed", report)
        self.assertIn("Time:", report)
        self.assertIn("Performance:", report)
        self.assertIn("Results:", report)
        self.assertIn("Successful: 10 items", report)
        self.assertIn("Failed: 0 items", report)
    
    def test_reset(self):
        """Test resetting progress tracker."""
        tracker = ProgressTracker("test_tracker", 100)
        tracker.update(10)
        
        tracker.reset()
        
        self.assertEqual(tracker.processed_items, 0)
        self.assertEqual(tracker.successful_items, 0)
        self.assertEqual(tracker.failed_items, 0)


class TestConnectionHealthMonitor(unittest.TestCase):
    """Test cases for ConnectionHealthMonitor class."""
    
    @patch('monitoring_utils.ConnectionManager')
    def test_init(self, mock_connection_manager):
        """Test initialization of ConnectionHealthMonitor."""
        # Reset singleton instance
        ConnectionHealthMonitor._instance = None
        
        monitor = ConnectionHealthMonitor.get_instance()
        self.assertIsNotNone(monitor)
        self.assertEqual(monitor.health_status["is_healthy"], False)
        
        # Test singleton pattern
        monitor2 = ConnectionHealthMonitor.get_instance()
        self.assertIs(monitor, monitor2)
    
    @patch('monitoring_utils.ConnectionManager')
    @patch('monitoring_utils.threading.Thread')
    def test_start_monitoring(self, mock_thread, mock_connection_manager):
        """Test starting monitoring."""
        # Reset singleton instance
        ConnectionHealthMonitor._instance = None
        
        monitor = ConnectionHealthMonitor.get_instance()
        monitor.start_monitoring()
        
        mock_thread.assert_called_once()
        mock_thread.return_value.start.assert_called_once()
    
    @patch('monitoring_utils.ConnectionManager')
    @patch('monitoring_utils.threading.Thread')
    def test_stop_monitoring(self, mock_thread, mock_connection_manager):
        """Test stopping monitoring."""
        # Reset singleton instance
        ConnectionHealthMonitor._instance = None
        
        monitor = ConnectionHealthMonitor.get_instance()
        monitor.health_check_thread = mock_thread.return_value
        mock_thread.return_value.is_alive.return_value = True
        
        monitor.stop_monitoring()
        
        monitor.stop_health_check.is_set.assert_called_once()
        mock_thread.return_value.join.assert_called_once()
    
    @patch('monitoring_utils.ConnectionManager')
    def test_check_connection_health(self, mock_connection_manager):
        """Test checking connection health."""
        # Reset singleton instance
        ConnectionHealthMonitor._instance = None
        
        # Mock connection manager
        mock_instance = MagicMock()
        mock_connection_manager.get_instance.return_value = mock_instance
        
        # Mock connection
        mock_connection = MagicMock()
        mock_instance.get_connection.return_value = mock_connection
        
        # Mock execute_query
        mock_connection.execute_query.return_value = [{"test": 1}]
        
        monitor = ConnectionHealthMonitor.get_instance()
        status = monitor.check_connection_health()
        
        self.assertTrue(status["is_healthy"])
        self.assertEqual(len(status["recent_successes"]), 1)
        self.assertEqual(len(status["recent_failures"]), 0)
        
        # Test failed connection
        mock_connection.execute_query.return_value = []
        
        status = monitor.check_connection_health()
        
        self.assertFalse(status["is_healthy"])
        self.assertEqual(len(status["recent_successes"]), 1)
        self.assertEqual(len(status["recent_failures"]), 1)
    
    @patch('monitoring_utils.ConnectionManager')
    def test_get_health_status(self, mock_connection_manager):
        """Test getting health status."""
        # Reset singleton instance
        ConnectionHealthMonitor._instance = None
        
        monitor = ConnectionHealthMonitor.get_instance()
        monitor.health_status["last_check_time"] = time.time()
        
        status = monitor.get_health_status()
        
        self.assertEqual(status, monitor.health_status)


class TestUtilityFunctions(unittest.TestCase):
    """Test cases for utility functions."""
    
    @patch('monitoring_utils.os.makedirs')
    @patch('monitoring_utils.os.path.exists')
    def test_create_monitoring_directory(self, mock_exists, mock_makedirs):
        """Test creating monitoring directory."""
        mock_exists.return_value = False
        
        create_monitoring_directory()
        
        mock_makedirs.assert_called_once_with("monitoring")
    
    @patch('monitoring_utils.ConnectionHealthMonitor')
    def test_monitor_connection_health(self, mock_monitor):
        """Test monitoring connection health."""
        mock_instance = MagicMock()
        mock_monitor.get_instance.return_value = mock_instance
        
        monitor_connection_health(30)
        
        self.assertEqual(os.environ['CONNECTION_HEALTH_CHECK_INTERVAL'], "30")
        mock_instance.start_monitoring.assert_called_once()
    
    @patch('monitoring_utils.ProgressTracker')
    def test_track_progress(self, mock_tracker):
        """Test tracking progress."""
        mock_tracker.return_value = MagicMock()
        
        tracker = track_progress("test_tracker", 100)
        
        mock_tracker.assert_called_once_with("test_tracker", 100, None)
        self.assertEqual(tracker, mock_tracker.return_value)
    
    @patch('monitoring_utils.PerformanceMetrics')
    def test_collect_performance_metrics(self, mock_metrics):
        """Test collecting performance metrics."""
        mock_metrics.return_value = MagicMock()
        
        metrics = collect_performance_metrics("test_metrics")
        
        mock_metrics.assert_called_once_with("test_metrics", None)
        self.assertEqual(metrics, mock_metrics.return_value)
    
    @patch('monitoring_utils.ConnectionHealthMonitor')
    def test_get_connection_health_status(self, mock_monitor):
        """Test getting connection health status."""
        mock_instance = MagicMock()
        mock_monitor.get_instance.return_value = mock_instance
        mock_instance.get_health_status.return_value = {"is_healthy": True}
        
        status = get_connection_health_status()
        
        mock_instance.get_health_status.assert_called_once()
        self.assertEqual(status, {"is_healthy": True})
    
    @patch('monitoring_utils.get_connection_health_status')
    @patch('monitoring_utils.os.path.exists')
    @patch('monitoring_utils.os.listdir')
    @patch('monitoring_utils.os.path.getsize')
    @patch('monitoring_utils.os.path.getmtime')
    def test_generate_monitoring_report(self, mock_getmtime, mock_getsize, mock_listdir, mock_exists, mock_get_status):
        """Test generating monitoring report."""
        mock_get_status.return_value = {
            "is_healthy": True,
            "active_connection_type": "pooler",
            "last_check_time": time.time(),
            "connection_types": {}
        }
        mock_exists.return_value = True
        mock_listdir.return_value = ["test_metrics.json"]
        mock_getsize.return_value = 1024
        mock_getmtime.return_value = time.time()
        
        report = generate_monitoring_report()
        
        self.assertIn("CryoProtect v2 Database Monitoring Report", report)
        self.assertIn("Connection Health:", report)
        self.assertIn("Status: Healthy", report)
        self.assertIn("Available Metrics Files:", report)
        self.assertIn("test_metrics.json", report)
    
    @patch('monitoring_utils.PerformanceMetrics')
    def test_timed_operation(self, mock_metrics):
        """Test timed operation context manager."""
        mock_metrics_instance = MagicMock()
        
        with timed_operation(mock_metrics_instance, "test_operation", 10):
            pass
        
        mock_metrics_instance.record_operation.assert_called_once()
        args, kwargs = mock_metrics_instance.record_operation.call_args
        self.assertEqual(kwargs["operation_type"], "test_operation")
        self.assertEqual(kwargs["success"], True)
        self.assertEqual(kwargs["items_processed"], 10)
        
        # Test with exception
        mock_metrics_instance.reset_mock()
        
        try:
            with timed_operation(mock_metrics_instance, "test_operation", 10):
                raise ValueError("Test exception")
        except ValueError:
            pass
        
        mock_metrics_instance.record_operation.assert_called_once()
        args, kwargs = mock_metrics_instance.record_operation.call_args
        self.assertEqual(kwargs["operation_type"], "test_operation")
        self.assertEqual(kwargs["success"], False)
        self.assertEqual(kwargs["items_processed"], 10)


if __name__ == '__main__':
    unittest.main()
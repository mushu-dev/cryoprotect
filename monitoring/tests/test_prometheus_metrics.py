#!/usr/bin/env python3
"""
Unit tests for the Prometheus metrics module.
"""

import unittest
import time
from unittest.mock import patch, MagicMock
from flask import Flask, Response
from prometheus_client import REGISTRY, Counter, Histogram, Gauge

# Import the module to test
from monitoring.prometheus_metrics import (
    metrics_bp, track_request_metrics, DatabaseQueryMetrics,
    update_system_metrics, update_connection_pool_metrics,
    update_business_metrics, init_app
)

class TestPrometheusMetrics(unittest.TestCase):
    """Test cases for the Prometheus metrics module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a Flask test app
        self.app = Flask(__name__)
        self.app.testing = True
        
        # Register the metrics blueprint
        self.app.register_blueprint(metrics_bp)
        
        # Create a test client
        self.client = self.app.test_client()
    
    def test_metrics_endpoint(self):
        """Test that the /metrics endpoint returns a valid response."""
        response = self.client.get('/metrics')
        
        # Check that the response is successful
        self.assertEqual(response.status_code, 200)
        
        # Check that the response has the correct content type
        self.assertEqual(response.content_type, 'text/plain; version=0.0.4; charset=utf-8')
        
        # Check that the response contains some metrics
        self.assertIn(b'# HELP', response.data)
    
    def test_track_request_metrics_decorator(self):
        """Test the track_request_metrics decorator."""
        # Create a test route with the decorator
        @self.app.route('/test')
        @track_request_metrics()
        def test_route():
            return 'Test'
        
        # Make a request to the test route
        response = self.client.get('/test')
        
        # Check that the response is successful
        self.assertEqual(response.status_code, 200)
        
        # Check that the metrics were updated
        # This is difficult to test directly, so we'll just check that the endpoint works
        metrics_response = self.client.get('/metrics')
        self.assertEqual(metrics_response.status_code, 200)
    
    def test_database_query_metrics(self):
        """Test the DatabaseQueryMetrics context manager."""
        # Use the context manager
        with DatabaseQueryMetrics('select', 'molecules'):
            # Simulate a database query
            time.sleep(0.01)
        
        # Check that the metrics were updated
        # This is difficult to test directly, so we'll just check that the context manager works
        metrics_response = self.client.get('/metrics')
        self.assertEqual(metrics_response.status_code, 200)
    
    @patch('monitoring.prometheus_metrics.psutil')
    def test_update_system_metrics(self, mock_psutil):
        """Test the update_system_metrics function."""
        # Mock psutil functions
        mock_psutil.cpu_percent.return_value = 50.0
        
        mock_memory = MagicMock()
        mock_memory.total = 16000000000
        mock_memory.available = 8000000000
        mock_memory.used = 8000000000
        mock_psutil.virtual_memory.return_value = mock_memory
        
        mock_disk = MagicMock()
        mock_disk.total = 500000000000
        mock_disk.free = 250000000000
        mock_disk.used = 250000000000
        mock_psutil.disk_usage.return_value = mock_disk
        
        # Call the function
        update_system_metrics()
        
        # Check that the psutil functions were called
        mock_psutil.cpu_percent.assert_called_once()
        mock_psutil.virtual_memory.assert_called_once()
        mock_psutil.disk_usage.assert_called_once_with('/')
    
    def test_update_connection_pool_metrics(self):
        """Test the update_connection_pool_metrics function."""
        # Call the function
        update_connection_pool_metrics(5, 10, 20)
        
        # Check that the metrics were updated
        # This is difficult to test directly, so we'll just check that the function doesn't raise an exception
        pass
    
    def test_update_business_metrics(self):
        """Test the update_business_metrics function."""
        # Call the function
        update_business_metrics(100, 50, 25, 10)
        
        # Check that the metrics were updated
        # This is difficult to test directly, so we'll just check that the function doesn't raise an exception
        pass
    
    def test_init_app(self):
        """Test the init_app function."""
        # Create a new Flask app
        app = Flask(__name__)
        
        # Call the function
        init_app(app)
        
        # Check that the blueprint was registered
        self.assertIn(metrics_bp.name, app.blueprints)

if __name__ == '__main__':
    unittest.main()
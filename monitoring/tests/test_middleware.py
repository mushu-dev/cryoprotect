#!/usr/bin/env python3
"""
Unit tests for the monitoring middleware module.
"""

import unittest
from unittest.mock import patch, MagicMock
from flask import Flask, g, request

# Import the module to test
from monitoring.middleware import (
    PrometheusMiddleware, track_database_metrics,
    update_connection_pool_status, update_entity_counts
)

class TestMonitoringMiddleware(unittest.TestCase):
    """Test cases for the monitoring middleware module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a Flask test app
        self.app = Flask(__name__)
        self.app.testing = True
        
        # Create a test client
        self.client = self.app.test_client()
    
    def test_prometheus_middleware_init(self):
        """Test the PrometheusMiddleware initialization."""
        # Initialize the middleware
        middleware = PrometheusMiddleware(self.app)
        
        # Check that the middleware was initialized
        self.assertEqual(middleware.app, self.app)
    
    def test_prometheus_middleware_init_app(self):
        """Test the PrometheusMiddleware init_app method."""
        # Initialize the middleware
        middleware = PrometheusMiddleware()
        middleware.init_app(self.app)
        
        # Check that the middleware was initialized
        self.assertEqual(middleware.app, self.app)
        
        # Create a test route
        @self.app.route('/test')
        def test_route():
            return 'Test'
        
        # Make a request to the test route
        response = self.client.get('/test')
        
        # Check that the response is successful
        self.assertEqual(response.status_code, 200)
    
    def test_track_database_metrics_decorator(self):
        """Test the track_database_metrics decorator."""
        # Create a test function with the decorator
        @track_database_metrics
        def test_function(operation='select', table='molecules'):
            return 'Test'
        
        # Call the function
        result = test_function()
        
        # Check that the function returns the expected result
        self.assertEqual(result, 'Test')
    
    @patch('monitoring.middleware.update_connection_pool_metrics')
    def test_update_connection_pool_status(self, mock_update_metrics):
        """Test the update_connection_pool_status function."""
        # Call the function
        update_connection_pool_status(5, 10, 20)
        
        # Check that the update_connection_pool_metrics function was called
        mock_update_metrics.assert_called_once_with(5, 10, 20)
    
    @patch('monitoring.middleware.update_business_metrics')
    def test_update_entity_counts(self, mock_update_metrics):
        """Test the update_entity_counts function."""
        # Create a mock Supabase client
        mock_supabase = MagicMock()
        
        # Mock the from_ method and its chain
        mock_from = MagicMock()
        mock_supabase.from_.return_value = mock_from
        
        mock_select = MagicMock()
        mock_from.select.return_value = mock_select
        
        mock_execute = MagicMock()
        mock_select.execute.return_value = mock_execute
        
        # Mock the count attribute
        mock_execute.count = 100
        
        # Call the function
        update_entity_counts(mock_supabase)
        
        # Check that the from_ method was called for each entity type
        self.assertEqual(mock_supabase.from_.call_count, 4)
        
        # Check that the update_business_metrics function was called
        mock_update_metrics.assert_called_once()
    
    @patch('monitoring.middleware.update_business_metrics')
    def test_update_entity_counts_exception(self, mock_update_metrics):
        """Test the update_entity_counts function with an exception."""
        # Create a mock Supabase client that raises an exception
        mock_supabase = MagicMock()
        mock_supabase.from_.side_effect = Exception('Test exception')
        
        # Call the function
        update_entity_counts(mock_supabase)
        
        # Check that the update_business_metrics function was not called
        mock_update_metrics.assert_not_called()

if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python3
"""
Simple test script for the API observability module.

This script tests the basic functionality of the API observability module
without depending on the full application.
"""

import unittest
import uuid
import time
import json
import logging
from unittest.mock import patch, MagicMock
from flask import Flask, request, g, jsonify

from api.observability import (
    ObservabilityMiddleware,
    trace_request,
    track_performance,
    report_error,
    get_correlation_id,
    get_request_context,
    init_observability,
    CORRELATION_ID_HEADER,
    REQUEST_ID_HEADER,
    TRACE_ID_HEADER,
    SPAN_ID_HEADER,
    TIMING_HEADER
)

class SimpleObservabilityTest(unittest.TestCase):
    """Simple test case for the API observability module."""
    
    def setUp(self):
        """Set up test case."""
        # Create a test logger
        self.logger = logging.getLogger('test_observability')
        
        # Create a test app
        self.app = Flask('test_app')
        self.app.config['TESTING'] = True
        
        # Initialize observability for the test app
        init_observability(self.app)
        
        # Set up a test client
        self.client = self.app.test_client()
        
        # Define a test endpoint
        @self.app.route('/test')
        def test_endpoint():
            return jsonify({'status': 'ok'})
        
        # Define a test endpoint that raises an error
        @self.app.route('/test-error')
        def test_error_endpoint():
            raise ValueError('Test error')
    
    def test_correlation_id(self):
        """Test correlation ID generation and propagation."""
        # Make a request with no correlation ID
        response = self.client.get('/test')
        
        # Check that a correlation ID was generated and included in the response
        self.assertIn(CORRELATION_ID_HEADER, response.headers)
        
        # Make a request with a custom correlation ID
        custom_id = str(uuid.uuid4())
        response = self.client.get('/test', headers={CORRELATION_ID_HEADER: custom_id})
        
        # Check that the custom correlation ID was used
        self.assertEqual(response.headers.get(CORRELATION_ID_HEADER), custom_id)
    
    def test_request_tracing(self):
        """Test request tracing headers."""
        # Make a request
        response = self.client.get('/test')
        
        # Check that tracing headers are present
        self.assertIn(CORRELATION_ID_HEADER, response.headers)
        self.assertIn(REQUEST_ID_HEADER, response.headers)
        self.assertIn(TRACE_ID_HEADER, response.headers)
        self.assertIn(SPAN_ID_HEADER, response.headers)
    
    def test_performance_timing(self):
        """Test performance timing header."""
        # Make a request
        response = self.client.get('/test')
        
        # Check that timing header is present
        self.assertIn(TIMING_HEADER, response.headers)
        
        # Check that timing value is a positive integer
        timing_value = int(response.headers.get(TIMING_HEADER))
        self.assertGreaterEqual(timing_value, 0)
    
    @patch('api.observability.log_with_context')
    def test_error_reporting(self, mock_log):
        """Test error reporting."""
        # Make a request that will trigger an error
        try:
            response = self.client.get('/test-error')
        except ValueError:
            pass  # Expected exception
        
        # Check that the error was logged
        mock_log.assert_any_call(
            unittest.mock.ANY,  # logger
            'error',  # log level
            unittest.mock.ANY,  # message
            context=unittest.mock.ANY,  # context dict
            exc_info=True
        )
    
    def test_trace_request_decorator(self):
        """Test the trace_request decorator."""
        # Define a test function with the decorator
        @trace_request
        def test_function(arg1, arg2=None):
            return f"{arg1}-{arg2}"
        
        # Call the function
        result = test_function('value1', arg2='value2')
        
        # Check the result
        self.assertEqual(result, 'value1-value2')
    
    def test_track_performance_decorator(self):
        """Test the track_performance decorator."""
        # Define a test function with the decorator
        @track_performance(threshold=1.0)
        def test_function():
            time.sleep(0.1)  # Small delay
            return 'result'
        
        # Call the function
        result = test_function()
        
        # Check the result
        self.assertEqual(result, 'result')

if __name__ == '__main__':
    unittest.main()
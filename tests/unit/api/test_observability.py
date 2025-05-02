"""
Unit tests for the API observability module.

This module tests the functionality of the API observability module, including:
- Request tracing with correlation IDs
- Performance timing middleware
- Detailed error reporting
- Integration with structured logging system
"""

import unittest
import uuid
import time
import json
import logging
from unittest.mock import patch, MagicMock, call
from flask import Flask, request, g, jsonify
from werkzeug.exceptions import BadRequest

from tests.base_test_case import BaseTestCase
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

class TestObservability(BaseTestCase):
    """Test cases for the API observability module."""
    
    def setUp(self):
        """Set up test case."""
        # Create a test logger
        self.logger = logging.getLogger('test_observability')
        
        # Create a test app
        self.test_app = Flask('test_app')
        self.test_app.config['TESTING'] = True
        
        # Set up a test client
        self.test_client = self.test_app.test_client()
        
        # Initialize observability for the test app
        init_observability(self.test_app)
    
    def test_correlation_id_generation(self):
        """Test correlation ID generation and propagation."""
        # Define a test endpoint
        @self.test_app.route('/test-correlation-id')
        def test_endpoint():
            # Get correlation ID
            correlation_id = get_correlation_id()
            return jsonify({'correlation_id': correlation_id})
        
        # Test with no correlation ID header
        response = self.test_client.get('/test-correlation-id')
        self.assertEqual(response.status_code, 200)
        
        # Check that a correlation ID was generated and included in the response
        self.assertIn(CORRELATION_ID_HEADER, response.headers)
        correlation_id = response.headers.get(CORRELATION_ID_HEADER)
        self.assertIsNotNone(correlation_id)
        
        # Verify the response body contains the same correlation ID
        response_data = json.loads(response.data)
        self.assertEqual(response_data['correlation_id'], correlation_id)
        
        # Test with a provided correlation ID header
        test_correlation_id = str(uuid.uuid4())
        response = self.test_client.get(
            '/test-correlation-id',
            headers={CORRELATION_ID_HEADER: test_correlation_id}
        )
        
        # Check that the provided correlation ID was used
        self.assertEqual(response.headers.get(CORRELATION_ID_HEADER), test_correlation_id)
        
        # Verify the response body contains the same correlation ID
        response_data = json.loads(response.data)
        self.assertEqual(response_data['correlation_id'], test_correlation_id)
    
    def test_request_tracing_headers(self):
        """Test that tracing headers are added to responses."""
        # Define a test endpoint
        @self.test_app.route('/test-tracing-headers')
        def test_endpoint():
            return jsonify({'status': 'ok'})
        
        # Make a request
        response = self.test_client.get('/test-tracing-headers')
        
        # Check that tracing headers are present
        self.assertIn(CORRELATION_ID_HEADER, response.headers)
        self.assertIn(REQUEST_ID_HEADER, response.headers)
        self.assertIn(TRACE_ID_HEADER, response.headers)
        self.assertIn(SPAN_ID_HEADER, response.headers)
    
    def test_performance_timing_header(self):
        """Test that performance timing header is added to responses."""
        # Define a test endpoint with a delay
        @self.test_app.route('/test-timing')
        def test_endpoint():
            time.sleep(0.1)  # Small delay to ensure timing is measurable
            return jsonify({'status': 'ok'})
        
        # Make a request
        response = self.test_client.get('/test-timing')
        
        # Check that timing header is present
        self.assertIn(TIMING_HEADER, response.headers)
        
        # Check that timing value is a positive integer
        timing_value = int(response.headers.get(TIMING_HEADER))
        self.assertGreater(timing_value, 0)
    
    @patch('api.observability.log_with_context')
    def test_slow_request_logging(self, mock_log_with_context):
        """Test that slow requests are logged with appropriate level."""
        # Define a test endpoint with a significant delay
        @self.test_app.route('/test-slow-request')
        def test_endpoint():
            time.sleep(0.2)  # Delay to trigger slow request logging
            return jsonify({'status': 'ok'})
        
        # Temporarily lower the warning threshold for testing
        original_warning_threshold = self.test_app.extensions.get('observability_thresholds', {}).get('warning', 1.0)
        self.test_app.extensions['observability_thresholds'] = {'warning': 0.1, 'critical': 0.5}
        
        try:
            # Make a request
            response = self.test_client.get('/test-slow-request')
            
            # Check that the slow request was logged
            mock_log_with_context.assert_any_call(
                unittest.mock.ANY,  # logger
                'warning',  # log level
                unittest.mock.ANY,  # message
                context=unittest.mock.ANY,  # context dict
            )
            
            # Get the context from the call
            context = None
            for call_args in mock_log_with_context.call_args_list:
                args, kwargs = call_args
                if args[1] == 'warning' and 'Slow request' in args[2]:
                    context = kwargs.get('context', {})
                    break
            
            # Check that the context contains performance information
            self.assertIsNotNone(context)
            self.assertIn('performance', context)
            self.assertEqual(context['performance']['threshold_exceeded'], 'warning')
            
        finally:
            # Restore original threshold
            if original_warning_threshold:
                self.test_app.extensions['observability_thresholds'] = {'warning': original_warning_threshold}
    
    @patch('api.observability.log_with_context')
    def test_error_reporting(self, mock_log_with_context):
        """Test error reporting functionality."""
        # Define a test endpoint that raises an error
        @self.test_app.route('/test-error')
        def test_endpoint():
            # Raise a test error
            raise BadRequest('Test error message')
        
        # Make a request that will trigger an error
        try:
            response = self.test_client.get('/test-error')
        except BadRequest:
            pass  # Expected exception
        
        # Check that the error was logged
        mock_log_with_context.assert_any_call(
            unittest.mock.ANY,  # logger
            'error',  # log level
            unittest.mock.ANY,  # message
            context=unittest.mock.ANY,  # context dict
            exc_info=True
        )
        
        # Get the context from the call
        context = None
        for call_args in mock_log_with_context.call_args_list:
            args, kwargs = call_args
            if args[1] == 'error' and 'Request failed' in args[2]:
                context = kwargs.get('context', {})
                break
        
        # Check that the context contains error information
        self.assertIsNotNone(context)
        self.assertIn('error', context)
        self.assertEqual(context['error']['type'], 'BadRequest')
        self.assertEqual(context['error']['message'], 'Test error message')
    
    @patch('api.observability.log_with_context')
    def test_trace_request_decorator(self, mock_log_with_context):
        """Test the trace_request decorator."""
        # Define a test function with the decorator
        @trace_request
        def test_function(arg1, arg2=None):
            return f"{arg1}-{arg2}"
        
        # Set up a request context
        with self.test_app.test_request_context('/test-path'):
            # Set correlation ID and trace ID
            g.correlation_id = 'test-correlation-id'
            g.trace_id = 'test-trace-id'
            g.span_id = 'test-span-id'
            
            # Call the function
            result = test_function('value1', arg2='value2')
            
            # Check the result
            self.assertEqual(result, 'value1-value2')
            
            # Check that function entry and exit were logged
            mock_log_with_context.assert_any_call(
                unittest.mock.ANY,  # logger
                'debug',  # log level
                'Function test_function started',
                context=unittest.mock.ANY
            )
            
            mock_log_with_context.assert_any_call(
                unittest.mock.ANY,  # logger
                'debug',  # log level
                'Function test_function completed',
                context=unittest.mock.ANY
            )
            
            # Get the context from the first call
            start_context = None
            for call_args in mock_log_with_context.call_args_list:
                args, kwargs = call_args
                if args[1] == 'debug' and 'started' in args[2]:
                    start_context = kwargs.get('context', {})
                    break
            
            # Check that the context contains tracing information
            self.assertIsNotNone(start_context)
            self.assertIn('trace', start_context)
            self.assertEqual(start_context['trace']['correlation_id'], 'test-correlation-id')
            self.assertEqual(start_context['trace']['trace_id'], 'test-trace-id')
            self.assertEqual(start_context['trace']['parent_span_id'], 'test-span-id')
    
    @patch('api.observability.log_with_context')
    def test_track_performance_decorator(self, mock_log_with_context):
        """Test the track_performance decorator."""
        # Define a test function with the decorator
        @track_performance(threshold=0.1)
        def slow_function():
            time.sleep(0.2)  # Delay to exceed threshold
            return 'result'
        
        # Set up a request context
        with self.test_app.test_request_context('/test-path'):
            # Set correlation ID and request ID
            g.correlation_id = 'test-correlation-id'
            g.request_id = 'test-request-id'
            
            # Call the function
            result = slow_function()
            
            # Check the result
            self.assertEqual(result, 'result')
            
            # Check that slow function was logged
            mock_log_with_context.assert_any_call(
                unittest.mock.ANY,  # logger
                'warning',  # log level
                unittest.mock.ANY,  # message containing "Slow function"
                context=unittest.mock.ANY
            )
            
            # Get the context from the call
            context = None
            for call_args in mock_log_with_context.call_args_list:
                args, kwargs = call_args
                if args[1] == 'warning' and 'Slow function' in args[2]:
                    context = kwargs.get('context', {})
                    break
            
            # Check that the context contains performance information
            self.assertIsNotNone(context)
            self.assertIn('performance', context)
            self.assertGreaterEqual(context['performance']['duration'], 0.1)
            self.assertEqual(context['performance']['threshold'], 0.1)
    
    def test_report_error_function(self):
        """Test the report_error function."""
        # Create a test error
        test_error = ValueError('Test error message')
        
        # Set up a mock for log_with_context
        with patch('api.observability.log_with_context') as mock_log_with_context:
            # Set up a request context
            with self.test_app.test_request_context('/test-path'):
                # Set correlation ID and request ID
                g.correlation_id = 'test-correlation-id'
                g.request_id = 'test-request-id'
                
                # Report the error
                report_error(test_error, {'custom_context': 'test-value'})
                
                # Check that the error was logged
                mock_log_with_context.assert_called_once_with(
                    unittest.mock.ANY,  # logger
                    'error',  # log level
                    'Error: ValueError: Test error message',
                    context=unittest.mock.ANY,
                    exc_info=True
                )
                
                # Get the context from the call
                context = mock_log_with_context.call_args[1]['context']
                
                # Check that the context contains error information
                self.assertIn('error', context)
                self.assertEqual(context['error']['type'], 'ValueError')
                self.assertEqual(context['error']['message'], 'Test error message')
                
                # Check that the context contains tracing information
                self.assertIn('trace', context)
                self.assertEqual(context['trace']['correlation_id'], 'test-correlation-id')
                self.assertEqual(context['trace']['request_id'], 'test-request-id')
                
                # Check that the context contains custom context
                self.assertEqual(context['custom_context'], 'test-value')
    
    def test_get_request_context(self):
        """Test the get_request_context function."""
        # Set up a request context
        with self.test_app.test_request_context('/test-path', method='POST'):
            # Set correlation ID and request ID
            g.correlation_id = 'test-correlation-id'
            g.request_id = 'test-request-id'
            g.trace_id = 'test-trace-id'
            g.span_id = 'test-span-id'
            g.user_id = 'test-user-id'
            
            # Get the request context
            context = get_request_context()
            
            # Check that the context contains tracing information
            self.assertIn('trace', context)
            self.assertEqual(context['trace']['correlation_id'], 'test-correlation-id')
            self.assertEqual(context['trace']['request_id'], 'test-request-id')
            self.assertEqual(context['trace']['trace_id'], 'test-trace-id')
            self.assertEqual(context['trace']['span_id'], 'test-span-id')
            
            # Check that the context contains request information
            self.assertIn('request', context)
            self.assertEqual(context['request']['method'], 'POST')
            self.assertEqual(context['request']['path'], '/test-path')
            
            # Check that the context contains user information
            self.assertEqual(context['user_id'], 'test-user-id')
            
            # Check that the context contains timestamp
            self.assertIn('timestamp', context)

if __name__ == '__main__':
    unittest.main()
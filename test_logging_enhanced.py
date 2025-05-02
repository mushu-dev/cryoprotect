#!/usr/bin/env python3
"""
Unit tests for the enhanced logging system.
"""

import os
import json
import logging
import unittest
import tempfile
import shutil
from unittest.mock import patch, MagicMock, ANY
from flask import Flask, request, g
from werkzeug.test import Client
from werkzeug.wrappers import Response

import logging_enhanced


class TestEnhancedLogging(unittest.TestCase):
    """Test cases for the enhanced logging system."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory for log files
        self.test_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(self.test_dir, "test.log")
        
        # Set environment variables for testing
        self.env_patcher = patch.dict('os.environ', {
            'LOG_LEVEL': 'DEBUG',
            'LOG_FILE': self.log_file,
            'LOG_TO_FILE': '1',
            'LOG_TO_CONSOLE': '0',
            'LOG_TO_ELK': '0',
            'LOG_JSON_FORMAT': '1',
            'LOG_ROTATION_SIZE': '1024',
            'LOG_BACKUP_COUNT': '3',
        })
        self.env_patcher.start()
        
        # Clear any existing handlers
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

    def tearDown(self):
        """Clean up after tests."""
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)
        
        # Stop environment variable patching
        self.env_patcher.stop()
        
        # Clear any existing handlers
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

    def test_setup_enhanced_logging(self):
        """Test that setup_enhanced_logging configures logging correctly."""
        # Set up enhanced logging
        logging_enhanced.setup_enhanced_logging()
        
        # Get the root logger
        root_logger = logging.getLogger()
        
        # Check log level
        self.assertEqual(root_logger.level, logging.DEBUG)
        
        # Check that handlers are set up correctly
        self.assertTrue(any(isinstance(h, logging.handlers.RotatingFileHandler) 
                           for h in root_logger.handlers))
        
        # Log a test message
        test_logger = logging.getLogger("test_logger")
        test_message = "Test log message"
        test_logger.info(test_message)
        
        # Check that the message was logged to the file
        with open(self.log_file, 'r') as f:
            log_content = f.read()
            self.assertIn(test_message, log_content)
            # Check that it's in JSON format
            try:
                log_lines = log_content.strip().split('\n')
                for line in log_lines:
                    json_log = json.loads(line)
                    if test_message in json_log.get('message', ''):
                        self.assertEqual(json_log.get('level'), 'INFO')
                        self.assertEqual(json_log.get('logger'), 'test_logger')
                        break
                else:
                    self.fail("Test message not found in JSON logs")
            except json.JSONDecodeError:
                self.fail("Log is not in valid JSON format")

    def test_contextual_filter(self):
        """Test that the ContextualFilter adds context to log records."""
        # Create a mock record
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="test_file.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        # Create a filter
        filter = logging_enhanced.ContextualFilter()
        
        # Apply the filter
        with patch('logging_enhanced.request', MagicMock()):
            with patch('logging_enhanced.g', MagicMock(
                correlation_id='test-correlation-id',
                request_id='test-request-id',
                user_id='test-user-id'
            )):
                logging_enhanced.request.remote_addr = '127.0.0.1'
                logging_enhanced.request.method = 'GET'
                logging_enhanced.request.path = '/test'
                logging_enhanced.request.headers = {'User-Agent': 'Test Agent'}
                
                result = filter.filter(record)
                
                # Check that the filter returns True
                self.assertTrue(result)
                
                # Check that context was added to the record
                self.assertEqual(record.correlation_id, 'test-correlation-id')
                self.assertEqual(record.request_id, 'test-request-id')
                self.assertEqual(record.remote_addr, '127.0.0.1')
                self.assertEqual(record.method, 'GET')
                self.assertEqual(record.path, '/test')
                self.assertEqual(record.user_agent, 'Test Agent')
                self.assertEqual(record.user_id, 'test-user-id')

    def test_custom_json_formatter(self):
        """Test that the CustomJsonFormatter formats logs correctly."""
        # Create a formatter
        formatter = logging_enhanced.CustomJsonFormatter('%(message)s')
        
        # Create a mock record
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="test_file.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        # Add custom attributes
        record.correlation_id = 'test-correlation-id'
        record.request_id = 'test-request-id'
        record.remote_addr = '127.0.0.1'
        record.method = 'GET'
        record.path = '/test'
        record.user_agent = 'Test Agent'
        record.user_id = 'test-user-id'
        record.hostname = 'test-host'
        record.process_id = 12345
        
        # Format the record
        formatted = formatter.format(record)
        
        # Parse the JSON
        log_entry = json.loads(formatted)
        
        # Check that all fields are present
        self.assertEqual(log_entry['message'], 'Test message')
        self.assertEqual(log_entry['level'], 'INFO')
        self.assertEqual(log_entry['logger'], 'test_logger')
        self.assertEqual(log_entry['correlation_id'], 'test-correlation-id')
        self.assertEqual(log_entry['request_id'], 'test-request-id')
        self.assertEqual(log_entry['remote_addr'], '127.0.0.1')
        self.assertEqual(log_entry['method'], 'GET')
        self.assertEqual(log_entry['path'], '/test')
        self.assertEqual(log_entry['user_agent'], 'Test Agent')
        self.assertEqual(log_entry['user_id'], 'test-user-id')
        self.assertEqual(log_entry['hostname'], 'test-host')
        self.assertEqual(log_entry['process_id'], 12345)
        self.assertIn('timestamp', log_entry)

    @patch('logging_enhanced.ElasticsearchHandler')
    def test_elasticsearch_handler(self, mock_handler):
        """Test that the ElasticsearchHandler sends logs to Elasticsearch."""
        # Set up environment for ELK
        with patch.dict('os.environ', {
            'LOG_TO_ELK': '1',
            'ELASTICSEARCH_HOST': 'test-es-host:9200',
            'ELASTICSEARCH_INDEX': 'test-index'
        }):
            # Set up enhanced logging
            logging_enhanced.setup_enhanced_logging()
            
            # Check that ElasticsearchHandler was created
            mock_handler.assert_called_once_with('test-es-host:9200', 'test-index')

    def test_request_logging_middleware(self):
        """Test that request logging middleware logs requests and responses."""
        # Create a Flask app
        app = Flask(__name__)
        
        # Set up enhanced logging for the app
        logging_enhanced.setup_enhanced_logging(app)
        
        # Add a test route
        @app.route('/test')
        def test_route():
            return 'Test response'
        
        # Create a test client
        client = app.test_client()
        
        # Mock the logger
        with patch.object(app.logger, 'info') as mock_info:
            # Make a request
            response = client.get('/test')
            
            # Check that the response is correct
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.data, b'Test response')
            
            # Check that the request was logged
            mock_info.assert_any_call(
                'Request started: GET /test',
                extra=ANY
            )
            
            # Check that the response was logged
            mock_info.assert_any_call(
                'Request completed: GET /test 200',
                extra=ANY
            )

    def test_log_with_context(self):
        """Test that log_with_context adds context to logs."""
        # Set up enhanced logging
        logging_enhanced.setup_enhanced_logging()
        
        # Get a logger
        logger = logging_enhanced.get_logger('test_logger')
        
        # Mock the logger's info method
        with patch.object(logger, 'info') as mock_info:
            # Log with context
            context = {'key1': 'value1', 'key2': 'value2'}
            logging_enhanced.log_with_context(logger, 'info', 'Test message', context)
            
            # Check that the logger was called with the right arguments
            mock_info.assert_called_once_with('Test message', extra=context, exc_info=None)


if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python3
"""
Test script for optimize_toxicity_queries.py

This script tests the functionality of optimize_toxicity_queries.py
without actually connecting to the database.
"""

import unittest
from unittest.mock import patch, MagicMock
import sys
import os
import logging

# Disable logging for tests
logging.disable(logging.CRITICAL)

# Import the module to test
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import optimize_toxicity_queries

class TestOptimizeToxicityQueries(unittest.TestCase):
    """Test cases for optimize_toxicity_queries.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a mock Supabase client
        self.mock_supabase = MagicMock()
        
        # Mock the response from Supabase
        mock_response = MagicMock()
        mock_response.data = [{"id": "test-id"}]
        self.mock_supabase.table.return_value.select.return_value.limit.return_value.execute.return_value = mock_response
        self.mock_supabase.postgrest.rpc.return_value.execute.return_value = mock_response

    def test_add_indexes(self):
        """Test add_indexes function"""
        result = optimize_toxicity_queries.add_indexes(self.mock_supabase)
        self.assertTrue(result)
        
        # Check that the correct number of indexes were created
        expected_calls = 5  # Number of indexes in the function
        actual_calls = self.mock_supabase.postgrest.rpc.call_count
        self.assertEqual(expected_calls, actual_calls)

    def test_create_materialized_views(self):
        """Test create_materialized_views function"""
        result = optimize_toxicity_queries.create_materialized_views(self.mock_supabase)
        self.assertTrue(result)
        
        # Check that the correct number of views were created
        expected_calls = 2  # Number of views in the function
        actual_calls = self.mock_supabase.postgrest.rpc.call_count
        self.assertEqual(expected_calls, actual_calls)

    def test_refresh_materialized_views(self):
        """Test refresh_materialized_views function"""
        result = optimize_toxicity_queries.refresh_materialized_views(self.mock_supabase)
        self.assertTrue(result)
        
        # Check that the correct number of views were refreshed
        expected_calls = 2  # Number of views in the function
        actual_calls = self.mock_supabase.postgrest.rpc.call_count
        self.assertEqual(expected_calls, actual_calls)

    def test_test_query_performance(self):
        """Test test_query_performance function"""
        result = optimize_toxicity_queries.test_query_performance(self.mock_supabase)
        self.assertTrue(result)
        
        # Check that queries were executed
        self.assertTrue(self.mock_supabase.postgrest.rpc.called)

    @patch('optimize_toxicity_queries.create_client')
    @patch('optimize_toxicity_queries.add_indexes')
    @patch('optimize_toxicity_queries.create_materialized_views')
    @patch('optimize_toxicity_queries.test_query_performance')
    def test_main_function(self, mock_test, mock_views, mock_indexes, mock_create_client):
        """Test main function with mocked dependencies"""
        # Set up mocks
        mock_create_client.return_value = self.mock_supabase
        mock_indexes.return_value = True
        mock_views.return_value = True
        mock_test.return_value = True
        
        # Mock Config and environment variables
        with patch('optimize_toxicity_queries.Config') as mock_config:
            mock_config_instance = MagicMock()
            mock_config_instance.SUPABASE_URL = "https://example.supabase.co"
            mock_config_instance.SUPABASE_KEY = "test-key"
            mock_config.return_value = mock_config_instance
            
            with patch.dict('os.environ', {
                'SUPABASE_URL': 'https://example.supabase.co',
                'SUPABASE_KEY': 'test-key'
            }):
                # Test with default arguments
                with patch('sys.argv', ['optimize_toxicity_queries.py']):
                    result = optimize_toxicity_queries.main()
                    self.assertTrue(result)
                    mock_indexes.assert_called_once()
                    mock_views.assert_called_once()
                    mock_test.assert_called_once()
                
                # Reset mocks
                mock_indexes.reset_mock()
                mock_views.reset_mock()
                mock_test.reset_mock()
                
                # Test with --skip-indexes
                with patch('sys.argv', ['optimize_toxicity_queries.py', '--skip-indexes']):
                    result = optimize_toxicity_queries.main()
                    self.assertTrue(result)
                    mock_indexes.assert_not_called()
                    mock_views.assert_called_once()
                    mock_test.assert_called_once()
                
                # Reset mocks
                mock_indexes.reset_mock()
                mock_views.reset_mock()
                mock_test.reset_mock()
                
                # Test with --skip-views
                with patch('sys.argv', ['optimize_toxicity_queries.py', '--skip-views']):
                    result = optimize_toxicity_queries.main()
                    self.assertTrue(result)
                    mock_indexes.assert_called_once()
                    mock_views.assert_not_called()
                    mock_test.assert_called_once()
                
                # Reset mocks
                mock_indexes.reset_mock()
                mock_views.reset_mock()
                mock_test.reset_mock()
                
                # Test with --skip-tests
                with patch('sys.argv', ['optimize_toxicity_queries.py', '--skip-tests']):
                    result = optimize_toxicity_queries.main()
                    self.assertTrue(result)
                    mock_indexes.assert_called_once()
                    mock_views.assert_called_once()
                    mock_test.assert_not_called()

if __name__ == '__main__':
    unittest.main()
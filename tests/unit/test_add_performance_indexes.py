#!/usr/bin/env python3
"""
Unit tests for the add_performance_indexes.py script.

These tests verify that the database performance optimization script correctly
identifies and creates necessary indexes using the direct PostgreSQL connection.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import json
from datetime import datetime

# Add the parent directory to the path so we can import the script
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Import the script to test
import add_performance_indexes

class TestAddPerformanceIndexes(unittest.TestCase):
    """Test cases for the add_performance_indexes.py script."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a mock for the SQL executor
        self.mock_sql_executor = MagicMock()
        self.mock_db = MagicMock()
        self.mock_sql_executor.get_db.return_value = self.mock_db
        
        # Create a mock for the transaction context manager
        self.mock_transaction = MagicMock()
        self.mock_db.transaction.return_value = self.mock_transaction
        
        # Create a mock for the execute_query function
        self.mock_execute_query = MagicMock()
        self.mock_sql_executor.execute_query = self.mock_execute_query
        
        # Sample performance metrics
        self.sample_metrics = {
            "test_duration": 10.5,
            "operations": {
                "get_molecules": {
                    "count": 5,
                    "min_ms": 50.0,
                    "max_ms": 150.0,
                    "avg_ms": 80.0,
                    "median_ms": 75.0,
                    "p95_ms": 140.0,
                    "p99_ms": 150.0
                },
                "get_mixture_by_id": {
                    "count": 5,
                    "min_ms": 100.0,
                    "max_ms": 600.0,
                    "avg_ms": 250.0,
                    "median_ms": 200.0,
                    "p95_ms": 550.0,
                    "p99_ms": 600.0
                }
            }
        }

    @patch('add_performance_indexes.sql_executor')
    def test_get_performance_indexes_sql(self, mock_sql_executor):
        """Test that the get_performance_indexes_sql function returns the expected SQL."""
        # Call the function
        result = add_performance_indexes.get_performance_indexes_sql()
        
        # Check that the result is a dictionary with the expected keys
        self.assertIsInstance(result, dict)
        self.assertIn('rls_indexes', result)
        self.assertIn('mixture_indexes', result)
        self.assertIn('prediction_indexes', result)
        self.assertIn('text_search_indexes', result)
        self.assertIn('additional_indexes', result)
        
        # Check that each value is a string containing SQL
        for key, value in result.items():
            self.assertIsInstance(value, str)
            self.assertIn('CREATE INDEX', value)

    @patch('add_performance_indexes.sql_executor')
    def test_verify_index_creation(self, mock_sql_executor):
        """Test that the verify_index_creation function works correctly."""
        # Set up the mock to return a result indicating the index exists
        mock_sql_executor.execute_query.return_value = [{'1': 1}]
        
        # Call the function
        result = add_performance_indexes.verify_index_creation('public', 'idx_test')
        
        # Check that the function returned True
        self.assertTrue(result)
        
        # Check that execute_query was called with the expected arguments
        mock_sql_executor.execute_query.assert_called_once()
        args, kwargs = mock_sql_executor.execute_query.call_args
        self.assertIn('SELECT 1 FROM pg_indexes', args[0])
        self.assertEqual(args[1], ['public', 'idx_test'])
        
        # Set up the mock to return an empty result
        mock_sql_executor.execute_query.reset_mock()
        mock_sql_executor.execute_query.return_value = []
        
        # Call the function again
        result = add_performance_indexes.verify_index_creation('public', 'idx_nonexistent')
        
        # Check that the function returned False
        self.assertFalse(result)

    @patch('add_performance_indexes.run_performance_test')
    @patch('add_performance_indexes.verify_index_creation')
    @patch('add_performance_indexes.sql_executor')
    def test_apply_performance_indexes(self, mock_sql_executor, mock_verify, mock_run_test):
        """Test that the apply_performance_indexes function works correctly."""
        # Set up the mocks
        metrics = add_performance_indexes.PerformanceMetrics()
        metrics.start_test()
        metrics.record_operation_time('get_molecules', 80.0)
        metrics.record_operation_time('get_mixture_by_id', 250.0)
        metrics.end_test()
        
        mock_run_test.return_value = metrics
        mock_verify.return_value = True
        
        # Call the function with test_only=True to avoid actual index creation
        result = add_performance_indexes.apply_performance_indexes(test_only=True)
        
        # Check that the function returned a dictionary with the expected keys
        self.assertIsInstance(result, dict)
        self.assertIn('before_metrics', result)
        self.assertIn('analysis_before', result)
        self.assertIn('recommendations_before', result)
        
        # Check that run_performance_test was called
        mock_run_test.assert_called_once()
        
        # Now test with actual index creation
        mock_run_test.reset_mock()
        mock_sql_executor.execute_query.reset_mock()
        
        # Call the function with test_only=False
        with patch('add_performance_indexes.generate_performance_report'):
            result = add_performance_indexes.apply_performance_indexes(test_only=False)
        
        # Check that run_performance_test was called twice (before and after)
        self.assertEqual(mock_run_test.call_count, 2)
        
        # Check that execute_query was called for each batch of indexes
        self.assertGreater(mock_sql_executor.execute_query.call_count, 0)

    def test_calculate_improvements(self):
        """Test that the calculate_improvements function works correctly."""
        # Create sample results
        results = {
            "before_metrics": {
                "operations": {
                    "get_molecules": {
                        "avg_ms": 100.0,
                        "p95_ms": 200.0
                    },
                    "get_mixture_by_id": {
                        "avg_ms": 300.0,
                        "p95_ms": 500.0
                    }
                }
            },
            "after_metrics": {
                "operations": {
                    "get_molecules": {
                        "avg_ms": 50.0,
                        "p95_ms": 100.0
                    },
                    "get_mixture_by_id": {
                        "avg_ms": 150.0,
                        "p95_ms": 250.0
                    }
                }
            }
        }
        
        # Call the function
        improvements = add_performance_indexes.calculate_improvements(results)
        
        # Check the results
        self.assertIsInstance(improvements, dict)
        self.assertIn('get_molecules', improvements)
        self.assertIn('get_mixture_by_id', improvements)
        
        # Check the calculated improvements
        self.assertEqual(improvements['get_molecules']['avg_ms_before'], 100.0)
        self.assertEqual(improvements['get_molecules']['avg_ms_after'], 50.0)
        self.assertEqual(improvements['get_molecules']['avg_improvement_pct'], 50.0)
        
        self.assertEqual(improvements['get_mixture_by_id']['p95_ms_before'], 500.0)
        self.assertEqual(improvements['get_mixture_by_id']['p95_ms_after'], 250.0)
        self.assertEqual(improvements['get_mixture_by_id']['p95_improvement_pct'], 50.0)

    def test_analyze_metrics(self):
        """Test that the analyze_metrics function correctly identifies bottlenecks."""
        # Call the function with sample metrics
        analysis = add_performance_indexes.analyze_metrics(self.sample_metrics)
        
        # Check the results
        self.assertIsInstance(analysis, dict)
        self.assertIn('bottlenecks', analysis)
        self.assertIn('slow_operations', analysis)
        
        # Check that get_mixture_by_id is identified as a slow operation (p95 > 500ms)
        slow_ops = [op['operation'] for op in analysis['slow_operations']]
        self.assertIn('get_mixture_by_id', slow_ops)

    def test_generate_recommendations(self):
        """Test that the generate_recommendations function works correctly."""
        # Call the function with sample metrics
        recommendations = add_performance_indexes.generate_recommendations(self.sample_metrics)
        
        # Check the results
        self.assertIsInstance(recommendations, list)
        self.assertGreater(len(recommendations), 0)
        
        # Check that there's a recommendation for get_mixture_by_id
        mixture_recs = [rec for rec in recommendations if rec.get('operation') == 'get_mixture_by_id']
        self.assertGreater(len(mixture_recs), 0)

    def test_performance_metrics_class(self):
        """Test the PerformanceMetrics class."""
        # Create an instance
        metrics = add_performance_indexes.PerformanceMetrics()
        
        # Start the test
        metrics.start_test()
        self.assertIsNotNone(metrics.start_time)
        
        # Record some operation times
        metrics.record_operation_time('op1', 100.0)
        metrics.record_operation_time('op1', 200.0)
        metrics.record_operation_time('op2', 300.0)
        
        # End the test
        metrics.end_test()
        self.assertIsNotNone(metrics.end_time)
        
        # Get the summary
        summary = metrics.get_summary()
        
        # Check the summary
        self.assertIsInstance(summary, dict)
        self.assertIn('test_duration', summary)
        self.assertIn('operations', summary)
        self.assertIn('op1', summary['operations'])
        self.assertIn('op2', summary['operations'])
        
        # Check the operation stats
        op1_stats = summary['operations']['op1']
        self.assertEqual(op1_stats['count'], 2)
        self.assertEqual(op1_stats['min_ms'], 100.0)
        self.assertEqual(op1_stats['max_ms'], 200.0)
        self.assertEqual(op1_stats['avg_ms'], 150.0)

if __name__ == '__main__':
    unittest.main()
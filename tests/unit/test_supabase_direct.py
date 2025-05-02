#!/usr/bin/env python3
"""
Unit tests for supabase_direct.py

These tests use mocking to test the SupabaseDirectConnection class
without requiring actual Supabase credentials.
"""

import os
import unittest
from unittest.mock import patch, MagicMock, call
import psycopg2
import sys
import logging

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from supabase_direct import SupabaseDirectConnection

# Disable logging for tests
logging.disable(logging.CRITICAL)

class TestSupabaseDirectConnection(unittest.TestCase):
    """Test cases for SupabaseDirectConnection class."""

    def setUp(self):
        """Set up test environment."""
        # Set required environment variables
        os.environ['SUPABASE_DB_HOST'] = 'test-host'
        os.environ['SUPABASE_DB_PASSWORD'] = 'test-password'
        
        # Reset the singleton instance
        SupabaseDirectConnection._instance = None
    
    def tearDown(self):
        """Clean up after tests."""
        # Remove environment variables
        if 'SUPABASE_DB_HOST' in os.environ:
            del os.environ['SUPABASE_DB_HOST']
        if 'SUPABASE_DB_PASSWORD' in os.environ:
            del os.environ['SUPABASE_DB_PASSWORD']
        
        # Reset the singleton instance
        SupabaseDirectConnection._instance = None
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_singleton_pattern(self, mock_pool):
        """Test that the class follows the singleton pattern."""
        # Create first instance
        instance1 = SupabaseDirectConnection.get_instance()
        
        # Create second instance
        instance2 = SupabaseDirectConnection.get_instance()
        
        # Verify both instances are the same object
        self.assertIs(instance1, instance2)
        
        # Verify pool was only created once
        mock_pool.assert_called_once()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_initialization_with_env_vars(self, mock_pool):
        """Test initialization with environment variables."""
        # Set custom environment variables
        os.environ['SUPABASE_DB_PORT'] = '5433'
        os.environ['SUPABASE_DB_NAME'] = 'test-db'
        os.environ['SUPABASE_DB_USER'] = 'test-user'
        os.environ['SUPABASE_DB_MIN_CONNECTIONS'] = '3'
        os.environ['SUPABASE_DB_MAX_CONNECTIONS'] = '15'
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Verify environment variables were used
        self.assertEqual(instance.db_host, 'test-host')
        self.assertEqual(instance.db_port, 5433)
        self.assertEqual(instance.db_name, 'test-db')
        self.assertEqual(instance.db_user, 'test-user')
        self.assertEqual(instance.db_password, 'test-password')
        self.assertEqual(instance.min_connections, 3)
        self.assertEqual(instance.max_connections, 15)
        
        # Verify pool was created with correct parameters
        mock_pool.assert_called_once_with(
            minconn=3,
            maxconn=15,
            host='test-host',
            port=5433,
            dbname='test-db',
            user='test-user',
            password='test-password',
            connect_timeout=30
        )
        
        # Clean up custom environment variables
        del os.environ['SUPABASE_DB_PORT']
        del os.environ['SUPABASE_DB_NAME']
        del os.environ['SUPABASE_DB_USER']
        del os.environ['SUPABASE_DB_MIN_CONNECTIONS']
        del os.environ['SUPABASE_DB_MAX_CONNECTIONS']
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_missing_required_env_vars(self, mock_pool):
        """Test initialization with missing required environment variables."""
        # Remove required environment variables
        del os.environ['SUPABASE_DB_HOST']
        del os.environ['SUPABASE_DB_PASSWORD']
        
        # Verify ValueError is raised
        with self.assertRaises(ValueError):
            SupabaseDirectConnection.get_instance()
        
        # Verify pool was not created
        mock_pool.assert_not_called()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_get_connection(self, mock_pool):
        """Test getting a connection from the pool."""
        # Set up mock pool
        mock_conn = MagicMock()
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Get connection
        conn = instance.get_connection()
        
        # Verify connection was retrieved from pool
        self.assertEqual(conn, mock_conn)
        mock_pool.return_value.getconn.assert_called_once()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_release_connection(self, mock_pool):
        """Test releasing a connection back to the pool."""
        # Set up mock pool
        mock_conn = MagicMock()
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Release connection
        instance.release_connection(mock_conn)
        
        # Verify connection was released to pool
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_close_all(self, mock_pool):
        """Test closing all connections in the pool."""
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Close all connections
        instance.close_all()
        
        # Verify all connections were closed
        mock_pool.return_value.closeall.assert_called_once()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_check_connection_health(self, mock_pool):
        """Test checking connection health."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.fetchone.return_value = (1,)
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Check connection health
        result = instance.check_connection_health()
        
        # Verify connection health was checked
        self.assertTrue(result)
        mock_cursor.execute.assert_called_once_with("SELECT 1")
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_query_select(self, mock_pool):
        """Test executing a SELECT query."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.description = True  # Indicates a SELECT query
        mock_cursor.fetchall.return_value = [{'test': 1}]
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute query
        result = instance.execute_query("SELECT 1 as test")
        
        # Verify query was executed
        self.assertEqual(result, [{'test': 1}])
        mock_cursor.execute.assert_called_once_with("SELECT 1 as test", None)
        mock_conn.commit.assert_called_once()
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_query_non_select(self, mock_pool):
        """Test executing a non-SELECT query."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.description = None  # Indicates a non-SELECT query
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute query
        result = instance.execute_query("INSERT INTO test VALUES (1)")
        
        # Verify query was executed
        self.assertIsNone(result)
        mock_cursor.execute.assert_called_once_with("INSERT INTO test VALUES (1)", None)
        mock_conn.commit.assert_called_once()
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_query_with_params(self, mock_pool):
        """Test executing a query with parameters."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.description = True  # Indicates a SELECT query
        mock_cursor.fetchall.return_value = [{'test': 'value'}]
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute query with parameters
        params = {'param': 'value'}
        result = instance.execute_query("SELECT %(param)s as test", params)
        
        # Verify query was executed with parameters
        self.assertEqual(result, [{'test': 'value'}])
        mock_cursor.execute.assert_called_once_with("SELECT %(param)s as test", params)
        mock_conn.commit.assert_called_once()
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_query_error(self, mock_pool):
        """Test error handling in execute_query."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.execute.side_effect = psycopg2.Error("Test error")
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute query with error
        with self.assertRaises(psycopg2.Error):
            instance.execute_query("SELECT 1")
        
        # Verify error was handled
        mock_conn.rollback.assert_called_once()
        self.assertEqual(instance.error_count, 1)
        self.assertEqual(instance.last_error, "Test error")
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_batch(self, mock_pool):
        """Test executing a batch of queries."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute batch
        queries = [
            "CREATE TABLE test (id INT)",
            "INSERT INTO test VALUES (1)",
            "INSERT INTO test VALUES (2)"
        ]
        result = instance.execute_batch(queries)
        
        # Verify batch was executed
        self.assertTrue(result)
        self.assertEqual(mock_cursor.execute.call_count, 3)
        mock_cursor.execute.assert_has_calls([
            call("CREATE TABLE test (id INT)"),
            call("INSERT INTO test VALUES (1)"),
            call("INSERT INTO test VALUES (2)")
        ])
        mock_conn.commit.assert_called_once()
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_batch_no_transaction(self, mock_pool):
        """Test executing a batch of queries without a transaction."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute batch without transaction
        queries = [
            "CREATE TABLE test (id INT)",
            "INSERT INTO test VALUES (1)"
        ]
        result = instance.execute_batch(queries, transaction=False)
        
        # Verify batch was executed without transaction
        self.assertTrue(result)
        self.assertEqual(mock_cursor.execute.call_count, 2)
        mock_conn.commit.assert_not_called()
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_batch_error_with_transaction(self, mock_pool):
        """Test error handling in execute_batch with transaction."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.execute.side_effect = [None, psycopg2.Error("Test error")]
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute batch with error
        queries = [
            "CREATE TABLE test (id INT)",
            "INSERT INTO test VALUES (1)"
        ]
        with self.assertRaises(psycopg2.Error):
            instance.execute_batch(queries)
        
        # Verify error was handled
        mock_conn.rollback.assert_called_once()
        self.assertEqual(instance.error_count, 1)
        self.assertEqual(instance.last_error, "Test error")
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_execute_batch_error_without_transaction(self, mock_pool):
        """Test error handling in execute_batch without transaction."""
        # Set up mock connection and cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.execute.side_effect = [None, psycopg2.Error("Test error")]
        mock_conn.cursor.return_value = mock_cursor
        mock_pool.return_value.getconn.return_value = mock_conn
        
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Execute batch with error but without transaction
        queries = [
            "CREATE TABLE test (id INT)",
            "INSERT INTO test VALUES (1)"
        ]
        result = instance.execute_batch(queries, transaction=False)
        
        # Verify error was handled but no exception was raised
        self.assertFalse(result)
        mock_conn.rollback.assert_not_called()
        self.assertEqual(instance.error_count, 1)
        self.assertEqual(instance.last_error, "Test error")
        mock_pool.return_value.putconn.assert_called_once_with(mock_conn)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_get_stats(self, mock_pool):
        """Test getting statistics."""
        # Create instance
        instance = SupabaseDirectConnection.get_instance()
        
        # Set some stats
        instance.query_count = 10
        instance.error_count = 2
        instance.last_error = "Test error"
        instance.last_query_time = 0.5
        
        # Get stats
        stats = instance.get_stats()
        
        # Verify stats
        self.assertEqual(stats['min_connections'], instance.min_connections)
        self.assertEqual(stats['max_connections'], instance.max_connections)
        self.assertEqual(stats['query_count'], 10)
        self.assertEqual(stats['error_count'], 2)
        self.assertEqual(stats['last_error'], "Test error")
        self.assertEqual(stats['last_query_time'], 0.5)

if __name__ == '__main__':
    unittest.main()
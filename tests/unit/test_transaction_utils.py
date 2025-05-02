#!/usr/bin/env python3
"""
Unit tests for transaction_utils.py

Tests the robust transaction management utilities including:
1. Safe transaction context manager
2. Transaction retry logic
3. Detection and handling of aborted transactions
4. Verification of transaction cleanup
"""

import unittest
from unittest.mock import MagicMock, patch, call
import threading
import time
import psycopg2

from transaction_utils import (
    safe_transaction,
    with_transaction_retry,
    execute_in_transaction,
    is_transaction_active
)

class MockConnection:
    """Mock database connection for testing transaction utilities."""
    
    def __init__(self, fail_on_execute=False, fail_on_commit=False, fail_on_rollback=False):
        self.autocommit = True
        self.fail_on_execute = fail_on_execute
        self.fail_on_commit = fail_on_commit
        self.fail_on_rollback = fail_on_rollback
        self.closed = False
        self.transaction_active = False
        self.cursor_mock = MagicMock()
        self.execute_query_called = False
        self.begin_transaction_called = False
        self.commit_transaction_called = False
        self.rollback_transaction_called = False
        self.commit_called = False
        self.rollback_called = False
        
    def cursor(self):
        """Return a mock cursor."""
        return self.cursor_mock
        
    def execute_query(self, query, params=None):
        """Mock execute_query method."""
        self.execute_query_called = True
        if self.fail_on_execute:
            raise Exception("Mock execute query failure")
        if "current transaction is aborted" in query:
            raise Exception("current transaction is aborted")
        return [{"result": 1}]
        
    def begin_transaction(self):
        """Mock begin_transaction method."""
        self.begin_transaction_called = True
        self.autocommit = False
        self.transaction_active = True
        return "mock_transaction"
        
    def commit_transaction(self, transaction):
        """Mock commit_transaction method."""
        self.commit_transaction_called = True
        if self.fail_on_commit:
            raise Exception("Mock commit failure")
        self.autocommit = True
        self.transaction_active = False
        
    def rollback_transaction(self, transaction):
        """Mock rollback_transaction method."""
        self.rollback_transaction_called = True
        if self.fail_on_rollback:
            raise Exception("Mock rollback failure")
        self.autocommit = True
        self.transaction_active = False
        
    def commit(self):
        """Mock commit method."""
        self.commit_called = True
        if self.fail_on_commit:
            raise Exception("Mock commit failure")
        self.autocommit = True
        self.transaction_active = False
        
    def rollback(self):
        """Mock rollback method."""
        self.rollback_called = True
        if self.fail_on_rollback:
            raise Exception("Mock rollback failure")
        self.autocommit = True
        self.transaction_active = False
        
    def close(self):
        """Mock close method."""
        self.closed = True


class TestTransactionUtils(unittest.TestCase):
    """Test cases for transaction_utils.py."""
    
    @patch('transaction_utils.get_db_connection')
    def test_safe_transaction_success(self, mock_get_db_connection):
        """Test safe_transaction context manager with successful transaction."""
        # Setup mock connection
        mock_conn = MockConnection()
        mock_get_db_connection.return_value = mock_conn
        
        # Use safe_transaction context manager
        with safe_transaction() as conn:
            self.assertEqual(conn, mock_conn)
            # Simulate some database operations
            conn.execute_query("SELECT 1")
        
        # Verify transaction was properly managed
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertTrue(mock_conn.commit_transaction_called)
        self.assertFalse(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_safe_transaction_exception(self, mock_get_db_connection):
        """Test safe_transaction context manager with exception during transaction."""
        # Setup mock connection
        mock_conn = MockConnection()
        mock_get_db_connection.return_value = mock_conn
        
        # Use safe_transaction context manager with an exception
        with self.assertRaises(ValueError):
            with safe_transaction() as conn:
                self.assertEqual(conn, mock_conn)
                # Simulate some database operations
                conn.execute_query("SELECT 1")
                # Raise an exception
                raise ValueError("Test exception")
        
        # Verify transaction was properly rolled back
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertFalse(mock_conn.commit_transaction_called)
        self.assertTrue(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_safe_transaction_aborted_transaction(self, mock_get_db_connection):
        """Test safe_transaction context manager with aborted transaction detection."""
        # Setup mock connection that fails on execute_query
        mock_conn = MockConnection(fail_on_execute=True)
        mock_get_db_connection.return_value = mock_conn
        
        # Use safe_transaction context manager
        with self.assertRaises(Exception):
            with safe_transaction() as conn:
                self.assertEqual(conn, mock_conn)
                # This will fail due to fail_on_execute=True
                conn.execute_query("SELECT 1")
        
        # Verify transaction was properly rolled back
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertFalse(mock_conn.commit_transaction_called)
        self.assertTrue(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_safe_transaction_commit_failure(self, mock_get_db_connection):
        """Test safe_transaction context manager with commit failure."""
        # Setup mock connection that fails on commit
        mock_conn = MockConnection(fail_on_commit=True)
        mock_get_db_connection.return_value = mock_conn
        
        # Use safe_transaction context manager
        with self.assertRaises(Exception):
            with safe_transaction() as conn:
                self.assertEqual(conn, mock_conn)
                # Simulate some database operations
                conn.execute_query("SELECT 1")
        
        # Verify transaction was properly rolled back after commit failure
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertTrue(mock_conn.commit_transaction_called)
        self.assertTrue(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_with_transaction_retry_success_first_attempt(self, mock_get_db_connection):
        """Test with_transaction_retry decorator with success on first attempt."""
        # Setup mock connection
        mock_conn = MockConnection()
        mock_get_db_connection.return_value = mock_conn
        
        # Define a function to decorate
        @with_transaction_retry(max_retries=3)
        def test_func(conn, arg1, arg2=None):
            conn.execute_query("SELECT 1")
            return f"{arg1}-{arg2}"
        
        # Call the decorated function
        result = test_func("test", arg2="value")
        
        # Verify result and transaction management
        self.assertEqual(result, "test-value")
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertTrue(mock_conn.commit_transaction_called)
        self.assertFalse(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_with_transaction_retry_success_after_retry(self, mock_get_db_connection):
        """Test with_transaction_retry decorator with success after retry."""
        # Setup mock connections - first fails, second succeeds
        mock_conn_fail = MockConnection(fail_on_execute=True)
        mock_conn_success = MockConnection()
        mock_get_db_connection.side_effect = [mock_conn_fail, mock_conn_success]
        
        # Define a function to decorate
        @with_transaction_retry(max_retries=3, retry_delay=0.01)
        def test_func(conn, arg):
            conn.execute_query("SELECT 1")
            return arg
        
        # Call the decorated function
        result = test_func("test_value")
        
        # Verify result and transaction management
        self.assertEqual(result, "test_value")
        self.assertTrue(mock_conn_fail.begin_transaction_called)
        self.assertFalse(mock_conn_fail.commit_transaction_called)
        self.assertTrue(mock_conn_fail.rollback_transaction_called)
        self.assertTrue(mock_conn_fail.closed)
        
        self.assertTrue(mock_conn_success.begin_transaction_called)
        self.assertTrue(mock_conn_success.commit_transaction_called)
        self.assertFalse(mock_conn_success.rollback_transaction_called)
        self.assertTrue(mock_conn_success.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_with_transaction_retry_all_attempts_fail(self, mock_get_db_connection):
        """Test with_transaction_retry decorator with all attempts failing."""
        # Setup mock connections - all fail
        mock_conn1 = MockConnection(fail_on_execute=True)
        mock_conn2 = MockConnection(fail_on_execute=True)
        mock_conn3 = MockConnection(fail_on_execute=True)
        mock_get_db_connection.side_effect = [mock_conn1, mock_conn2, mock_conn3]
        
        # Define a function to decorate
        @with_transaction_retry(max_retries=3, retry_delay=0.01)
        def test_func(conn, arg):
            conn.execute_query("SELECT 1")
            return arg
        
        # Call the decorated function
        with self.assertRaises(Exception):
            test_func("test_value")
        
        # Verify all connections were properly managed
        for conn in [mock_conn1, mock_conn2, mock_conn3]:
            self.assertTrue(conn.begin_transaction_called)
            self.assertFalse(conn.commit_transaction_called)
            self.assertTrue(conn.rollback_transaction_called)
            self.assertTrue(conn.closed)
    
    @patch('transaction_utils.get_db_connection')
    def test_execute_in_transaction(self, mock_get_db_connection):
        """Test execute_in_transaction function."""
        # Setup mock connection
        mock_conn = MockConnection()
        mock_get_db_connection.return_value = mock_conn
        
        # Define a function to execute in transaction
        def test_func(conn, arg1, arg2=None):
            conn.execute_query("SELECT 1")
            return f"{arg1}-{arg2}"
        
        # Call execute_in_transaction
        result = execute_in_transaction(test_func, "test", arg2="value")
        
        # Verify result and transaction management
        self.assertEqual(result, "test-value")
        self.assertTrue(mock_conn.begin_transaction_called)
        self.assertTrue(mock_conn.commit_transaction_called)
        self.assertFalse(mock_conn.rollback_transaction_called)
        self.assertTrue(mock_conn.closed)
    
    def test_is_transaction_active_true(self):
        """Test is_transaction_active function with active transaction."""
        # Setup mock connection with active transaction
        mock_conn = MockConnection()
        mock_conn.autocommit = False
        mock_conn.transaction_active = True
        
        # Check if transaction is active
        self.assertTrue(is_transaction_active(mock_conn))
    
    def test_is_transaction_active_false(self):
        """Test is_transaction_active function with no active transaction."""
        # Setup mock connection with no active transaction
        mock_conn = MockConnection()
        mock_conn.autocommit = True
        mock_conn.transaction_active = False
        
        # Check if transaction is active
        self.assertFalse(is_transaction_active(mock_conn))
    
    def test_is_transaction_active_error(self):
        """Test is_transaction_active function with error during check."""
        # Setup mock connection that fails on execute
        mock_conn = MockConnection(fail_on_execute=True)
        
        # Check if transaction is active - should return False for unknown state
        self.assertFalse(is_transaction_active(mock_conn))


if __name__ == '__main__':
    unittest.main()
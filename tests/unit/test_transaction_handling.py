#!/usr/bin/env python3
"""
Unit tests for transaction handling in sql_executor.py.

These tests verify that the transaction handling works correctly with different
adapter types, including PoolerAdapter.
"""

import unittest
import logging
import sys
import os

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the modules to test
import sql_executor

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MockAdapter:
    """Base mock adapter with transaction tracking."""
    
    def __init__(self):
        self.transactions = []
        self.commits = []
        self.rollbacks = []
    
    def begin_transaction(self):
        logger.info(f"{self.__class__.__name__}: begin_transaction called")
        tx_id = f"tx-{len(self.transactions) + 1}"
        self.transactions.append(tx_id)
        return tx_id
    
    def commit_transaction(self, transaction):
        logger.info(f"{self.__class__.__name__}: commit_transaction called with {transaction}")
        self.commits.append(transaction)
        return True
    
    def rollback_transaction(self, transaction):
        logger.info(f"{self.__class__.__name__}: rollback_transaction called with {transaction}")
        self.rollbacks.append(transaction)
        return True
    
    def execute_query(self, query, params=None):
        logger.info(f"{self.__class__.__name__}: execute_query called with {query}")
        return []


class MockPoolerAdapter(MockAdapter):
    """Mock pooler adapter that implements the ConnectionAdapter interface."""
    pass


class MockFailingAdapter(MockAdapter):
    """Mock adapter that fails during rollback."""
    
    def rollback_transaction(self, transaction):
        logger.info(f"{self.__class__.__name__}: rollback_transaction called with {transaction}")
        self.rollbacks.append(transaction)
        raise RuntimeError("Simulated rollback failure")


class TestTransactionHandling(unittest.TestCase):
    """Test cases for transaction handling in sql_executor.py."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Save original get_db function
        self.original_get_db = sql_executor.get_db
    
    def tearDown(self):
        """Tear down test fixtures."""
        # Restore original get_db function
        sql_executor.get_db = self.original_get_db
    
    def test_with_transaction_decorator_success(self):
        """Test that the with_transaction decorator commits on success."""
        # Create a mock adapter
        mock_adapter = MockAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Define a function to decorate
        @sql_executor.with_transaction
        def test_func(transaction, arg1, arg2=None):
            self.assertEqual(arg1, "value1")
            self.assertEqual(arg2, "value2")
            return "result"
        
        # Call the decorated function
        result = test_func("value1", arg2="value2")
        
        # Verify the result
        self.assertEqual(result, "result")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 1)
        self.assertEqual(len(mock_adapter.rollbacks), 0)
        self.assertEqual(mock_adapter.commits[0], mock_adapter.transactions[0])
    
    def test_with_transaction_decorator_exception(self):
        """Test that the with_transaction decorator rolls back on exception."""
        # Create a mock adapter
        mock_adapter = MockAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Define a function to decorate that raises an exception
        @sql_executor.with_transaction
        def test_func(transaction, arg1):
            self.assertEqual(arg1, "value1")
            raise ValueError("Test exception")
        
        # Call the decorated function and expect an exception
        with self.assertRaises(ValueError):
            test_func("value1")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 0)
        self.assertEqual(len(mock_adapter.rollbacks), 1)
        self.assertEqual(mock_adapter.rollbacks[0], mock_adapter.transactions[0])
    
    def test_with_transaction_decorator_pooler_adapter(self):
        """Test that the with_transaction decorator works with PoolerAdapter."""
        # Create a mock pooler adapter
        mock_adapter = MockPoolerAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Define a function to decorate
        @sql_executor.with_transaction
        def test_func(transaction, arg1):
            self.assertEqual(arg1, "value1")
            return "result"
        
        # Call the decorated function
        result = test_func("value1")
        
        # Verify the result
        self.assertEqual(result, "result")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 1)
        self.assertEqual(len(mock_adapter.rollbacks), 0)
        self.assertEqual(mock_adapter.commits[0], mock_adapter.transactions[0])
    
    def test_with_transaction_decorator_rollback_failure(self):
        """Test that the with_transaction decorator handles rollback failures."""
        # Create a mock failing adapter
        mock_adapter = MockFailingAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Define a function to decorate that raises an exception
        @sql_executor.with_transaction
        def test_func(transaction, arg1):
            self.assertEqual(arg1, "value1")
            raise ValueError("Test exception")
        
        # Call the decorated function and expect the original exception
        with self.assertRaises(ValueError):
            test_func("value1")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 0)
        self.assertEqual(len(mock_adapter.rollbacks), 1)
        self.assertEqual(mock_adapter.rollbacks[0], mock_adapter.transactions[0])
    
    def test_transaction_context_manager_success(self):
        """Test that the transaction context manager commits on success."""
        # Create a mock adapter
        mock_adapter = MockAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Use the transaction context manager
        with sql_executor.transaction() as transaction:
            self.assertEqual(transaction, mock_adapter.transactions[0])
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 1)
        self.assertEqual(len(mock_adapter.rollbacks), 0)
        self.assertEqual(mock_adapter.commits[0], mock_adapter.transactions[0])
    
    def test_transaction_context_manager_exception(self):
        """Test that the transaction context manager rolls back on exception."""
        # Create a mock adapter
        mock_adapter = MockAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Use the transaction context manager with an exception
        with self.assertRaises(ValueError):
            with sql_executor.transaction() as transaction:
                self.assertEqual(transaction, mock_adapter.transactions[0])
                raise ValueError("Test exception")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 0)
        self.assertEqual(len(mock_adapter.rollbacks), 1)
        self.assertEqual(mock_adapter.rollbacks[0], mock_adapter.transactions[0])
    
    def test_transaction_context_manager_pooler_adapter(self):
        """Test that the transaction context manager works with PoolerAdapter."""
        # Create a mock pooler adapter
        mock_adapter = MockPoolerAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Use the transaction context manager
        with sql_executor.transaction() as transaction:
            self.assertEqual(transaction, mock_adapter.transactions[0])
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 1)
        self.assertEqual(len(mock_adapter.rollbacks), 0)
        self.assertEqual(mock_adapter.commits[0], mock_adapter.transactions[0])
    
    def test_transaction_context_manager_rollback_failure(self):
        """Test that the transaction context manager handles rollback failures."""
        # Create a mock failing adapter
        mock_adapter = MockFailingAdapter()
        
        # Replace get_db function
        sql_executor.get_db = lambda: mock_adapter
        
        # Use the transaction context manager with an exception
        with self.assertRaises(ValueError):
            with sql_executor.transaction() as transaction:
                self.assertEqual(transaction, mock_adapter.transactions[0])
                raise ValueError("Test exception")
        
        # Verify transaction operations
        self.assertEqual(len(mock_adapter.transactions), 1)
        self.assertEqual(len(mock_adapter.commits), 0)
        self.assertEqual(len(mock_adapter.rollbacks), 1)
        self.assertEqual(mock_adapter.rollbacks[0], mock_adapter.transactions[0])


if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python3
"""
Unit tests for db_connection_utils.py

This module contains tests for the resilient database connection utilities,
including the CircuitBreaker pattern, ConnectionManager, and context managers.
"""

import os
import time
import unittest
from unittest.mock import patch, MagicMock, call
import threading
import psycopg2
from psycopg2.pool import ThreadedConnectionPool

# Import the module to test
import db_connection_utils
from db_connection_utils import (
    CircuitBreaker,
    ConnectionManager,
    ConnectionWrapper,
    MCPConnection,
    get_db_connection,
    safe_transaction
)


class TestCircuitBreaker(unittest.TestCase):
    """Tests for the CircuitBreaker class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.circuit_breaker = CircuitBreaker(
            name="test_breaker",
            failure_threshold=2,
            recovery_timeout=0.1  # Short timeout for testing
        )
    
    def test_initial_state(self):
        """Test the initial state of the circuit breaker."""
        self.assertEqual(self.circuit_breaker.state, "CLOSED")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        self.assertEqual(self.circuit_breaker.name, "test_breaker")
        self.assertEqual(self.circuit_breaker.failure_threshold, 2)
        self.assertEqual(self.circuit_breaker.recovery_timeout, 0.1)
    
    def test_successful_execution(self):
        """Test successful function execution."""
        # Define a test function
        def test_func(x, y):
            return x + y
        
        # Execute the function through the circuit breaker
        result = self.circuit_breaker.execute(test_func, 2, 3)
        
        # Check the result and state
        self.assertEqual(result, 5)
        self.assertEqual(self.circuit_breaker.state, "CLOSED")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        self.assertEqual(self.circuit_breaker.stats["success_count"], 1)
    
    def test_failed_execution(self):
        """Test failed function execution."""
        # Define a test function that raises an exception
        def failing_func():
            raise ValueError("Test error")
        
        # Execute the function through the circuit breaker
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Check the state after one failure
        self.assertEqual(self.circuit_breaker.state, "CLOSED")
        self.assertEqual(self.circuit_breaker.failure_count, 1)
        self.assertEqual(self.circuit_breaker.stats["failure_count"], 1)
        
        # Execute the function again to trigger the circuit breaker
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Check the state after reaching the failure threshold
        self.assertEqual(self.circuit_breaker.state, "OPEN")
        self.assertEqual(self.circuit_breaker.failure_count, 2)
        self.assertEqual(self.circuit_breaker.stats["failure_count"], 2)
        self.assertEqual(self.circuit_breaker.stats["open_count"], 1)
    
    def test_open_circuit(self):
        """Test behavior when circuit is open."""
        # First, open the circuit
        def failing_func():
            raise ValueError("Test error")
        
        # Execute the function twice to open the circuit
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Now try to execute a function when the circuit is open
        with self.assertRaises(ConnectionError):
            self.circuit_breaker.execute(lambda: True)
    
    def test_half_open_circuit(self):
        """Test transition to half-open state after recovery timeout."""
        # First, open the circuit
        def failing_func():
            raise ValueError("Test error")
        
        # Execute the function twice to open the circuit
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Wait for recovery timeout
        time.sleep(0.2)  # Longer than recovery_timeout
        
        # Execute a successful function to close the circuit
        result = self.circuit_breaker.execute(lambda: True)
        
        # Check the result and state
        self.assertTrue(result)
        self.assertEqual(self.circuit_breaker.state, "CLOSED")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
    
    def test_half_open_to_open(self):
        """Test transition from half-open to open state on failure."""
        # First, open the circuit
        def failing_func():
            raise ValueError("Test error")
        
        # Execute the function twice to open the circuit
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Wait for recovery timeout
        time.sleep(0.2)  # Longer than recovery_timeout
        
        # Execute a failing function to reopen the circuit
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Check the state
        self.assertEqual(self.circuit_breaker.state, "OPEN")
        self.assertEqual(self.circuit_breaker.stats["open_count"], 2)
    
    def test_reset(self):
        """Test manual reset of the circuit breaker."""
        # First, open the circuit
        def failing_func():
            raise ValueError("Test error")
        
        # Execute the function twice to open the circuit
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        with self.assertRaises(ValueError):
            self.circuit_breaker.execute(failing_func)
        
        # Reset the circuit breaker
        self.circuit_breaker.reset()
        
        # Check the state after reset
        self.assertEqual(self.circuit_breaker.state, "CLOSED")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        
        # Execute a function after reset
        result = self.circuit_breaker.execute(lambda: True)
        self.assertTrue(result)


class TestConnectionManager(unittest.TestCase):
    """Tests for the ConnectionManager class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Reset the singleton instance
        ConnectionManager._instance = None
        
        # Mock environment variables
        self.env_patcher = patch.dict(os.environ, {
            'POOLER_ENABLED': 'true',
            'SUPABASE_DB_HOST': 'test-host',
            'SUPABASE_DB_PORT': '6543',
            'SUPABASE_DB_NAME': 'test-db',
            'SUPABASE_DB_USER': 'test-user',
            'SUPABASE_DB_PASSWORD': 'test-password',
            'SUPABASE_DB_IP_ADDRESS': '127.0.0.1',
            'DB_CONNECTION_ORDER': 'pooler,direct,ip,mcp'
        })
        self.env_patcher.start()
        
        # Create the connection manager
        self.manager = ConnectionManager.get_instance()
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.env_patcher.stop()
    
    def test_singleton_pattern(self):
        """Test that ConnectionManager is a singleton."""
        manager1 = ConnectionManager.get_instance()
        manager2 = ConnectionManager.get_instance()
        
        self.assertIs(manager1, manager2)
        
        # Ensure we can't create a new instance directly
        with self.assertRaises(RuntimeError):
            ConnectionManager()
    
    def test_load_config(self):
        """Test loading configuration from environment variables."""
        config = self.manager.config
        
        # Check pooler config
        self.assertTrue(config["pooler"]["enabled"])
        self.assertEqual(config["pooler"]["host"], "test-host")
        self.assertEqual(config["pooler"]["port"], 6543)
        self.assertEqual(config["pooler"]["dbname"], "test-db")
        self.assertEqual(config["pooler"]["user"], "test-user")
        self.assertEqual(config["pooler"]["password"], "test-password")
        
        # Check connection order
        self.assertEqual(self.manager.connection_order, ["pooler", "direct", "ip", "mcp"])
    
    @patch('socket.gethostbyname')
    def test_resolve_hostname(self, mock_gethostbyname):
        """Test hostname resolution."""
        mock_gethostbyname.return_value = "192.168.1.1"
        
        # Test successful resolution
        ip = self.manager._resolve_hostname("test-host")
        self.assertEqual(ip, "192.168.1.1")
        mock_gethostbyname.assert_called_with("test-host")
        
        # Test failed resolution
        mock_gethostbyname.side_effect = socket.gaierror("Test error")
        ip = self.manager._resolve_hostname("invalid-host")
        self.assertIsNone(ip)
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_create_connection_pool(self, mock_pool_class):
        """Test creating a connection pool."""
        # Mock the pool creation
        mock_pool = MagicMock()
        mock_pool_class.return_value = mock_pool
        
        # Test creating a pooler connection pool
        result = self.manager._create_connection_pool("pooler")
        
        # Check the result
        self.assertTrue(result)
        self.assertIn("pooler", self.manager.pools)
        self.assertEqual(self.manager.pools["pooler"]["pool"], mock_pool)
        
        # Check the pool creation arguments
        mock_pool_class.assert_called_with(
            minconn=1,
            maxconn=10,
            host="test-host",
            port=6543,
            dbname="test-db",
            user="test-user",
            password="test-password",
            connect_timeout=10
        )
    
    @patch('db_connection_utils.ConnectionManager._create_connection_pool')
    @patch('db_connection_utils.ConnectionManager._get_mcp_connection')
    def test_get_connection(self, mock_get_mcp, mock_create_pool):
        """Test getting a connection."""
        # Mock the pool creation
        mock_create_pool.return_value = True
        
        # Mock the connection pool
        mock_pool = MagicMock()
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.fetchone.return_value = (1,)
        mock_conn.cursor.return_value.__enter__.return_value = mock_cursor
        mock_pool.getconn.return_value = mock_conn
        
        # Set up the pools
        self.manager.pools["pooler"] = {
            "pool": mock_pool,
            "created_at": time.time(),
            "last_used": time.time(),
            "stats": {
                "connections_acquired": 0,
                "connections_released": 0,
                "errors": 0
            }
        }
        
        # Test getting a connection
        connection = self.manager.get_connection()
        
        # Check the result
        self.assertIsNotNone(connection)
        self.assertEqual(self.manager.active_pool, "pooler")
        self.assertEqual(self.manager.stats["successful_connections"], 1)
        
        # Check that the connection was acquired from the pool
        mock_pool.getconn.assert_called_once()
        mock_conn.cursor.assert_called_once()
        mock_cursor.execute.assert_called_with("SELECT 1")
    
    def test_get_stats(self):
        """Test getting connection statistics."""
        # Set up some test data
        self.manager.pools["pooler"] = {
            "pool": MagicMock(),
            "created_at": time.time(),
            "last_used": time.time(),
            "stats": {
                "connections_acquired": 5,
                "connections_released": 4,
                "errors": 1
            }
        }
        
        self.manager.active_pool = "pooler"
        self.manager.stats["total_connections"] = 10
        self.manager.stats["successful_connections"] = 8
        self.manager.stats["failed_connections"] = 2
        
        # Get the stats
        stats = self.manager.get_stats()
        
        # Check the stats
        self.assertEqual(stats["total_connections"], 10)
        self.assertEqual(stats["successful_connections"], 8)
        self.assertEqual(stats["failed_connections"], 2)
        self.assertEqual(stats["active_pool"], "pooler")
        self.assertIn("pooler", stats["pools"])
        self.assertEqual(stats["pools"]["pooler"]["stats"]["connections_acquired"], 5)
        self.assertEqual(stats["pools"]["pooler"]["stats"]["connections_released"], 4)
        self.assertEqual(stats["pools"]["pooler"]["stats"]["errors"], 1)


class TestConnectionWrapper(unittest.TestCase):
    """Tests for the ConnectionWrapper class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.mock_manager = MagicMock()
        self.mock_connection = MagicMock()
        self.wrapper = ConnectionWrapper(
            connection=self.mock_connection,
            manager=self.mock_manager,
            pool_type="pooler"
        )
    
    def test_close(self):
        """Test closing the connection."""
        # Close the connection
        self.wrapper.close()
        
        # Check that the connection was released
        self.mock_manager.release_connection.assert_called_with(
            self.mock_connection, "pooler"
        )
        
        # Check that the connection is marked as closed
        self.assertTrue(self.wrapper.closed)
        
        # Closing again should not release the connection again
        self.mock_manager.release_connection.reset_mock()
        self.wrapper.close()
        self.mock_manager.release_connection.assert_not_called()
    
    def test_cursor(self):
        """Test getting a cursor."""
        # Get a cursor
        self.wrapper.cursor(cursor_factory="test_factory")
        
        # Check that the cursor was requested from the connection
        self.mock_connection.cursor.assert_called_with(cursor_factory="test_factory")
    
    def test_commit(self):
        """Test committing a transaction."""
        # Commit the transaction
        self.wrapper.commit()
        
        # Check that the commit was called on the connection
        self.mock_connection.commit.assert_called_once()
    
    def test_rollback(self):
        """Test rolling back a transaction."""
        # Roll back the transaction
        self.wrapper.rollback()
        
        # Check that the rollback was called on the connection
        self.mock_connection.rollback.assert_called_once()
    
    def test_execute_query_select(self):
        """Test executing a SELECT query."""
        # Mock the cursor
        mock_cursor = MagicMock()
        mock_cursor.fetchall.return_value = [{"id": 1, "name": "test"}]
        self.mock_connection.cursor.return_value.__enter__.return_value = mock_cursor
        
        # Execute a SELECT query
        result = self.wrapper.execute_query("SELECT * FROM test")
        
        # Check the result
        self.assertEqual(result, [{"id": 1, "name": "test"}])
        
        # Check that the query was executed
        mock_cursor.execute.assert_called_with("SELECT * FROM test", None)
        mock_cursor.fetchall.assert_called_once()
        
        # Check that commit was not called (for SELECT)
        self.mock_connection.commit.assert_not_called()
    
    def test_execute_query_insert(self):
        """Test executing an INSERT query."""
        # Mock the cursor
        mock_cursor = MagicMock()
        mock_cursor.rowcount = 1
        self.mock_connection.cursor.return_value.__enter__.return_value = mock_cursor
        
        # Execute an INSERT query
        result = self.wrapper.execute_query("INSERT INTO test VALUES (1, 'test')")
        
        # Check the result
        self.assertEqual(result, 1)
        
        # Check that the query was executed
        mock_cursor.execute.assert_called_with("INSERT INTO test VALUES (1, 'test')", None)
        mock_cursor.fetchall.assert_not_called()
        
        # Check that commit was called (for INSERT)
        self.mock_connection.commit.assert_called_once()
    
    def test_execute_query_error(self):
        """Test error handling in execute_query."""
        # Mock the cursor to raise an exception
        mock_cursor = MagicMock()
        mock_cursor.execute.side_effect = psycopg2.Error("Test error")
        self.mock_connection.cursor.return_value.__enter__.return_value = mock_cursor
        
        # Execute a query that will fail
        with self.assertRaises(psycopg2.Error):
            self.wrapper.execute_query("SELECT * FROM test")
        
        # Check that rollback was called
        self.mock_connection.rollback.assert_called_once()


class TestContextManagers(unittest.TestCase):
    """Tests for the context manager functions."""
    
    @patch('db_connection_utils.ConnectionManager.get_instance')
    def test_get_db_connection(self, mock_get_instance):
        """Test the get_db_connection context manager."""
        # Mock the connection manager
        mock_manager = MagicMock()
        mock_get_instance.return_value = mock_manager
        
        # Mock the connection
        mock_connection = MagicMock()
        mock_manager.get_connection.return_value = mock_connection
        
        # Use the context manager
        with get_db_connection() as conn:
            # Check that we got the connection
            self.assertEqual(conn, mock_connection)
            
            # Do something with the connection
            conn.execute_query("SELECT 1")
        
        # Check that the connection was closed
        mock_connection.close.assert_called_once()
    
    @patch('db_connection_utils.ConnectionManager.get_instance')
    def test_get_db_connection_error(self, mock_get_instance):
        """Test error handling in get_db_connection."""
        # Mock the connection manager to return None
        mock_manager = MagicMock()
        mock_manager.get_connection.return_value = None
        mock_get_instance.return_value = mock_manager
        
        # Use the context manager
        with self.assertRaises(ConnectionError):
            with get_db_connection() as conn:
                pass
    
    @patch('db_connection_utils.ConnectionManager.get_instance')
    def test_safe_transaction(self, mock_get_instance):
        """Test the safe_transaction context manager."""
        # Mock the connection manager
        mock_manager = MagicMock()
        mock_get_instance.return_value = mock_manager
        
        # Mock the connection
        mock_connection = MagicMock()
        mock_manager.get_connection.return_value = mock_connection
        
        # Use the context manager
        with safe_transaction() as conn:
            # Check that we got the connection
            self.assertEqual(conn, mock_connection)
            
            # Do something with the connection
            conn.execute_query("INSERT INTO test VALUES (1, 'test')")
        
        # Check that the transaction was committed
        mock_connection.commit.assert_called_once()
        
        # Check that the connection was closed
        mock_connection.close.assert_called_once()
    
    @patch('db_connection_utils.ConnectionManager.get_instance')
    def test_safe_transaction_error(self, mock_get_instance):
        """Test error handling in safe_transaction."""
        # Mock the connection manager
        mock_manager = MagicMock()
        mock_get_instance.return_value = mock_manager
        
        # Mock the connection
        mock_connection = MagicMock()
        mock_manager.get_connection.return_value = mock_connection
        
        # Use the context manager with an error
        with self.assertRaises(ValueError):
            with safe_transaction() as conn:
                # Do something that raises an error
                raise ValueError("Test error")
        
        # Check that the transaction was rolled back
        mock_connection.rollback.assert_called_once()
        
        # Check that the connection was closed
        mock_connection.close.assert_called_once()


if __name__ == '__main__':
    unittest.main()
"""
Unit tests for the enhanced connection pool.

This module tests the functionality of the enhanced connection pool with
advanced features like health checks, dynamic sizing, and connection validation.
"""

import unittest
import logging
import time
import asyncio
from unittest.mock import MagicMock, AsyncMock, patch, call
import threading
from typing import List, Dict, Any, Tuple

# Import module under test
from unified_importer.core.enhanced_connection_pool import (
    EnhancedConnectionPool, EnhancedConnectionInfo, ConnectionState, CircuitBreaker
)


class TestCircuitBreaker(unittest.TestCase):
    """Test cases for the CircuitBreaker class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger('test_circuit_breaker')
        self.logger.setLevel(logging.DEBUG)
        
        # Create circuit breaker with small thresholds for testing
        self.circuit_breaker = CircuitBreaker(
            failure_threshold=3,
            reset_timeout=1.0,
            half_open_timeout=0.5,
            logger=self.logger
        )
    
    def test_initial_state(self):
        """Test that the circuit breaker starts in the closed state."""
        self.assertEqual(self.circuit_breaker.state, "closed")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        self.assertTrue(self.circuit_breaker.allow_request())
    
    def test_failure_tracking(self):
        """Test that failures are properly tracked."""
        # Record failures
        self.circuit_breaker.record_failure()
        self.assertEqual(self.circuit_breaker.failure_count, 1)
        self.assertEqual(self.circuit_breaker.state, "closed")
        
        self.circuit_breaker.record_failure()
        self.assertEqual(self.circuit_breaker.failure_count, 2)
        self.assertEqual(self.circuit_breaker.state, "closed")
        
        # This should open the circuit
        self.circuit_breaker.record_failure()
        self.assertEqual(self.circuit_breaker.failure_count, 3)
        self.assertEqual(self.circuit_breaker.state, "open")
    
    def test_circuit_transitions(self):
        """Test circuit state transitions."""
        # Open the circuit
        for _ in range(3):
            self.circuit_breaker.record_failure()
        
        self.assertEqual(self.circuit_breaker.state, "open")
        self.assertFalse(self.circuit_breaker.allow_request())
        
        # Wait for reset timeout
        time.sleep(1.1)
        
        # Circuit should be half-open now
        self.assertTrue(self.circuit_breaker.allow_request())
        self.assertEqual(self.circuit_breaker.state, "half-open")
        
        # Success in half-open state should close the circuit
        self.circuit_breaker.record_success()
        self.assertEqual(self.circuit_breaker.state, "closed")
        self.assertEqual(self.circuit_breaker.failure_count, 0)
        
        # Failure in half-open state should reopen the circuit
        for _ in range(3):
            self.circuit_breaker.record_failure()
        
        time.sleep(1.1)  # Wait for reset timeout
        self.assertEqual(self.circuit_breaker.state, "half-open")
        
        self.circuit_breaker.record_failure()
        self.assertEqual(self.circuit_breaker.state, "open")


class TestEnhancedConnectionPool(unittest.TestCase):
    """Test cases for the EnhancedConnectionPool class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger('test_enhanced_connection_pool')
        self.logger.setLevel(logging.DEBUG)
        
        # Create mock connection params
        self.connection_params = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass'
        }
        
        # Patch the _create_connection method to return a mock connection
        self.patcher = patch.object(
            EnhancedConnectionPool, '_create_connection', 
            return_value=(MagicMock(), EnhancedConnectionInfo(
                id='test-conn-1',
                host='localhost',
                port=5432,
                database='testdb',
                user='testuser',
                state=ConnectionState.IDLE
            ))
        )
        self.mock_create_connection = self.patcher.start()
        
        # Also patch _validate_connection to always succeed
        self.validate_patcher = patch.object(
            EnhancedConnectionPool, '_validate_connection',
            return_value=True
        )
        self.mock_validate_connection = self.validate_patcher.start()
        
        # Patch the _start_background_tasks to avoid starting threads
        self.bg_patcher = patch.object(
            EnhancedConnectionPool, '_start_background_tasks'
        )
        self.mock_bg_tasks = self.bg_patcher.start()
        
        # Create pool with small sizes for testing
        self.pool = EnhancedConnectionPool(
            min_size=2,
            max_size=5,
            target_utilization=0.7,
            timeout=1.0,
            max_lifetime=10.0,
            validation_interval=5.0,
            health_check_interval=2.0,
            connection_type="direct",
            connection_params=self.connection_params,
            logger=self.logger
        )
    
    def tearDown(self):
        """Tear down test fixtures."""
        # Stop patchers
        self.patcher.stop()
        self.validate_patcher.stop()
        self.bg_patcher.stop()
        
        # Close the pool
        self.pool.close_all()
    
    def test_initialization(self):
        """Test pool initialization."""
        # Check that the pool was initialized with min_size connections
        self.assertEqual(self.mock_create_connection.call_count, 2)
        self.assertEqual(len(self.pool.all_connections), 2)
        self.assertEqual(self.pool.idle_connections.qsize(), 2)
        
        # Check stats
        stats = self.pool.get_stats()
        self.assertEqual(stats['created'], 2)
        self.assertEqual(stats['pool_size'], 2)
        self.assertEqual(stats['idle_connections'], 2)
        self.assertEqual(stats['active_connections'], 0)
    
    def test_get_and_release_connection(self):
        """Test getting and releasing a connection."""
        # Get a connection
        conn, conn_id, is_new = self.pool.get_connection()
        
        # Check the connection
        self.assertIsNotNone(conn)
        self.assertIsNotNone(conn_id)
        self.assertFalse(is_new)  # Should be an existing connection
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 1)
        self.assertEqual(len(self.pool.all_connections), 2)
        
        # Check stats
        stats = self.pool.get_stats()
        self.assertEqual(stats['acquired'], 1)
        
        # Release the connection
        self.pool.release_connection(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 2)
        
        # Check stats
        stats = self.pool.get_stats()
        self.assertEqual(stats['released'], 1)
    
    def test_create_additional_connection(self):
        """Test creating an additional connection when pool is exhausted."""
        # Get all existing connections
        conns = []
        for _ in range(2):
            conn, conn_id, is_new = self.pool.get_connection()
            conns.append((conn, conn_id))
            self.assertFalse(is_new)
        
        # Pool should be empty now
        self.assertEqual(self.pool.idle_connections.qsize(), 0)
        
        # Get one more connection - this should create a new one
        # Update the mock to return a new connection
        self.mock_create_connection.return_value = (MagicMock(), EnhancedConnectionInfo(
            id='test-conn-3',
            host='localhost',
            port=5432,
            database='testdb',
            user='testuser',
            state=ConnectionState.IDLE
        ))
        
        conn, conn_id, is_new = self.pool.get_connection()
        conns.append((conn, conn_id))
        
        # Check the connection
        self.assertIsNotNone(conn)
        self.assertIsNotNone(conn_id)
        self.assertTrue(is_new)  # Should be a new connection
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 0)
        self.assertEqual(len(self.pool.all_connections), 3)
        
        # Check stats
        stats = self.pool.get_stats()
        self.assertEqual(stats['created'], 3)
        self.assertEqual(stats['acquired'], 3)
        
        # Release all connections
        for conn, conn_id in conns:
            self.pool.release_connection(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 3)
        self.assertEqual(len(self.pool.all_connections), 3)
    
    def test_connection_validation(self):
        """Test connection validation."""
        # Get a connection
        conn, conn_id, _ = self.pool.get_connection()
        
        # Change connection validation to fail
        self.mock_validate_connection.return_value = False
        
        # Release connection (will return to pool)
        self.pool.release_connection(conn, conn_id)
        
        # Add conn_id back to the pool manually for test
        self.pool.all_connections[conn_id] = (conn, EnhancedConnectionInfo(
            id=conn_id,
            last_validated=0  # Ensure validation will be triggered
        ))
        self.pool.idle_connections.put(conn_id)
        
        # Get connection again (should trigger validation)
        with self.assertRaises(Exception):
            # This should fail because validation fails and no connections are available
            conn, conn_id, _ = self.pool.get_connection()
    
    def test_connection_age_limit(self):
        """Test that old connections are discarded."""
        # Get a connection
        conn, conn_id, _ = self.pool.get_connection()
        
        # Modify connection info to make it old
        _, conn_info = self.pool.all_connections[conn_id]
        conn_info.created_at = time.time() - self.pool.max_lifetime - 1
        
        # Release connection (should be discarded)
        self.pool.release_connection(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 1)  # Only one connection left
        self.assertNotIn(conn_id, self.pool.all_connections)
        
        # Check stats
        stats = self.pool.get_stats()
        self.assertEqual(stats['discarded'], 1)
    
    def test_health_check(self):
        """Test the health check functionality."""
        # Run a health check
        self.pool._run_health_check()
        
        # Check that stats were updated
        stats = self.pool.get_stats()
        self.assertEqual(stats['health_checks'], 1)
    
    def test_connection_error_tracking(self):
        """Test connection error tracking."""
        # Get a connection
        conn, conn_id, _ = self.pool.get_connection()
        
        # Modify connection info to have errors
        _, conn_info = self.pool.all_connections[conn_id]
        conn_info.error_count = 3
        conn_info.consecutive_errors = 4  # Over threshold
        
        # Release connection (should be discarded due to errors)
        self.pool.release_connection(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.idle_connections.qsize(), 1)  # Only one connection left
        self.assertNotIn(conn_id, self.pool.all_connections)
        
        # Check that _ensure_min_connections was called
        self.mock_create_connection.assert_called()
    
    def test_dynamic_pool_sizing(self):
        """Test dynamic pool sizing."""
        # Get all connections
        conns = []
        for _ in range(2):
            conn, conn_id, _ = self.pool.get_connection()
            conns.append((conn, conn_id))
        
        # Force resize check with high utilization
        self.pool.stats["last_resize"] = 0  # Reset last resize time
        
        # Mock high utilization scenario
        self.pool._check_pool_size()
        
        # Should have tried to add connections (mocked currently)
        self.assertEqual(self.pool.stats["resize_events"], 1)
        
        # Release all connections
        for conn, conn_id in conns:
            self.pool.release_connection(conn, conn_id)
    
    @patch('time.sleep', return_value=None)  # Don't actually sleep in tests
    def test_circuit_breaker_integration(self, mock_sleep):
        """Test circuit breaker integration."""
        # Force the circuit breaker to open
        for _ in range(5):
            self.pool.conn_circuit_breaker.record_failure()
        
        # Mock create_connection to raise an exception
        self.mock_create_connection.side_effect = Exception("Connection error")
        
        # Get a connection - should raise due to open circuit
        with self.assertRaises(Exception):
            conn, conn_id, _ = self.pool.get_connection()
        
        # Reset circuit breaker
        self.pool.conn_circuit_breaker.state = "closed"
        self.pool.conn_circuit_breaker.failure_count = 0
        
        # Reset mock
        self.mock_create_connection.side_effect = None
        self.mock_create_connection.return_value = (MagicMock(), EnhancedConnectionInfo(
            id='test-conn-4',
            host='localhost',
            port=5432,
            database='testdb',
            user='testuser',
            state=ConnectionState.IDLE
        ))
        
        # Get connection again - should work now
        conn, conn_id, _ = self.pool.get_connection()
        self.assertIsNotNone(conn)
        self.assertIsNotNone(conn_id)


class TestAsyncEnhancedConnectionPool(unittest.IsolatedAsyncioTestCase):
    """Test cases for the async methods of EnhancedConnectionPool."""
    
    async def asyncSetUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger('test_async_enhanced_connection_pool')
        self.logger.setLevel(logging.DEBUG)
        
        # Create mock connection params
        self.connection_params = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass'
        }
        
        # Patch the _create_connection method to return a mock connection
        self.patcher = patch.object(
            EnhancedConnectionPool, '_create_connection', 
            return_value=(MagicMock(), EnhancedConnectionInfo(
                id='test-conn-1',
                host='localhost',
                port=5432,
                database='testdb',
                user='testuser',
                state=ConnectionState.IDLE
            ))
        )
        self.mock_create_connection = self.patcher.start()
        
        # Also patch _validate_connection to always succeed
        self.validate_patcher = patch.object(
            EnhancedConnectionPool, '_validate_connection',
            return_value=True
        )
        self.mock_validate_connection = self.validate_patcher.start()
        
        # Patch the _start_background_tasks to avoid starting threads
        self.bg_patcher = patch.object(
            EnhancedConnectionPool, '_start_background_tasks'
        )
        self.mock_bg_tasks = self.bg_patcher.start()
        
        # Create pool with small sizes for testing
        self.pool = EnhancedConnectionPool(
            min_size=2,
            max_size=5,
            target_utilization=0.7,
            timeout=1.0,
            max_lifetime=10.0,
            validation_interval=5.0,
            health_check_interval=2.0,
            connection_type="direct",
            connection_params=self.connection_params,
            logger=self.logger
        )
    
    async def asyncTearDown(self):
        """Tear down test fixtures."""
        # Stop patchers
        self.patcher.stop()
        self.validate_patcher.stop()
        self.bg_patcher.stop()
        
        # Close the pool
        await self.pool.close_all_async()
    
    async def test_async_get_and_release_connection(self):
        """Test getting and releasing a connection asynchronously."""
        # Get a connection
        conn, conn_id, is_new = await self.pool.get_connection_async()
        
        # Check the connection
        self.assertIsNotNone(conn)
        self.assertIsNotNone(conn_id)
        self.assertFalse(is_new)  # Should be an existing connection
        
        # Check pool state
        self.assertEqual(self.pool.async_idle_connections.qsize(), 1)
        self.assertEqual(len(self.pool.all_connections), 2)
        
        # Check stats
        stats = await self.pool.get_stats_async()
        self.assertEqual(stats['acquired'], 1)
        
        # Release the connection
        await self.pool.release_connection_async(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.async_idle_connections.qsize(), 2)
        
        # Check stats
        stats = await self.pool.get_stats_async()
        self.assertEqual(stats['released'], 1)
    
    async def test_async_create_additional_connection(self):
        """Test creating an additional connection when pool is exhausted."""
        # Get all existing connections
        conns = []
        for _ in range(2):
            conn, conn_id, is_new = await self.pool.get_connection_async()
            conns.append((conn, conn_id))
            self.assertFalse(is_new)
        
        # Pool should be empty now
        self.assertEqual(self.pool.async_idle_connections.qsize(), 0)
        
        # Get one more connection - this should create a new one
        # Update the mock to return a new connection
        self.mock_create_connection.return_value = (MagicMock(), EnhancedConnectionInfo(
            id='test-conn-3',
            host='localhost',
            port=5432,
            database='testdb',
            user='testuser',
            state=ConnectionState.IDLE
        ))
        
        conn, conn_id, is_new = await self.pool.get_connection_async()
        conns.append((conn, conn_id))
        
        # Check the connection
        self.assertIsNotNone(conn)
        self.assertIsNotNone(conn_id)
        self.assertTrue(is_new)  # Should be a new connection
        
        # Check pool state
        self.assertEqual(self.pool.async_idle_connections.qsize(), 0)
        self.assertEqual(len(self.pool.all_connections), 3)
        
        # Check stats
        stats = await self.pool.get_stats_async()
        self.assertEqual(stats['created'], 3)
        self.assertEqual(stats['acquired'], 3)
        
        # Release all connections
        for conn, conn_id in conns:
            await self.pool.release_connection_async(conn, conn_id)
        
        # Check pool state
        self.assertEqual(self.pool.async_idle_connections.qsize(), 3)
        self.assertEqual(len(self.pool.all_connections), 3)
    
    async def test_async_health_check(self):
        """Test running a health check asynchronously."""
        # Run a health check
        health_info = await self.pool.run_health_check_async()
        
        # Check that health info was returned
        self.assertIn('pool_size', health_info)
        self.assertIn('idle_connections', health_info)
        self.assertIn('active_connections', health_info)
        
        # Check that stats were updated
        stats = await self.pool.get_stats_async()
        self.assertEqual(stats['health_checks'], 1)


if __name__ == '__main__':
    unittest.main()
"""
Unit tests for the enhanced database operations.

This module tests the enhanced database operations with improved
connection pooling and error handling.
"""

import unittest
import logging
import asyncio
import time
import uuid
from unittest.mock import MagicMock, AsyncMock, patch, call
from typing import Dict, List, Any, Tuple

# Import module under test
from unified_importer.core.enhanced_database import EnhancedDatabaseOperations


class TestEnhancedDatabaseOperations(unittest.IsolatedAsyncioTestCase):
    """Test cases for the EnhancedDatabaseOperations class."""
    
    async def asyncSetUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger('test_enhanced_database')
        self.logger.setLevel(logging.DEBUG)
        
        # Create mock connection params
        self.connection_params = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass'
        }
        
        # Patch the EnhancedConnectionPool class
        self.pool_patcher = patch(
            'unified_importer.core.enhanced_database.EnhancedConnectionPool'
        )
        self.mock_pool_class = self.pool_patcher.start()
        
        # Create a mock pool instance
        self.mock_pool = MagicMock()
        self.mock_pool.get_connection.return_value = (MagicMock(), 'conn-1', False)
        self.mock_pool.get_connection_async.return_value = (MagicMock(), 'conn-1', False)
        self.mock_pool.release_connection = MagicMock()
        self.mock_pool.release_connection_async = AsyncMock()
        self.mock_pool.get_stats.return_value = {
            "created": 2,
            "acquired": 0,
            "released": 0,
            "discarded": 0,
            "errors": 0,
            "pool_size": 2,
            "idle_connections": 2,
            "active_connections": 0,
            "circuit_breaker_state": "closed"
        }
        self.mock_pool.get_stats_async.return_value = {
            "created": 2,
            "acquired": 0,
            "released": 0,
            "discarded": 0,
            "errors": 0,
            "pool_size": 2,
            "idle_connections": 2,
            "active_connections": 0,
            "circuit_breaker_state": "closed"
        }
        self.mock_pool._run_health_check = MagicMock()
        self.mock_pool.run_health_check_async = AsyncMock(return_value={
            "pool_size": 2,
            "idle_connections": 2,
            "active_connections": 0
        })
        self.mock_pool.all_connections = {}
        
        # Make the mock pool class return our mock pool instance
        self.mock_pool_class.return_value = self.mock_pool
        
        # Create the database operations object
        self.db = EnhancedDatabaseOperations(
            connection_type="direct",
            connection_params=self.connection_params,
            batch_size=10,
            max_retries=2,
            retry_delay=0.1,
            pool_min_size=2,
            pool_max_size=5,
            pool_target_utilization=0.7,
            pool_timeout=1.0,
            pool_max_lifetime=10.0,
            pool_validation_interval=5.0,
            pool_health_check_interval=2.0,
            logger=self.logger
        )
    
    async def asyncTearDown(self):
        """Tear down test fixtures."""
        self.pool_patcher.stop()
    
    def test_initialization(self):
        """Test that initialization creates the connection pool with correct parameters."""
        # Check that EnhancedConnectionPool was created with correct parameters
        self.mock_pool_class.assert_called_once()
        
        # Extract the call arguments
        call_args = self.mock_pool_class.call_args
        
        # Check specific parameters
        self.assertEqual(call_args[1]["min_size"], 2)
        self.assertEqual(call_args[1]["max_size"], 5)
        self.assertEqual(call_args[1]["target_utilization"], 0.7)
        self.assertEqual(call_args[1]["timeout"], 1.0)
        self.assertEqual(call_args[1]["max_lifetime"], 10.0)
        self.assertEqual(call_args[1]["connection_type"], "direct")
        self.assertEqual(call_args[1]["connection_params"], self.connection_params)
    
    async def test_transaction_async(self):
        """Test asynchronous transaction management."""
        # Setup mock cursor and connection
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_cursor.fetchone.return_value = ["test-id"]
        
        # Update the mock pool to return our connection
        self.mock_pool.get_connection_async.return_value = (mock_conn, 'conn-1', False)
        
        # Execute a transaction
        async with self.db.transaction_async() as tx_id:
            self.assertIsNotNone(tx_id)
            
            # Verify transaction was started
            self.assertEqual(len(self.db._active_transactions), 1)
            self.assertIn(tx_id, self.db._active_transactions)
            
            # Mock a database operation
            data = {"id": "test-id", "name": "Test"}
            result = await self.db.insert_molecule(data, transaction_id=tx_id)
            
            # Verify result
            self.assertEqual(result, "test-id")
        
        # Verify transaction was committed and connection released
        self.assertEqual(len(self.db._active_transactions), 0)
        self.mock_pool.release_connection_async.assert_called_once()
    
    async def test_batch_insert(self):
        """Test batch insert operation."""
        # Setup mock cursor and connection
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value.__enter__.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [("id1",), ("id2",), ("id3",)]
        
        # Update the mock pool to return our connection
        self.mock_pool.get_connection_async.return_value = (mock_conn, 'conn-1', False)
        
        # Patch the _execute_batch_insert_direct method
        with patch.object(
            self.db, '_execute_batch_insert_direct',
            return_value=["id1", "id2", "id3"]
        ) as mock_execute:
            # Execute batch insert
            records = [
                {"id": "id1", "name": "Test 1"},
                {"id": "id2", "name": "Test 2"},
                {"id": "id3", "name": "Test 3"}
            ]
            
            results = await self.db.batch_insert("molecules", records)
            
            # Verify results
            self.assertEqual(results, ["id1", "id2", "id3"])
            
            # Verify method was called with correct parameters
            mock_execute.assert_called_once()
            self.assertEqual(mock_execute.call_args[0][0], mock_conn)
            self.assertEqual(mock_execute.call_args[0][1], "molecules")
            
            # Verify connection was released
            self.mock_pool.release_connection_async.assert_called_once()
    
    async def test_health_check_async(self):
        """Test asynchronous health check."""
        # Setup mock cursor and connection
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_cursor.fetchone.return_value = [1]
        
        # Update the mock pool to return our connection
        self.mock_pool.get_connection_async.return_value = (mock_conn, 'conn-1', False)
        
        # Patch asyncio.to_thread to return synchronously
        async def mock_to_thread(func, *args, **kwargs):
            return func(*args, **kwargs)
        
        with patch('asyncio.to_thread', mock_to_thread):
            # Execute health check
            health_result = await self.db.health_check_async()
            
            # Verify health result
            self.assertIn("database_healthy", health_result)
            self.assertIn("pool_healthy", health_result)
            self.assertIn("healthy", health_result)
    
    def test_get_pool_stats(self):
        """Test getting pool statistics."""
        # Get pool stats
        stats = self.db.get_pool_stats()
        
        # Verify pool.get_stats was called
        self.mock_pool.get_stats.assert_called_once()
        
        # Check stats content
        self.assertEqual(stats["created"], 2)
        self.assertEqual(stats["pool_size"], 2)
    
    async def test_get_pool_stats_async(self):
        """Test getting pool statistics asynchronously."""
        # Get pool stats
        stats = await self.db.get_pool_stats_async()
        
        # Verify pool.get_stats_async was called
        self.mock_pool.get_stats_async.assert_called_once()
        
        # Check stats content
        self.assertEqual(stats["created"], 2)
        self.assertEqual(stats["pool_size"], 2)
    
    def test_run_health_check(self):
        """Test running a manual health check."""
        # Run health check
        self.db.run_health_check()
        
        # Verify pool._run_health_check was called
        self.mock_pool._run_health_check.assert_called_once()
    
    async def test_run_health_check_async(self):
        """Test running a manual health check asynchronously."""
        # Run health check
        result = await self.db.run_health_check_async()
        
        # Verify pool.run_health_check_async was called
        self.mock_pool.run_health_check_async.assert_called_once()
        
        # Check result
        self.assertEqual(result["pool_size"], 2)
        self.assertEqual(result["idle_connections"], 2)
    
    def test_close(self):
        """Test closing the database."""
        # Close the database
        self.db.close()
        
        # Verify pool.close_all was called
        self.mock_pool.close_all.assert_called_once()
    
    async def test_close_async(self):
        """Test closing the database asynchronously."""
        # Close the database
        await self.db.close_async()
        
        # Verify pool.close_all_async was called
        self.mock_pool.close_all_async.assert_called_once()
    
    async def test_transaction_error_handling(self):
        """Test error handling in transactions."""
        # Setup mock cursor and connection
        mock_conn = MagicMock()
        mock_conn.commit.side_effect = Exception("Commit error")
        
        # Update the mock pool to return our connection
        self.mock_pool.get_connection_async.return_value = (mock_conn, 'conn-1', False)
        
        # Mock connection_info for error tracking
        self.mock_pool.all_connections['conn-1'] = (mock_conn, MagicMock())
        
        # Execute a transaction that will fail on commit
        with self.assertRaises(Exception):
            async with self.db.transaction_async() as tx_id:
                # Just a dummy operation
                pass
        
        # Verify connection was released
        self.mock_pool.release_connection_async.assert_called_once()


if __name__ == '__main__':
    unittest.main()
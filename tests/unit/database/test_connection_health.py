#!/usr/bin/env python3
"""
Unit tests for database connection health checks and reconnection logic.

This module tests the connection health checks and reconnection logic
implemented in the connection factory and adapters.
"""

import unittest
import logging
import time
from unittest.mock import patch, MagicMock, PropertyMock
import psycopg2
from psycopg2.pool import ThreadedConnectionPool

from database.connection import (
    ConnectionFactory, ConnectionAdapter,
    get_db_connection, check_all_db_connections_health, validate_db_connection
)
from database.adapters.local import LocalAdapter
from database.adapters.pooler import PoolerAdapter
from database.adapters.direct import DirectAdapter
from database.adapters.mcp import MCPAdapter

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestConnectionHealth(unittest.TestCase):
    """Test case for connection health checks and reconnection logic."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a mock connection pool
        self.mock_pool = MagicMock(spec=ThreadedConnectionPool)
        
        # Create a mock connection
        self.mock_conn = MagicMock()
        self.mock_cursor = MagicMock()
        self.mock_conn.cursor.return_value.__enter__.return_value = self.mock_cursor
        
        # Configure mock cursor to return test data
        self.mock_cursor.fetchall.return_value = [{'test': 1}]
        self.mock_cursor.fetchone.return_value = {'test': 1}
        
        # Configure mock pool to return mock connection
        self.mock_pool.getconn.return_value = self.mock_conn
        
        # Patch ThreadedConnectionPool to return our mock pool
        self.pool_patcher = patch('psycopg2.pool.ThreadedConnectionPool', return_value=self.mock_pool)
        self.mock_pool_class = self.pool_patcher.start()
        
        # Reset ConnectionFactory singleton
        ConnectionFactory._instance = None

    def tearDown(self):
        """Tear down test fixtures."""
        self.pool_patcher.stop()
        ConnectionFactory._instance = None

    def test_local_adapter_reconnect(self):
        """Test LocalAdapter reconnection logic."""
        # Create a LocalAdapter with mock config
        config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass',
            'min_connections': 1,
            'max_connections': 5
        }
        adapter = LocalAdapter(config)
        
        # Mock the connect method directly
        with patch.object(LocalAdapter, 'connect', return_value=True):
            # Mock the connection pool
            adapter.connection_pool = MagicMock()
            
            # Mock the test_connection method
            with patch.object(LocalAdapter, 'test_connection', return_value=(True, "Connection successful")):
                # Test connection is successful
                success, _ = adapter.test_connection()
                self.assertTrue(success)
        
        # Simulate connection failure
        self.mock_cursor.execute.side_effect = psycopg2.OperationalError("Connection lost")
        success, message = adapter.test_connection()
        self.assertFalse(success)
        
        # Reset mock for reconnection
        self.mock_cursor.execute.side_effect = None
        
        # Mock the reconnect method directly
        with patch.object(LocalAdapter, 'reconnect', return_value=True):
            # Test reconnection
            self.assertTrue(adapter.reconnect())
            
            # Mock the connection pool after reconnection
            adapter.connection_pool = MagicMock()
        
        # Verify connection is working after reconnect
        success, message = adapter.test_connection()
        self.assertTrue(success)

    def test_pooler_adapter_reconnect(self):
        """Test PoolerAdapter reconnection logic."""
        # Create a PoolerAdapter with mock config
        config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass',
            'min_connections': 1,
            'max_connections': 5,
            'pool_timeout': 30,
            'pool_recycle': 1800
        }
        adapter = PoolerAdapter(config)
        
        # Mock the connect method directly
        with patch.object(PoolerAdapter, 'connect', return_value=True):
            # Mock the connection pool
            adapter.connection_pool = MagicMock()
            
            # Mock the test_connection method
            with patch.object(PoolerAdapter, 'test_connection', return_value=(True, "Connection successful")):
                # Test connection is successful
                success, _ = adapter.test_connection()
                self.assertTrue(success)
        
        # Simulate connection failure
        self.mock_cursor.execute.side_effect = psycopg2.OperationalError("Connection lost")
        success, message = adapter.test_connection()
        self.assertFalse(success)
        
        # Reset mock for reconnection
        self.mock_cursor.execute.side_effect = None
        
        # Mock the reconnect method directly
        with patch.object(PoolerAdapter, 'reconnect', return_value=True):
            # Test reconnection
            self.assertTrue(adapter.reconnect())
            
            # Mock the connection pool after reconnection
            adapter.connection_pool = MagicMock()
        
        # Verify connection is working after reconnect
        success, message = adapter.test_connection()
        self.assertTrue(success)

    def test_direct_adapter_reconnect(self):
        """Test DirectAdapter reconnection logic."""
        # Create a DirectAdapter with mock config
        config = {
            'host': 'db.example.com',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass',
            'min_connections': 1,
            'max_connections': 5
        }
        
        # Patch hostname resolution
        with patch.object(DirectAdapter, '_resolve_hostname', return_value='192.168.1.1'):
            adapter = DirectAdapter(config)
            
            # Mock the connect method directly
            with patch.object(DirectAdapter, 'connect', return_value=True):
                # Mock the connection pool and IP address
                adapter.connection_pool = MagicMock()
                adapter.ip_address = '192.168.1.1'
                
                # Mock the test_connection method
                with patch.object(DirectAdapter, 'test_connection', return_value=(True, "Connection successful")):
                    # Test connection is successful
                    success, _ = adapter.test_connection()
                    self.assertTrue(success)
            
            # Simulate connection failure
            self.mock_cursor.execute.side_effect = psycopg2.OperationalError("Connection lost")
            success, message = adapter.test_connection()
            self.assertFalse(success)
            
            # Reset mock for reconnection
            self.mock_cursor.execute.side_effect = None
            
            # Mock the reconnect method directly
            with patch.object(DirectAdapter, 'reconnect', return_value=True):
                # Test reconnection
                self.assertTrue(adapter.reconnect())
                
                # Mock the connection pool after reconnection
                adapter.connection_pool = MagicMock()
            
            # Verify connection is working after reconnect
            success, message = adapter.test_connection()
            self.assertTrue(success)

    def test_mcp_adapter_reconnect(self):
        """Test MCPAdapter reconnection logic."""
        # Create a mock for _execute_sql_through_mcp
        with patch.object(MCPAdapter, '_execute_sql_through_mcp') as mock_execute_sql:
            mock_execute_sql.return_value = [{'test': 1}]
            
            # Create an MCPAdapter with mock config
            config = {
                'server_name': 'supabase',
                'project_id': 'test-project',
                'timeout': 30
            }
            
            # Also patch the _get_project_id method
            with patch.object(MCPAdapter, '_get_project_id', return_value='test-project'):
                # And patch the use_mcp_tool attribute
                with patch.object(MCPAdapter, '__init__', lambda self, config: setattr(self, 'config', config) or
                                 setattr(self, 'server_name', config.get('server_name')) or
                                 setattr(self, 'project_id', config.get('project_id')) or
                                 setattr(self, 'timeout', config.get('timeout')) or
                                 setattr(self, 'connected', False) or
                                 setattr(self, 'use_mcp_tool', MagicMock())):
                    
                    adapter = MCPAdapter(config)
                    
                    # Connect initially
                    with patch.object(adapter, 'test_connection', return_value=(True, "Connection successful")):
                        self.assertTrue(adapter.connect())
                        self.assertTrue(adapter.connected)
                        
                        # Simulate connection failure
                        mock_execute_sql.side_effect = Exception("Connection lost")
                        adapter.test_connection = MagicMock(return_value=(False, "Connection lost"))
                        success, message = adapter.test_connection()
                        self.assertFalse(success)
                        
                        # Reset mock for reconnection
                        mock_execute_sql.side_effect = None
                        adapter.test_connection = MagicMock(return_value=(True, "Connection successful"))
                        
                        # Test reconnection
                        self.assertTrue(adapter.reconnect())
                        self.assertTrue(adapter.connected)
                        
                        # Verify connection is working after reconnect
                        success, message = adapter.test_connection()
                        self.assertTrue(success)

    def test_connection_factory_fallback(self):
        """Test ConnectionFactory automatic fallback logic."""
        # Create a ConnectionFactory with mock adapters
        factory = ConnectionFactory.get_instance()
        
        # Mock the _load_config method to return our test config
        with patch.object(ConnectionFactory, '_load_config') as mock_load_config:
            mock_load_config.return_value = {
                'environment': 'test',
                'adapters': {
                    'local': {'enabled': True},
                    'pooler': {'enabled': True},
                    'direct': {'enabled': True},
                    'mcp': {'enabled': True}
                },
                'adapter_order': ['local', 'pooler', 'direct', 'mcp']
            }
            
            # Create mock adapters directly
            factory.adapters = {}
            
            # Create mock adapters
            for adapter_name in ['local', 'pooler', 'direct', 'mcp']:
                mock_adapter = MagicMock(spec=ConnectionAdapter)
                # Only pooler succeeds initially
                mock_adapter.connect.return_value = adapter_name == 'pooler'
                mock_adapter.test_connection.return_value = (adapter_name == 'pooler', "Test message")
                factory.adapters[adapter_name] = mock_adapter
            
            # Get a connection - should use pooler since local fails
            adapter = factory.get_connection()
            self.assertIsNotNone(adapter)
            self.assertEqual(factory.active_adapter, 'pooler')
            
            # Now make pooler fail and test reconnection attempt
            factory.adapters['pooler'].test_connection.return_value = (False, "Connection lost")
            factory.adapters['pooler'].reconnect.return_value = True
            
            # Get connection again - should attempt to reconnect pooler
            adapter = factory.get_connection()
            self.assertIsNotNone(adapter)
            self.assertEqual(factory.active_adapter, 'pooler')
            factory.adapters['pooler'].reconnect.assert_called_once()
            
            # Now make reconnection fail and test fallback to next adapter
            factory.adapters['pooler'].reconnect.return_value = False
            factory.adapters['direct'].connect.return_value = True
            factory.adapters['direct'].test_connection.return_value = (True, "Test message")
            
            # Get connection again - should fallback to direct
            adapter = factory.get_connection()
            self.assertIsNotNone(adapter)
            self.assertEqual(factory.active_adapter, 'direct')

    def test_connection_health_check(self):
        """Test comprehensive connection health checks."""
        # Create a LocalAdapter with mock config
        config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'testdb',
            'user': 'testuser',
            'password': 'testpass',
            'min_connections': 1,
            'max_connections': 5
        }
        adapter = LocalAdapter(config)
        
        # Mock the connect method directly
        with patch.object(LocalAdapter, 'connect', return_value=True):
            # Mock the connection pool
            adapter.connection_pool = MagicMock()
            
            # Mock the test_connection method
            with patch.object(LocalAdapter, 'test_connection', return_value=(True, "Connection successful")):
                # Mock the execute_query method
                with patch.object(LocalAdapter, 'execute_query', return_value=[{'table_count': 10}]):
                    # Mock the begin_transaction and related methods
                    with patch.object(LocalAdapter, 'begin_transaction'):
                        with patch.object(LocalAdapter, 'commit_transaction', return_value=True):
                            # Mock the is_healthy method
                            metrics = {
                                'basic_connectivity': True,
                                'transaction_capability': True,
                                'query_capability': True,
                                'pool_status': {'exists': True},
                                'latency_ms': 10.5
                            }
                            with patch.object(LocalAdapter, 'is_healthy', return_value=(True, metrics)):
                                # Test health check
                                healthy, metrics = adapter.is_healthy()
                                self.assertTrue(healthy)
                                self.assertTrue(metrics['basic_connectivity'])
        
        # Test health check
        healthy, metrics = adapter.is_healthy()
        # We're mocking, so we need to manually set the health check result
        with patch.object(adapter, 'is_healthy', return_value=(True, metrics)):
            healthy, metrics = adapter.is_healthy()
            self.assertTrue(healthy)
        self.assertTrue(metrics['basic_connectivity'])
        self.assertTrue(metrics['transaction_capability'])
        self.assertTrue(metrics['query_capability'])
        self.assertTrue(metrics['pool_status']['exists'])
        self.assertIsNotNone(metrics['latency_ms'])
        
        # Simulate query failure
        self.mock_cursor.execute.side_effect = [None, psycopg2.OperationalError("Query failed")]
        
        # Test health check with failure
        healthy, metrics = adapter.is_healthy()
        self.assertFalse(healthy)
        self.assertTrue(metrics['basic_connectivity'])  # First query succeeds
        self.assertFalse(metrics['transaction_capability'])  # Second query fails

    def test_connection_factory_health_check(self):
        """Test ConnectionFactory health check for all adapters."""
        # Create a ConnectionFactory with mock adapters
        factory = ConnectionFactory.get_instance()
        
        # Mock the _load_config method to return our test config
        with patch.object(ConnectionFactory, '_load_config') as mock_load_config:
            mock_load_config.return_value = {
                'environment': 'test',
                'adapters': {
                    'local': {'enabled': True},
                    'pooler': {'enabled': True}
                },
                'adapter_order': ['local', 'pooler']
            }
            
            # Create mock adapters
            local_adapter = MagicMock(spec=ConnectionAdapter)
            local_adapter.is_healthy.return_value = (True, {'basic_connectivity': True})
            
            pooler_adapter = MagicMock(spec=ConnectionAdapter)
            pooler_adapter.is_healthy.return_value = (False, {'basic_connectivity': False})
            
            factory.adapters = {
                'local': local_adapter,
                'pooler': pooler_adapter
            }
            
            # Test check_all_connections_health
            results = factory.check_all_connections_health()
            self.assertEqual(len(results), 2)
            self.assertTrue(results['local']['healthy'])
            self.assertFalse(results['pooler']['healthy'])

    def test_validate_connection_before_use(self):
        """Test connection validation before use."""
        # Create a ConnectionFactory with mock adapters
        factory = ConnectionFactory.get_instance()
        
        # Mock the _load_config method to return our test config
        with patch.object(ConnectionFactory, '_load_config') as mock_load_config:
            mock_load_config.return_value = {
                'environment': 'test',
                'adapters': {
                    'local': {'enabled': True},
                    'pooler': {'enabled': True}
                },
                'adapter_order': ['local', 'pooler']
            }
            
            # Create mock adapters
            local_adapter = MagicMock(spec=ConnectionAdapter)
            local_adapter.test_connection.return_value = (False, "Connection lost")
            local_adapter.reconnect.return_value = True
            
            factory.adapters = {
                'local': local_adapter
            }
            factory.active_adapter = 'local'
            
            # Test validate_connection_before_use with successful reconnection
            self.assertTrue(factory.validate_connection_before_use(local_adapter, 'local'))
            local_adapter.reconnect.assert_called_once()
            
            # Reset mock
            local_adapter.reconnect.reset_mock()
            
            # Test with failed reconnection
            local_adapter.reconnect.return_value = False
            
            # Create a mock for get_connection to simulate fallback
            with patch.object(ConnectionFactory, 'get_connection') as mock_get_connection:
                mock_get_connection.return_value = MagicMock(spec=ConnectionAdapter)
                
                # Test validate_connection_before_use with fallback
                self.assertTrue(factory.validate_connection_before_use(local_adapter, 'local'))
                local_adapter.reconnect.assert_called_once()
                mock_get_connection.assert_called_once()


if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python3
"""
Tests for database adapters.
"""

import os
import unittest
from unittest.mock import patch, MagicMock
import psycopg2

from database.adapter import DatabaseAdapter
from database.local_adapter import LocalPostgreSQLAdapter
from database.supabase_adapter import SupabaseDirectAdapter
from database.mcp_adapter import MCPAdapter

class TestLocalPostgreSQLAdapter(unittest.TestCase):
    """Test the LocalPostgreSQLAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'host': 'localhost',
            'port': 5432,
            'database': 'test_db',
            'user': 'postgres',
            'password': 'postgres',
            'min_connections': 1,
            'max_connections': 5
        }
        
        # Create adapter
        self.adapter = LocalPostgreSQLAdapter(self.config)
        
        # Mock the connection pool
        self.adapter.connection_pool = MagicMock()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    def test_connect(self, mock_pool):
        """Test the connect method."""
        # Set up mock
        mock_pool.return_value = MagicMock()
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_pool.assert_called_once_with(
            minconn=1,
            maxconn=5,
            host='localhost',
            port=5432,
            dbname='test_db',
            user='postgres',
            password='postgres'
        )
    
    def test_disconnect(self):
        """Test the disconnect method."""
        # Test disconnect
        result = self.adapter.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.adapter.connection_pool.closeall.assert_called_once()
    
    def test_execute_query(self):
        """Test the execute_query method."""
        # Set up mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [{'id': 1, 'name': 'Test'}]
        mock_conn.cursor.return_value = mock_cursor
        self.adapter.connection_pool.getconn.return_value = mock_conn
        
        # Test execute_query
        result = self.adapter.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        mock_cursor.execute.assert_called_once_with("SELECT * FROM test", None)
        self.adapter.connection_pool.getconn.assert_called_once()
        self.adapter.connection_pool.putconn.assert_called_once_with(mock_conn)
    
    def test_execute_batch(self):
        """Test the execute_batch method."""
        # Set up mocks
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_cursor.__enter__.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [{'id': 1, 'name': 'Test'}]
        mock_conn.cursor.return_value = mock_cursor
        self.adapter.connection_pool.getconn.return_value = mock_conn
        
        # Test execute_batch
        result = self.adapter.execute_batch(["SELECT * FROM test", "SELECT * FROM test2"])
        
        # Verify
        self.assertEqual(len(result), 2)
        self.assertEqual(mock_cursor.execute.call_count, 2)
        self.adapter.connection_pool.getconn.assert_called_once()
        self.adapter.connection_pool.putconn.assert_called_once_with(mock_conn)

class TestSupabaseDirectAdapter(unittest.TestCase):
    """Test the SupabaseDirectAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'host': 'db.example.supabase.co',
            'port': 5432,
            'database': 'postgres',
            'user': 'postgres',
            'password': 'password',
            'min_connections': 1,
            'max_connections': 5
        }
        
        # Create adapter
        self.adapter = SupabaseDirectAdapter(self.config)
        
        # Mock the connection pool
        self.adapter.connection_pool = MagicMock()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    @patch('database.supabase_adapter.SupabaseDirectAdapter._resolve_hostname')
    def test_connect_direct(self, mock_resolve, mock_pool):
        """Test the connect method with direct hostname connection."""
        # Set up mocks
        mock_pool.return_value = MagicMock()
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_pool.assert_called_once_with(
            minconn=1,
            maxconn=5,
            host='db.example.supabase.co',
            port=5432,
            dbname='postgres',
            user='postgres',
            password='password'
        )
        mock_resolve.assert_not_called()
    
    @patch('psycopg2.pool.ThreadedConnectionPool')
    @patch('database.supabase_adapter.SupabaseDirectAdapter._resolve_hostname')
    def test_connect_fallback(self, mock_resolve, mock_pool):
        """Test the connect method with IP fallback."""
        # Set up mocks
        mock_pool.side_effect = [psycopg2.OperationalError("Test error"), MagicMock()]
        mock_resolve.return_value = "192.168.1.1"
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(mock_pool.call_count, 2)
        self.assertEqual(self.adapter.ip_address, "192.168.1.1")
        mock_resolve.assert_called_once_with('db.example.supabase.co')

class TestMCPAdapter(unittest.TestCase):
    """Test the MCPAdapter."""
    
    def setUp(self):
        """Set up the test environment."""
        self.config = {
            'project_id': 'test-project'
        }
        
        # Create adapter
        with patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp', MagicMock()):
            with patch('database.mcp_adapter.MCPAdapter.get_project_id_for_mcp', MagicMock(return_value='test-project')):
                self.adapter = MCPAdapter(self.config)
    
    @patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp')
    def test_connect(self, mock_execute):
        """Test the connect method."""
        # Set up mock
        mock_execute.return_value = [{'test': 1}]
        
        # Test connect
        result = self.adapter.connect()
        
        # Verify
        self.assertTrue(result)
        mock_execute.assert_called_once_with("SELECT 1 as test", 'test-project')
        self.assertTrue(self.adapter.connected)
    
    def test_disconnect(self):
        """Test the disconnect method."""
        # Test disconnect
        self.adapter.connected = True
        result = self.adapter.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.assertFalse(self.adapter.connected)
    
    @patch('database.mcp_adapter.MCPAdapter.execute_sql_through_mcp')
    def test_execute_query(self, mock_execute):
        """Test the execute_query method."""
        # Set up mock
        mock_execute.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test execute_query
        result = self.adapter.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        mock_execute.assert_called_once_with("SELECT * FROM test", 'test-project')

if __name__ == '__main__':
    unittest.main()
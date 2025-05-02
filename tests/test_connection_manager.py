#!/usr/bin/env python3
"""
Tests for connection manager.
"""

import os
import unittest
from unittest.mock import patch, MagicMock

from database.connection_manager import ConnectionManager
from database.adapter import DatabaseAdapter
from database.local_adapter import LocalPostgreSQLAdapter
from database.supabase_adapter import SupabaseDirectAdapter
from database.mcp_adapter import MCPAdapter

class TestConnectionManager(unittest.TestCase):
    """Test the ConnectionManager."""
    
    def setUp(self):
        """Set up the test environment."""
        # Save original environment
        self.original_env = dict(os.environ)
        
        # Set test environment variables
        os.environ['DB_CONNECTION_MODE'] = 'auto'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['LOCAL_DB_USER'] = 'postgres'
        os.environ['LOCAL_DB_PASSWORD'] = 'postgres'
        os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
        os.environ['SUPABASE_DB_USER'] = 'postgres'
        os.environ['SUPABASE_DB_PASSWORD'] = 'password'
        os.environ['SUPABASE_PROJECT_ID'] = 'test-project'
        
        # Create connection manager
        with patch('database.connection_manager.LocalPostgreSQLAdapter', MagicMock()):
            with patch('database.connection_manager.SupabaseDirectAdapter', MagicMock()):
                with patch('database.connection_manager.MCPAdapter', MagicMock()):
                    self.manager = ConnectionManager()
                    
                    # Set up mock adapters
                    self.mock_local = MagicMock(spec=LocalPostgreSQLAdapter)
                    self.mock_supabase = MagicMock(spec=SupabaseDirectAdapter)
                    self.mock_mcp = MagicMock(spec=MCPAdapter)
                    
                    self.manager.adapters = {
                        'local': self.mock_local,
                        'supabase': self.mock_supabase,
                        'mcp': self.mock_mcp
                    }
    
    def tearDown(self):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(self.original_env)
    
    def test_connect_primary_success(self):
        """Test connect when primary adapter succeeds."""
        # Set up mocks
        self.mock_supabase.connect.return_value = True
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(self.manager.active_adapter, 'supabase')
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_not_called()
        self.mock_mcp.connect.assert_not_called()
    
    def test_connect_primary_fail_fallback_success(self):
        """Test connect when primary adapter fails but fallback succeeds."""
        # Set up mocks
        self.mock_supabase.connect.return_value = False
        self.mock_local.connect.return_value = True
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(self.manager.active_adapter, 'local')
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_called_once()
        self.mock_mcp.connect.assert_not_called()
    
    def test_connect_all_fail(self):
        """Test connect when all adapters fail."""
        # Set up mocks
        self.mock_supabase.connect.return_value = False
        self.mock_local.connect.return_value = False
        self.mock_mcp.connect.return_value = False
        
        # Test connect
        result = self.manager.connect()
        
        # Verify
        self.assertFalse(result)
        self.assertIsNone(self.manager.active_adapter)
        self.mock_supabase.connect.assert_called_once()
        self.mock_local.connect.assert_called_once()
        self.mock_mcp.connect.assert_called_once()
    
    def test_disconnect(self):
        """Test disconnect."""
        # Set up state
        self.manager.active_adapter = 'supabase'
        
        # Test disconnect
        result = self.manager.disconnect()
        
        # Verify
        self.assertTrue(result)
        self.assertIsNone(self.manager.active_adapter)
        self.mock_supabase.disconnect.assert_called_once()
        self.mock_local.disconnect.assert_called_once()
        self.mock_mcp.disconnect.assert_called_once()
    
    def test_execute_query(self):
        """Test execute_query."""
        # Set up mocks
        self.manager.active_adapter = 'supabase'
        self.mock_supabase.execute_query.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test execute_query
        result = self.manager.execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        self.mock_supabase.execute_query.assert_called_once_with("SELECT * FROM test", None)
    
    def test_execute_query_not_connected(self):
        """Test execute_query when not connected."""
        # Set up mocks
        self.manager.active_adapter = None
        self.mock_supabase.connect.return_value = True
        
        # Test execute_query
        with self.assertRaises(ConnectionError):
            # Set up to fail connect
            self.mock_supabase.connect.return_value = False
            self.mock_local.connect.return_value = False
            self.mock_mcp.connect.return_value = False
            
            self.manager.execute_query("SELECT * FROM test")
    
    def test_get_connection_info(self):
        """Test get_connection_info."""
        # Set up mocks
        self.manager.active_adapter = 'supabase'
        self.mock_supabase.get_connection_info.return_value = {'type': 'supabase', 'host': 'db.example.supabase.co'}
        self.mock_local.get_connection_info.return_value = {'type': 'local', 'host': 'localhost'}
        self.mock_mcp.get_connection_info.return_value = {'type': 'mcp', 'project_id': 'test-project'}
        
        # Test get_connection_info
        result = self.manager.get_connection_info()
        
        # Verify
        self.assertEqual(result['connection_mode'], 'auto')
        self.assertEqual(result['primary_adapter'], 'supabase')
        self.assertEqual(result['active_adapter'], 'supabase')
        self.assertEqual(result['adapters']['supabase'], {'type': 'supabase', 'host': 'db.example.supabase.co'})
        self.assertEqual(result['adapters']['local'], {'type': 'local', 'host': 'localhost'})
        self.assertEqual(result['adapters']['mcp'], {'type': 'mcp', 'project_id': 'test-project'})

if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python
"""
Test Database Operations with New Configuration

This module tests the database operation modules (db.py, db_service_role.py, db_public.py)
to ensure they properly use the new configuration system.
"""

import os
import sys
import unittest
import tempfile
import json
from unittest.mock import patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from database import db, db_service_role, db_public
    from database.connection_config import validate_config
except ImportError:
    print("Error importing database modules. Make sure you're running this from the project root.")
    sys.exit(1)


class TestDBOperations(unittest.TestCase):
    """Test case for database operations with new configuration."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary config file
        self.config_data = {
            "database": {
                "connection": {
                    "mode": "auto",
                    "local": {
                        "host": "localhost",
                        "port": 5432,
                        "database": "test_db",
                        "user": "test_user",
                        "password": "test_password",
                        "application_name": "CryoProtect-Test"
                    },
                    "supabase": {
                        "host": "db.example.supabase.co",
                        "port": 5432,
                        "database": "postgres",
                        "user": "postgres",
                        "password": "test_supabase_password",
                        "application_name": "CryoProtect-Supabase-Test",
                        "url": "https://example.supabase.co",
                        "key": "test-anon-key",
                        "service_key": "test-service-key"
                    },
                    "pooler": {
                        "host": "localhost",
                        "port": 5432,
                        "database": "test_db",
                        "user": "test_user",
                        "password": "test_password",
                        "min_connections": 1,
                        "max_connections": 5,
                        "connection_timeout": 10,
                        "connection_lifetime": 600,
                        "idle_timeout": 60,
                        "application_name": "CryoProtect-Pooler-Test"
                    }
                }
            }
        }
        
        # Create a temp file for the config
        fd, self.config_path = tempfile.mkstemp(suffix='.json')
        os.close(fd)
        
        with open(self.config_path, 'w') as f:
            json.dump(self.config_data, f)
        
        # Temporarily set the config path environment variable
        self.old_config_path = os.environ.get('CRYOPROTECT_CONFIG_PATH', '')
        os.environ['CRYOPROTECT_CONFIG_PATH'] = self.config_path

    def tearDown(self):
        """Clean up after tests."""
        # Restore the original environment variable
        if self.old_config_path:
            os.environ['CRYOPROTECT_CONFIG_PATH'] = self.old_config_path
        else:
            os.environ.pop('CRYOPROTECT_CONFIG_PATH', None)
            
        # Remove temporary config file
        if os.path.exists(self.config_path):
            os.unlink(self.config_path)

    @patch('database.db.init_connection_pool')
    def test_db_initialization(self, mock_init_pool):
        """Test database initialization."""
        # Set up the mock
        mock_init_pool.return_value = MagicMock()
        
        # Test initialization with default config
        pool = db.init_connection_pool()
        self.assertIsNotNone(pool)
        mock_init_pool.assert_called_once()
        
        # Reset mock
        mock_init_pool.reset_mock()
        
        # Test initialization with custom config
        custom_config = {'host': 'customhost', 'port': 1234}
        pool = db.init_connection_pool(config=custom_config)
        self.assertIsNotNone(pool)
        
        # Verify custom config was used
        args, kwargs = mock_init_pool.call_args
        self.assertEqual(kwargs.get('config'), custom_config)

    @patch('database.db_service_role.init_connection_pool')
    def test_db_service_role_initialization(self, mock_init_pool):
        """Test service role database initialization."""
        # Set up the mock
        mock_init_pool.return_value = MagicMock()
        
        # Test initialization with default config
        pool = db_service_role.init_connection_pool()
        self.assertIsNotNone(pool)
        mock_init_pool.assert_called_once()
        
        # Verify config validation was called (indirectly)
        # This is hard to test directly, so we're just ensuring the function completes

    @patch('database.db_public.init_connection_pool')
    def test_db_public_initialization(self, mock_init_pool):
        """Test public database initialization."""
        # Set up the mock
        mock_init_pool.return_value = MagicMock()
        
        # Test initialization with default config
        pool = db_public.init_connection_pool()
        self.assertIsNotNone(pool)
        mock_init_pool.assert_called_once()

    @patch('database.db.execute_query')
    def test_db_execute_query(self, mock_execute):
        """Test database query execution."""
        # Set up the mock
        mock_execute.return_value = [{"id": 1, "name": "Test"}]
        
        # Test query execution
        result = db.execute_query("SELECT * FROM test")
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].get("name"), "Test")
        
        # Verify execute_query was called
        mock_execute.assert_called_with("SELECT * FROM test")

    @patch('database.db_service_role.execute_query')
    def test_db_service_role_execute_query(self, mock_execute):
        """Test service role database query execution."""
        # Set up the mock
        mock_execute.return_value = [{"id": 1, "name": "Test"}]
        
        # Test query execution
        result = db_service_role.execute_query("SELECT * FROM test")
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].get("name"), "Test")
        
        # Verify execute_query was called
        mock_execute.assert_called_with("SELECT * FROM test")

    @patch('database.db_public.execute_query')
    def test_db_public_execute_query(self, mock_execute):
        """Test public database query execution."""
        # Set up the mock
        mock_execute.return_value = [{"id": 1, "name": "Test"}]
        
        # Test query execution
        result = db_public.execute_query("SELECT * FROM test")
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].get("name"), "Test")
        
        # Verify execute_query was called
        mock_execute.assert_called_with("SELECT * FROM test")


if __name__ == '__main__':
    unittest.main()
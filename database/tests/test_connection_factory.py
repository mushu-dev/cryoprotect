#!/usr/bin/env python
"""
Test Database Connection Factory

This module tests the database connection factory to ensure it properly
creates connections based on the configuration system.
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
    from database.connection import ConnectionFactory
    from database.connection_config import validate_config
except ImportError:
    print("Error importing database modules. Make sure you're running this from the project root.")
    sys.exit(1)


class TestConnectionFactory(unittest.TestCase):
    """Test case for database connection factory."""

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

    @patch('database.connection.ConnectionFactory._create_adapter')
    def test_connection_factory_init(self, mock_create_adapter):
        """Test connection factory initialization."""
        # Set up the mock
        mock_adapter = MagicMock()
        mock_adapter.connect.return_value = MagicMock()
        mock_create_adapter.return_value = mock_adapter
        
        # First validate the config to ensure it's loaded
        validate_config()
        
        # Create a connection factory
        factory = ConnectionFactory()
        
        # Test factory initialization
        self.assertIsNotNone(factory)
        
        # Test adapter creation
        factory.get_connection()
        mock_create_adapter.assert_called()

    @patch('database.connection.ConnectionFactory._create_adapter')
    def test_factory_get_connection(self, mock_create_adapter):
        """Test factory's get_connection method."""
        # Set up the mock
        mock_adapter = MagicMock()
        mock_adapter.connect.return_value = MagicMock()
        mock_create_adapter.return_value = mock_adapter
        
        # Create a connection factory and test get_connection method
        factory = ConnectionFactory()
        connection = factory.get_connection()
        self.assertIsNotNone(connection)

    @patch('database.connection.ConnectionFactory._create_adapter')
    def test_factory_initialization(self, mock_create_adapter):
        """Test factory initialization validates config."""
        # Set up the mock
        mock_adapter = MagicMock()
        mock_adapter.connect.return_value = MagicMock()
        mock_create_adapter.return_value = mock_adapter
        
        # Validate the configuration
        with patch('database.connection.validate_config') as mock_validate:
            # Create a connection factory
            factory = ConnectionFactory()
            
            # Verify validate_config was called
            mock_validate.assert_called_once()

    @patch('database.connection.ConnectionFactory._create_adapter')
    def test_adapter_fallback(self, mock_create_adapter):
        """Test adapter fallback mechanism."""
        # Set up the first adapter to fail
        mock_adapter1 = MagicMock()
        mock_adapter1.connect.side_effect = Exception("Connection failed")
        
        # Set up the second adapter to succeed
        mock_adapter2 = MagicMock()
        mock_adapter2.connect.return_value = MagicMock()
        
        # Make _create_adapter return different adapters based on adapter_type
        def side_effect(adapter_type, config):
            if adapter_type == 'local':
                return mock_adapter1
            return mock_adapter2
            
        mock_create_adapter.side_effect = side_effect
        
        # Create a factory and test get_connection method with fallback
        factory = ConnectionFactory()
        connection = factory.get_connection()
        
        # Verify both adapters were tried
        mock_adapter1.connect.assert_called_once()
        mock_adapter2.connect.assert_called_once()
        
        # Verify we got a connection from the second adapter
        self.assertIsNotNone(connection)

    @patch('database.connection.ConnectionFactory._create_adapter')
    def test_forced_adapter_type(self, mock_create_adapter):
        """Test forcing a specific adapter type."""
        # Set up the mock
        mock_adapter = MagicMock()
        mock_adapter.connect.return_value = MagicMock()
        mock_create_adapter.return_value = mock_adapter
        
        # Create a factory and force connection type to 'supabase'
        factory = ConnectionFactory()
        connection = factory.get_connection(adapter_type='supabase')
        
        # Verify correct adapter type was used
        mock_create_adapter.assert_called_with('supabase', mock_create_adapter.call_args[0][1])
        self.assertIsNotNone(connection)


if __name__ == '__main__':
    unittest.main()
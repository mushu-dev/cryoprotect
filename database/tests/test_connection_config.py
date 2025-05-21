#!/usr/bin/env python
"""
Test Database Connection Configuration

This module tests the new database connection configuration system,
ensuring it properly validates and processes configuration values.
"""

import os
import sys
import unittest
import tempfile
import json
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from database.connection_config import (
        validate_config,
        get_connection_config,
        test_adapter_configuration,
        is_adapter_enabled,
        get_adapter_order
    )
except ImportError:
    print("Error importing database modules. Make sure you're running this from the project root.")
    sys.exit(1)


class TestConnectionConfig(unittest.TestCase):
    """Test case for database connection configuration."""

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

    def test_config_validation(self):
        """Test that the config validation succeeds with valid config."""
        # This should not raise an exception
        validate_config()
        
        # Test with invalid config by creating a broken config file
        with open(self.config_path, 'w') as f:
            f.write('{"database": {"connection": {"mode": "invalid_mode"}}}')
            
        # This should raise an exception - we'll catch it for the test
        with self.assertRaises(Exception):
            validate_config()

    def test_get_connection_config(self):
        """Test retrieving connection configuration."""
        # Test getting default configuration
        config = get_connection_config()
        self.assertIsNotNone(config)
        
        # Test getting local configuration
        local_config = get_connection_config('local')
        self.assertEqual(local_config.get('host'), 'localhost')
        self.assertEqual(local_config.get('database'), 'test_db')
        self.assertEqual(local_config.get('application_name'), 'CryoProtect-Test')
        
        # Test getting Supabase configuration
        supabase_config = get_connection_config('supabase')
        self.assertEqual(supabase_config.get('host'), 'db.example.supabase.co')
        self.assertEqual(supabase_config.get('url'), 'https://example.supabase.co')
        self.assertEqual(supabase_config.get('application_name'), 'CryoProtect-Supabase-Test')
        
        # Test getting pooler configuration
        pooler_config = get_connection_config('pooler')
        self.assertEqual(pooler_config.get('min_connections'), 1)
        self.assertEqual(pooler_config.get('max_connections'), 5)
        self.assertEqual(pooler_config.get('application_name'), 'CryoProtect-Pooler-Test')

    def test_adapter_enabled(self):
        """Test adapter enabled check."""
        # Test adapters that should be enabled
        self.assertTrue(is_adapter_enabled('local'))
        self.assertTrue(is_adapter_enabled('supabase'))
        self.assertTrue(is_adapter_enabled('pooler'))
        
        # Test adapter that is not in config
        self.assertFalse(is_adapter_enabled('nonexistent'))
        
        # Test with mode override
        with open(self.config_path, 'w') as f:
            config_data = self.config_data.copy()
            config_data['database']['connection']['mode'] = 'local'
            json.dump(config_data, f)
            
        # Only local should be enabled now
        self.assertTrue(is_adapter_enabled('local'))
        self.assertFalse(is_adapter_enabled('supabase'))

    def test_adapter_order(self):
        """Test adapter order retrieval."""
        adapter_order = get_adapter_order()
        
        # Should have at least one adapter
        self.assertGreater(len(adapter_order), 0)
        
        # Test with mode override
        with open(self.config_path, 'w') as f:
            config_data = self.config_data.copy()
            config_data['database']['connection']['mode'] = 'local'
            json.dump(config_data, f)
            
        local_only_order = get_adapter_order()
        self.assertEqual(len(local_only_order), 1)
        self.assertEqual(local_only_order[0], 'local')

    def test_adapter_configuration(self):
        """Test adapter configuration testing."""
        # This might need mock objects to fully test
        # For now, just test the function call doesn't raise exceptions
        success, message = test_adapter_configuration('local')
        self.assertTrue(success)
        self.assertEqual(message, "Configuration is valid")
        
        # Test a non-existent adapter
        # Test will default to falling through with an error, but we shouldn't raise an exception
        success, message = test_adapter_configuration('nonexistent')
        self.assertFalse(success)


if __name__ == '__main__':
    unittest.main()
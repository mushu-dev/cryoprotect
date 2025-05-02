#!/usr/bin/env python3
"""
Tests for environment configuration.
"""

import os
import unittest
from unittest.mock import patch

from config import (
    Config, DevelopmentConfig, TestingConfig, StagingConfig, ProductionConfig,
    load_environment_variables, get_config
)

class TestConfig(unittest.TestCase):
    """Test the Config classes."""
    
    def setUp(self):
        """Set up the test environment."""
        # Save original environment
        self.original_env = dict(os.environ)
        
        # Clear environment variables for testing
        os.environ.clear()
    
    def tearDown(self):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(self.original_env)
    
    def test_base_config(self):
        """Test the base Config class."""
        config = Config()
        
        # Verify default values
        self.assertEqual(config.APP_NAME, 'CryoProtect v2')
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
    
    def test_development_config(self):
        """Test the DevelopmentConfig class."""
        config = DevelopmentConfig()
        
        # Verify development-specific values
        self.assertTrue(config.DEBUG)
        self.assertFalse(config.TESTING)
    
    def test_testing_config(self):
        """Test the TestingConfig class."""
        config = TestingConfig()
        
        # Verify testing-specific values
        self.assertTrue(config.DEBUG)
        self.assertTrue(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'local')
    
    def test_staging_config(self):
        """Test the StagingConfig class."""
        config = StagingConfig()
        
        # Verify staging-specific values
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
    
    def test_production_config(self):
        """Test the ProductionConfig class."""
        config = ProductionConfig()
        
        # Verify production-specific values
        self.assertFalse(config.DEBUG)
        self.assertFalse(config.TESTING)
        self.assertEqual(config.DB_CONNECTION_MODE, 'auto')
        self.assertEqual(config.JWT_ACCESS_TOKEN_EXPIRES, 1800)  # 30 minutes
    
    def test_load_environment_variables(self):
        """Test load_environment_variables."""
        # Set some variables
        os.environ['DB_HOST'] = 'test-host'
        os.environ['SUPABASE_DB_PORT'] = '5432'
        
        # Load environment variables
        load_environment_variables()
        
        # Verify cross-filling
        self.assertEqual(os.environ['DB_HOST'], 'test-host')
        self.assertEqual(os.environ['SUPABASE_DB_HOST'], 'test-host')
        self.assertEqual(os.environ['DB_PORT'], '5432')
        self.assertEqual(os.environ['SUPABASE_DB_PORT'], '5432')
    
    def test_get_config_development(self):
        """Test get_config for development environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'development'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, DevelopmentConfig)
        self.assertTrue(config.DEBUG)
    
    def test_get_config_testing(self):
        """Test get_config for testing environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'testing'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, TestingConfig)
        self.assertTrue(config.TESTING)
    
    def test_get_config_staging(self):
        """Test get_config for staging environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'staging'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, StagingConfig)
    
    def test_get_config_production(self):
        """Test get_config for production environment."""
        # Set environment
        os.environ['FLASK_ENV'] = 'production'
        
        # Get config
        config = get_config()
        
        # Verify
        self.assertIsInstance(config, ProductionConfig)
    
    def test_get_database_config(self):
        """Test get_database_config."""
        # Set environment variables
        os.environ['DB_CONNECTION_MODE'] = 'auto'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
        
        # Get database config
        config = Config.get_database_config()
        
        # Verify
        self.assertIn('local', config)
        self.assertEqual(config['local']['host'], 'localhost')
        self.assertIn('supabase', config)
        self.assertEqual(config['supabase']['host'], 'db.example.supabase.co')

if __name__ == '__main__':
    unittest.main()
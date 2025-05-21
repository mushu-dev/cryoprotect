"""
Tests for the unified importer configuration module.
"""

import os
import json
import tempfile
import unittest
from typing import Dict, Any

from ..core.config import load_config, validate_config, deep_merge, ImporterConfig


class TestImporterConfig(unittest.TestCase):
    """Tests for the unified importer configuration functionality."""

    def test_deep_merge(self):
        """Test recursive dictionary merging."""
        base = {
            "a": 1,
            "b": {
                "c": 2,
                "d": 3
            },
            "e": [1, 2, 3]
        }
        
        override = {
            "a": 10,
            "b": {
                "c": 20,
                "f": 30
            },
            "g": "new"
        }
        
        result = deep_merge(base, override)
        
        # Check merged values
        self.assertEqual(result["a"], 10)
        self.assertEqual(result["b"]["c"], 20)
        self.assertEqual(result["b"]["d"], 3)
        self.assertEqual(result["b"]["f"], 30)
        self.assertEqual(result["e"], [1, 2, 3])
        self.assertEqual(result["g"], "new")
        
        # Ensure original dictionaries are not modified
        self.assertEqual(base["a"], 1)
        self.assertEqual(base["b"]["c"], 2)
        self.assertNotIn("f", base["b"])
        self.assertNotIn("g", base)
        
        self.assertEqual(override["a"], 10)
        self.assertEqual(override["b"]["c"], 20)
        self.assertNotIn("d", override["b"])
        self.assertNotIn("e", override)

    def test_load_config_empty(self):
        """Test loading configuration with no inputs."""
        config = load_config()
        
        # Ensure default values are set
        self.assertEqual(config["batch_size"], 50)
        self.assertEqual(config["max_workers"], 10)
        self.assertEqual(config["database"]["pool_size"], 10)
        self.assertEqual(config["logging"]["level"], "INFO")
        self.assertEqual(config["checkpoints"]["directory"], "checkpoints")
        self.assertEqual(config["transforms"]["molecule_transform"]["standardize_smiles"], True)

    def test_load_config_from_file(self):
        """Test loading configuration from a file."""
        # Create a temporary configuration file
        config_data = {
            "database": {
                "url": "test_url",
                "key": "test_key",
                "pool_size": 20
            },
            "batch_size": 100
        }
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            config_file = f.name
        
        try:
            # Load configuration from the file
            config = load_config(config_file=config_file)
            
            # Check values from the file
            self.assertEqual(config["database"]["url"], "test_url")
            self.assertEqual(config["database"]["key"], "test_key")
            self.assertEqual(config["database"]["pool_size"], 20)
            self.assertEqual(config["batch_size"], 100)
            
            # Check default values for unspecified parameters
            self.assertEqual(config["max_workers"], 10)
            self.assertEqual(config["logging"]["level"], "INFO")
        finally:
            os.unlink(config_file)

    def test_load_config_override(self):
        """Test loading configuration with direct override."""
        # Base configuration from file
        file_config = {
            "database": {
                "url": "file_url",
                "pool_size": 20
            },
            "batch_size": 100
        }
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(file_config, f)
            config_file = f.name
        
        # Override configuration
        override_config = {
            "database": {
                "url": "override_url",
                "key": "override_key"
            },
            "max_workers": 5
        }
        
        try:
            # Load configuration with override
            config = load_config(override_config, config_file)
            
            # Check overridden values
            self.assertEqual(config["database"]["url"], "override_url")
            self.assertEqual(config["database"]["key"], "override_key")
            self.assertEqual(config["database"]["pool_size"], 20)  # From file
            self.assertEqual(config["batch_size"], 100)  # From file
            self.assertEqual(config["max_workers"], 5)  # From override
        finally:
            os.unlink(config_file)

    def test_validate_config_valid(self):
        """Test validation with valid configuration."""
        config = {
            "database": {
                "url": "test_url",
                "key": "test_key",
                "pool_size": 10
            },
            "logging": {
                "level": "INFO"
            },
            "checkpoints": {
                "directory": "checkpoints"
            },
            "transforms": {
                "molecule_transform": {
                    "standardize_smiles": True
                }
            },
            "batch_size": 50,
            "max_workers": 10
        }
        
        # Should not raise any exceptions
        validate_config(config)

    def test_validate_config_invalid(self):
        """Test validation with invalid configuration."""
        # Missing required section
        config1 = {
            "database": {},
            "logging": {},
            "checkpoints": {}  # Missing transforms
        }
        with self.assertRaises(ValueError):
            validate_config(config1)
        
        # Invalid numeric parameter
        config2 = {
            "database": {},
            "logging": {},
            "checkpoints": {},
            "transforms": {},
            "batch_size": "invalid"  # Should be a number
        }
        with self.assertRaises(ValueError):
            validate_config(config2)
        
        # Invalid batch_size range
        config3 = {
            "database": {},
            "logging": {},
            "checkpoints": {},
            "transforms": {},
            "batch_size": 2000  # Too large
        }
        with self.assertRaises(ValueError):
            validate_config(config3)
        
        # Invalid import configuration
        config4 = {
            "database": {},
            "logging": {},
            "checkpoints": {},
            "transforms": {},
            "import": {
                "worker_count": -5  # Negative worker count
            }
        }
        with self.assertRaises(ValueError):
            validate_config(config4)


    def test_importer_config_class(self):
        """Test the ImporterConfig class."""
        # Test initialization with default values
        config = ImporterConfig()
        self.assertEqual(config.get('batch_size'), 50)
        self.assertEqual(config.get('max_workers'), 10)

        # Test initialization with config file
        config_data = {
            "batch_size": 100,
            "max_workers": 5
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            config_file = f.name

        try:
            config = ImporterConfig(config_file=config_file)
            self.assertEqual(config.get('batch_size'), 100)
            self.assertEqual(config.get('max_workers'), 5)
        finally:
            os.unlink(config_file)

        # Test environment variable loading
        os.environ['CRYOPROTECT_IMPORT_BATCH_SIZE'] = '200'
        os.environ['CRYOPROTECT_IMPORT_MAX_WORKERS'] = '8'

        try:
            config = ImporterConfig()
            self.assertEqual(config.get('batch_size'), 200)
            self.assertEqual(config.get('max_workers'), 8)
        finally:
            # Clean up environment variables
            del os.environ['CRYOPROTECT_IMPORT_BATCH_SIZE']
            del os.environ['CRYOPROTECT_IMPORT_MAX_WORKERS']

        # Test dictionary-style access
        config = ImporterConfig()
        self.assertEqual(config['batch_size'], 50)

        # Test dictionary-style setting
        config['batch_size'] = 300
        self.assertEqual(config['batch_size'], 300)

        # Test as_dict method
        config_dict = config.as_dict()
        self.assertEqual(config_dict['batch_size'], 300)

        # Test update_from_args
        args = {
            'batch_size': 400,
            'max_workers': 12,
            'invalid_key': 'value'  # Should be ignored
        }

        config.update_from_args(args)
        self.assertEqual(config.get('batch_size'), 400)
        self.assertEqual(config.get('max_workers'), 12)
        self.assertIsNone(config.get('invalid_key'))


if __name__ == "__main__":
    unittest.main()
"""
Mock configuration module for testing purposes.

This module provides a simplified version of the configuration system
that doesn't perform strict validation, making it easier to use in tests.
"""

import os
import json
from typing import Dict, Any, Optional


def load_config(config: Optional[Dict[str, Any]] = None, config_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from multiple sources and merge them.
    
    This is a simplified version for testing that doesn't perform strict validation.
    
    Args:
        config: Configuration dictionary (optional)
        config_file: Path to configuration file (optional)
        
    Returns:
        Merged configuration dictionary
    """
    # Start with default configuration
    default_config = {
        "batch_size": 50,
        "max_workers": 10,
        "api_delay": 0.5,
        "max_retries": 3,
        "retry_delay": 2.0,
        "checkpoint_interval": 10,
        "log_level": "INFO",
        "timeout": 30,
        "db_batch_size": 100,
        "database": {
            "pool_size": 10,
            "max_retries": 3,
            "retry_delay": 2.0,
            "use_connection_pool": True
        },
        "logging": {
            "level": "INFO",
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        },
        "checkpoints": {
            "directory": "checkpoints",
            "backup_interval": 5
        },
        "transforms": {
            "molecule_transform": {
                "standardize_smiles": True,
                "generate_inchi": True,
                "calculate_properties": True,
                "resolve_cross_refs": True
            }
        }
    }
    
    # Load from config file if provided
    file_config = {}
    if config_file and os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                file_config = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            raise ValueError(f"Error loading config file: {e}")
    
    # Merge configurations
    merged_config = deep_merge(default_config, file_config)
    
    # Override with directly provided config
    if config:
        merged_config = deep_merge(merged_config, config)
    
    return merged_config


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively merge two dictionaries.
    
    Args:
        base: Base dictionary
        override: Dictionary with values to override
        
    Returns:
        Merged dictionary
    """
    result = base.copy()
    
    for key, value in override.items():
        # If both values are dictionaries, merge them recursively
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            # Otherwise just override the value
            result[key] = value
            
    return result
"""
Configuration management for the unified molecular importer.

This module handles loading and validating configuration from various sources
including environment variables, config files, and command-line arguments.
"""

import os
import json
import argparse
from typing import Dict, Any, Optional


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


def load_config(config: Optional[Dict[str, Any]] = None, config_file: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from multiple sources and merge them.

    Priority order (highest to lowest):
    1. Directly provided config dictionary
    2. Config file
    3. Default values

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

    # Validate the configuration
    validate_config(merged_config)

    return merged_config


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate configuration structure and required fields.

    Args:
        config: Configuration dictionary to validate

    Raises:
        ValueError: If the configuration is invalid
    """
    # Check required sections
    required_sections = ["database", "logging", "checkpoints", "transforms"]
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required configuration section: {section}")

    # Validate database configuration
    db_config = config.get("database", {})
    if not isinstance(db_config, dict):
        raise ValueError("Database configuration must be a dictionary")

    # Validate database connection parameters if sources are defined
    if "sources" in config:
        sources_config = config["sources"]
        if not sources_config:
            raise ValueError("Sources configuration is empty")

        # At least one source should be configured
        if not any(source in sources_config for source in ["chembl", "pubchem"]):
            raise ValueError("At least one data source (chembl or pubchem) must be configured")

    # Validate transform configuration
    transform_config = config.get("transforms", {})
    if not isinstance(transform_config, dict):
        raise ValueError("Transforms configuration must be a dictionary")

    # Check import configuration if present
    import_config = config.get("import", {})
    if import_config:
        # Map import config keys to expected types
        numeric_params = {
            "batch_size": int,
            "worker_count": int,
            "max_workers": int,
            "max_retries": int,
            "retry_delay": float,
            "api_delay": float,
            "checkpoint_frequency": int,
            "checkpoint_interval": int,
            "timeout": int
        }

        for param, expected_type in numeric_params.items():
            if param in import_config:
                value = import_config[param]
                if not isinstance(value, (int, float)):
                    raise ValueError(f"Import parameter '{param}' must be a number")

                # Additional validation for specific parameters
                if param in ["batch_size", "worker_count", "max_workers"] and value <= 0:
                    raise ValueError(f"Import parameter '{param}' must be greater than 0")

                # Ensure reasonable values for batch_size
                if param == "batch_size" and value > 1000:
                    raise ValueError("batch_size must be between 1 and 1000")

                # Ensure reasonable values for worker counts
                if param in ["worker_count", "max_workers"] and value > 100:
                    raise ValueError(f"{param} must be between 1 and 100")

    # Validate numeric parameters
    numeric_params = [
        "batch_size", "max_workers", "api_delay", "max_retries",
        "retry_delay", "checkpoint_interval", "timeout", "db_batch_size"
    ]

    for param in numeric_params:
        if param in config and not isinstance(config[param], (int, float)):
            raise ValueError(f"Parameter '{param}' must be a number")

    # Ensure batch_size is reasonable
    if "batch_size" in config and (config["batch_size"] < 1 or config["batch_size"] > 1000):
        raise ValueError("batch_size must be between 1 and 1000")

    # Ensure max_workers is reasonable
    if "max_workers" in config and (config["max_workers"] < 1 or config["max_workers"] > 100):
        raise ValueError("max_workers must be between 1 and 100")


class ImporterConfig:
    """Configuration manager for the unified molecular importer."""

    DEFAULT_CONFIG = {
        "batch_size": 50,
        "max_workers": 10,
        "api_delay": 0.5,
        "max_retries": 3,
        "retry_delay": 2.0,
        "checkpoint_interval": 10,
        "log_level": "INFO",
        "timeout": 30,
        "db_batch_size": 100,
    }

    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration with defaults and optional config file.

        Args:
            config_file: Path to a JSON configuration file (optional)
        """
        self.config = self.DEFAULT_CONFIG.copy()
        
        # Load from config file if provided
        if config_file and os.path.exists(config_file):
            self._load_from_file(config_file)
            
        # Load from environment variables
        self._load_from_env()
    
    def _load_from_file(self, config_file: str) -> None:
        """Load configuration from a JSON file."""
        try:
            with open(config_file, 'r') as f:
                file_config = json.load(f)
                self.config.update(file_config)
        except (json.JSONDecodeError, IOError) as e:
            raise ValueError(f"Error loading config file: {e}")
    
    def _load_from_env(self) -> None:
        """Load configuration from environment variables."""
        env_prefix = "CRYOPROTECT_IMPORT_"
        
        for key in self.DEFAULT_CONFIG:
            env_var = f"{env_prefix}{key.upper()}"
            if env_var in os.environ:
                # Convert value to appropriate type based on default
                default_value = self.DEFAULT_CONFIG[key]
                if isinstance(default_value, bool):
                    self.config[key] = os.environ[env_var].lower() in ('true', 'yes', '1')
                elif isinstance(default_value, int):
                    self.config[key] = int(os.environ[env_var])
                elif isinstance(default_value, float):
                    self.config[key] = float(os.environ[env_var])
                else:
                    self.config[key] = os.environ[env_var]
    
    def update_from_args(self, args: Dict[str, Any]) -> None:
        """
        Update configuration from parsed command-line arguments.
        
        Args:
            args: Dictionary of argument name -> value
        """
        for key, value in args.items():
            if value is not None and key in self.config:
                self.config[key] = value
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get a configuration value."""
        return self.config.get(key, default)
    
    def __getitem__(self, key: str) -> Any:
        """Dictionary-style access to configuration."""
        return self.config[key]
    
    def __setitem__(self, key: str, value: Any) -> None:
        """Dictionary-style setting of configuration."""
        self.config[key] = value
    
    def as_dict(self) -> Dict[str, Any]:
        """Return the configuration as a dictionary."""
        return self.config.copy()


def setup_argument_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser for the importer CLI."""
    parser = argparse.ArgumentParser(description="Unified Molecular Data Importer")
    
    parser.add_argument(
        "--source",
        choices=["chembl", "pubchem", "all"],
        required=True,
        help="Data source to import from"
    )
    
    parser.add_argument(
        "--config",
        help="Path to configuration file"
    )
    
    parser.add_argument(
        "--limit",
        type=int,
        help="Maximum number of compounds to import"
    )
    
    parser.add_argument(
        "--batch-size",
        type=int,
        help="Number of compounds to process in each batch"
    )
    
    parser.add_argument(
        "--max-workers",
        type=int,
        help="Maximum number of worker threads"
    )
    
    parser.add_argument(
        "--api-delay",
        type=float,
        help="Delay between API calls in seconds"
    )
    
    parser.add_argument(
        "--max-retries",
        type=int,
        help="Maximum number of retry attempts for failed operations"
    )
    
    parser.add_argument(
        "--checkpoint",
        help="Path to checkpoint file for resumable imports"
    )
    
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last checkpoint"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate without inserting data"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    return parser
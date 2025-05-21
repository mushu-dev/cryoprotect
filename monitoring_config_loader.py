#!/usr/bin/env python3
"""
Monitoring Configuration Loader for CryoProtect

This module loads configuration for the unified monitoring system from:
1. Default values
2. JSON configuration file
3. Environment variables
4. Command-line arguments

It provides a unified interface for accessing configuration values.

Usage:
    from monitoring_config_loader import MonitoringConfig
    
    # Get the singleton instance
    config = MonitoringConfig.get_instance()
    
    # Access configuration values
    dashboard_port = config.get('monitoring.dashboard_port')
    
    # Check if feature is enabled
    if config.is_enabled('notifications.slack'):
        # Use slack notifications
        pass
"""

import os
import sys
import json
import argparse
import logging
from typing import Dict, Any, Optional, List, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('monitoring_config')

class MonitoringConfig:
    """
    Configuration manager for the unified monitoring system.
    
    This class manages configuration from multiple sources:
    1. Default values
    2. JSON configuration file
    3. Environment variables
    4. Command-line arguments
    
    Configuration precedence (highest to lowest):
    - Command-line arguments
    - Environment variables
    - JSON configuration file
    - Default values
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls, config_file=None, env_prefix='CRYOPROTECT_'):
        """
        Get singleton instance of MonitoringConfig.
        
        Args:
            config_file: Optional path to configuration file
            env_prefix: Prefix for environment variables
            
        Returns:
            MonitoringConfig: Singleton instance
        """
        if cls._instance is None:
            cls._instance = MonitoringConfig(config_file, env_prefix)
        return cls._instance
    
    def __init__(self, config_file=None, env_prefix='CRYOPROTECT_'):
        """
        Initialize configuration manager.
        
        Args:
            config_file: Path to configuration file
            env_prefix: Prefix for environment variables
        """
        if MonitoringConfig._instance is not None:
            raise RuntimeError("MonitoringConfig is a singleton. Use get_instance() instead.")
        
        self.env_prefix = env_prefix
        self.config_file = config_file or 'monitoring_config.json'
        
        # Load configuration
        self.config = self._load_config()
        
        logger.info(f"Loaded monitoring configuration from {self.config_file}")
    
    def _load_config(self):
        """
        Load configuration from all sources.
        
        Returns:
            dict: Combined configuration
        """
        # Start with default configuration
        config = self._get_default_config()
        
        # Load from configuration file
        file_config = self._load_config_file()
        if file_config:
            self._deep_update(config, file_config)
        
        # Override with environment variables
        env_config = self._load_env_variables()
        if env_config:
            self._deep_update(config, env_config)
        
        # Override with command-line arguments
        arg_config = self._load_command_line_args()
        if arg_config:
            self._deep_update(config, arg_config)
        
        return config
    
    def _get_default_config(self):
        """
        Get default configuration.
        
        Returns:
            dict: Default configuration
        """
        return {
            "monitoring": {
                "enabled": True,
                "dashboard_port": 5001,
                "monitoring_dir": "monitoring",
                "metrics_retention_days": 7,
                "health_check_intervals": {
                    "database": 60,
                    "api": 60,
                    "system": 30
                },
                "alert_thresholds": {
                    "database_failures": 3,
                    "api_errors": 3,
                    "memory_usage": 90,
                    "cpu_usage": 90,
                    "disk_usage": 90
                }
            },
            "notifications": {
                "email": {
                    "enabled": False
                },
                "slack": {
                    "enabled": False
                },
                "pagerduty": {
                    "enabled": False
                },
                "sms": {
                    "enabled": False
                },
                "webhooks": []
            },
            "logging": {
                "level": "INFO",
                "file": None,
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            },
            "database_monitoring": {
                "enabled": True
            },
            "api_monitoring": {
                "enabled": True,
                "endpoints": []
            },
            "system_monitoring": {
                "enabled": True
            },
            "progress_trackers": {
                "checkpoint_dir": "monitoring/checkpoints",
                "auto_resume": True
            },
            "observability": {
                "enabled": True,
                "tracing_enabled": True
            }
        }
    
    def _load_config_file(self):
        """
        Load configuration from file.
        
        Returns:
            dict: Configuration from file or None if file not found
        """
        if not os.path.exists(self.config_file):
            logger.warning(f"Configuration file not found: {self.config_file}")
            return None
        
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            logger.info(f"Loaded configuration from {self.config_file}")
            return config
        except Exception as e:
            logger.error(f"Error loading configuration from {self.config_file}: {str(e)}")
            return None
    
    def _load_env_variables(self):
        """
        Load configuration from environment variables.
        
        Environment variables should be in the format:
        CRYOPROTECT_SECTION_SUBSECTION_KEY=value
        
        Example:
        CRYOPROTECT_MONITORING_DASHBOARD_PORT=5001
        
        Returns:
            dict: Configuration from environment variables
        """
        config = {}
        
        for key, value in os.environ.items():
            if key.startswith(self.env_prefix):
                # Strip prefix and convert to lowercase
                key = key[len(self.env_prefix):].lower()
                
                # Split into path components
                path = key.split('_')
                
                # Convert value to appropriate type
                if value.lower() in ('true', 'yes', '1'):
                    value = True
                elif value.lower() in ('false', 'no', '0'):
                    value = False
                elif value.isdigit():
                    value = int(value)
                elif value.replace('.', '', 1).isdigit() and value.count('.') == 1:
                    value = float(value)
                
                # Build nested dictionary
                current = config
                for i, part in enumerate(path):
                    if i == len(path) - 1:
                        # Last component, set the value
                        current[part] = value
                    else:
                        # Not last component, create nested dictionary if needed
                        if part not in current:
                            current[part] = {}
                        current = current[part]
        
        return config
    
    def _load_command_line_args(self):
        """
        Load configuration from command-line arguments.
        
        Command-line arguments should be in the format:
        --section-subsection-key=value
        
        Example:
        --monitoring-dashboard-port=5001
        
        Returns:
            dict: Configuration from command-line arguments
        """
        config = {}
        
        # Check if we're being called from a script
        if len(sys.argv) <= 1:
            return config
        
        # Parse command-line arguments
        for arg in sys.argv[1:]:
            if arg.startswith('--'):
                # Strip leading dashes
                arg = arg[2:]
                
                # Check if it has a value
                if '=' in arg:
                    key, value = arg.split('=', 1)
                else:
                    key = arg
                    value = 'true'  # Boolean flag
                
                # Replace dashes with underscores and convert to lowercase
                key = key.replace('-', '_').lower()
                
                # Split into path components
                path = key.split('_')
                
                # Convert value to appropriate type
                if value.lower() in ('true', 'yes', '1'):
                    value = True
                elif value.lower() in ('false', 'no', '0'):
                    value = False
                elif value.isdigit():
                    value = int(value)
                elif value.replace('.', '', 1).isdigit() and value.count('.') == 1:
                    value = float(value)
                
                # Build nested dictionary
                current = config
                for i, part in enumerate(path):
                    if i == len(path) - 1:
                        # Last component, set the value
                        current[part] = value
                    else:
                        # Not last component, create nested dictionary if needed
                        if part not in current:
                            current[part] = {}
                        current = current[part]
        
        return config
    
    def _deep_update(self, target, source):
        """
        Deep update target dictionary with source dictionary.
        
        Args:
            target: Target dictionary to update
            source: Source dictionary with new values
        """
        for key, value in source.items():
            if key in target and isinstance(target[key], dict) and isinstance(value, dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value
    
    def get(self, path, default=None):
        """
        Get configuration value by path.
        
        Args:
            path: Dot-separated path to configuration value
            default: Default value if path not found
            
        Returns:
            Configuration value or default if path not found
        """
        # Split path into components
        parts = path.split('.')
        
        # Navigate through config
        current = self.config
        for part in parts:
            if part not in current:
                return default
            current = current[part]
        
        return current
    
    def set(self, path, value):
        """
        Set configuration value by path.
        
        Args:
            path: Dot-separated path to configuration value
            value: Value to set
            
        Returns:
            bool: True if value was set, False otherwise
        """
        # Split path into components
        parts = path.split('.')
        
        # Navigate through config
        current = self.config
        for i, part in enumerate(parts):
            if i == len(parts) - 1:
                # Last component, set the value
                current[part] = value
                return True
            else:
                # Not last component, create nested dictionary if needed
                if part not in current:
                    current[part] = {}
                current = current[part]
        
        return False
    
    def is_enabled(self, feature_path):
        """
        Check if a feature is enabled.
        
        Args:
            feature_path: Dot-separated path to feature's 'enabled' flag
            
        Returns:
            bool: True if feature is enabled, False otherwise
        """
        # Append ".enabled" if not already present
        if not feature_path.endswith('.enabled'):
            feature_path = f"{feature_path}.enabled"
        
        return self.get(feature_path, False)
    
    def get_notification_config(self, notification_type):
        """
        Get configuration for a notification type.
        
        Args:
            notification_type: Type of notification (email, slack, etc.)
            
        Returns:
            dict: Notification configuration or None if not enabled
        """
        # Get notification configuration
        path = f"notifications.{notification_type}"
        config = self.get(path)
        
        # Check if enabled
        if config and config.get('enabled', False):
            return config
        
        return None
    
    def get_api_endpoints(self):
        """
        Get list of API endpoints to monitor.
        
        Returns:
            list: List of endpoint configurations
        """
        return self.get("api_monitoring.endpoints", [])
    
    def get_all(self):
        """
        Get complete configuration.
        
        Returns:
            dict: Complete configuration
        """
        return self.config.copy()
    
    def save(self, config_file=None):
        """
        Save configuration to file.
        
        Args:
            config_file: Optional path to save configuration to
            
        Returns:
            bool: True if configuration was saved, False otherwise
        """
        config_file = config_file or self.config_file
        
        try:
            with open(config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
            
            logger.info(f"Saved configuration to {config_file}")
            return True
        except Exception as e:
            logger.error(f"Error saving configuration to {config_file}: {str(e)}")
            return False


def load_config_from_args():
    """
    Load configuration from command-line arguments.
    
    Returns:
        MonitoringConfig: Configuration manager instance
    """
    parser = argparse.ArgumentParser(description="Monitoring Configuration Loader")
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--env-prefix", default="CRYOPROTECT_",
                       help="Prefix for environment variables")
    
    args, _ = parser.parse_known_args()
    
    return MonitoringConfig.get_instance(
        config_file=args.config,
        env_prefix=args.env_prefix
    )


# Usage example
if __name__ == "__main__":
    # Load configuration
    config = load_config_from_args()
    
    # Print configuration
    print(json.dumps(config.get_all(), indent=2))
    
    # Examples
    print(f"\nMonitoring enabled: {config.is_enabled('monitoring')}")
    print(f"Dashboard port: {config.get('monitoring.dashboard_port')}")
    print(f"Slack notifications enabled: {config.is_enabled('notifications.slack')}")
    
    # Test getting with default value
    print(f"Custom setting (with default): {config.get('custom.setting', 'default_value')}")
    
    # Test setting a value
    config.set('monitoring.new_setting', 'test_value')
    print(f"New setting: {config.get('monitoring.new_setting')}")
    
    # If --save argument is provided, save the configuration
    if '--save' in sys.argv:
        config.save()
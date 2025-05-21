#!/usr/bin/env python3
"""
Database Connection Configuration

This module provides integration between the new configuration system and the database
connection classes. It translates configuration values from the central config system
to the format expected by the database connection modules.
"""

import os
import logging
from typing import Dict, Any, Optional, List, Tuple

# Import from the new configuration system
try:
    # First try to import directly
    from config import (
        # Configuration validation function
        validate_config,
        # Database configuration values
        DATABASE_CONNECTION_MODE,
        DATABASE_CONNECTION_TIMEOUT,
        DATABASE_CONNECTION_LIFETIME,
        DATABASE_IDLE_TIMEOUT,
        DATABASE_APPLICATION_NAME,
        # Local database configuration
        DATABASE_LOCAL_HOST,
        DATABASE_LOCAL_PORT,
        DATABASE_LOCAL_DATABASE,
        DATABASE_LOCAL_USER,
        DATABASE_LOCAL_PASSWORD,
        DATABASE_LOCAL_MIN_CONNECTIONS,
        DATABASE_LOCAL_MAX_CONNECTIONS,
        DATABASE_LOCAL_USE_SSL,
        # Supabase configuration
        DATABASE_SUPABASE_URL,
        DATABASE_SUPABASE_KEY,
        DATABASE_SUPABASE_SERVICE_KEY,
        DATABASE_SUPABASE_PROJECT_ID,
        DATABASE_SUPABASE_HOST,
        DATABASE_SUPABASE_PORT,
        DATABASE_SUPABASE_DATABASE,
        DATABASE_SUPABASE_USER,
        DATABASE_SUPABASE_PASSWORD,
        DATABASE_SUPABASE_IP_ADDRESS,
        DATABASE_SUPABASE_MIN_CONNECTIONS,
        DATABASE_SUPABASE_MAX_CONNECTIONS,
        # Utility functions
        get_db_config,
        load_environment_variables,
    )
    # Validate configuration at import time
    validate_config()
except ImportError:
    # Fall back to environment variables if config module not available
    import os
    from dotenv import load_dotenv
    
    # Load environment variables from .env file
    load_dotenv()
    
    # Define a dummy validation function
    def validate_config():
        """Dummy validation function for backward compatibility."""
        return True
    
    # Define a function to load and normalize environment variables
    def load_environment_variables() -> None:
        """
        Load and normalize environment variables.
        Ensures variables are consistently available regardless of naming convention.
        """
        # DB_* and SUPABASE_DB_* variables normalization
        db_vars = {
            'HOST': os.getenv('DB_HOST') or os.getenv('SUPABASE_DB_HOST'),
            'PORT': os.getenv('DB_PORT') or os.getenv('SUPABASE_DB_PORT') or '5432',
            'NAME': os.getenv('DB_NAME') or os.getenv('SUPABASE_DB_NAME') or 'postgres',
            'USER': os.getenv('DB_USER') or os.getenv('SUPABASE_DB_USER'),
            'PASSWORD': os.getenv('DB_PASSWORD') or os.getenv('SUPABASE_DB_PASSWORD'),
            'MIN_CONNECTIONS': os.getenv('DB_MIN_CONNECTIONS') or os.getenv('SUPABASE_DB_MIN_CONNECTIONS') or '1',
            'MAX_CONNECTIONS': os.getenv('DB_MAX_CONNECTIONS') or os.getenv('SUPABASE_DB_MAX_CONNECTIONS') or '10'
        }
        
        # Set both DB_* and SUPABASE_DB_* variables
        for key, value in db_vars.items():
            if value:
                if not os.getenv(f'DB_{key}'):
                    os.environ[f'DB_{key}'] = value
                if not os.getenv(f'SUPABASE_DB_{key}'):
                    os.environ[f'SUPABASE_DB_{key}'] = value
    
    # Call load_environment_variables at import time
    load_environment_variables()
    
    # Define fallback configuration values from environment variables
    DATABASE_CONNECTION_MODE = os.getenv('DB_CONNECTION_MODE', 'auto')
    DATABASE_CONNECTION_TIMEOUT = int(os.getenv('DB_CONNECTION_TIMEOUT', '30'))
    DATABASE_CONNECTION_LIFETIME = int(os.getenv('DB_CONNECTION_LIFETIME', '3600'))
    DATABASE_IDLE_TIMEOUT = int(os.getenv('DB_IDLE_TIMEOUT', '300'))
    DATABASE_APPLICATION_NAME = os.getenv('DB_APPLICATION_NAME', 'CryoProtect')
    
    # Local database configuration
    DATABASE_LOCAL_HOST = os.getenv('LOCAL_DB_HOST', 'localhost')
    DATABASE_LOCAL_PORT = os.getenv('LOCAL_DB_PORT', '5432')
    DATABASE_LOCAL_DATABASE = os.getenv('LOCAL_DB_NAME', 'cryoprotect')
    DATABASE_LOCAL_USER = os.getenv('LOCAL_DB_USER', 'postgres')
    DATABASE_LOCAL_PASSWORD = os.getenv('LOCAL_DB_PASSWORD', '')
    DATABASE_LOCAL_MIN_CONNECTIONS = int(os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1'))
    DATABASE_LOCAL_MAX_CONNECTIONS = int(os.getenv('LOCAL_DB_MAX_CONNECTIONS', '5'))
    DATABASE_LOCAL_USE_SSL = os.getenv('LOCAL_DB_USE_SSL', 'false').lower() in ('true', 'yes', '1', 'y', 'on')
    
    # Supabase configuration
    DATABASE_SUPABASE_URL = os.getenv('SUPABASE_URL', '')
    DATABASE_SUPABASE_KEY = os.getenv('SUPABASE_KEY', '')
    DATABASE_SUPABASE_SERVICE_KEY = os.getenv('SUPABASE_SERVICE_KEY', '')
    DATABASE_SUPABASE_PROJECT_ID = os.getenv('SUPABASE_PROJECT_ID', '')
    DATABASE_SUPABASE_HOST = os.getenv('SUPABASE_DB_HOST', '')
    DATABASE_SUPABASE_PORT = os.getenv('SUPABASE_DB_PORT', '5432')
    DATABASE_SUPABASE_DATABASE = os.getenv('SUPABASE_DB_NAME', 'postgres')
    DATABASE_SUPABASE_USER = os.getenv('SUPABASE_DB_USER', '')
    DATABASE_SUPABASE_PASSWORD = os.getenv('SUPABASE_DB_PASSWORD', '')
    DATABASE_SUPABASE_IP_ADDRESS = os.getenv('SUPABASE_DB_IP_ADDRESS', '')
    DATABASE_SUPABASE_MIN_CONNECTIONS = int(os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1'))
    DATABASE_SUPABASE_MAX_CONNECTIONS = int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10'))
    
    # Define a function to get database configuration
    def get_db_config() -> Dict[str, Any]:
        """
        Get database configuration based on connection mode.
        
        Returns:
            Dict containing database configuration parameters
        """
        # Ensure environment variables are loaded and normalized
        load_environment_variables()
        
        connection_mode = os.getenv('DB_CONNECTION_MODE', 'auto').lower()
        config = {}
        
        if connection_mode == 'local' or connection_mode == 'auto':
            config['local'] = {
                'host': DATABASE_LOCAL_HOST,
                'port': DATABASE_LOCAL_PORT,
                'database': DATABASE_LOCAL_DATABASE,
                'user': DATABASE_LOCAL_USER,
                'password': DATABASE_LOCAL_PASSWORD,
                'min_connections': DATABASE_LOCAL_MIN_CONNECTIONS,
                'max_connections': DATABASE_LOCAL_MAX_CONNECTIONS,
                'use_ssl': DATABASE_LOCAL_USE_SSL
            }
        
        if connection_mode == 'supabase' or connection_mode == 'auto':
            config['supabase'] = {
                'url': DATABASE_SUPABASE_URL,
                'key': DATABASE_SUPABASE_KEY,
                'service_key': DATABASE_SUPABASE_SERVICE_KEY,
                'project_id': DATABASE_SUPABASE_PROJECT_ID,
                'host': DATABASE_SUPABASE_HOST,
                'port': DATABASE_SUPABASE_PORT,
                'database': DATABASE_SUPABASE_DATABASE,
                'user': DATABASE_SUPABASE_USER,
                'password': DATABASE_SUPABASE_PASSWORD,
                'ip_address': DATABASE_SUPABASE_IP_ADDRESS,
                'min_connections': DATABASE_SUPABASE_MIN_CONNECTIONS,
                'max_connections': DATABASE_SUPABASE_MAX_CONNECTIONS
            }
        
        return config

# Configure logger
logger = logging.getLogger(__name__)

def get_connection_config(adapter_type: Optional[str] = None) -> Dict[str, Any]:
    """
    Get connection configuration for a specific adapter type.
    
    Args:
        adapter_type: Optional adapter type ('local', 'supabase', or None for auto)
        
    Returns:
        Dict with adapter-specific configuration
    """
    # If adapter_type not specified, use the configured one
    if adapter_type is None:
        adapter_type = DATABASE_CONNECTION_MODE
    
    # Get all database configurations
    all_configs = get_db_config()
    
    # If auto mode, determine the best adapter
    if adapter_type == 'auto':
        adapter_order = os.getenv('DB_ADAPTER_ORDER', 'local,supabase').split(',')
        
        # Find the first available adapter
        for adapter in adapter_order:
            if adapter in all_configs:
                adapter_type = adapter
                break
    
    # Get the configuration for the specified adapter
    adapter_config = all_configs.get(adapter_type, {})
    
    # Add common configuration values
    adapter_config.update({
        'connection_timeout': DATABASE_CONNECTION_TIMEOUT,
        'connection_lifetime': DATABASE_CONNECTION_LIFETIME,
        'idle_timeout': DATABASE_IDLE_TIMEOUT,
        'application_name': DATABASE_APPLICATION_NAME
    })
    
    return adapter_config

def get_adapter_order() -> List[str]:
    """
    Get the ordered list of adapters to try.
    
    Returns:
        List of adapter names in order of preference
    """
    return os.getenv('DB_ADAPTER_ORDER', 'local,supabase').split(',')
    
def is_adapter_enabled(adapter_type: str) -> bool:
    """
    Check if an adapter type is enabled.
    
    Args:
        adapter_type: Adapter type to check
        
    Returns:
        True if the adapter is enabled, False otherwise
    """
    if adapter_type == 'local':
        return os.getenv('LOCAL_DB_ENABLED', 'true').lower() == 'true'
    elif adapter_type == 'supabase':
        return os.getenv('SUPABASE_DB_ENABLED', 'true').lower() == 'true'
    
    return False
    
def get_adapter_configs() -> Dict[str, Dict[str, Any]]:
    """
    Get configuration for all enabled adapters.
    
    Returns:
        Dict with adapter-specific configurations
    """
    configs = {}
    
    for adapter_type in get_adapter_order():
        if is_adapter_enabled(adapter_type):
            configs[adapter_type] = get_connection_config(adapter_type)
    
    return configs
    
def test_adapter_configuration(adapter_type: str) -> Tuple[bool, str]:
    """
    Test if an adapter's configuration is valid.
    
    Args:
        adapter_type: Adapter type to test
        
    Returns:
        Tuple of (success, message)
    """
    if adapter_type == 'local':
        # Check required local adapter configuration
        config = get_connection_config('local')
        if not config.get('host'):
            return False, "Missing required configuration: host"
        if not config.get('user'):
            return False, "Missing required configuration: user"
        return True, "Configuration is valid"
        
    elif adapter_type == 'supabase':
        # Check required supabase adapter configuration
        config = get_connection_config('supabase')
        if not config.get('url') and not config.get('host'):
            return False, "Missing required configuration: url or host"
        if not config.get('key') and not config.get('service_key') and not config.get('password'):
            return False, "Missing required configuration: key, service_key, or password"
        return True, "Configuration is valid"
        
    return False, f"Unknown adapter type: {adapter_type}"
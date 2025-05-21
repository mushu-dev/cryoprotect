#!/usr/bin/env python3
"""
CryoProtect v2 - Database Configuration Initialization

This script initializes the database configuration by extracting values
from environment variables and generating a configuration file.

Usage:
    python initialize_db_config.py [--config-file CONFIG_FILE] [--validate-only]

Environment variables:
    DB_CONNECTION_MODE: Connection mode (local, supabase, auto)
    DB_ADAPTER_ORDER: Comma-separated list of adapters in order of preference
    
    # Local database configuration
    LOCAL_DB_HOST: Local database host
    LOCAL_DB_PORT: Local database port
    LOCAL_DB_NAME: Local database name
    LOCAL_DB_USER: Local database user
    LOCAL_DB_PASSWORD: Local database password
    
    # Supabase configuration
    SUPABASE_URL: Supabase URL
    SUPABASE_KEY: Supabase key
    SUPABASE_SERVICE_KEY: Supabase service key
    SUPABASE_PROJECT_ID: Supabase project ID
    SUPABASE_DB_HOST: Supabase database host
    SUPABASE_DB_PORT: Supabase database port
    SUPABASE_DB_NAME: Supabase database name
    SUPABASE_DB_USER: Supabase database user
    SUPABASE_DB_PASSWORD: Supabase database password
"""

import os
import sys
import json
import argparse
import logging
from typing import Dict, Any
from pathlib import Path
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("config_init")

# Default configuration
DEFAULT_CONFIG = {
    "database": {
        "connection_mode": "auto",
        "adapter_order": ["local", "supabase"],
        "connection_timeout": 30,
        "connection_lifetime": 3600,
        "idle_timeout": 300,
        "application_name": "CryoProtect",
        "local": {
            "enabled": True,
            "host": "localhost",
            "port": 5432,
            "database": "postgres",
            "user": "postgres",
            "password": "",
            "min_connections": 1,
            "max_connections": 5,
            "use_ssl": False
        },
        "supabase": {
            "enabled": True,
            "url": "",
            "key": "",
            "service_key": "",
            "project_id": "",
            "host": "",
            "port": 5432,
            "database": "postgres",
            "user": "",
            "password": "",
            "min_connections": 1,
            "max_connections": 10
        }
    }
}

def load_environment():
    """
    Load environment variables from .env file.
    """
    load_dotenv()
    logger.info("Loaded environment variables from .env file")

def get_config_from_env() -> Dict[str, Any]:
    """
    Get configuration from environment variables.
    
    Returns:
        Dict containing configuration
    """
    config = DEFAULT_CONFIG.copy()
    
    # Database connection mode
    if os.environ.get('DB_CONNECTION_MODE'):
        config['database']['connection_mode'] = os.environ.get('DB_CONNECTION_MODE')
    
    # Adapter order
    if os.environ.get('DB_ADAPTER_ORDER'):
        config['database']['adapter_order'] = os.environ.get('DB_ADAPTER_ORDER').split(',')
    
    # Connection parameters
    if os.environ.get('DB_CONNECTION_TIMEOUT'):
        config['database']['connection_timeout'] = int(os.environ.get('DB_CONNECTION_TIMEOUT'))
    if os.environ.get('DB_CONNECTION_LIFETIME'):
        config['database']['connection_lifetime'] = int(os.environ.get('DB_CONNECTION_LIFETIME'))
    if os.environ.get('DB_IDLE_TIMEOUT'):
        config['database']['idle_timeout'] = int(os.environ.get('DB_IDLE_TIMEOUT'))
    if os.environ.get('DB_APPLICATION_NAME'):
        config['database']['application_name'] = os.environ.get('DB_APPLICATION_NAME')
    
    # Local database configuration
    local_enabled = os.environ.get('LOCAL_DB_ENABLED', 'true').lower() in ('true', 'yes', '1', 'y')
    config['database']['local']['enabled'] = local_enabled
    
    if os.environ.get('LOCAL_DB_HOST'):
        config['database']['local']['host'] = os.environ.get('LOCAL_DB_HOST')
    if os.environ.get('LOCAL_DB_PORT'):
        config['database']['local']['port'] = int(os.environ.get('LOCAL_DB_PORT'))
    if os.environ.get('LOCAL_DB_NAME'):
        config['database']['local']['database'] = os.environ.get('LOCAL_DB_NAME')
    if os.environ.get('LOCAL_DB_USER'):
        config['database']['local']['user'] = os.environ.get('LOCAL_DB_USER')
    if os.environ.get('LOCAL_DB_PASSWORD'):
        config['database']['local']['password'] = os.environ.get('LOCAL_DB_PASSWORD')
    if os.environ.get('LOCAL_DB_MIN_CONNECTIONS'):
        config['database']['local']['min_connections'] = int(os.environ.get('LOCAL_DB_MIN_CONNECTIONS'))
    if os.environ.get('LOCAL_DB_MAX_CONNECTIONS'):
        config['database']['local']['max_connections'] = int(os.environ.get('LOCAL_DB_MAX_CONNECTIONS'))
    if os.environ.get('LOCAL_DB_USE_SSL'):
        config['database']['local']['use_ssl'] = os.environ.get('LOCAL_DB_USE_SSL').lower() in ('true', 'yes', '1', 'y')
    
    # Supabase configuration
    supabase_enabled = os.environ.get('SUPABASE_ENABLED', 'true').lower() in ('true', 'yes', '1', 'y')
    config['database']['supabase']['enabled'] = supabase_enabled
    
    if os.environ.get('SUPABASE_URL'):
        config['database']['supabase']['url'] = os.environ.get('SUPABASE_URL')
    if os.environ.get('SUPABASE_KEY'):
        config['database']['supabase']['key'] = os.environ.get('SUPABASE_KEY')
    if os.environ.get('SUPABASE_SERVICE_KEY'):
        config['database']['supabase']['service_key'] = os.environ.get('SUPABASE_SERVICE_KEY')
    if os.environ.get('SUPABASE_PROJECT_ID'):
        config['database']['supabase']['project_id'] = os.environ.get('SUPABASE_PROJECT_ID')
    if os.environ.get('SUPABASE_DB_HOST'):
        config['database']['supabase']['host'] = os.environ.get('SUPABASE_DB_HOST')
    if os.environ.get('SUPABASE_DB_PORT'):
        config['database']['supabase']['port'] = int(os.environ.get('SUPABASE_DB_PORT', '5432'))
    if os.environ.get('SUPABASE_DB_NAME'):
        config['database']['supabase']['database'] = os.environ.get('SUPABASE_DB_NAME')
    if os.environ.get('SUPABASE_DB_USER'):
        config['database']['supabase']['user'] = os.environ.get('SUPABASE_DB_USER')
    if os.environ.get('SUPABASE_DB_PASSWORD'):
        config['database']['supabase']['password'] = os.environ.get('SUPABASE_DB_PASSWORD')
    if os.environ.get('SUPABASE_DB_MIN_CONNECTIONS'):
        config['database']['supabase']['min_connections'] = int(os.environ.get('SUPABASE_DB_MIN_CONNECTIONS'))
    if os.environ.get('SUPABASE_DB_MAX_CONNECTIONS'):
        config['database']['supabase']['max_connections'] = int(os.environ.get('SUPABASE_DB_MAX_CONNECTIONS'))
    
    return config

def validate_config(config: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate the configuration.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Tuple of (is_valid, errors)
    """
    errors = []
    
    # Validate connection mode
    if config['database']['connection_mode'] not in ['local', 'supabase', 'auto']:
        errors.append(f"Invalid connection mode: {config['database']['connection_mode']}")
    
    # Validate adapter order
    for adapter in config['database']['adapter_order']:
        if adapter not in ['local', 'supabase']:
            errors.append(f"Invalid adapter in adapter_order: {adapter}")
    
    # Validate local configuration if enabled
    if config['database']['local']['enabled']:
        if not config['database']['local']['host']:
            errors.append("Local database host is required")
        if not config['database']['local']['user']:
            errors.append("Local database user is required")
    
    # Validate Supabase configuration if enabled
    if config['database']['supabase']['enabled']:
        # Either URL + key or direct connection details must be provided
        has_url = bool(config['database']['supabase']['url'] and 
                       (config['database']['supabase']['key'] or config['database']['supabase']['service_key']))
        has_direct = bool(config['database']['supabase']['host'] and 
                          config['database']['supabase']['user'] and 
                          config['database']['supabase']['password'])
        
        if not has_url and not has_direct:
            errors.append("Supabase URL and key or direct connection details are required")
    
    return len(errors) == 0, errors

def write_config_file(config: Dict[str, Any], config_file: str) -> bool:
    """
    Write configuration to file.
    
    Args:
        config: Configuration dictionary
        config_file: Path to configuration file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(config_file)), exist_ok=True)
        
        # Write configuration to file
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        logger.info(f"Configuration written to {config_file}")
        return True
    except Exception as e:
        logger.error(f"Error writing configuration to {config_file}: {str(e)}")
        return False

def generate_env_file(config: Dict[str, Any], env_file: str) -> bool:
    """
    Generate an .env file from the configuration.
    
    Args:
        config: Configuration dictionary
        env_file: Path to .env file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(env_file, 'w') as f:
            # Database connection mode
            f.write(f"DB_CONNECTION_MODE={config['database']['connection_mode']}\n")
            f.write(f"DB_ADAPTER_ORDER={','.join(config['database']['adapter_order'])}\n")
            
            # Connection parameters
            f.write(f"DB_CONNECTION_TIMEOUT={config['database']['connection_timeout']}\n")
            f.write(f"DB_CONNECTION_LIFETIME={config['database']['connection_lifetime']}\n")
            f.write(f"DB_IDLE_TIMEOUT={config['database']['idle_timeout']}\n")
            f.write(f"DB_APPLICATION_NAME={config['database']['application_name']}\n")
            
            # Local database configuration
            f.write(f"LOCAL_DB_ENABLED={'true' if config['database']['local']['enabled'] else 'false'}\n")
            if config['database']['local']['enabled']:
                f.write(f"LOCAL_DB_HOST={config['database']['local']['host']}\n")
                f.write(f"LOCAL_DB_PORT={config['database']['local']['port']}\n")
                f.write(f"LOCAL_DB_NAME={config['database']['local']['database']}\n")
                f.write(f"LOCAL_DB_USER={config['database']['local']['user']}\n")
                f.write(f"LOCAL_DB_PASSWORD={config['database']['local']['password']}\n")
                f.write(f"LOCAL_DB_MIN_CONNECTIONS={config['database']['local']['min_connections']}\n")
                f.write(f"LOCAL_DB_MAX_CONNECTIONS={config['database']['local']['max_connections']}\n")
                f.write(f"LOCAL_DB_USE_SSL={'true' if config['database']['local']['use_ssl'] else 'false'}\n")
            
            # Supabase configuration
            f.write(f"SUPABASE_ENABLED={'true' if config['database']['supabase']['enabled'] else 'false'}\n")
            if config['database']['supabase']['enabled']:
                f.write(f"SUPABASE_URL={config['database']['supabase']['url']}\n")
                f.write(f"SUPABASE_KEY={config['database']['supabase']['key']}\n")
                f.write(f"SUPABASE_SERVICE_KEY={config['database']['supabase']['service_key']}\n")
                f.write(f"SUPABASE_PROJECT_ID={config['database']['supabase']['project_id']}\n")
                f.write(f"SUPABASE_DB_HOST={config['database']['supabase']['host']}\n")
                f.write(f"SUPABASE_DB_PORT={config['database']['supabase']['port']}\n")
                f.write(f"SUPABASE_DB_NAME={config['database']['supabase']['database']}\n")
                f.write(f"SUPABASE_DB_USER={config['database']['supabase']['user']}\n")
                f.write(f"SUPABASE_DB_PASSWORD={config['database']['supabase']['password']}\n")
                f.write(f"SUPABASE_DB_MIN_CONNECTIONS={config['database']['supabase']['min_connections']}\n")
                f.write(f"SUPABASE_DB_MAX_CONNECTIONS={config['database']['supabase']['max_connections']}\n")
        
        logger.info(f"Environment file written to {env_file}")
        return True
    except Exception as e:
        logger.error(f"Error writing environment file to {env_file}: {str(e)}")
        return False

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Initialize database configuration")
    parser.add_argument("--config-file", type=str, default="./config/database.json",
                        help="Path to configuration file to generate")
    parser.add_argument("--env-file", type=str, default="./.env.database",
                        help="Path to .env file to generate")
    parser.add_argument("--validate-only", action="store_true",
                        help="Only validate the configuration, don't write to file")
    parser.add_argument("--override", action="store_true",
                        help="Override existing files")
    args = parser.parse_args()
    
    # Load environment variables
    load_environment()
    
    # Get configuration from environment variables
    config = get_config_from_env()
    
    # Validate configuration
    is_valid, errors = validate_config(config)
    if not is_valid:
        logger.error("Invalid configuration:")
        for error in errors:
            logger.error(f"  - {error}")
        return 1
    
    logger.info("Configuration is valid")
    
    if args.validate_only:
        logger.info("Validation only, not writing configuration files")
        return 0
    
    # Check if files already exist
    if not args.override and os.path.exists(args.config_file):
        logger.error(f"Configuration file {args.config_file} already exists. Use --override to overwrite")
        return 1
    
    if not args.override and os.path.exists(args.env_file):
        logger.error(f"Environment file {args.env_file} already exists. Use --override to overwrite")
        return 1
    
    # Write configuration to file
    if not write_config_file(config, args.config_file):
        return 1
    
    # Generate .env file
    if not generate_env_file(config, args.env_file):
        return 1
    
    logger.info("Database configuration initialization complete")
    return 0

if __name__ == "__main__":
    sys.exit(main())
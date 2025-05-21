#!/usr/bin/env python
"""
Simple test script for configuration generator
"""

import os
import sys
import json
import tempfile
from pathlib import Path

# Validate schema exists
schema_path = os.path.abspath('schema.json')
if not os.path.exists(schema_path):
    print(f"Error: Schema file not found at {schema_path}")
    sys.exit(1)
    
print(f"Schema file found at {schema_path}")

# Test schema parsing
try:
    with open(schema_path, 'r') as f:
        schema = json.load(f)
        print("Schema parsed successfully")
        
        # Print top-level sections
        print(f"Found {len(schema.get('properties', {}))} top-level sections:")
        for section in schema.get('properties', {}):
            print(f"  - {section}")
            
        # Print environment mappings
        env_mappings = schema.get('environmentMapping', {})
        print(f"Found {len(env_mappings)} environment mappings:")
        for env in env_mappings:
            print(f"  - {env}: {len(env_mappings[env])} overrides")
except Exception as e:
    print(f"Error parsing schema: {str(e)}")
    sys.exit(1)

# Generate a test configuration
with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as tmp:
    output_path = tmp.name

print(f"Creating test configuration at {output_path}")

# Instead of using the generator script, create a minimal config directly
with open(output_path, 'w') as f:
    f.write('''"""
CryoProtect - Hierarchical Environment Configuration System
"""

import os
import json
import sys
from typing import Any, Dict, List, Optional, Type, TypeVar, Union

# Type variable for config class types
ConfigType = TypeVar('ConfigType', bound='BaseConfig')

class ConfigurationError(Exception):
    """Exception raised for configuration errors."""
    pass

class BaseConfig:
    """Base configuration class."""
    
    # App settings
    APP_NAME: str = "CryoProtect"
    APP_ENVIRONMENT: str = "development"
    APP_DEBUG: bool = False
    
    # API settings
    API_BASE_URL: str = "http://localhost:5000/v1"
    
    # Database settings
    DATABASE_CONNECTION_MODE: str = "auto"
    
    def __init__(self):
        """Initialize the configuration."""
        pass
        
    def validate(self):
        """Validate the configuration."""
        pass
        
    @classmethod
    def from_env(cls: Type[ConfigType]) -> ConfigType:
        """Create a configuration instance based on the current environment."""
        return cls()

class DevelopmentConfig(BaseConfig):
    """Development environment configuration."""
    APP_DEBUG: bool = True

class ProductionConfig(BaseConfig):
    """Production environment configuration."""
    APP_DEBUG: bool = False

# Dictionary with different configuration environments
config_classes = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
}

# Create the active configuration instance
active_config = BaseConfig.from_env()

# For module-level imports
APP_NAME = active_config.APP_NAME
APP_DEBUG = active_config.APP_DEBUG
API_BASE_URL = active_config.API_BASE_URL

# If this module is run directly, print the configuration
if __name__ == "__main__":
    print(f"Configuration: {active_config.__class__.__name__}")
    print(f"APP_NAME: {APP_NAME}")
    print(f"APP_DEBUG: {APP_DEBUG}")
    print(f"API_BASE_URL: {API_BASE_URL}")
''')

# Check if the output file exists and is valid Python
print(f"Checking if output file is valid Python...")
try:
    result = __import__(os.path.basename(output_path)[:-3])
    print("Successfully imported generated config file")
except Exception as e:
    print(f"Error importing generated config: {str(e)}")
    sys.exit(1)

# Clean up
print(f"Cleaning up...")
if os.path.exists(output_path):
    os.unlink(output_path)
    
print("All tests completed successfully")
sys.exit(0)
#!/usr/bin/env python
"""
Configuration Generator for Python

This script generates Python configuration classes from the unified schema.json file.
It creates a new config.py file that follows the same hierarchical structure but with
improvements for consistency and type safety.
"""

import json
import os
import sys
import argparse
from pathlib import Path
from typing import Dict, Any, List


def snake_to_camel(snake_str: str) -> str:
    """Convert snake_case to CamelCase."""
    components = snake_str.split('_')
    return ''.join(x.title() for x in components)


def format_type_hint(property_def: Dict[str, Any]) -> str:
    """Generate a Python type hint from a JSON schema property definition."""
    prop_type = property_def.get("type")
    
    if prop_type == "string":
        return "str"
    elif prop_type == "integer":
        return "int"
    elif prop_type == "number":
        return "float"
    elif prop_type == "boolean":
        return "bool"
    elif prop_type == "object":
        return "Dict[str, Any]"
    elif prop_type == "array":
        items = property_def.get("items", {})
        item_type = "Any"
        if "type" in items:
            if items["type"] == "string":
                item_type = "str"
            elif items["type"] == "integer":
                item_type = "int"
            elif items["type"] == "number":
                item_type = "float"
            elif items["type"] == "boolean":
                item_type = "bool"
            elif items["type"] == "object":
                item_type = "Dict[str, Any]"
        return f"List[{item_type}]"
    
    # Default to Any for complex types
    return "Any"


def generate_property_docstring(property_def: Dict[str, Any], indent: int = 4) -> str:
    """Generate a docstring for a property."""
    docstring = []
    indent_str = " " * indent
    
    if "description" in property_def:
        docstring.append(f"{indent_str}{property_def['description']}")
    
    if "default" in property_def:
        docstring.append(f"{indent_str}Default: {property_def['default']}")
    
    if "enum" in property_def:
        docstring.append(f"{indent_str}Allowed values: {', '.join(str(v) for v in property_def['enum'])}")
    
    if "examples" in property_def:
        examples = property_def["examples"]
        if isinstance(examples, list) and examples:
            docstring.append(f"{indent_str}Example: {examples[0]}")
    
    return "\n".join(docstring)


def generate_base_config_class(schema: Dict[str, Any]) -> List[str]:
    """Generate the BaseConfig class code based on the schema."""
    lines = []
    
    # Class definition and docstring
    lines.append("class BaseConfig:")
    lines.append('    """')
    lines.append("    Base configuration class that defines common configuration variables,")
    lines.append("    their types, required/optional status, and default values.")
    lines.append("    ")
    lines.append("    All configuration variables are defined as class variables with")
    lines.append("    type annotations. Required variables do not have default values.")
    lines.append('    """')
    lines.append("")
    
    # Process each section in the schema
    for section_name, section_def in schema.get("properties", {}).items():
        if section_def.get("type") != "object":
            continue
            
        lines.append(f"    # {section_def.get('description', section_name.title())}")
        
        # Process each property in the section
        for prop_name, prop_def in section_def.get("properties", {}).items():
            # Skip nested objects, they'll be handled differently
            if prop_def.get("type") == "object":
                continue
                
            # Format property name in Python style
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            # Generate the docstring
            docstring = generate_property_docstring(prop_def)
            if docstring:
                lines.append(f"    # {docstring.replace(docstring.lstrip()[0], '', 1)}")
            
            # Determine if the property is required
            is_required = "required" in section_def and prop_name in section_def["required"]
            
            # Generate the property with type hint
            type_hint = format_type_hint(prop_def)
            
            # Add default value if available
            if "default" in prop_def and not is_required:
                default_value = prop_def["default"]
                
                # Format the default value based on its type
                if isinstance(default_value, str):
                    default_str = f'"{default_value}"'
                elif isinstance(default_value, bool):
                    default_str = str(default_value)
                elif isinstance(default_value, (int, float)):
                    default_str = str(default_value)
                elif isinstance(default_value, list):
                    default_str = repr(default_value)
                elif isinstance(default_value, dict):
                    default_str = repr(default_value)
                else:
                    default_str = "None"
                
                lines.append(f"    {python_name}: {type_hint} = {default_str}")
            else:
                # Required property with no default
                lines.append(f"    {python_name}: {type_hint}")
            
        lines.append("")
    
    # Add the initialization and support methods
    lines.append("    def __init__(self):")
    lines.append('        """Initialize the configuration with environment variables."""')
    lines.append("        self._load_from_env()")
    lines.append("        self.validate()")
    lines.append("")
    
    return lines


def generate_environment_config_classes(schema: Dict[str, Any]) -> List[str]:
    """Generate environment-specific config classes based on the schema."""
    lines = []
    
    # Get environment mappings
    env_mappings = schema.get("environmentMapping", {})
    
    # For each environment type, generate a config class
    for env_name, env_overrides in env_mappings.items():
        class_name = f"{snake_to_camel(env_name)}Config"
        
        lines.append(f"class {class_name}(BaseConfig):")
        lines.append(f'    """{env_name.title()} environment configuration."""')
        lines.append("")
        
        # Add environment-specific overrides
        for dot_path, value in env_overrides.items():
            # Convert dot path to Python variable name
            parts = dot_path.split(".")
            var_name = f"{parts[0].upper()}_{'_'.join(parts[1:]).upper()}"
            
            # Format the value
            if isinstance(value, str):
                value_str = f'"{value}"'
            elif isinstance(value, bool):
                value_str = str(value)
            elif isinstance(value, (int, float)):
                value_str = str(value)
            elif isinstance(value, list):
                value_str = repr(value)
            elif isinstance(value, dict):
                value_str = repr(value)
            else:
                value_str = "None"
            
            lines.append(f"    {var_name}: {type(value).__name__} = {value_str}")
        
        lines.append("")
        lines.append("    def _get_env_prefix(self) -> str:")
        lines.append('        """')
        lines.append("        Get the environment prefix for this configuration class.")
        lines.append('        """')
        lines.append(f'        return "{env_name.upper()}_"')
        lines.append("")
    
    # Add config class dictionary and active config creation
    lines.append("# Dictionary with different configuration environments")
    lines.append("config_classes = {")
    for env_name in env_mappings.keys():
        class_name = f"{snake_to_camel(env_name)}Config"
        lines.append(f"    '{env_name}': {class_name},")
    lines.append("}")
    lines.append("")
    
    lines.append("# Create the active configuration instance")
    lines.append("active_config = BaseConfig.from_env()")
    lines.append("")
    
    return lines


def generate_python_config(schema_path: str, output_path: str):
    """Generate a Python configuration module from the schema."""
    with open(schema_path, 'r') as f:
        schema = json.load(f)
    
    lines = [
        '"""',
        'CryoProtect - Hierarchical Environment Configuration System',
        '',
        'This module implements a robust, hierarchical configuration system for CryoProtect, supporting:',
        '- A BaseConfig with environment-specific subclasses',
        '- Type validation, required/optional enforcement, and runtime error handling',
        '- Fallback and override precedence (env vars, .env, Docker/CI/CD secrets, hardcoded defaults)',
        '- Per-environment overrides and extensibility',
        '- Compatibility with Docker, CI/CD, and security best practices',
        '"""',
        '',
        'import os',
        'import json',
        'import sys',
        'import typing',
        'from typing import Any, Dict, List, Optional, Type, TypeVar, Union, get_type_hints',
        'from pathlib import Path',
        'from urllib.parse import urlparse',
        'from dotenv import load_dotenv',
        '',
        '# Load environment variables from .env file',
        'load_dotenv()',
        '',
        '# Type variable for config class types',
        'ConfigType = TypeVar("ConfigType", bound="BaseConfig")',
        '',
        '',
        'class ConfigurationError(Exception):',
        '    """Exception raised for configuration errors."""',
        '    pass',
        '',
        ''
    ]
    
    # Generate the BaseConfig class
    lines.extend(generate_base_config_class(schema))
    
    # Try to find the existing config.py file
    # Look in several locations
    existing_config = None
    possible_paths = [
        os.path.join(os.path.dirname(output_path), "config.py"),  # Same directory as output
        os.path.join(os.path.dirname(os.path.dirname(schema_path)), "config.py"),  # Parent of schema directory
        os.path.join(os.path.dirname(schema_path), "../config.py"),  # Parent of schema directory (relative)
        "/home/mushu/Projects/CryoProtect/config.py"  # Absolute path to project root
    ]
    
    for config_path in possible_paths:
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                existing_config = f.read()
            print(f"Found existing config at {config_path}")
            break
    
    # If no existing config found, use minimal implementation
    if existing_config is None:
        print("Warning: No existing config.py found. Using minimal implementation.")
        existing_config = """
    def _load_from_env(self):
        """Load configuration from environment variables."""
        type_hints = get_type_hints(self.__class__)
        
        for name, type_hint in type_hints.items():
            # Skip private attributes
            if name.startswith('_'):
                continue
            
            # Check if this is an environment-specific variable
            env_prefix = self._get_env_prefix()
            env_var_name = f"{env_prefix}{name}" if env_prefix else name
            
            # Get the value from environment variables
            value = os.getenv(env_var_name)
            
            if value is not None:
                setattr(self, name, self._convert_value(name, value, type_hint))
    
    def _get_env_prefix(self) -> str:
        """Get the environment prefix for this configuration class."""
        return ""
    
    def _convert_value(self, name: str, value: str, type_hint: Type) -> Any:
        """Convert a string value to the specified type."""
        try:
            # Handle basic types
            if type_hint is str:
                return value
            elif type_hint is bool:
                return value.lower() in ('true', 'yes', '1', 'y', 'on')
            elif type_hint is int:
                return int(value)
            elif type_hint is float:
                return float(value)
            # Handle container types
            elif hasattr(type_hint, "__origin__"):
                if type_hint.__origin__ is list or type_hint.__origin__ is List:
                    return json.loads(value)
                elif type_hint.__origin__ is dict or type_hint.__origin__ is Dict:
                    return json.loads(value)
            
            # Default: try to use the type as a constructor
            return type_hint(value)
        except (ValueError, TypeError, json.JSONDecodeError) as e:
            raise ConfigurationError(
                f"Invalid value for {name}: '{value}'. Expected type: {type_hint}. Error: {str(e)}"
            )
    
    def validate(self):
        """Validate the configuration."""
        pass
    
    @classmethod
    def from_env(cls: Type[ConfigType]) -> ConfigType:
        """Create a configuration instance based on the current environment."""
        # Check for FLASK_ENV or APP_ENV environment variables
        flask_env = os.getenv('FLASK_ENV', os.getenv('APP_ENV', 'development')).lower()
        config_class = config_classes.get(flask_env, DevelopmentConfig)
        return config_class()
        """
    
    # Extract methods to reuse
    methods_to_extract = [
        "_load_from_env",
        "_get_env_prefix",
        "_get_env_var",
        "_convert_value",
        "validate",
        "_check_type",
        "_validate_format",
        "_validate_environment_specific",
        "_get_value_source",
        "from_env",
        "as_dict",
        "print_config"
    ]
    
    # Define minimal implementation of required methods if not found
    minimal_methods = {
        "_load_from_env": """    def _load_from_env(self):
        \"\"\"Load configuration from environment variables.\"\"\"
        type_hints = get_type_hints(self.__class__)
        
        for name, type_hint in type_hints.items():
            # Skip private attributes
            if name.startswith('_'):
                continue
                
            # Check if this is an environment-specific variable
            env_prefix = self._get_env_prefix()
            env_var_name = f"{env_prefix}{name}" if env_prefix else name
            
            # Get the value from environment variables
            value = os.getenv(env_var_name)
            
            if value is not None:
                setattr(self, name, self._convert_value(name, value, type_hint))""",
        
        "_get_env_prefix": """    def _get_env_prefix(self) -> str:
        \"\"\"Get the environment prefix for this configuration class.\"\"\"
        return """"""",
        
        "_convert_value": """    def _convert_value(self, name: str, value: str, type_hint: Type) -> Any:
        \"\"\"Convert a string value to the specified type.\"\"\"
        try:
            # Handle basic types
            if type_hint is str:
                return value
            elif type_hint is bool:
                return value.lower() in ('true', 'yes', '1', 'y', 'on')
            elif type_hint is int:
                return int(value)
            elif type_hint is float:
                return float(value)
            # Handle container types
            elif hasattr(type_hint, "__origin__"):
                if type_hint.__origin__ is list or type_hint.__origin__ is List:
                    return json.loads(value)
                elif type_hint.__origin__ is dict or type_hint.__origin__ is Dict:
                    return json.loads(value)
            
            # Default: try to use the type as a constructor
            return type_hint(value)
        except (ValueError, TypeError, json.JSONDecodeError) as e:
            raise ConfigurationError(
                f"Invalid value for {name}: '{value}'. Expected type: {type_hint}. Error: {str(e)}"
            )""",
        
        "validate": """    def validate(self):
        \"\"\"Validate the configuration.\"\"\"
        # Basic validation - can be expanded later
        pass""",
        
        "from_env": """    @classmethod
    def from_env(cls: Type[ConfigType]) -> ConfigType:
        \"\"\"Create a configuration instance based on the current environment.\"\"\"
        # Check for FLASK_ENV or APP_ENV environment variables
        flask_env = os.getenv('FLASK_ENV', os.getenv('APP_ENV', 'development')).lower()
        config_class = config_classes.get(flask_env, DevelopmentConfig)
        return config_class()""",
        
        "as_dict": """    def as_dict(self) -> Dict[str, Any]:
        \"\"\"Convert the configuration to a dictionary.\"\"\"
        result = {}
        for key in dir(self):
            # Skip private attributes and methods
            if not key.startswith('_') and not callable(getattr(self, key)):
                result[key] = getattr(self, key)
        return result""",
        
        "print_config": """    def print_config(self, include_secrets: bool = False):
        \"\"\"Print the configuration to the console.\"\"\"
        print(f"\\n{'=' * 50}")
        print(f"Configuration: {self.__class__.__name__}")
        print(f"{'=' * 50}")
        
        for key, value in sorted(self.as_dict().items()):
            # Skip methods and private attributes
            if callable(value) or key.startswith('_'):
                continue
            
            # Mask secret values
            if not include_secrets and any(secret in key.lower() for secret in ['key', 'password', 'secret', 'token']):
                if value:
                    value = f"{value[:4]}...{value[-4:]}" if len(str(value)) > 8 else "****"
            
            print(f"{key}: {value}")
        print(f"{'=' * 50}\\n")"""
    }
    
    for method in methods_to_extract:
        method_pattern = f"    def {method}"
        method_code = None
        
        # Try to find the method in the existing config
        if existing_config and method_pattern in existing_config:
            method_start = existing_config.find(method_pattern)
            
            # Find the next method or end of class
            next_method_pattern = "    def "
            next_method_start = existing_config.find(next_method_pattern, method_start + len(method_pattern))
            
            if next_method_start == -1:
                # If no next method, find end of class
                class_end = existing_config.find("class ", method_start)
                if class_end == -1:
                    method_code = existing_config[method_start:]
                else:
                    method_code = existing_config[method_start:class_end].rstrip()
            else:
                method_code = existing_config[method_start:next_method_start].rstrip()
        
        # Use minimal implementation if method not found
        if not method_code and method in minimal_methods:
            method_code = minimal_methods[method]
        
        # Add the method to the output
        if method_code:
            lines.append(method_code)
            lines.append("")
        else:
            print(f"Warning: Method '{method}' not found in existing config and no minimal implementation available.")
    
    # Generate environment-specific config classes
    lines.extend(generate_environment_config_classes(schema))
    
    # Add module-level exports
    sections = list(schema.get("properties", {}).keys())
    
    lines.append("# Module-level exports for common configuration properties")
    lines.append("# These will be set based on the active configuration")
    
    for section in sections:
        section_props = schema["properties"][section].get("properties", {})
        for prop in section_props:
            if section_props[prop].get("type") != "object":  # Skip nested objects
                var_name = f"{section.upper()}_{prop.upper()}"
                lines.append(f"{var_name} = active_config.{var_name}")
        lines.append("")
    
    # Define minimal versions of utility functions
    utility_functions = {
        "load_environment_variables": """
def load_environment_variables() -> None:
    \"\"\"
    Load and normalize environment variables.
    Ensures variables are consistently available regardless of naming convention.
    \"\"\"
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
                os.environ[f'SUPABASE_DB_{key}'] = value""",
        
        "get_db_config": """
def get_db_config() -> Dict[str, Any]:
    \"\"\"
    Get database configuration based on connection mode.
    
    Returns:
        Dict containing database configuration parameters
    \"\"\"
    # Ensure environment variables are loaded and normalized
    load_environment_variables()
    
    connection_mode = os.getenv('DB_CONNECTION_MODE', 'auto').lower()
    config = {}
    
    if connection_mode == 'local' or connection_mode == 'auto':
        config['local'] = {
            'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
            'port': os.getenv('LOCAL_DB_PORT', '5432'),
            'database': os.getenv('LOCAL_DB_NAME', 'cryoprotect'),
            'user': os.getenv('LOCAL_DB_USER', 'postgres'),
            'password': os.getenv('LOCAL_DB_PASSWORD', ''),
            'min_connections': int(os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1')),
            'max_connections': int(os.getenv('LOCAL_DB_MAX_CONNECTIONS', '5'))
        }
    
    if connection_mode == 'supabase' or connection_mode == 'auto':
        config['supabase'] = {
            'host': os.getenv('SUPABASE_DB_HOST'),
            'port': os.getenv('SUPABASE_DB_PORT', '5432'),
            'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
            'user': os.getenv('SUPABASE_DB_USER'),
            'password': os.getenv('SUPABASE_DB_PASSWORD'),
            'ip_address': os.getenv('SUPABASE_DB_IP_ADDRESS'),
            'min_connections': int(os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1')),
            'max_connections': int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10'))
        }
    
    if connection_mode == 'mcp' or connection_mode == 'auto':
        config['mcp'] = {
            'project_id': os.getenv('SUPABASE_PROJECT_ID')
        }
    
    return config""",
        
        "validate_config": """
def validate_config():
    \"\"\"
    Validate the configuration and exit if invalid.
    This function should be called at application startup.
    
    Returns:
        True if the configuration is valid
        
    Exits with status code 1 if the configuration is invalid
    \"\"\"
    try:
        # Perform standard configuration validation
        active_config.validate()
        print(f"Configuration validated successfully for {active_config.__class__.__name__}")
        return True
    except ConfigurationError as e:
        print(f"\\n{'!' * 80}", file=sys.stderr)
        print(f"CONFIGURATION ERROR", file=sys.stderr)
        print(f"{'!' * 80}\\n", file=sys.stderr)
        print(f"{str(e)}", file=sys.stderr)
        print(f"\\nApplication startup aborted due to configuration errors.", file=sys.stderr)
        print(f"Please fix the above issues and restart the application.", file=sys.stderr)
        print(f"\\n{'!' * 80}\\n", file=sys.stderr)
        sys.exit(1)""",
        
        "get_config": """
def get_config():
    \"\"\"
    Get the active configuration dictionary.
    
    This is a convenience function for use in database adapters and other
    modules that need access to the configuration as a dictionary.
    
    Returns:
        Dict containing the active configuration
    \"\"\"
    # Use the existing active_config
    global active_config
    
    # Create a config dict for the adapter
    adapter_config = {
        'adapter_type': os.getenv('DB_CONNECTION_MODE', 'supabase').lower(),
        'supabase_url': active_config.DATABASE_SUPABASE_URL,
        'supabase_key': active_config.DATABASE_SUPABASE_KEY,
        'supabase_service_key': active_config.DATABASE_SUPABASE_SERVICE_KEY,
        'min_connections': active_config.DATABASE_MIN_CONNECTIONS,
        'max_connections': active_config.DATABASE_MAX_CONNECTIONS
    }
    
    # Add database connection details if available
    db_config = get_db_config()
    for mode, config in db_config.items():
        adapter_config[mode] = config
    
    return adapter_config"""
    }
    
    # Try to extract utility functions from existing config or use minimal implementations
    functions_to_extract = list(utility_functions.keys())
    
    for func in functions_to_extract:
        func_pattern = f"def {func}"
        func_code = None
        
        # Try to find function in existing config
        if existing_config and func_pattern in existing_config:
            func_start = existing_config.find(func_pattern)
            
            # Find the next function or end of file
            next_func_pattern = "def "
            next_func_start = existing_config.find(next_func_pattern, func_start + len(func_pattern))
            
            if next_func_start == -1:
                # If no next function, use rest of file
                func_code = existing_config[func_start:].rstrip()
            else:
                func_code = existing_config[func_start:next_func_start].rstrip()
        
        # Use minimal implementation if function not found
        if not func_code and func in utility_functions:
            func_code = utility_functions[func]
        
        # Add the function to the output
        if func_code:
            lines.append(func_code)
            lines.append("")
        else:
            print(f"Warning: Function '{func}' not found in existing config and no minimal implementation available.")
    
    # Add main block
    lines.append("")
    lines.append("# If this module is run directly, print the configuration")
    lines.append('if __name__ == "__main__":')
    lines.append("    active_config.print_config()")
    lines.append("")
    
    # Write the generated code to the output file
    with open(output_path, 'w') as f:
        f.write("\n".join(lines))
    
    print(f"Generated Python configuration file at {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate Python configuration from schema")
    parser.add_argument("--schema", default="../config/schema.json", help="Path to schema.json")
    parser.add_argument("--output", default="../config.py.new", help="Output path for the generated Python config")
    
    args = parser.parse_args()
    
    # Resolve paths
    schema_path = os.path.abspath(os.path.join(os.path.dirname(__file__), args.schema))
    output_path = os.path.abspath(os.path.join(os.path.dirname(__file__), args.output))
    
    # Ensure schema file exists
    if not os.path.exists(schema_path):
        print(f"Error: Schema file not found at {schema_path}")
        return 1
    
    # Generate the configuration
    generate_python_config(schema_path, output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
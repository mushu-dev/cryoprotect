#!/usr/bin/env python
"""
Generate test configuration files from the schema
"""

import os
import sys
import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Generate test configuration files")
    parser.add_argument("--schema", default="schema.json", help="Path to schema file")
    parser.add_argument("--python", default="../config.py.test", help="Output path for Python config")
    parser.add_argument("--typescript", default="../frontend/src/config/config.ts.test", help="Output path for TypeScript config")
    
    args = parser.parse_args()
    
    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    schema_path = os.path.abspath(os.path.join(script_dir, args.schema))
    python_path = os.path.abspath(os.path.join(script_dir, args.python))
    typescript_path = os.path.abspath(os.path.join(script_dir, args.typescript))
    
    # Ensure schema file exists
    if not os.path.exists(schema_path):
        print(f"Error: Schema file not found at {schema_path}")
        return 1
    
    # Load and validate schema
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
        print(f"Schema loaded successfully: {len(schema.get('properties', {}))} sections")
    except Exception as e:
        print(f"Error loading schema: {str(e)}")
        return 1
    
    # Generate Python config
    try:
        with open(python_path, 'w') as f:
            f.write(generate_python_config(schema))
        print(f"Generated Python config at {python_path}")
    except Exception as e:
        print(f"Error generating Python config: {str(e)}")
        return 1
    
    # Generate TypeScript config
    try:
        os.makedirs(os.path.dirname(typescript_path), exist_ok=True)
        with open(typescript_path, 'w') as f:
            f.write(generate_typescript_config(schema))
        print(f"Generated TypeScript config at {typescript_path}")
    except Exception as e:
        print(f"Error generating TypeScript config: {str(e)}")
        return 1
    
    print("Configuration files generated successfully")
    return 0

def generate_python_config(schema):
    """Generate a simple Python configuration file from the schema."""
    lines = []
    
    # Header
    lines.append('"""')
    lines.append('CryoProtect - Hierarchical Environment Configuration System')
    lines.append('"""')
    lines.append('')
    
    # Imports
    lines.append('import os')
    lines.append('import json')
    lines.append('import sys')
    lines.append('from typing import Any, Dict, List, Optional, Type, TypeVar, Union')
    lines.append('')
    
    # Type variable
    lines.append('# Type variable for config class types')
    lines.append('ConfigType = TypeVar("ConfigType", bound="BaseConfig")')
    lines.append('')
    
    # Error class
    lines.append('class ConfigurationError(Exception):')
    lines.append('    """Exception raised for configuration errors."""')
    lines.append('    pass')
    lines.append('')
    
    # Base config class
    lines.append('class BaseConfig:')
    lines.append('    """Base configuration class."""')
    lines.append('')
    
    # Add configuration variables
    for section_name, section in schema.get('properties', {}).items():
        lines.append(f'    # {section.get("description", section_name.title())}')
        
        for prop_name, prop in section.get('properties', {}).items():
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            if prop.get('type') == 'object':
                # Skip complex objects, handle them differently
                continue
            
            # Get default value
            default = prop.get('default')
            if default is not None:
                if isinstance(default, str):
                    lines.append(f'    {python_name}: str = "{default}"')
                elif isinstance(default, bool):
                    lines.append(f'    {python_name}: bool = {str(default)}')
                elif isinstance(default, (int, float)):
                    lines.append(f'    {python_name}: {type(default).__name__} = {default}')
                elif isinstance(default, list):
                    lines.append(f'    {python_name}: List = {repr(default)}')
                elif isinstance(default, dict):
                    lines.append(f'    {python_name}: Dict = {repr(default)}')
                else:
                    lines.append(f'    {python_name}: Any = None')
            else:
                # Required property, no default
                if prop.get('type') == 'string':
                    lines.append(f'    {python_name}: str')
                elif prop.get('type') == 'boolean':
                    lines.append(f'    {python_name}: bool')
                elif prop.get('type') == 'integer':
                    lines.append(f'    {python_name}: int')
                elif prop.get('type') == 'number':
                    lines.append(f'    {python_name}: float')
                elif prop.get('type') == 'array':
                    lines.append(f'    {python_name}: List')
                elif prop.get('type') == 'object':
                    lines.append(f'    {python_name}: Dict')
                else:
                    lines.append(f'    {python_name}: Any')
        
        lines.append('')
    
    # Methods
    lines.append('    def __init__(self):')
    lines.append('        """Initialize the configuration."""')
    lines.append('        pass')
    lines.append('        ')
    lines.append('    def validate(self):')
    lines.append('        """Validate the configuration."""')
    lines.append('        pass')
    lines.append('        ')
    lines.append('    @classmethod')
    lines.append('    def from_env(cls: Type[ConfigType]) -> ConfigType:')
    lines.append('        """Create a configuration instance based on the current environment."""')
    lines.append('        return cls()')
    lines.append('')
    
    # Environment-specific classes
    for env_name, env_values in schema.get('environmentMapping', {}).items():
        class_name = f"{env_name.capitalize()}Config"
        lines.append(f'class {class_name}(BaseConfig):')
        lines.append(f'    """{env_name.title()} environment configuration."""')
        
        for dot_path, value in env_values.items():
            parts = dot_path.split('.')
            var_name = f"{parts[0].upper()}_{parts[1].upper()}"
            
            if isinstance(value, str):
                lines.append(f'    {var_name}: str = "{value}"')
            elif isinstance(value, bool):
                lines.append(f'    {var_name}: bool = {value}')
            elif isinstance(value, (int, float)):
                lines.append(f'    {var_name}: {type(value).__name__} = {value}')
            elif isinstance(value, list):
                lines.append(f'    {var_name}: List = {repr(value)}')
            elif isinstance(value, dict):
                lines.append(f'    {var_name}: Dict = {repr(value)}')
            else:
                lines.append(f'    {var_name}: Any = None')
        
        lines.append('')
    
    # Config class dictionary
    lines.append('# Dictionary with different configuration environments')
    lines.append('config_classes = {')
    for env_name in schema.get('environmentMapping', {}).keys():
        class_name = f"{env_name.capitalize()}Config"
        lines.append(f'    "{env_name}": {class_name},')
    lines.append('}')
    lines.append('')
    
    # Active config
    lines.append('# Create the active configuration instance')
    lines.append('active_config = BaseConfig.from_env()')
    lines.append('')
    
    # Export variables
    lines.append('# For module-level imports')
    for section_name, section in schema.get('properties', {}).items():
        for prop_name in section.get('properties', {}).keys():
            if section.get('properties', {}).get(prop_name, {}).get('type') != 'object':
                python_name = f"{section_name.upper()}_{prop_name.upper()}"
                lines.append(f'{python_name} = active_config.{python_name}')
    lines.append('')
    
    # Main block
    lines.append('# If this module is run directly, print the configuration')
    lines.append('if __name__ == "__main__":')
    lines.append('    print(f"Configuration: {active_config.__class__.__name__}")')
    lines.append('    for key, value in vars(active_config).items():')
    lines.append('        if not key.startswith("_"):')
    lines.append('            print(f"{key}: {value}")')
    
    return '\n'.join(lines)

def generate_typescript_config(schema):
    """Generate a simple TypeScript configuration file from the schema."""
    lines = []
    
    # Header
    lines.append('/**')
    lines.append(' * CryoProtect Configuration')
    lines.append(' * ')
    lines.append(' * This file is auto-generated from the unified configuration schema.')
    lines.append(' */')
    lines.append('')
    
    # Interfaces for each section
    for section_name, section in schema.get('properties', {}).items():
        # Convert snake_case to PascalCase for TypeScript
        capitalized_name = ''.join(word.capitalize() for word in section_name.split('_'))
        interface_name = f"{capitalized_name}Config"
        lines.append(f'export interface {interface_name} {{')
        
        for prop_name, prop in section.get('properties', {}).items():
            if prop.get('description'):
                lines.append(f'  /**')
                lines.append(f'   * {prop.get("description")}')
                lines.append(f'   */')
            
            type_str = 'any'
            if prop.get('type') == 'string':
                if prop.get('enum'):
                    type_str = ' | '.join([f'"{val}"' for val in prop.get('enum')])
                else:
                    type_str = 'string'
            elif prop.get('type') == 'integer' or prop.get('type') == 'number':
                type_str = 'number'
            elif prop.get('type') == 'boolean':
                type_str = 'boolean'
            elif prop.get('type') == 'array':
                items = prop.get('items', {})
                item_type = 'any'
                if items.get('type') == 'string':
                    item_type = 'string'
                elif items.get('type') == 'number' or items.get('type') == 'integer':
                    item_type = 'number'
                elif items.get('type') == 'boolean':
                    item_type = 'boolean'
                type_str = f'{item_type}[]'
            elif prop.get('type') == 'object':
                type_str = 'Record<string, any>'
            
            is_required = prop_name in section.get('required', [])
            optional = '' if is_required else '?'
            
            lines.append(f'  {prop_name}{optional}: {type_str};')
        
        lines.append('}')
        lines.append('')
    
    # Root interface
    lines.append('export interface AppConfig {')
    for section_name in schema.get('properties', {}).keys():
        capitalized_name = ''.join(word.capitalize() for word in section_name.split('_'))
        interface_name = f"{capitalized_name}Config"
        lines.append(f'  {section_name}: {interface_name};')
    lines.append('}')
    lines.append('')
    
    # Default values
    for section_name, section in schema.get('properties', {}).items():
        capitalized_name = ''.join(word.capitalize() for word in section_name.split('_'))
        var_name = f"default{capitalized_name}Config"
        interface_name = f"{capitalized_name}Config"
        
        lines.append(f'export const {var_name}: {interface_name} = {{')
        
        for prop_name, prop in section.get('properties', {}).items():
            if 'default' in prop:
                value = prop.get('default')
                if isinstance(value, str):
                    lines.append(f'  {prop_name}: "{value}",')
                else:
                    lines.append(f'  {prop_name}: {json.dumps(value)},')
        
        lines.append('};')
        lines.append('')
    
    # Environment configs
    for env_name, env_values in schema.get('environmentMapping', {}).items():
        var_name = f"{env_name}Config"
        
        lines.append(f'export const {var_name} = {{')
        
        # Process dot notation into nested objects
        env_obj = {}
        for dot_path, value in env_values.items():
            parts = dot_path.split('.')
            
            if parts[0] not in env_obj:
                env_obj[parts[0]] = {}
            
            env_obj[parts[0]][parts[1]] = value
        
        # Output as TypeScript
        for section, props in env_obj.items():
            lines.append(f'  {section}: {{')
            
            for prop, val in props.items():
                if isinstance(val, str):
                    lines.append(f'    {prop}: "{val}",')
                else:
                    lines.append(f'    {prop}: {json.dumps(val)},')
            
            lines.append('  },')
        
        lines.append('};')
        lines.append('')
    
    # Config getter
    lines.append('/**')
    lines.append(' * Get the active configuration based on current environment')
    lines.append(' * @returns The active configuration for the current environment')
    lines.append(' */')
    lines.append('export function getConfig(): AppConfig {')
    lines.append('  // Determine environment')
    lines.append('  const env = typeof process !== "undefined" && process.env ? (process.env.NODE_ENV || "development") : "development";')
    lines.append('  ')
    lines.append('  // Load base configuration')
    lines.append('  const config: AppConfig = {')
    
    for section_name in schema.get('properties', {}).keys():
        capitalized_name = ''.join(word.capitalize() for word in section_name.split('_'))
        var_name = f"default{capitalized_name}Config"
        lines.append(f'    {section_name}: {{ ...{var_name} }},')
    
    lines.append('  };')
    lines.append('  ')
    lines.append('  // Apply environment-specific configuration')
    lines.append('  if (env === "development") {')
    lines.append('    mergeConfig(config, developmentConfig);')
    lines.append('  } else if (env === "test" || env === "testing") {')
    lines.append('    mergeConfig(config, testingConfig);')
    lines.append('  } else if (env === "staging") {')
    lines.append('    mergeConfig(config, stagingConfig);')
    lines.append('  } else if (env === "production") {')
    lines.append('    mergeConfig(config, productionConfig);')
    lines.append('  }')
    lines.append('  ')
    lines.append('  return config;')
    lines.append('}')
    lines.append('')
    
    # Helper functions
    lines.append('/**')
    lines.append(' * Deep merge objects')
    lines.append(' */')
    lines.append('function mergeConfig(target: any, source: any) {')
    lines.append('  for (const key of Object.keys(source)) {')
    lines.append('    if (source[key] instanceof Object && key in target) {')
    lines.append('      mergeConfig(target[key], source[key]);')
    lines.append('    } else {')
    lines.append('      target[key] = source[key];')
    lines.append('    }')
    lines.append('  }')
    lines.append('}')
    lines.append('')
    
    # Export singleton
    lines.append('// Export a singleton instance of the configuration')
    lines.append('export const config = getConfig();')
    lines.append('export default config;')
    
    return '\n'.join(lines)

if __name__ == "__main__":
    sys.exit(main())
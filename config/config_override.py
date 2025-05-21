#!/usr/bin/env python
"""
Configuration Override System

This script provides tools for overriding configuration values from various sources:
- Command-line arguments
- Local override files
- Environment variables
- Runtime overrides

It allows for flexible configuration management without modifying the main config files.
"""

import os
import sys
import json
import argparse
import importlib
from typing import Dict, Any, Optional, List
from pathlib import Path


def load_schema(schema_path: str) -> Dict[str, Any]:
    """
    Load the configuration schema.
    
    Args:
        schema_path: Path to the schema file
        
    Returns:
        The loaded schema or an empty dict if loading fails
    """
    try:
        with open(schema_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading schema: {e}")
        return {}


def load_override_file(override_path: str) -> Dict[str, Any]:
    """
    Load configuration overrides from a JSON file.
    
    Args:
        override_path: Path to the override file
        
    Returns:
        Dictionary of overrides or an empty dict if loading fails
    """
    try:
        with open(override_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading override file: {e}")
        return {}


def get_config_path_value(config: Any, path: str, default: Any = None) -> Any:
    """
    Get a value from the configuration using a dot-notation path.
    
    Args:
        config: The configuration object
        path: Dot-notation path (e.g., "app.debug")
        default: Default value to return if the path doesn't exist
        
    Returns:
        The value at the specified path or the default
    """
    parts = path.split('.')
    
    # Handle attribute-style access (Python modules)
    if hasattr(config, parts[0].upper() + '_' + parts[1].upper()):
        # Python convention: APP_DEBUG for app.debug
        return getattr(config, parts[0].upper() + '_' + parts[1].upper(), default)
    
    # Handle dict-style access (JSON, TypeScript)
    current = config
    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return default
    
    return current


def set_config_path_value(config: Any, path: str, value: Any) -> bool:
    """
    Set a value in the configuration using a dot-notation path.
    
    Args:
        config: The configuration object
        path: Dot-notation path (e.g., "app.debug")
        value: The value to set
        
    Returns:
        True if successful, False otherwise
    """
    parts = path.split('.')
    
    # Handle attribute-style access (Python modules)
    if hasattr(config, parts[0].upper() + '_' + parts[1].upper()):
        setattr(config, parts[0].upper() + '_' + parts[1].upper(), value)
        return True
    
    # Handle dict-style access (JSON, TypeScript)
    if not isinstance(config, dict):
        return False
        
    current = config
    for i, part in enumerate(parts):
        if i == len(parts) - 1:
            # Last part, set the value
            current[part] = value
        else:
            # Not the last part, navigate down
            if part not in current:
                current[part] = {}
            current = current[part]
    
    return True


def convert_value(value: str, target_type: str) -> Any:
    """
    Convert a string value to the target type.
    
    Args:
        value: String value to convert
        target_type: Target type ("string", "boolean", "integer", "number", "array", "object")
        
    Returns:
        The converted value
    """
    if target_type == "string":
        return value
    elif target_type == "boolean":
        return value.lower() in ('true', 'yes', '1', 'y', 'on')
    elif target_type == "integer":
        return int(value)
    elif target_type == "number":
        return float(value)
    elif target_type == "array":
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            # Try comma-separated list
            return [item.strip() for item in value.split(',')]
    elif target_type == "object":
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            # Try key-value pairs
            result = {}
            for pair in value.split(','):
                if ':' in pair:
                    k, v = pair.split(':', 1)
                    result[k.strip()] = v.strip()
            return result
    
    # Default: return as-is
    return value


def apply_override_file(config: Any, override_path: str, schema: Dict[str, Any]) -> List[str]:
    """
    Apply overrides from a file to the configuration.
    
    Args:
        config: The configuration object
        override_path: Path to the override file
        schema: The configuration schema
        
    Returns:
        List of messages about the applied overrides
    """
    messages = []
    overrides = load_override_file(override_path)
    
    if not overrides:
        messages.append(f"No overrides found in {override_path}")
        return messages
    
    # Apply each override
    for path, value in overrides.items():
        parts = path.split('.')
        if len(parts) != 2:
            messages.append(f"Invalid path: {path} (must be in format 'section.property')")
            continue
        
        section, prop = parts
        
        # Validate against schema
        if section not in schema.get("properties", {}) or prop not in schema.get("properties", {}).get(section, {}).get("properties", {}):
            messages.append(f"Unknown configuration path: {path}")
            continue
        
        # Get the type from schema
        prop_type = schema.get("properties", {}).get(section, {}).get("properties", {}).get(prop, {}).get("type", "string")
        
        # Convert the value if needed
        if isinstance(value, str) and prop_type != "string":
            try:
                value = convert_value(value, prop_type)
            except Exception as e:
                messages.append(f"Error converting value for {path}: {e}")
                continue
        
        # Apply the override
        if set_config_path_value(config, path, value):
            messages.append(f"Applied override: {path} = {value}")
        else:
            messages.append(f"Failed to apply override: {path} = {value}")
    
    return messages


def apply_command_line_overrides(config: Any, args: Dict[str, Any], schema: Dict[str, Any]) -> List[str]:
    """
    Apply overrides from command-line arguments to the configuration.
    
    Args:
        config: The configuration object
        args: Dictionary of command-line arguments
        schema: The configuration schema
        
    Returns:
        List of messages about the applied overrides
    """
    messages = []
    
    # Apply each override
    for key, value in args.items():
        if not value or key == "config":
            continue
            
        # Normalize key (convert --app-debug to app.debug)
        path = key.lstrip('-').replace('-', '.')
        
        # Validate against schema
        parts = path.split('.')
        if len(parts) != 2:
            messages.append(f"Invalid path: {path} (must be in format 'section.property')")
            continue
        
        section, prop = parts
        
        if section not in schema.get("properties", {}) or prop not in schema.get("properties", {}).get(section, {}).get("properties", {}):
            messages.append(f"Unknown configuration path: {path}")
            continue
        
        # Get the type from schema
        prop_type = schema.get("properties", {}).get(section, {}).get("properties", {}).get(prop, {}).get("type", "string")
        
        # Convert the value if needed
        if isinstance(value, str) and prop_type != "string":
            try:
                value = convert_value(value, prop_type)
            except Exception as e:
                messages.append(f"Error converting value for {path}: {e}")
                continue
        
        # Apply the override
        if set_config_path_value(config, path, value):
            messages.append(f"Applied command-line override: {path} = {value}")
        else:
            messages.append(f"Failed to apply command-line override: {path} = {value}")
    
    return messages


def apply_environment_overrides(config: Any, schema: Dict[str, Any]) -> List[str]:
    """
    Apply overrides from environment variables to the configuration.
    
    Args:
        config: The configuration object
        schema: The configuration schema
        
    Returns:
        List of messages about the applied overrides
    """
    messages = []
    
    # Get all environment variables
    for key, value in os.environ.items():
        # Skip non-override variables
        if not key.startswith('OVERRIDE_'):
            continue
        
        # Extract the path (convert OVERRIDE_APP_DEBUG to app.debug)
        parts = key[9:].lower().split('_')
        if len(parts) < 2:
            continue
            
        section = parts[0]
        prop = '_'.join(parts[1:])
        path = f"{section}.{prop}"
        
        # Validate against schema
        if section not in schema.get("properties", {}) or prop not in schema.get("properties", {}).get(section, {}).get("properties", {}):
            messages.append(f"Unknown configuration path from environment: {path}")
            continue
        
        # Get the type from schema
        prop_type = schema.get("properties", {}).get(section, {}).get("properties", {}).get(prop, {}).get("type", "string")
        
        # Convert the value
        try:
            typed_value = convert_value(value, prop_type)
        except Exception as e:
            messages.append(f"Error converting value for {path}: {e}")
            continue
        
        # Apply the override
        if set_config_path_value(config, path, typed_value):
            messages.append(f"Applied environment override: {path} = {typed_value}")
        else:
            messages.append(f"Failed to apply environment override: {path} = {typed_value}")
    
    return messages


def build_override_argument_parser(schema: Dict[str, Any]) -> argparse.ArgumentParser:
    """
    Build an argument parser with options based on the schema.
    
    Args:
        schema: The configuration schema
        
    Returns:
        An ArgumentParser instance
    """
    parser = argparse.ArgumentParser(description="Application with configuration overrides")
    parser.add_argument("--config", help="Path to config override file")
    
    # Add arguments for each configuration property
    for section_name, section in schema.get("properties", {}).items():
        for prop_name, prop in section.get("properties", {}).items():
            # Skip complex objects, they can't be easily overridden from command line
            if prop.get("type") == "object":
                continue
                
            arg_name = f"--{section_name}-{prop_name}"
            help_text = prop.get("description", "")
            
            if "default" in prop:
                help_text += f" (default: {prop['default']})"
            
            parser.add_argument(arg_name, help=help_text)
    
    return parser


def load_python_config(config_path: str) -> Optional[Any]:
    """
    Load a Python configuration module.
    
    Args:
        config_path: Path to the Python configuration file
        
    Returns:
        The loaded module or None if loading fails
    """
    try:
        spec = importlib.util.spec_from_file_location("config_module", config_path)
        if spec is None or spec.loader is None:
            print(f"Failed to load Python config: Invalid module spec")
            return None
            
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        print(f"Failed to load Python config: {e}")
        return None


def generate_override_file(schema_path: str, output_path: str, env: str = None) -> bool:
    """
    Generate a template override file based on the schema.
    
    Args:
        schema_path: Path to the schema file
        output_path: Path to write the override file
        env: Optional environment to include those defaults
        
    Returns:
        True if successful, False otherwise
    """
    schema = load_schema(schema_path)
    if not schema:
        return False
    
    # Build the overrides template
    overrides = {}
    
    # If environment is specified, include those values
    env_overrides = {}
    if env and env in schema.get("environmentMapping", {}):
        env_overrides = schema.get("environmentMapping", {}).get(env, {})
    
    # Add commented examples for each property
    for section_name, section in schema.get("properties", {}).items():
        overrides[section_name] = {}
        
        for prop_name, prop in section.get("properties", {}).items():
            # Skip complex objects
            if prop.get("type") == "object":
                continue
                
            path = f"{section_name}.{prop_name}"
            
            # Use environment override if available
            if path in env_overrides:
                overrides[section_name][prop_name] = env_overrides[path]
            # Otherwise use the default value if available
            elif "default" in prop:
                overrides[section_name][prop_name] = prop["default"]
            else:
                # Use a placeholder based on type
                if prop.get("type") == "string":
                    overrides[section_name][prop_name] = ""
                elif prop.get("type") == "boolean":
                    overrides[section_name][prop_name] = False
                elif prop.get("type") == "integer":
                    overrides[section_name][prop_name] = 0
                elif prop.get("type") == "number":
                    overrides[section_name][prop_name] = 0.0
                elif prop.get("type") == "array":
                    overrides[section_name][prop_name] = []
                else:
                    overrides[section_name][prop_name] = None
    
    # Write the override file
    try:
        with open(output_path, 'w') as f:
            json.dump(overrides, f, indent=2)
        return True
    except Exception as e:
        print(f"Error writing override file: {e}")
        return False


def print_config_info(config: Any, schema: Dict[str, Any]) -> None:
    """
    Print information about the configuration.
    
    Args:
        config: The configuration object
        schema: The configuration schema
    """
    # Determine the type of configuration
    if hasattr(config, "BaseConfig"):
        config_type = "Python module"
        # Check if it has an active_config
        if hasattr(config, "active_config"):
            config_class = config.active_config.__class__.__name__
        else:
            config_class = "Unknown"
    elif isinstance(config, dict):
        config_type = "Dictionary"
        config_class = "N/A"
    else:
        config_type = str(type(config))
        config_class = "N/A"
    
    print(f"Configuration Type: {config_type}")
    print(f"Active Configuration: {config_class}")
    print(f"Schema Sections: {', '.join(schema.get('properties', {}).keys())}")
    
    # Print some example values
    print("\nSample Configuration Values:")
    for section_name, section in schema.get("properties", {}).items():
        for prop_name, prop in section.get("properties", {}).items():
            # Skip complex objects
            if prop.get("type") == "object":
                continue
                
            path = f"{section_name}.{prop_name}"
            value = get_config_path_value(config, path, "<not set>")
            print(f"  {path}: {value}")
            
            # Only show a few values
            if len(section.get("properties", {})) > 5:
                remaining = len(section.get("properties", {})) - 5
                print(f"  ... and {remaining} more properties in {section_name}")
                break


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Configuration Override System")
    parser.add_argument("--schema", default="schema.json", help="Path to schema file")
    parser.add_argument("--config", default="../config.py", help="Path to configuration file")
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Generate subcommand
    generate_parser = subparsers.add_parser("generate", help="Generate an override file")
    generate_parser.add_argument("--output", required=True, help="Path to output file")
    generate_parser.add_argument("--env", choices=["development", "testing", "staging", "production"], help="Environment to include")
    
    # Apply subcommand
    apply_parser = subparsers.add_parser("apply", help="Apply overrides to configuration")
    apply_parser.add_argument("--override", required=True, help="Path to override file")
    apply_parser.add_argument("--output", help="Path to output JSON (defaults to stdout)")
    
    # Info subcommand
    info_parser = subparsers.add_parser("info", help="Show information about the configuration")
    
    args = parser.parse_args()
    
    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    schema_path = os.path.abspath(os.path.join(script_dir, args.schema))
    config_path = os.path.abspath(os.path.join(script_dir, args.config))
    
    # Load schema
    schema = load_schema(schema_path)
    if not schema:
        return 1
    
    if args.command == "generate":
        output_path = os.path.abspath(args.output)
        if generate_override_file(schema_path, output_path, args.env):
            print(f"Generated override file: {output_path}")
            return 0
        else:
            return 1
    elif args.command == "apply" or args.command == "info":
        # Load configuration
        config = None
        if config_path.endswith('.py'):
            config = load_python_config(config_path)
        elif config_path.endswith('.json'):
            try:
                with open(config_path, 'r') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"Failed to load JSON config: {e}")
        
        if config is None:
            return 1
        
        if args.command == "info":
            print_config_info(config, schema)
            return 0
        else:  # apply
            override_path = os.path.abspath(args.override)
            messages = apply_override_file(config, override_path, schema)
            
            for message in messages:
                print(message)
            
            # Output the result
            if args.output:
                output_path = os.path.abspath(args.output)
                try:
                    # Convert to dictionary if it's a module
                    if hasattr(config, "active_config") and hasattr(config.active_config, "as_dict"):
                        config_dict = config.active_config.as_dict()
                    elif not isinstance(config, dict):
                        config_dict = {}
                        for section_name, section in schema.get("properties", {}).items():
                            config_dict[section_name] = {}
                            for prop_name in section.get("properties", {}).keys():
                                path = f"{section_name}.{prop_name}"
                                value = get_config_path_value(config, path)
                                if value is not None:
                                    config_dict[section_name][prop_name] = value
                    else:
                        config_dict = config
                    
                    with open(output_path, 'w') as f:
                        json.dump(config_dict, f, indent=2)
                    print(f"Wrote result to {output_path}")
                except Exception as e:
                    print(f"Failed to write output: {e}")
                    return 1
            
            return 0
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
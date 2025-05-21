#!/usr/bin/env python
"""
Configuration Validation System

This script validates the configuration files against the schema and checks for consistency
across different environments and between backend and frontend.
"""

import os
import sys
import json
import argparse
import importlib.util
from pathlib import Path
from typing import Dict, Any, List, Optional, Set, Tuple

# ANSI colors for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def colorize(text: str, color_code: str) -> str:
    """Add color to terminal output."""
    return f"{color_code}{text}{Colors.END}"


def validate_schema(schema_path: str) -> Tuple[bool, Dict[str, Any], List[str]]:
    """
    Validate the schema itself for correctness and completeness.
    
    Args:
        schema_path: Path to the schema file
        
    Returns:
        Tuple of (is_valid, schema, error_messages)
    """
    errors = []
    schema = None
    
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
    except Exception as e:
        errors.append(f"Failed to load schema: {str(e)}")
        return False, {}, errors
    
    # Check for required top-level properties
    if "properties" not in schema:
        errors.append("Schema missing 'properties' section")
        return False, schema, errors
    
    if "environmentMapping" not in schema:
        errors.append("Schema missing 'environmentMapping' section")
    
    # Validate all properties have proper type definitions
    for section_name, section in schema.get("properties", {}).items():
        if "type" not in section:
            errors.append(f"Section '{section_name}' missing 'type' property")
            continue
            
        if section["type"] != "object":
            errors.append(f"Section '{section_name}' has invalid type: {section['type']}. Must be 'object'")
            continue
            
        if "properties" not in section:
            errors.append(f"Section '{section_name}' missing 'properties' object")
            continue
            
        # Check each property in the section
        for prop_name, prop in section["properties"].items():
            if "type" not in prop:
                errors.append(f"Property '{section_name}.{prop_name}' missing 'type' property")
            
            # If property is required, ensure it has a default or is in environment overrides
            if "required" in section and prop_name in section["required"]:
                if "default" not in prop:
                    # Check if the property is specified in environment overrides
                    path = f"{section_name}.{prop_name}"
                    has_env_override = False
                    
                    for env, overrides in schema.get("environmentMapping", {}).items():
                        if path in overrides:
                            has_env_override = True
                            break
                    
                    if not has_env_override:
                        errors.append(f"Required property '{path}' has no default value and no environment overrides")
    
    # Check environment mappings point to valid properties
    for env_name, env_mapping in schema.get("environmentMapping", {}).items():
        for path, value in env_mapping.items():
            parts = path.split(".")
            if len(parts) != 2:
                errors.append(f"Invalid path '{path}' in environment mapping for '{env_name}'. Must be in format 'section.property'")
                continue
                
            section_name, prop_name = parts
            
            if section_name not in schema.get("properties", {}):
                errors.append(f"Path '{path}' in environment mapping for '{env_name}' references non-existent section '{section_name}'")
                continue
                
            if prop_name not in schema.get("properties", {}).get(section_name, {}).get("properties", {}):
                errors.append(f"Path '{path}' in environment mapping for '{env_name}' references non-existent property '{prop_name}'")
                continue
                
            # Check value type matches property type
            prop_type = schema.get("properties", {}).get(section_name, {}).get("properties", {}).get(prop_name, {}).get("type")
            
            if prop_type == "string" and not isinstance(value, str):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be a string, but got {type(value).__name__}")
            elif prop_type == "number" and not isinstance(value, (int, float)):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be a number, but got {type(value).__name__}")
            elif prop_type == "integer" and not isinstance(value, int):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be an integer, but got {type(value).__name__}")
            elif prop_type == "boolean" and not isinstance(value, bool):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be a boolean, but got {type(value).__name__}")
            elif prop_type == "array" and not isinstance(value, list):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be an array, but got {type(value).__name__}")
            elif prop_type == "object" and not isinstance(value, dict):
                errors.append(f"Value for '{path}' in environment mapping for '{env_name}' should be an object, but got {type(value).__name__}")
    
    is_valid = len(errors) == 0
    return is_valid, schema, errors


def check_python_config(python_config_path: str, schema: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Check if the Python configuration matches the schema.
    
    Args:
        python_config_path: Path to the Python configuration file
        schema: The schema to validate against
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    warnings = []
    
    if not os.path.exists(python_config_path):
        errors.append(f"Python config file not found: {python_config_path}")
        return False, errors
    
    # Load the Python configuration as a module
    try:
        spec = importlib.util.spec_from_file_location("config", python_config_path)
        if spec is None or spec.loader is None:
            errors.append(f"Failed to load Python config: Invalid module spec")
            return False, errors
            
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception as e:
        errors.append(f"Failed to load Python config: {str(e)}")
        return False, errors
    
    # Get the BaseConfig class
    if not hasattr(module, "BaseConfig"):
        errors.append("Python config missing BaseConfig class")
        return False, errors
    
    base_config = module.BaseConfig
    
    # Check if all properties from the schema are present in the BaseConfig
    for section_name, section in schema.get("properties", {}).items():
        for prop_name, prop in section.get("properties", {}).items():
            if prop.get("type") == "object":
                # Skip complex objects, they might be handled differently
                continue
                
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            if not hasattr(base_config, python_name):
                errors.append(f"Property '{section_name}.{prop_name}' missing from Python config as '{python_name}'")
    
    # Check if environment-specific configurations are properly defined
    for env_name, env_mapping in schema.get("environmentMapping", {}).items():
        class_name = f"{env_name.capitalize()}Config"
        
        if not hasattr(module, class_name):
            errors.append(f"Python config missing environment class '{class_name}'")
            continue
        
        env_class = getattr(module, class_name)
        
        for path, value in env_mapping.items():
            parts = path.split(".")
            if len(parts) != 2:
                continue  # Already validated in schema check
                
            section_name, prop_name = parts
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            if not hasattr(env_class, python_name):
                warnings.append(f"Environment override '{path}' for '{env_name}' missing from Python config class '{class_name}'")
    
    # Check for module-level exports
    for section_name, section in schema.get("properties", {}).items():
        for prop_name, prop in section.get("properties", {}).items():
            if prop.get("type") == "object":
                # Skip complex objects, they might be handled differently
                continue
                
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            if not hasattr(module, python_name):
                warnings.append(f"Module-level export for '{section_name}.{prop_name}' missing from Python config as '{python_name}'")
    
    # Output warnings
    if warnings:
        print("\nWarnings:")
        for warning in warnings:
            print(colorize(f"  - {warning}", Colors.WARNING))
    
    is_valid = len(errors) == 0
    return is_valid, errors


def check_typescript_config(typescript_config_path: str, schema: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Check if the TypeScript configuration matches the schema.
    
    Args:
        typescript_config_path: Path to the TypeScript configuration file
        schema: The schema to validate against
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    if not os.path.exists(typescript_config_path):
        errors.append(f"TypeScript config file not found: {typescript_config_path}")
        return False, errors
    
    # Load the TypeScript configuration content
    try:
        with open(typescript_config_path, 'r') as f:
            content = f.read()
    except Exception as e:
        errors.append(f"Failed to read TypeScript config: {str(e)}")
        return False, errors
    
    # Check for interfaces
    for section_name, section in schema.get("properties", {}).items():
        # Convert section name to PascalCase
        pascal_section = ''.join(word.capitalize() for word in section_name.split('_'))
        interface_name = f"{pascal_section}Config"
        
        if f"interface {interface_name}" not in content and f"export interface {interface_name}" not in content:
            errors.append(f"TypeScript config missing interface '{interface_name}' for section '{section_name}'")
            continue
        
        # Check for all properties
        for prop_name, prop in section.get("properties", {}).items():
            prop_pattern = f"{prop_name}:"
            prop_nullable_pattern = f"{prop_name}?:"
            
            if prop_pattern not in content and prop_nullable_pattern not in content:
                errors.append(f"TypeScript config missing property '{prop_name}' in interface '{interface_name}'")
    
    # Check for environment configurations
    for env_name in schema.get("environmentMapping", {}).keys():
        var_name = f"{env_name}Config"
        
        if f"const {var_name}" not in content and f"export const {var_name}" not in content:
            errors.append(f"TypeScript config missing environment configuration '{var_name}'")
    
    # Check for the root interface and config getter
    if "interface AppConfig" not in content and "export interface AppConfig" not in content:
        errors.append("TypeScript config missing root AppConfig interface")
    
    if "function getConfig" not in content and "export function getConfig" not in content:
        errors.append("TypeScript config missing getConfig function")
    
    if "const config =" not in content and "export const config =" not in content:
        errors.append("TypeScript config missing exported config instance")
    
    is_valid = len(errors) == 0
    return is_valid, errors


def check_for_inconsistencies(schema: Dict[str, Any], python_config_path: str, typescript_config_path: str) -> List[str]:
    """
    Check for inconsistencies between Python and TypeScript configurations.
    
    Args:
        schema: The schema to validate against
        python_config_path: Path to the Python configuration file
        typescript_config_path: Path to the TypeScript configuration file
        
    Returns:
        List of inconsistency warnings
    """
    warnings = []
    
    # Skip if either file doesn't exist
    if not os.path.exists(python_config_path) or not os.path.exists(typescript_config_path):
        warnings.append("Skipping cross-language consistency check due to missing files")
        return warnings
    
    # Load Python config
    try:
        spec = importlib.util.spec_from_file_location("config", python_config_path)
        if spec is None or spec.loader is None:
            warnings.append("Failed to load Python config for consistency check")
            return warnings
            
        py_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(py_module)
    except Exception:
        warnings.append("Failed to load Python config for consistency check")
        return warnings
    
    # Load TypeScript config content
    try:
        with open(typescript_config_path, 'r') as f:
            ts_content = f.read()
    except Exception:
        warnings.append("Failed to read TypeScript config for consistency check")
        return warnings
    
    # Check for mismatches in default values
    for section_name, section in schema.get("properties", {}).items():
        pascal_section = ''.join(word.capitalize() for word in section_name.split('_'))
        for prop_name, prop in section.get("properties", {}).items():
            if "default" in prop and prop.get("type") != "object":
                # Check Python default
                python_name = f"{section_name.upper()}_{prop_name.upper()}"
                py_default = None
                
                if hasattr(py_module.BaseConfig, python_name):
                    py_default = getattr(py_module.BaseConfig, python_name)
                
                # Check TypeScript default
                ts_var_name = f"default{pascal_section}Config"
                ts_default_pattern = f"{prop_name}: "
                
                if py_default is not None:
                    if isinstance(py_default, str):
                        ts_check = f'{prop_name}: "{py_default}"'
                    elif isinstance(py_default, bool):
                        ts_check = f'{prop_name}: {str(py_default).lower()}'
                    else:
                        ts_check = f'{prop_name}: {py_default}'
                    
                    if ts_check not in ts_content and ts_var_name in ts_content:
                        warnings.append(f"Default value mismatch for '{section_name}.{prop_name}': Python={py_default}")
    
    # Check for environment-specific overrides
    for env_name, env_mapping in schema.get("environmentMapping", {}).items():
        py_class_name = f"{env_name.capitalize()}Config"
        ts_var_name = f"{env_name}Config"
        
        if not hasattr(py_module, py_class_name) or ts_var_name not in ts_content:
            continue
        
        py_class = getattr(py_module, py_class_name)
        
        for path, value in env_mapping.items():
            parts = path.split(".")
            if len(parts) != 2:
                continue
                
            section_name, prop_name = parts
            python_name = f"{section_name.upper()}_{prop_name.upper()}"
            
            if hasattr(py_class, python_name):
                py_value = getattr(py_class, python_name)
                
                # Check if TypeScript has matching value
                if isinstance(py_value, str):
                    ts_check = f'{prop_name}: "{py_value}"'
                elif isinstance(py_value, bool):
                    ts_check = f'{prop_name}: {str(py_value).lower()}'
                else:
                    ts_check = f'{prop_name}: {py_value}'
                
                if ts_check not in ts_content and f"{section_name}: {{" in ts_content:
                    warnings.append(f"Environment override mismatch for '{path}' in '{env_name}': Python={py_value}")
    
    return warnings


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Validate configuration against schema")
    parser.add_argument("--schema", default="schema.json", help="Path to schema file")
    parser.add_argument("--python", default="../config.py", help="Path to Python config file")
    parser.add_argument("--typescript", default="../frontend/src/config/config.ts", help="Path to TypeScript config file")
    parser.add_argument("--check-only", action="store_true", help="Only check, don't try to fix issues")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show verbose output")
    
    args = parser.parse_args()
    
    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    schema_path = os.path.abspath(os.path.join(script_dir, args.schema))
    python_path = os.path.abspath(os.path.join(script_dir, args.python))
    typescript_path = os.path.abspath(os.path.join(script_dir, args.typescript))
    
    print(colorize("Configuration Validation", Colors.HEADER + Colors.BOLD))
    print(f"Schema: {schema_path}")
    print(f"Python: {python_path}")
    print(f"TypeScript: {typescript_path}")
    print("")
    
    # Step 1: Validate the schema
    print(colorize("Step 1: Validating Schema", Colors.BOLD))
    schema_valid, schema, schema_errors = validate_schema(schema_path)
    
    if schema_valid:
        print(colorize("  ✓ Schema is valid", Colors.GREEN))
    else:
        print(colorize("  ✗ Schema validation failed", Colors.FAIL))
        for error in schema_errors:
            print(colorize(f"    - {error}", Colors.FAIL))
        return 1
    
    # Step 2: Check Python configuration
    print(colorize("\nStep 2: Checking Python Configuration", Colors.BOLD))
    python_valid, python_errors = check_python_config(python_path, schema)
    
    if python_valid:
        print(colorize("  ✓ Python configuration is valid", Colors.GREEN))
    else:
        print(colorize("  ✗ Python configuration has errors", Colors.FAIL))
        for error in python_errors:
            print(colorize(f"    - {error}", Colors.FAIL))
    
    # Step 3: Check TypeScript configuration
    print(colorize("\nStep 3: Checking TypeScript Configuration", Colors.BOLD))
    typescript_valid, typescript_errors = check_typescript_config(typescript_path, schema)
    
    if typescript_valid:
        print(colorize("  ✓ TypeScript configuration is valid", Colors.GREEN))
    else:
        print(colorize("  ✗ TypeScript configuration has errors", Colors.FAIL))
        for error in typescript_errors:
            print(colorize(f"    - {error}", Colors.FAIL))
    
    # Step 4: Check for inconsistencies
    print(colorize("\nStep 4: Checking Cross-Language Consistency", Colors.BOLD))
    inconsistencies = check_for_inconsistencies(schema, python_path, typescript_path)
    
    if not inconsistencies:
        print(colorize("  ✓ No inconsistencies found", Colors.GREEN))
    else:
        print(colorize("  ! Inconsistencies found", Colors.WARNING))
        for warning in inconsistencies:
            print(colorize(f"    - {warning}", Colors.WARNING))
    
    # Summary
    print(colorize("\nValidation Summary:", Colors.BOLD))
    if schema_valid and python_valid and typescript_valid and not inconsistencies:
        print(colorize("  ✓ All checks passed!", Colors.GREEN + Colors.BOLD))
        return 0
    else:
        issues = []
        if not schema_valid:
            issues.append("Schema issues")
        if not python_valid:
            issues.append("Python configuration issues")
        if not typescript_valid:
            issues.append("TypeScript configuration issues")
        if inconsistencies:
            issues.append("Cross-language inconsistencies")
            
        print(colorize(f"  ✗ Validation failed: {', '.join(issues)}", Colors.FAIL))
        
        # If not check-only, suggest fixes
        if not args.check_only:
            print(colorize("\nSuggested Actions:", Colors.BOLD))
            print("  1. Fix any schema issues first")
            print("  2. Regenerate configuration files using 'generate_configs.py'")
            print("  3. Run this validation script again to verify fixes")
        
        return 1


if __name__ == "__main__":
    sys.exit(main())
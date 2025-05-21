#!/usr/bin/env python
"""
Configuration Generator

This script generates configuration files for different languages from a unified schema.
Currently supports Python and TypeScript.
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path


def validate_schema(schema_path):
    """Validate the schema for correctness and completeness."""
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
        
        # Check for required top-level sections
        required_sections = ["app", "api", "database", "auth"]
        for section in required_sections:
            if section not in schema.get("properties", {}):
                print(f"Warning: Required section '{section}' is missing from schema")
        
        # Check for environmentMapping
        if "environmentMapping" not in schema:
            print("Warning: Missing 'environmentMapping' section in schema")
        else:
            # Check for all required environments
            required_envs = ["development", "testing", "staging", "production"]
            for env in required_envs:
                if env not in schema["environmentMapping"]:
                    print(f"Warning: Environment '{env}' is missing from environmentMapping")
        
        print("Schema validation completed successfully")
        return True
    except Exception as e:
        print(f"Error validating schema: {str(e)}")
        return False


def generate_python_config(schema_path, python_path):
    """Generate Python configuration file."""
    script_path = os.path.join(os.path.dirname(schema_path), "generate_python_config.py")
    
    # Run the Python generator script
    try:
        result = subprocess.run(
            [sys.executable, script_path, "--schema", schema_path, "--output", python_path],
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error generating Python config: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False


def generate_typescript_config(schema_path, typescript_path):
    """Generate TypeScript configuration file."""
    script_path = os.path.join(os.path.dirname(schema_path), "generate_typescript_config.js")
    
    # Ensure script is executable
    os.chmod(script_path, 0o755)
    
    # Run the TypeScript generator script
    try:
        # Check if Node.js is available
        try:
            subprocess.run(["node", "--version"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Warning: Node.js not found, skipping TypeScript config generation")
            return False
        
        # Run the generator
        result = subprocess.run(
            ["node", script_path],
            check=True,
            capture_output=True,
            text=True,
            env={**os.environ, "SCHEMA_PATH": schema_path, "OUTPUT_PATH": typescript_path}
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error generating TypeScript config: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Generate configuration files from unified schema")
    parser.add_argument("--schema", default="schema.json", help="Path to schema.json")
    parser.add_argument("--python", default="../config.py.new", help="Output path for Python config")
    parser.add_argument("--typescript", default="../frontend/src/config/config.ts", help="Output path for TypeScript config")
    parser.add_argument("--validate-only", action="store_true", help="Only validate the schema, don't generate files")
    
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
    
    # Validate schema
    if not validate_schema(schema_path):
        print("Schema validation failed")
        return 1
    
    if args.validate_only:
        print("Schema validated successfully. Skipping generation as --validate-only was specified.")
        return 0
    
    # Generate configs
    success = True
    
    # Python
    print(f"Generating Python config at {python_path}...")
    if not generate_python_config(schema_path, python_path):
        print("Failed to generate Python configuration")
        success = False
    
    # TypeScript
    print(f"Generating TypeScript config at {typescript_path}...")
    if not generate_typescript_config(schema_path, typescript_path):
        print("Failed to generate TypeScript configuration")
        success = False
    
    if success:
        print("\nAll configuration files generated successfully")
        return 0
    else:
        print("\nOne or more configuration files failed to generate")
        return 1


if __name__ == "__main__":
    sys.exit(main())
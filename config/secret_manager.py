#!/usr/bin/env python
"""
Secret Manager for Configuration

This script provides tools for managing secrets in the configuration system.
It supports generating, encrypting, and storing secrets in various backends:
- .env files (development)
- Environment variables (any environment)
- Docker secrets (production)
- Secure credential storage (where available)
"""

import os
import sys
import json
import base64
import getpass
import secrets
import argparse
import hashlib
from typing import Dict, Any, Optional
from pathlib import Path

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


def generate_random_secret(length: int = 32) -> str:
    """Generate a cryptographically secure random string."""
    return secrets.token_urlsafe(length)


def hash_secret(secret: str, salt: Optional[str] = None) -> tuple:
    """
    Hash a secret with a salt using a secure algorithm.
    
    Args:
        secret: The secret to hash
        salt: Optional salt, will be generated if not provided
        
    Returns:
        Tuple of (hash, salt)
    """
    if salt is None:
        salt = secrets.token_hex(16)
    
    # Use a secure hashing algorithm (SHA-256)
    key = hashlib.pbkdf2_hmac(
        'sha256',
        secret.encode('utf-8'),
        salt.encode('utf-8'),
        100000  # Number of iterations
    )
    
    key_hex = key.hex()
    return key_hex, salt


def identify_secrets_from_schema(schema_path: str) -> Dict[str, Any]:
    """
    Identify secret properties in the schema.
    
    Args:
        schema_path: Path to the schema file
        
    Returns:
        Dictionary of secret properties with their paths and metadata
    """
    secrets_dict = {}
    
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
    except Exception as e:
        print(colorize(f"Error loading schema: {str(e)}", Colors.FAIL))
        return {}
    
    # Look for properties marked as secrets (format: "password" or has "secret" in the name)
    for section_name, section in schema.get("properties", {}).items():
        for prop_name, prop in section.get("properties", {}).items():
            # Check if it's a secret
            is_secret = False
            
            # Check property format
            if prop.get("format") == "password":
                is_secret = True
            
            # Check property name
            if any(secret_term in prop_name.lower() for secret_term in ["secret", "password", "key", "token", "credential"]):
                is_secret = True
            
            if is_secret:
                env_var = f"{section_name.upper()}_{prop_name.upper()}"
                
                # Skip if it has a default value that is non-sensitive
                if "default" in prop and "example" in prop:
                    # If default matches example, it's likely a demo/development value
                    if prop["default"] == prop["example"]:
                        continue
                
                secrets_dict[f"{section_name}.{prop_name}"] = {
                    "env_var": env_var,
                    "description": prop.get("description", ""),
                    "required": prop_name in section.get("required", []),
                    "type": prop.get("type", "string")
                }
    
    return secrets_dict


def read_current_secrets_from_env() -> Dict[str, str]:
    """
    Read current secret values from environment variables.
    
    Returns:
        Dictionary of secret values
    """
    env_secrets = {}
    
    for key, value in os.environ.items():
        # Look for keys that might contain secrets
        if any(secret_term in key.lower() for secret_term in ["secret", "password", "key", "token", "credential"]):
            env_secrets[key] = value
    
    return env_secrets


def write_secrets_to_env_file(env_file: str, secrets_dict: Dict[str, str], overwrite: bool = False) -> bool:
    """
    Write secrets to a .env file.
    
    Args:
        env_file: Path to the .env file
        secrets_dict: Dictionary of secrets to write
        overwrite: Whether to overwrite existing values
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(env_file)), exist_ok=True)
        
        # Read existing .env file if it exists
        existing_env = {}
        if os.path.exists(env_file):
            with open(env_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    if '=' in line:
                        key, value = line.split('=', 1)
                        existing_env[key.strip()] = value.strip()
        
        # Merge new secrets with existing
        merged_env = {**existing_env}
        for key, value in secrets_dict.items():
            if key in merged_env and not overwrite:
                print(f"Skipping existing key: {key}")
            else:
                merged_env[key] = value
        
        # Write to .env file
        with open(env_file, 'w') as f:
            for key, value in sorted(merged_env.items()):
                f.write(f"{key}={value}\n")
        
        return True
    except Exception as e:
        print(colorize(f"Error writing .env file: {str(e)}", Colors.FAIL))
        return False


def create_docker_secrets(secrets_dict: Dict[str, str], output_dir: str) -> bool:
    """
    Create Docker secret files for use with Docker Compose or Kubernetes.
    
    Args:
        secrets_dict: Dictionary of secrets to create
        output_dir: Directory to write secret files
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create the output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Write each secret to a file
        for key, value in secrets_dict.items():
            secret_file = os.path.join(output_dir, key)
            with open(secret_file, 'w') as f:
                f.write(value)
            
            # Set restrictive permissions
            os.chmod(secret_file, 0o600)
        
        return True
    except Exception as e:
        print(colorize(f"Error creating Docker secrets: {str(e)}", Colors.FAIL))
        return False


def generate_secrets_for_schema(schema_path: str, env_file: str, docker_secrets_dir: Optional[str] = None) -> bool:
    """
    Generate secrets for all properties marked as secrets in the schema.
    
    Args:
        schema_path: Path to the schema file
        env_file: Path to the .env file
        docker_secrets_dir: Optional directory for Docker secrets
        
    Returns:
        True if successful, False otherwise
    """
    # Identify secrets from schema
    secrets = identify_secrets_from_schema(schema_path)
    
    if not secrets:
        print("No secrets identified in schema")
        return True
    
    print(f"Identified {len(secrets)} secrets in schema:")
    for path, info in secrets.items():
        required = "Required" if info.get("required") else "Optional"
        print(f"  - {path} ({required}): {info.get('description', '')}")
    
    # Read existing secrets from environment
    env_secrets = read_current_secrets_from_env()
    
    # Generate new secrets
    generated_secrets = {}
    for path, info in secrets.items():
        env_var = info["env_var"]
        
        # Check if already exists in environment
        if env_var in env_secrets:
            print(f"Using existing value for {env_var}")
            value = env_secrets[env_var]
        else:
            # Prompt for sensitive secrets, generate others automatically
            if any(sensitive in path.lower() for sensitive in ["password", "secret"]):
                prompt = f"Enter value for {path} ({info.get('description', '')}): "
                value = getpass.getpass(prompt)
                if not value and info.get("required"):
                    print(colorize(f"Error: {path} is required", Colors.FAIL))
                    return False
            else:
                value = generate_random_secret()
                print(f"Generated random value for {path}")
        
        if value:
            generated_secrets[env_var] = value
    
    # Write to .env file
    if env_file:
        print(f"\nWriting secrets to {env_file}")
        if not write_secrets_to_env_file(env_file, generated_secrets, overwrite=False):
            return False
    
    # Write Docker secrets if requested
    if docker_secrets_dir:
        print(f"\nWriting Docker secrets to {docker_secrets_dir}")
        if not create_docker_secrets(generated_secrets, docker_secrets_dir):
            return False
    
    print(colorize("\nSecret generation complete!", Colors.GREEN))
    return True


def list_secrets(schema_path: str) -> None:
    """
    List all secrets defined in the schema.
    
    Args:
        schema_path: Path to the schema file
    """
    secrets = identify_secrets_from_schema(schema_path)
    
    if not secrets:
        print("No secrets identified in schema")
        return
    
    print(colorize("Secrets in schema:", Colors.BOLD))
    print(f"{'Path':<30} {'Environment Variable':<30} {'Required':<10} {'Description'}")
    print(f"{'-' * 30} {'-' * 30} {'-' * 10} {'-' * 40}")
    
    for path, info in sorted(secrets.items()):
        required = "Yes" if info.get("required") else "No"
        description = info.get("description", "")
        if len(description) > 40:
            description = description[:37] + "..."
        print(f"{path:<30} {info['env_var']:<30} {required:<10} {description}")


def rotate_secret(env_var: str, env_file: str, docker_secrets_dir: Optional[str] = None) -> bool:
    """
    Rotate a specific secret.
    
    Args:
        env_var: Environment variable name
        env_file: Path to the .env file
        docker_secrets_dir: Optional directory for Docker secrets
        
    Returns:
        True if successful, False otherwise
    """
    # Check if the secret exists
    existing_value = os.environ.get(env_var)
    
    # Prompt for new value
    print(f"Rotating secret for {env_var}")
    if existing_value:
        print("  Current value is set in environment")
    
    new_value = getpass.getpass("  Enter new value (leave empty to generate random): ")
    if not new_value:
        new_value = generate_random_secret()
        print("  Generated random value")
    
    # Write to .env file
    if env_file:
        print(f"  Writing to {env_file}")
        if not write_secrets_to_env_file(env_file, {env_var: new_value}, overwrite=True):
            return False
    
    # Write Docker secret if requested
    if docker_secrets_dir:
        print(f"  Writing Docker secret to {docker_secrets_dir}")
        if not create_docker_secrets({env_var: new_value}, docker_secrets_dir):
            return False
    
    print(colorize("  Secret rotation complete!", Colors.GREEN))
    return True


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Manage configuration secrets")
    parser.add_argument("--schema", default="schema.json", help="Path to schema file")
    parser.add_argument("--env-file", default="../.env", help="Path to .env file")
    parser.add_argument("--docker-secrets", help="Path to Docker secrets directory")
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Generate subcommand
    generate_parser = subparsers.add_parser("generate", help="Generate secrets")
    
    # List subcommand
    list_parser = subparsers.add_parser("list", help="List secrets")
    
    # Rotate subcommand
    rotate_parser = subparsers.add_parser("rotate", help="Rotate a secret")
    rotate_parser.add_argument("env_var", help="Environment variable to rotate")
    
    args = parser.parse_args()
    
    # Resolve paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    schema_path = os.path.abspath(os.path.join(script_dir, args.schema))
    env_file = os.path.abspath(os.path.join(script_dir, args.env_file))
    docker_secrets_dir = os.path.abspath(os.path.join(script_dir, args.docker_secrets)) if args.docker_secrets else None
    
    if args.command == "generate":
        return 0 if generate_secrets_for_schema(schema_path, env_file, docker_secrets_dir) else 1
    elif args.command == "list":
        list_secrets(schema_path)
        return 0
    elif args.command == "rotate":
        return 0 if rotate_secret(args.env_var, env_file, docker_secrets_dir) else 1
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
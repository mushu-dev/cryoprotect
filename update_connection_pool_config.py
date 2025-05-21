#!/usr/bin/env python3
"""
Update Connection Pool Configuration

This script updates the connection pool configuration parameters in the config.py file
based on the recommended settings from the stress test results.

Usage:
    python update_connection_pool_config.py [--minimal] [--medium] [--heavy]

Options:
    --minimal   Apply minimal recommended configuration
    --medium    Apply medium workload configuration (default)
    --heavy     Apply heavy workload configuration
"""

import os
import sys
import re
import argparse
import shutil
from datetime import datetime


def backup_config_file(config_path):
    """Create a backup of the config file."""
    backup_dir = "backups"
    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
        
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.basename(config_path)
    backup_path = os.path.join(backup_dir, f"{filename}.{timestamp}")
    
    shutil.copy2(config_path, backup_path)
    print(f"Created backup at: {backup_path}")
    return backup_path


def get_minimal_config():
    """Get minimal configuration for light workloads."""
    return {
        'SUPABASE_MIN_CONNECTIONS': 2,
        'SUPABASE_MAX_CONNECTIONS': 10,
        'SUPABASE_CONNECTION_TIMEOUT': 20,
        'SUPABASE_CONNECTION_LIFETIME': 3600,  # 1 hour
        'SUPABASE_IDLE_TIMEOUT': 300,  # 5 minutes
        'SUPABASE_HEALTH_CHECK_INTERVAL': 60,  # 1 minute
        'SUPABASE_RETRY_ATTEMPTS': 3
    }


def get_medium_config():
    """Get recommended configuration for medium workloads."""
    return {
        'SUPABASE_MIN_CONNECTIONS': 3,
        'SUPABASE_MAX_CONNECTIONS': 20,
        'SUPABASE_CONNECTION_TIMEOUT': 15,
        'SUPABASE_CONNECTION_LIFETIME': 1800,  # 30 minutes
        'SUPABASE_IDLE_TIMEOUT': 240,  # 4 minutes
        'SUPABASE_HEALTH_CHECK_INTERVAL': 45,  # 45 seconds
        'SUPABASE_RETRY_ATTEMPTS': 3,
        'SUPABASE_INITIAL_RETRY_DELAY': 0.2,
        'SUPABASE_MAX_RETRY_DELAY': 5,
        'SUPABASE_RETRY_JITTER_FACTOR': 0.1,
        'SUPABASE_VALIDATION_QUERY': "'SELECT 1'",
        'SUPABASE_VALIDATION_TIMEOUT': 5
    }


def get_heavy_config():
    """Get recommended configuration for heavy workloads."""
    return {
        'SUPABASE_MIN_CONNECTIONS': 5,
        'SUPABASE_MAX_CONNECTIONS': 30,
        'SUPABASE_CONNECTION_TIMEOUT': 10,
        'SUPABASE_CONNECTION_LIFETIME': 1800,  # 30 minutes
        'SUPABASE_IDLE_TIMEOUT': 180,  # 3 minutes
        'SUPABASE_HEALTH_CHECK_INTERVAL': 30,  # 30 seconds
        'SUPABASE_RETRY_ATTEMPTS': 5,
        'SUPABASE_INITIAL_RETRY_DELAY': 0.1,
        'SUPABASE_MAX_RETRY_DELAY': 5,
        'SUPABASE_RETRY_JITTER_FACTOR': 0.15,
        'SUPABASE_VALIDATION_QUERY': "'SELECT 1'",
        'SUPABASE_VALIDATION_TIMEOUT': 5,
        'SUPABASE_CIRCUIT_BREAKER_THRESHOLD': 5,
        'SUPABASE_CIRCUIT_BREAKER_TIMEOUT': 30,
        'SUPABASE_CIRCUIT_BREAKER_RESET': 2
    }


def update_config_file(config_path, config_values):
    """
    Update the config file with new connection pool parameters.
    
    Args:
        config_path: Path to the config.py file
        config_values: Dictionary of configuration values to update
    """
    # Read the config file
    with open(config_path, 'r') as f:
        content = f.read()
    
    # Create a backup
    backup_path = backup_config_file(config_path)
    
    # Check if the connection pool section exists
    connection_pool_section = re.search(r'# Connection pool settings.*?(?=\n\n|$)', content, re.DOTALL)
    if connection_pool_section:
        # Replace existing connection pool section
        section_content = connection_pool_section.group(0)
        new_section = "# Connection pool settings\n"
        for key, value in config_values.items():
            new_section += f"{key}: {type(value).__name__} = {value}\n"
        
        content = content.replace(section_content, new_section)
    else:
        # Add new connection pool section at the end of the file or before a specific section
        new_section = "\n# Connection pool settings\n"
        for key, value in config_values.items():
            new_section += f"{key}: {type(value).__name__} = {value}\n"
        
        # Try to find a good place to insert the section
        # First try to find the BaseConfig class
        base_config_match = re.search(r'class BaseConfig\(.*?\):', content)
        if base_config_match:
            # Find the end of the class definition
            class_end_pattern = r'class BaseConfig\(.*?\):.*?(?=\n\s*\n|\Z)'
            class_match = re.search(class_end_pattern, content, re.DOTALL)
            if class_match:
                class_content = class_match.group(0)
                indented_section = "\n    " + new_section.replace("\n", "\n    ")
                # Remove trailing whitespace from indented section
                indented_section = re.sub(r'\s+$', '', indented_section, flags=re.MULTILINE)
                updated_class = class_content + indented_section
                content = content.replace(class_content, updated_class)
            else:
                # Just append to the end of the file
                content += new_section
        else:
            # Just append to the end of the file
            content += new_section
    
    # Write updated content back to the file
    with open(config_path, 'w') as f:
        f.write(content)
    
    print(f"Updated config file: {config_path}")
    print(f"Added/updated the following configuration parameters:")
    for key, value in config_values.items():
        print(f"  {key} = {value}")


def main():
    """Main function to update the config file."""
    parser = argparse.ArgumentParser(description="Update connection pool configuration parameters.")
    parser.add_argument("--minimal", action="store_true", help="Apply minimal configuration for light workloads")
    parser.add_argument("--medium", action="store_true", help="Apply medium workload configuration (default)")
    parser.add_argument("--heavy", action="store_true", help="Apply heavy workload configuration")
    parser.add_argument("--config", default="config.py", help="Path to the config.py file (default: config.py)")
    args = parser.parse_args()
    
    # Determine which configuration to apply
    if args.minimal:
        config_values = get_minimal_config()
        config_name = "minimal (light workload)"
    elif args.heavy:
        config_values = get_heavy_config()
        config_name = "heavy workload"
    else:
        # Default to medium config
        config_values = get_medium_config()
        config_name = "medium workload"
    
    # Check if config file exists
    config_path = args.config
    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)
    
    print(f"Applying {config_name} configuration to {config_path}")
    
    # Update the config file
    update_config_file(config_path, config_values)
    
    print(f"\nConnection pool configuration updated successfully!")
    print(f"You may need to restart your application for changes to take effect.")
    

if __name__ == "__main__":
    main()
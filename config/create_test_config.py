#\!/usr/bin/env python3
"""
Create Test Configuration File

This script creates a test configuration file for the database connection
using the known working Supabase credentials that were verified through MCP.
It also creates a shell script to set the environment variable.
"""

import os
import json
import subprocess
from pathlib import Path

# Supabase connection parameters from MCP
SUPABASE_URL = "https://tsdlmynydfuypiugmkev.supabase.co"
SUPABASE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRzZGxteW55ZGZ1eXBpdWdta2V2Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDQ3ODAxODMsImV4cCI6MjA2MDM1NjE4M30.T6udmD-3lhlTy4bY2p0y2lX5-11yvyn425PWrlnIPLU"
DATABASE_HOST = "db.tsdlmynydfuypiugmkev.supabase.co"
DATABASE_PORT = 5432
DATABASE_NAME = "postgres"
DATABASE_USER = "postgres"

def create_test_config():
    """Create a test configuration file for the database connection."""
    current_dir = Path(__file__).parent.absolute()
    project_root = current_dir.parent.absolute()
    
    # Prompt for the database password
    password = input("Enter the database password: ")
    
    # Create configuration structure
    config = {
        "database": {
            "connection": {
                "mode": "supabase",
                "local": {
                    "host": "localhost",
                    "port": 5432,
                    "database": "cryoprotect",
                    "user": "postgres",
                    "password": "postgres",
                    "application_name": "CryoProtect-Local",
                    "min_connections": 1,
                    "max_connections": 5
                },
                "supabase": {
                    "host": DATABASE_HOST,
                    "port": DATABASE_PORT,
                    "database": DATABASE_NAME,
                    "user": DATABASE_USER,
                    "password": password,
                    "application_name": "CryoProtect-Test",
                    "url": SUPABASE_URL,
                    "key": SUPABASE_KEY,
                    "use_ssl": True,
                    "min_connections": 1,
                    "max_connections": 5
                }
            }
        }
    }
    
    # Create config directory if it doesn't exist
    config_dir = current_dir
    config_dir.mkdir(exist_ok=True)
    
    # Write config to file
    config_path = config_dir / "test_config.json"
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    # Create a shell script to set the environment variable
    env_script_path = project_root / "set_test_config_env.sh"
    with open(env_script_path, 'w') as f:
        f.write(f"""#\!/bin/bash
# This script sets the environment variable for the test configuration
export CRYOPROTECT_CONFIG_PATH="{config_path}"
echo "Environment variable CRYOPROTECT_CONFIG_PATH has been set to {config_path}"
echo "To persist this variable, add the following line to your ~/.bashrc file:"
echo "export CRYOPROTECT_CONFIG_PATH=\"{config_path}\""
""")
    
    # Make the script executable
    env_script_path.chmod(0o755)
    
    print(f"\nTest configuration created at {config_path}")
    print(f"\nA script has been created to set the environment variable:")
    print(f"{env_script_path}")
    print("\nTo use this configuration, run the following command:")
    print(f"source {env_script_path}")
    
    # Ask if user wants to set the environment variable now
    set_now = input("\nDo you want to set the environment variable now? (y/n): ")
    if set_now.lower() in ('y', 'yes'):
        # Set the environment variable in the current process
        os.environ["CRYOPROTECT_CONFIG_PATH"] = str(config_path)
        print(f"Environment variable CRYOPROTECT_CONFIG_PATH set to {config_path}")
        print("Note: This will only affect the current process. To use this in other terminals,")
        print(f"run 'source {env_script_path}' in each terminal.")

if __name__ == "__main__":
    create_test_config()

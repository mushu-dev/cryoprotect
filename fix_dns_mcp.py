#!/usr/bin/env python3
"""
CryoProtect v2 - Fix DNS with MCP

This script attempts to resolve Supabase connection issues by:
1. Using MCP tools to get project info
2. Testing connection to Supabase
3. Updating the .env file with correct URL and other connection details
"""

import os
import json
import subprocess
import socket
from dotenv import load_dotenv, set_key

# Load environment variables
load_dotenv()

# Get Supabase configuration
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")

print("Current Supabase Configuration:")
print(f"  URL: {SUPABASE_URL}")
print(f"  Project ID: {SUPABASE_PROJECT_ID}")
print(f"  Key Type: {'service_role' if 'service_role' in SUPABASE_KEY else 'anon'}")
print()

# Function to use MCP to get project info
def get_project_info():
    print("Attempting to get project info via MCP...")
    
    try:
        # Create temporary args file
        args = {"id": SUPABASE_PROJECT_ID}
        args_str = json.dumps(args)
        
        # Try getting project details using MCP
        cmd = [
            "npx", 
            "-y", 
            "@supabase/mcp-server-supabase@latest", 
            "get_project", 
            "--access-token", 
            "sbp_2ef753d5e351cd40412c70c7ed9852c59f18a559",
            "--args", 
            args_str
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error getting project info: {result.stderr}")
            return None
        
        try:
            project_info = json.loads(result.stdout)
            print("Project information retrieved successfully!")
            return project_info
        except json.JSONDecodeError:
            print(f"Failed to parse project info: {result.stdout}")
            return None
            
    except Exception as e:
        print(f"Exception while getting project info: {str(e)}")
        return None

# Function to check connection with project info
def check_connection_with_info(project_info):
    print("\nChecking connection with retrieved project info...")
    
    if not project_info:
        print("No project info available to check connection.")
        return False
    
    # Get and validate project URL
    project_domain = project_info.get("domain")
    if not project_domain:
        print("Project domain not found in project info.")
        return False
    
    # Check if domain follows expected format
    if not project_domain.endswith(".supabase.co"):
        print(f"Unexpected domain format: {project_domain}")
        return False
    
    # Try to resolve the hostname
    try:
        print(f"Resolving hostname {project_domain}...")
        ip_address = socket.gethostbyname(project_domain)
        print(f"Hostname resolved successfully: {project_domain} -> {ip_address}")
    except socket.gaierror as e:
        print(f"Failed to resolve hostname {project_domain}: {str(e)}")
        return False
    
    # Try to resolve the database hostname
    db_hostname = f"db.{project_domain}"
    try:
        print(f"Resolving database hostname {db_hostname}...")
        db_ip_address = socket.gethostbyname(db_hostname)
        print(f"Database hostname resolved successfully: {db_hostname} -> {db_ip_address}")
    except socket.gaierror as e:
        print(f"Failed to resolve database hostname {db_hostname}: {str(e)}")
        return False
    
    return True

# Function to update .env file with correct connection details
def update_env_file(project_info):
    print("\nUpdating .env file with correct connection details...")
    
    if not project_info:
        print("No project info available to update .env file.")
        return False
    
    # Validate project domain
    project_domain = project_info.get("domain")
    if not project_domain or not project_domain.endswith(".supabase.co"):
        print(f"Invalid project domain: {project_domain}")
        return False
    
    # Prepare new values
    new_url = f"https://{project_domain}"
    db_hostname = f"db.{project_domain}"
    
    # Backup .env file
    try:
        with open(".env", "r") as f:
            env_content = f.read()
            
        with open(".env.backup", "w") as f:
            f.write(env_content)
            
        print("Created backup of .env file as .env.backup")
    except Exception as e:
        print(f"Failed to backup .env file: {str(e)}")
        return False
    
    # Update .env file
    try:
        dotenv_path = ".env"
        set_key(dotenv_path, "SUPABASE_URL", new_url)
        set_key(dotenv_path, "SUPABASE_DB_HOST", db_hostname)
        print(f"Updated SUPABASE_URL to {new_url}")
        print(f"Updated SUPABASE_DB_HOST to {db_hostname}")
        
        # Also add a comment explaining the change
        with open(".env", "r") as f:
            lines = f.readlines()
        
        # Add comment if not already there
        comment_line = "# Updated by fix_dns_mcp.py script\n"
        if not any(comment_line == line for line in lines):
            with open(".env", "w") as f:
                f.write(comment_line)
                f.writelines(lines)
        
        return True
    except Exception as e:
        print(f"Failed to update .env file: {str(e)}")
        return False

# Main function
def main():
    # Print header
    print("=" * 80)
    print("CryoProtect v2 - Fix DNS with MCP")
    print("=" * 80)
    print("This script will use MCP to fix Supabase connection issues.")
    print()
    
    # Step 1: Get project info
    project_info = get_project_info()
    
    if not project_info:
        print("\nFailed to get project info. Trying to use existing Project ID...")
        # Create a minimal project info structure
        if SUPABASE_PROJECT_ID:
            project_info = {
                "domain": f"{SUPABASE_PROJECT_ID}.supabase.co"
            }
            print(f"Created project info with domain: {project_info['domain']}")
        else:
            print("\nCannot proceed without a valid Project ID. Please set SUPABASE_PROJECT_ID in .env file.")
            return
    
    # Print project info
    print("\nProject Information:")
    print(f"  Name: {project_info.get('name', 'Unknown')}")
    print(f"  Domain: {project_info.get('domain', 'Unknown')}")
    print(f"  Region: {project_info.get('region', 'Unknown')}")
    print(f"  Status: {project_info.get('status', 'Unknown')}")
    
    # Step 2: Check connection with project info
    connection_ok = check_connection_with_info(project_info)
    
    # Step 3: Update .env file if successful
    if connection_ok:
        print("\nConnection check successful!")
        update_success = update_env_file(project_info)
        
        if update_success:
            print("\nSuccessfully updated .env file with correct connection details.")
            print("\nNext steps:")
            print("1. Restart your application to apply changes")
            print("2. Run 'python supabase_connection_diagnostic.py' to verify the connection")
        else:
            print("\nFailed to update .env file.")
    else:
        print("\nConnection check failed. Cannot update .env file.")
        print("\nYou may need to:")
        print("1. Check if your DNS server can resolve *.supabase.co domains")
        print("2. Try using an alternative DNS server (8.8.8.8 or 1.1.1.1)")
        print("3. Manually update your .env file with these values:")
        print(f"   SUPABASE_URL=https://{project_info.get('domain', SUPABASE_PROJECT_ID + '.supabase.co')}")
        print(f"   SUPABASE_DB_HOST=db.{project_info.get('domain', SUPABASE_PROJECT_ID + '.supabase.co')}")

if __name__ == "__main__":
    main()

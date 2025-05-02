#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Supabase MCP Connection Issues

This script checks and fixes issues with the Supabase MCP tools:
1. Verifies Node.js and npm installation
2. Installs or updates the Supabase MCP package
3. Tests the MCP connection
4. Updates the MCP token if needed
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# ANSI color codes for better readability
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print(f"{Colors.HEADER}{Colors.BOLD}CryoProtect v2 - Fix Supabase MCP Connection Issues{Colors.ENDC}")
print(f"{Colors.CYAN}This tool will diagnose and fix issues with the Supabase MCP tools.{Colors.ENDC}")
print()

# Extract Supabase configuration
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")
if not SUPABASE_PROJECT_ID:
    print(f"{Colors.WARNING}WARNING: SUPABASE_PROJECT_ID not found in .env file.{Colors.ENDC}")
    SUPABASE_URL = os.getenv("SUPABASE_URL", "")
    if SUPABASE_URL:
        import urllib.parse
        parsed_url = urllib.parse.urlparse(SUPABASE_URL)
        hostname = parsed_url.netloc
        if hostname and hostname.endswith('.supabase.co'):
            SUPABASE_PROJECT_ID = hostname.split('.')[0]
            print(f"{Colors.GREEN}✓ Extracted project ID from URL: {SUPABASE_PROJECT_ID}{Colors.ENDC}")
        else:
            print(f"{Colors.FAIL}ERROR: Could not extract project ID from SUPABASE_URL{Colors.ENDC}")
            sys.exit(1)
    else:
        print(f"{Colors.FAIL}ERROR: Neither SUPABASE_PROJECT_ID nor SUPABASE_URL found in .env file{Colors.ENDC}")
        sys.exit(1)

# Step 1: Check Node.js and npm installation
def check_node_environment():
    print(f"{Colors.BOLD}[1] Checking Node.js and npm installation...{Colors.ENDC}")
    
    # Check Node.js
    try:
        result = subprocess.run(['node', '--version'], capture_output=True, text=True, check=False)
        if result.returncode == 0:
            node_version = result.stdout.strip()
            print(f"  {Colors.GREEN}✓ Node.js installed: {node_version}{Colors.ENDC}")
        else:
            print(f"  {Colors.FAIL}✗ Node.js not installed or not in PATH{Colors.ENDC}")
            print(f"  {Colors.BOLD}Please install Node.js from https://nodejs.org/{Colors.ENDC}")
            return False
    except FileNotFoundError:
        print(f"  {Colors.FAIL}✗ Node.js not installed or not in PATH{Colors.ENDC}")
        print(f"  {Colors.BOLD}Please install Node.js from https://nodejs.org/{Colors.ENDC}")
        return False
    
    # Check npm
    try:
        result = subprocess.run(['npm', '--version'], capture_output=True, text=True, check=False)
        if result.returncode == 0:
            npm_version = result.stdout.strip()
            print(f"  {Colors.GREEN}✓ npm installed: {npm_version}{Colors.ENDC}")
            return True
        else:
            print(f"  {Colors.FAIL}✗ npm not installed or not in PATH{Colors.ENDC}")
            print(f"  {Colors.BOLD}npm should be included with Node.js installation{Colors.ENDC}")
            return False
    except FileNotFoundError:
        print(f"  {Colors.FAIL}✗ npm not installed or not in PATH{Colors.ENDC}")
        print(f"  {Colors.BOLD}npm should be included with Node.js installation{Colors.ENDC}")
        return False

# Step 2: Install or update Supabase MCP package
def install_or_update_mcp():
    print(f"\n{Colors.BOLD}[2] Installing/updating Supabase MCP package...{Colors.ENDC}")
    
    try:
        # Check if the package is already installed
        result = subprocess.run(
            ['npx', '-y', '@supabase/mcp-server-supabase@latest', '--version'],
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode == 0:
            version = result.stdout.strip()
            print(f"  {Colors.GREEN}✓ @supabase/mcp-server-supabase is installed: {version}{Colors.ENDC}")
            
            # Check if it needs to be updated
            print(f"  Checking for updates...")
            install_result = subprocess.run(
                ['npm', 'install', '-g', '@supabase/mcp-server-supabase@latest'],
                capture_output=True,
                text=True,
                check=False
            )
            
            if install_result.returncode == 0:
                print(f"  {Colors.GREEN}✓ @supabase/mcp-server-supabase is up to date{Colors.ENDC}")
            else:
                print(f"  {Colors.WARNING}⚠ Could not update package globally. Trying local installation...{Colors.ENDC}")
                local_install_result = subprocess.run(
                    ['npm', 'install', '@supabase/mcp-server-supabase@latest'],
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                if local_install_result.returncode == 0:
                    print(f"  {Colors.GREEN}✓ @supabase/mcp-server-supabase installed locally{Colors.ENDC}")
                else:
                    print(f"  {Colors.FAIL}✗ Failed to install package: {local_install_result.stderr}{Colors.ENDC}")
                    return False
        else:
            print(f"  {Colors.WARNING}⚠ @supabase/mcp-server-supabase is not installed. Installing now...{Colors.ENDC}")
            
            install_result = subprocess.run(
                ['npm', 'install', '-g', '@supabase/mcp-server-supabase@latest'],
                capture_output=True,
                text=True,
                check=False
            )
            
            if install_result.returncode == 0:
                print(f"  {Colors.GREEN}✓ @supabase/mcp-server-supabase installed globally{Colors.ENDC}")
            else:
                print(f"  {Colors.WARNING}⚠ Could not install package globally. Trying local installation...{Colors.ENDC}")
                local_install_result = subprocess.run(
                    ['npm', 'install', '@supabase/mcp-server-supabase@latest'],
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                if local_install_result.returncode == 0:
                    print(f"  {Colors.GREEN}✓ @supabase/mcp-server-supabase installed locally{Colors.ENDC}")
                else:
                    print(f"  {Colors.FAIL}✗ Failed to install package: {local_install_result.stderr}{Colors.ENDC}")
                    return False
        
        return True
    except Exception as e:
        print(f"  {Colors.FAIL}✗ Error installing/updating MCP package: {str(e)}{Colors.ENDC}")
        return False

# Step 3: Check and update MCP token
def check_and_update_mcp_token():
    print(f"\n{Colors.BOLD}[3] Checking and updating MCP token...{Colors.ENDC}")
    
    # Find and read supabase_mcp_tools.py
    mcp_tools_path = './supabase_mcp_tools.py'
    if not os.path.exists(mcp_tools_path):
        print(f"  {Colors.FAIL}✗ Could not find supabase_mcp_tools.py{Colors.ENDC}")
        return False
    
    try:
        with open(mcp_tools_path, 'r') as f:
            content = f.read()
        
        # Check if access token is embedded in the script
        import re
        token_match = re.search(r'--access-token\s+"([^"]+)"', content)
        current_token = token_match.group(1) if token_match else ""
        
        if not current_token:
            print(f"  {Colors.WARNING}⚠ Could not find access token in supabase_mcp_tools.py{Colors.ENDC}")
        else:
            print(f"  {Colors.GREEN}✓ Found existing access token in supabase_mcp_tools.py{Colors.ENDC}")
        
        # Ask user for new token
        print()
        print(f"{Colors.CYAN}To generate a new Supabase access token:{Colors.ENDC}")
        print(f"1. Visit https://supabase.com/dashboard/account/tokens")
        print(f"2. Create a new token with an appropriate name and expiration")
        print(f"3. Copy the token and paste it below")
        print()
        
        token_input = input(f"{Colors.BOLD}Enter new Supabase access token (leave empty to keep existing): {Colors.ENDC}")
        
        if not token_input:
            if current_token:
                print(f"  {Colors.GREEN}✓ Keeping existing token{Colors.ENDC}")
                new_token = current_token
            else:
                print(f"  {Colors.FAIL}✗ No token provided and no existing token found{Colors.ENDC}")
                return False
        else:
            new_token = token_input
            print(f"  {Colors.GREEN}✓ Using new token{Colors.ENDC}")
        
        # Update the token in the script
        if token_match:
            new_content = re.sub(r'--access-token\s+"([^"]+)"', f'--access-token "{new_token}"', content)
        else:
            # If no token found, we need to add it to the appropriate places
            # This is a more complex update that requires understanding the structure of the file
            print(f"  {Colors.WARNING}⚠ Could not find existing token pattern to replace{Colors.ENDC}")
            print(f"  {Colors.BOLD}Manual update required:{Colors.ENDC}")
            print(f"  Please edit {mcp_tools_path} and ensure the access token is set to: {new_token}")
            return False
        
        # Write updated content back to the file
        with open(mcp_tools_path, 'w') as f:
            f.write(new_content)
        
        print(f"  {Colors.GREEN}✓ Updated access token in supabase_mcp_tools.py{Colors.ENDC}")
        return True
        
    except Exception as e:
        print(f"  {Colors.FAIL}✗ Error updating MCP token: {str(e)}{Colors.ENDC}")
        return False

# Step 4: Test MCP connection
def test_mcp_connection():
    print(f"\n{Colors.BOLD}[4] Testing MCP connection...{Colors.ENDC}")
    
    try:
        # Import the supabase_mcp_tools module
        import sys
        sys.path.append('.')
        from supabase_mcp_tools import execute_sql_on_supabase
        
        try:
            print(f"  Testing SQL execution via MCP...")
            result = execute_sql_on_supabase(SUPABASE_PROJECT_ID, "SELECT 1 AS test")
            
            if result:
                print(f"  {Colors.GREEN}✓ MCP connection successful!{Colors.ENDC}")
                print(f"  Result: {json.dumps(result, indent=2)}")
                return True
            else:
                print(f"  {Colors.FAIL}✗ MCP connection failed: No result returned{Colors.ENDC}")
                return False
        except Exception as e:
            print(f"  {Colors.FAIL}✗ MCP connection failed: {str(e)}{Colors.ENDC}")
            
            # Check for common error messages
            error_str = str(e)
            if "ENOTFOUND" in error_str or "getaddrinfo" in error_str:
                print(f"  {Colors.WARNING}⚠ This appears to be a DNS resolution issue.{Colors.ENDC}")
                print(f"  {Colors.BOLD}Recommended fix:{Colors.ENDC}")
                print(f"  Run fix_supabase_dns.py to resolve DNS issues")
            elif "Unauthorized" in error_str or "401" in error_str:
                print(f"  {Colors.WARNING}⚠ This appears to be an authentication issue.{Colors.ENDC}")
                print(f"  {Colors.BOLD}Recommended fix:{Colors.ENDC}")
                print(f"  Generate a new Supabase access token and update it in supabase_mcp_tools.py")
            elif "not found" in error_str.lower() or "command not found" in error_str.lower():
                print(f"  {Colors.WARNING}⚠ The MCP package could not be found.{Colors.ENDC}")
                print(f"  {Colors.BOLD}Recommended fix:{Colors.ENDC}")
                print(f"  Run this script again to reinstall the MCP package")
            
            return False
    except ImportError:
        print(f"  {Colors.FAIL}✗ Could not import supabase_mcp_tools module{Colors.ENDC}")
        return False

# Main function
def main():
    # Check Node.js and npm
    node_ok = check_node_environment()
    if not node_ok:
        print(f"\n{Colors.FAIL}Node.js or npm is not installed. Please install them first.{Colors.ENDC}")
        return
    
    # Install or update MCP package
    mcp_ok = install_or_update_mcp()
    if not mcp_ok:
        print(f"\n{Colors.FAIL}Failed to install or update MCP package.{Colors.ENDC}")
        return
    
    # Check and update MCP token
    token_ok = check_and_update_mcp_token()
    
    # Test connection
    connection_ok = test_mcp_connection()
    
    # Summary
    print(f"\n{Colors.HEADER}{Colors.BOLD}Fix Summary:{Colors.ENDC}")
    print(f"  Node.js/npm:       {'✓' if node_ok else '✗'}")
    print(f"  MCP Package:       {'✓' if mcp_ok else '✗'}")
    print(f"  MCP Token:         {'✓' if token_ok else '?'}")
    print(f"  MCP Connection:    {'✓' if connection_ok else '✗'}")
    
    print(f"\n{Colors.BOLD}Overall Status:{Colors.ENDC}")
    if node_ok and mcp_ok and connection_ok:
        print(f"{Colors.GREEN}✓ All checks passed! Your Supabase MCP tools are working correctly.{Colors.ENDC}")
    else:
        print(f"{Colors.FAIL}✗ Some checks failed. See above for specific issues and recommended fixes.{Colors.ENDC}")
    
    # Next steps
    print(f"\n{Colors.BOLD}Next Steps:{Colors.ENDC}")
    if connection_ok:
        print(f"1. You can now use the Supabase MCP tools to interact with your database")
        print(f"2. Test your application's connection to Supabase")
    else:
        print(f"1. Fix any remaining issues mentioned above")
        print(f"2. If DNS issues persist, run fix_supabase_dns.py")
        print(f"3. Try alternative solution: use direct database connection instead of MCP")

if __name__ == "__main__":
    main()
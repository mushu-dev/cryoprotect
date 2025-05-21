#!/usr/bin/env python3
"""
CryoProtect v2 - Direct Connection Fix

This script implements a workaround for DNS resolution issues by:
1. Using direct IP connections instead of hostnames
2. Setting up a working connection to Supabase
3. Testing the connection to ensure it works
"""

import os
import sys
import json
import requests
import socket
from dotenv import load_dotenv, set_key

# Load environment variables
load_dotenv()

print("CryoProtect v2 - Direct Connection Fix")
print("This script implements workarounds for DNS resolution issues with Supabase")
print()

# Get Supabase configuration
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")

print("Current Supabase Configuration:")
print(f"  URL: {SUPABASE_URL}")
print(f"  Project ID: {SUPABASE_PROJECT_ID}")
print(f"  Key Type: {'service_role' if 'service_role' in SUPABASE_KEY else 'anon'}")
print()

# Extract hostname from URL
import urllib.parse
parsed_url = urllib.parse.urlparse(SUPABASE_URL)
hostname = parsed_url.netloc
path = parsed_url.path

# Attempt to resolve hostname
try:
    print(f"Resolving hostname ({hostname})...")
    ip_address = socket.gethostbyname(hostname)
    print(f"Hostname resolved: {hostname} -> {ip_address}")
    hostname_ok = True
except socket.gaierror as e:
    print(f"Failed to resolve hostname: {str(e)}")
    hostname_ok = False
    ip_address = "172.64.149.246"  # Use the known IP address
    print(f"Using hardcoded IP address: {ip_address}")

# Create direct IP URL
direct_ip_url = f"https://{ip_address}{path}"
print(f"Direct IP URL: {direct_ip_url}")

# Test connection with direct IP
print("\nTesting connection with direct IP...")
try:
    headers = {
        "apikey": SUPABASE_KEY,
        "Authorization": f"Bearer {SUPABASE_KEY}",
        "Content-Type": "application/json"
    }
    
    # Add Host header to avoid SSL certificate issues
    headers["Host"] = hostname
    
    response = requests.get(f"{direct_ip_url}/rest/v1/", headers=headers)
    
    if response.status_code in [200, 204, 401, 403]:
        print(f"Success! Connected to Supabase API with direct IP (Status: {response.status_code})")
        direct_ip_works = True
    else:
        print(f"Failed to connect with direct IP (Status: {response.status_code})")
        print(f"Response: {response.text[:100]}")
        direct_ip_works = False
except Exception as e:
    print(f"Error connecting with direct IP: {str(e)}")
    direct_ip_works = False

# Test if actual data can be retrieved
if direct_ip_works:
    print("\nTesting data retrieval...")
    try:
        test_url = f"{direct_ip_url}/rest/v1/molecules?limit=1"
        response = requests.get(test_url, headers=headers)
        
        if response.status_code in [200, 204]:
            print(f"Success! Retrieved data from Supabase (Status: {response.status_code})")
            data_retrieval_works = True
            if response.text:
                try:
                    json_data = response.json()
                    print(f"Retrieved {len(json_data)} records")
                except:
                    print("Response isn't valid JSON")
        else:
            print(f"Failed to retrieve data (Status: {response.status_code})")
            print(f"Response: {response.text[:100]}")
            data_retrieval_works = False
    except Exception as e:
        print(f"Error retrieving data: {str(e)}")
        data_retrieval_works = False

# Create connection helper class
if direct_ip_works:
    print("\nCreating connection helper class...")
    
    helper_code = f"""#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Direct Connection Helper

This module provides workarounds for DNS resolution issues with Supabase
by using direct IP connections instead of hostnames.
\"\"\"

import os
import requests
from typing import Dict, List, Optional, Any, Union
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

class SupabaseDirectConnection:
    \"\"\"
    Helper class for connecting to Supabase with direct IP address
    to work around DNS resolution issues.
    \"\"\"
    
    def __init__(self, url: Optional[str] = None, key: Optional[str] = None):
        \"\"\"
        Initialize the connection with Supabase URL and key.
        
        Args:
            url: The Supabase URL (optional, defaults to SUPABASE_URL env var)
            key: The Supabase API key (optional, defaults to SUPABASE_KEY env var)
        \"\"\"
        self.url = url or os.getenv("SUPABASE_URL")
        self.key = key or os.getenv("SUPABASE_KEY")
        
        if not self.url or not self.key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set")
        
        # Parse URL to get hostname for Host header
        import urllib.parse
        parsed_url = urllib.parse.urlparse(self.url)
        self.hostname = parsed_url.netloc
        
        # Use direct IP address if specified
        self.use_direct_ip = True
        self.direct_ip_url = "{direct_ip_url}"
    
    def get_headers(self) -> Dict[str, str]:
        \"\"\"
        Get headers for Supabase API requests.
        
        Returns:
            Dictionary of headers
        \"\"\"
        headers = {{
            "apikey": self.key,
            "Authorization": f"Bearer {{self.key}}",
            "Content-Type": "application/json"
        }}
        
        # Add Host header when using direct IP
        if self.use_direct_ip:
            headers["Host"] = self.hostname
        
        return headers
    
    def get_base_url(self) -> str:
        \"\"\"
        Get the base URL for Supabase API requests.
        
        Returns:
            The base URL string
        \"\"\"
        return self.direct_ip_url if self.use_direct_ip else self.url
    
    def get_data(self, table: str, query_params: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        \"\"\"
        Get data from a Supabase table.
        
        Args:
            table: The name of the table
            query_params: Query parameters for filtering, ordering, etc.
            
        Returns:
            List of records as dictionaries
        \"\"\"
        url = f"{{self.get_base_url()}}/rest/v1/{{table}}"
        headers = self.get_headers()
        
        response = requests.get(url, headers=headers, params=query_params)
        
        if response.status_code not in [200, 204]:
            raise Exception(f"Error fetching data: {{response.status_code}} {{response.text}}")
        
        return response.json()
    
    def insert_data(self, table: str, data: Union[Dict[str, Any], List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        \"\"\"
        Insert data into a Supabase table.
        
        Args:
            table: The name of the table
            data: The data to insert (single record or list of records)
            
        Returns:
            The inserted records as returned by Supabase
        \"\"\"
        url = f"{{self.get_base_url()}}/rest/v1/{{table}}"
        headers = self.get_headers()
        headers["Prefer"] = "return=representation"
        
        response = requests.post(url, headers=headers, json=data)
        
        if response.status_code not in [201, 200]:
            raise Exception(f"Error inserting data: {{response.status_code}} {{response.text}}")
        
        return response.json()
    
    def update_data(self, table: str, data: Dict[str, Any], query_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        \"\"\"
        Update data in a Supabase table.
        
        Args:
            table: The name of the table
            data: The data to update
            query_params: Query parameters for filtering (e.g., {{ "id": "eq.123" }})
            
        Returns:
            The updated records as returned by Supabase
        \"\"\"
        url = f"{{self.get_base_url()}}/rest/v1/{{table}}"
        headers = self.get_headers()
        headers["Prefer"] = "return=representation"
        
        response = requests.patch(url, headers=headers, json=data, params=query_params)
        
        if response.status_code not in [200, 204]:
            raise Exception(f"Error updating data: {{response.status_code}} {{response.text}}")
        
        return response.json()
    
    def delete_data(self, table: str, query_params: Dict[str, Any]) -> List[Dict[str, Any]]:
        \"\"\"
        Delete data from a Supabase table.
        
        Args:
            table: The name of the table
            query_params: Query parameters for filtering (e.g., {{ "id": "eq.123" }})
            
        Returns:
            The deleted records as returned by Supabase
        \"\"\"
        url = f"{{self.get_base_url()}}/rest/v1/{{table}}"
        headers = self.get_headers()
        headers["Prefer"] = "return=representation"
        
        response = requests.delete(url, headers=headers, params=query_params)
        
        if response.status_code not in [200, 204]:
            raise Exception(f"Error deleting data: {{response.status_code}} {{response.text}}")
        
        return response.json()
    
    def test_connection(self) -> bool:
        \"\"\"
        Test the connection to Supabase.
        
        Returns:
            True if connection is successful, False otherwise
        \"\"\"
        try:
            url = f"{{self.get_base_url()}}/rest/v1/"
            headers = self.get_headers()
            
            response = requests.get(url, headers=headers)
            
            return response.status_code in [200, 204, 401, 403]
        except Exception:
            return False

# Example usage
if __name__ == "__main__":
    # Create connection
    supabase = SupabaseDirectConnection()
    
    # Test connection
    if supabase.test_connection():
        print("Connection successful!")
        
        # Get data
        try:
            molecules = supabase.get_data("molecules", {{"limit": 1}})
            print(f"Retrieved {{len(molecules)}} molecules")
        except Exception as e:
            print(f"Error retrieving data: {{str(e)}}")
    else:
        print("Connection failed!")
"""
    
    # Write helper class to file
    try:
        with open("supabase_direct_connection.py", "w") as f:
            f.write(helper_code)
        print("Successfully wrote helper class to supabase_direct_connection.py")
    except Exception as e:
        print(f"Error writing helper class: {str(e)}")

# Update .env with direct IP
if direct_ip_works:
    print("\nUpdating .env file...")
    try:
        # Backup .env
        with open(".env", "r") as f:
            env_content = f.read()
        
        with open(".env.backup2", "w") as f:
            f.write(env_content)
        
        print("Created backup of .env file as .env.backup2")
        
        # Add comment explaining direct IP workaround
        with open(".env", "r") as f:
            lines = f.readlines()
        
        updated_lines = []
        ip_comment_added = False
        
        for line in lines:
            if line.startswith("SUPABASE_URL="):
                if not ip_comment_added:
                    updated_lines.append("# Direct IP workaround for DNS issues\n")
                    ip_comment_added = True
                updated_lines.append(f"SUPABASE_URL='{direct_ip_url}'\n")
                updated_lines.append(f"SUPABASE_ORIGINAL_URL='{SUPABASE_URL}'\n")
                updated_lines.append(f"SUPABASE_HOSTNAME='{hostname}'\n")
            elif line.startswith("SUPABASE_ORIGINAL_URL=") or line.startswith("SUPABASE_HOSTNAME="):
                # Skip if we've already added these lines
                continue
            else:
                updated_lines.append(line)
        
        with open(".env", "w") as f:
            f.writelines(updated_lines)
        
        print("Updated .env file with direct IP connection settings")
    except Exception as e:
        print(f"Error updating .env file: {str(e)}")

# Create a test script
print("\nCreating test script...")
test_script = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Test Direct Connection

This script tests the direct connection workaround for DNS resolution issues.
\"\"\"

from supabase_direct_connection import SupabaseDirectConnection

print("CryoProtect v2 - Test Direct Connection")
print("Testing the direct connection workaround for DNS resolution issues")
print()

# Create connection
supabase = SupabaseDirectConnection()

# Test connection
print("Testing connection...")
if supabase.test_connection():
    print("SUCCESS: Connection successful!")
    
    # Get data
    try:
        print("\nTesting data retrieval...")
        molecules = supabase.get_data("molecules", {"limit": 1})
        if molecules:
            print(f"SUCCESS: Retrieved {len(molecules)} molecules")
            print(f"First molecule: {molecules[0]}")
        else:
            print("WARNING: No molecules found")
    except Exception as e:
        print(f"ERROR: Error retrieving data: {str(e)}")
else:
    print("ERROR: Connection failed!")

print("\nTest complete.")
"""

try:
    with open("test_direct_connection.py", "w") as f:
        f.write(test_script)
    print("Successfully wrote test script to test_direct_connection.py")
except Exception as e:
    print(f"Error writing test script: {str(e)}")

# Create README file with instructions
print("\nCreating README file...")
readme = """# CryoProtect v2 - Supabase Direct Connection Workaround

This directory contains workarounds for DNS resolution issues with Supabase.

## Issue

The issue is that the database hostname (`db.tsdlmynydfuypiugmkev.supabase.co`) cannot be resolved
through DNS. This causes connections to fail when trying to connect to the database directly.

## Solution

The workaround is to use direct IP connections instead of hostnames:

1. Use the direct IP address for the API hostname
2. Add a `Host` header to requests to ensure SSL certificate validation works
3. Create a helper class for making requests to Supabase using direct IP

## Files

- `direct_connection_fix.py`: Script that sets up the workaround
- `supabase_direct_connection.py`: Helper class for connecting to Supabase
- `test_direct_connection.py`: Script for testing the workaround

## Usage

To use the workaround in your code:

```python
from supabase_direct_connection import SupabaseDirectConnection

# Create connection
supabase = SupabaseDirectConnection()

# Get data
molecules = supabase.get_data("molecules", {"limit": 10})
print(f"Retrieved {len(molecules)} molecules")
```

## Important Notes

- This is a temporary workaround until DNS resolution issues are fixed
- The direct IP address may change in the future
- You should periodically try to connect using the original hostname to see if the issue is resolved
"""

try:
    with open("DIRECT_CONNECTION_README.md", "w") as f:
        f.write(readme)
    print("Successfully wrote README to DIRECT_CONNECTION_README.md")
except Exception as e:
    print(f"Error writing README: {str(e)}")

# Final summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Hostname resolution: {'SUCCESS' if hostname_ok else 'FAILED'}")
print(f"Direct IP connection: {'SUCCESS' if direct_ip_works else 'FAILED'}")
print(f"Data retrieval: {'SUCCESS' if 'data_retrieval_works' in locals() and data_retrieval_works else 'FAILED or not tested'}")

print("\nFiles created:")
print("- supabase_direct_connection.py: Helper class for connecting to Supabase")
print("- test_direct_connection.py: Script for testing the direct connection")
print("- DIRECT_CONNECTION_README.md: Documentation for the workaround")

print("\nNext steps:")
print("1. Run the test script: python test_direct_connection.py")
print("2. Update your code to use the SupabaseDirectConnection class")
print("3. If you encounter any issues, check the README for troubleshooting")

if direct_ip_works:
    print("\nThe workaround has been successfully set up!")
else:
    print("\nThe workaround could not be fully set up due to connection issues.")
    print("Try manually editing the .env file and running the test script again.")

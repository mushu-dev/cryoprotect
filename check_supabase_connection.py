#!/usr/bin/env python3
"""
CryoProtect v2 - Check Supabase Connection

This script checks if the application can connect to Supabase without authentication,
which helps diagnose if the issue is with authentication or the connection itself.
"""

import os
import json
import requests
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

print(f"Checking Supabase connection with:")
print(f"  URL: {SUPABASE_URL}")
print(f"  Using Key Type: {'anon' if 'anon' in SUPABASE_KEY else 'service_role'}")
print()

# Function to check Supabase connection
def check_connection():
    print("Testing Supabase connection...")
    
    # Try to get health status
    health_url = f"{SUPABASE_URL}/rest/v1/"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json"
    }
    
    try:
        response = requests.get(health_url, headers=headers)
        
        if response.status_code in [200, 204]:
            print("[SUCCESS] Connection successful!")
            return True
        else:
            print(f"[FAILED] Connection failed: Status code {response.status_code}")
            print(f"Response: {response.text}")
            return False
    except Exception as e:
        print(f"[ERROR] Connection error: {str(e)}")
        return False

# Function to check if tables are accessible
def check_tables():
    print("\nChecking if tables are accessible...")
    
    tables_to_check = ["molecules", "property_types", "molecular_properties"]
    
    for table in tables_to_check:
        print(f"  Checking table '{table}'...")
        
        table_url = f"{SUPABASE_URL}/rest/v1/{table}?limit=1"
        headers = {
            "apikey": SUPABASE_KEY,
            "Content-Type": "application/json"
        }
        
        try:
            response = requests.get(table_url, headers=headers)
            
            if response.status_code in [200, 204]:
                print(f"  [SUCCESS] Table '{table}' is accessible")
            else:
                print(f"  [FAILED] Table '{table}' is not accessible: Status code {response.status_code}")
                print(f"  Response: {response.text}")
        except Exception as e:
            print(f"  [ERROR] Error checking table '{table}': {str(e)}")

# Main execution
if __name__ == "__main__":
    connection_ok = check_connection()
    
    if connection_ok:
        check_tables()
    
    print("\nConnection check complete.")
    print("\nNext steps:")
    print("1. If connection is successful but authentication fails, run fix_supabase_auth.py")
    print("2. If connection fails, check your Supabase URL and key in the .env file")
    print("3. If tables are not accessible, check your Supabase schema")

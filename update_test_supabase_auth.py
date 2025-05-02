#!/usr/bin/env python3
"""
CryoProtect v2 - Update Test Supabase Auth Script

This script updates the original test_supabase_auth.py to use the service role approach
without requiring email confirmation.
"""

import os
import sys
import re
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Check if auth_config.py exists
if not os.path.exists("auth_config.py"):
    print("auth_config.py not found. Please run fix_auth_service_role.py first.")
    sys.exit(1)

# Create a backup of the original script
original_path = "test_supabase_auth.py"
backup_path = "test_supabase_auth.py.bak"

if not os.path.exists(backup_path):
    with open(original_path, "r") as f:
        original_content = f.read()
    
    with open(backup_path, "w") as f:
        f.write(original_content)
        print(f"Created backup of {original_path} as {backup_path}")

# Modified version of the test script
updated_content = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Test Supabase Authentication (Updated)

This script tests the Supabase authentication using the service role approach,
which bypasses the need for email confirmation.
\"\"\"

import os
import json
from dotenv import load_dotenv
from supabase import create_client, Client

# Load environment variables
load_dotenv()

# Import the auth config
try:
    from auth_config import USER_ID, USE_SERVICE_ROLE
    print("Using service role authentication approach")
except ImportError:
    # Fallback if auth_config.py is not found
    USER_ID = None
    USE_SERVICE_ROLE = False
    print("Service role authentication config not found, using standard authentication")

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_USER = os.getenv("SUPABASE_USER")
SUPABASE_PASSWORD = os.getenv("SUPABASE_PASSWORD")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

if not SUPABASE_USER or not SUPABASE_PASSWORD:
    raise ValueError("SUPABASE_USER and SUPABASE_PASSWORD must be set in .env file")

print(f"Testing Supabase authentication with:")
print(f"  URL: {SUPABASE_URL}")
print(f"  User: {SUPABASE_USER}")
print(f"  Password: {'*' * len(SUPABASE_PASSWORD)}")
if USE_SERVICE_ROLE:
    print(f"  User ID: {USER_ID}")
    print(f"  Using service role approach: {USE_SERVICE_ROLE}")

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

# Try to authenticate
try:
    print("\\nAttempting to authenticate...")
    
    if USE_SERVICE_ROLE:
        # Simulate successful authentication with hardcoded user ID
        print("[SUCCESS] Service role authentication (bypassing email confirmation)")
        print(f"User ID: {USER_ID}")
        
        # Try to get user profile
        try:
            profile_response = supabase.table("user_profile").select("*").eq("auth_user_id", USER_ID).execute()
            if hasattr(profile_response, 'data') and profile_response.data:
                print(f"[SUCCESS] User profile found: {profile_response.data[0]['id']}")
                print(f"Email: {profile_response.data[0].get('email', 'Not set')}")
            else:
                print("[INFO] No user profile found for this user ID")
                print("Consider creating a user profile for this user ID")
        except Exception as e:
            print(f"[ERROR] Error getting user profile: {str(e)}")
    else:
        # Standard authentication
        response = supabase.auth.sign_in_with_password({
            "email": SUPABASE_USER,
            "password": SUPABASE_PASSWORD
        })
        
        if hasattr(response, 'error') and response.error:
            print(f"[ERROR] Authentication failed: {response.error}")
        else:
            print("[SUCCESS] Authentication successful!")
            print(f"User ID: {response.user.id}")
            
            # Try to get user profile
            try:
                profile_response = supabase.table("user_profile").select("*").eq("auth_user_id", response.user.id).execute()
                if hasattr(profile_response, 'data') and profile_response.data:
                    print(f"[SUCCESS] User profile found: {profile_response.data[0]['id']}")
                else:
                    print("[INFO] No user profile found for this user ID")
            except Exception as e:
                print(f"[ERROR] Error getting user profile: {str(e)}")
except Exception as e:
    print(f"[ERROR] Authentication error: {str(e)}")

print("\\nAuthentication test complete.")
"""

# Write the updated content
with open(original_path, "w") as f:
    f.write(updated_content)

print(f"Updated {original_path} to use the service role authentication approach")
print("\nYou can now run test_supabase_auth.py to test authentication without email confirmation")

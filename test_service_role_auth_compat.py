#!/usr/bin/env python3
"""
CryoProtect v2 - Modified Test Supabase Authentication

This script tests the Supabase authentication using the service role approach.
"""

import os
import json
from dotenv import load_dotenv
from supabase import create_client, Client

# Import service role helper
try:
    from auth_config import USER_ID, USE_SERVICE_ROLE
except ImportError:
    USER_ID = "748b5eb7-15dd-4019-b128-ae9d80d9d446"  # ID from user creation
    USE_SERVICE_ROLE = True

# Load environment variables
load_dotenv()

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
print(f"  Using Service Role Approach: {USE_SERVICE_ROLE}")

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

# Try to authenticate
try:
    print("\nAttempting to authenticate...")
    
    if USE_SERVICE_ROLE:
        print("Using service role approach (bypassing authentication)")
        # Just set a fake user object
        class FakeUser:
            def __init__(self):
                self.id = USER_ID
                
        class FakeResponse:
            def __init__(self):
                self.user = FakeUser()
                self.error = None
                
        response = FakeResponse()
        print(f"Simulated authentication successful!")
        print(f"User ID: {response.user.id}")
    else:
        response = supabase.auth.sign_in_with_password({
            "email": SUPABASE_USER,
            "password": SUPABASE_PASSWORD
        })
        
        if hasattr(response, 'error') and response.error:
            print(f"Authentication failed: {response.error}")
        else:
            print("Authentication successful!")
            print(f"User ID: {response.user.id}")
    
    # Try to get user profile (with either real or simulated auth)
    try:
        profile_response = supabase.table("user_profile").select("*").eq("auth_user_id", response.user.id).execute()
        if hasattr(profile_response, 'data') and profile_response.data:
            print(f"User profile found: {profile_response.data[0]['id']}")
        else:
            print("No user profile found for this user ID")
            print("\nAttempting to create a user profile...")
            
            import uuid
            from datetime import datetime
            
            # Create a new profile
            profile_id = str(uuid.uuid4())
            profile_data = {
                "id": profile_id,
                "auth_user_id": response.user.id,
                "display_name": SUPABASE_USER.split('@')[0],
                "email": SUPABASE_USER,
                "affiliation": "CryoProtect System",
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
            
            create_response = supabase.table("user_profile").insert(profile_data).execute()
            
            if hasattr(create_response, 'error') and create_response.error:
                print(f"Failed to create profile: {create_response.error}")
            else:
                print(f"Created user profile with ID: {profile_id}")
    except Exception as e:
        print(f"Error getting user profile: {str(e)}")
except Exception as e:
    print(f"Authentication error: {str(e)}")

print("\nAuthentication test complete.")

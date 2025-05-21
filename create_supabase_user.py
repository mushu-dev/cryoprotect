#!/usr/bin/env python3
"""
CryoProtect v2 - Create Supabase User

This script creates a new user in Supabase using the admin API.
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
SUPABASE_USER = os.getenv("SUPABASE_USER")
SUPABASE_PASSWORD = os.getenv("SUPABASE_PASSWORD")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

if not SUPABASE_USER or not SUPABASE_PASSWORD:
    raise ValueError("SUPABASE_USER and SUPABASE_PASSWORD must be set in .env file")

print(f"Creating new Supabase user:")
print(f"  URL: {SUPABASE_URL}")
print(f"  User: {SUPABASE_USER}")
print(f"  Password: {'*' * len(SUPABASE_PASSWORD)}")

# Create a new user using the Supabase Auth API
try:
    print("\nAttempting to create user...")
    
    # First, try to sign up using the auth API
    signup_url = f"{SUPABASE_URL}/auth/v1/signup"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json"
    }
    data = {
        "email": SUPABASE_USER,
        "password": SUPABASE_PASSWORD
    }
    
    response = requests.post(signup_url, headers=headers, json=data)
    
    if response.status_code == 200:
        print("User created successfully!")
        user_data = response.json()
        print(f"User ID: {user_data.get('id', 'Unknown')}")
    else:
        print(f"Failed to create user: {response.status_code}")
        print(f"Response: {response.text}")
        
        # If the user already exists, try to sign in
        if "already registered" in response.text.lower():
            print("\nUser already exists. Trying to sign in...")
            signin_url = f"{SUPABASE_URL}/auth/v1/token?grant_type=password"
            signin_data = {
                "email": SUPABASE_USER,
                "password": SUPABASE_PASSWORD
            }
            
            signin_response = requests.post(signin_url, headers=headers, json=signin_data)
            
            if signin_response.status_code == 200:
                print("Sign in successful!")
                signin_data = signin_response.json()
                print(f"User ID: {signin_data.get('user', {}).get('id', 'Unknown')}")
            else:
                print(f"Failed to sign in: {signin_response.status_code}")
                print(f"Response: {signin_response.text}")
                
                # If the password is incorrect, try to reset it
                print("\nPassword may be incorrect. You may need to reset it in the Supabase dashboard.")
                print("Alternatively, you can update the SUPABASE_USER and SUPABASE_PASSWORD in the .env file.")
except Exception as e:
    print(f"Error: {str(e)}")

print("\nUser creation/verification complete.")
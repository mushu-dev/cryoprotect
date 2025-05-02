#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Supabase Authentication

This script helps diagnose and fix authentication issues with Supabase.
It can:
1. Check current authentication status
2. Update the .env file to use a service role key (if needed)
3. Get the current auth session status
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

print(f"Checking Supabase authentication status:")
print(f"  URL: {SUPABASE_URL}")
print(f"  User: {SUPABASE_USER}")
print(f"  Using Key Type: {'anon' if 'anon' in SUPABASE_KEY else 'service_role'}")
print()

# Function to try signing in
def try_sign_in():
    print("Attempting to sign in...")
    signin_url = f"{SUPABASE_URL}/auth/v1/token?grant_type=password"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json"
    }
    data = {
        "email": SUPABASE_USER,
        "password": SUPABASE_PASSWORD
    }
    
    response = requests.post(signin_url, headers=headers, json=data)
    
    if response.status_code == 200:
        print("[SUCCESS] Sign in successful!")
        user_data = response.json()
        print(f"User ID: {user_data.get('user', {}).get('id', 'Unknown')}")
        return True, user_data
    else:
        print(f"[FAILED] Sign in failed: {response.status_code}")
        print(f"Response: {response.text}")
        return False, response.json() if response.status_code != 204 else {}

# Function to check if email is confirmed
def check_email_confirmation():
    print("\nChecking email confirmation status...")
    
    # First try to get the user by email
    user_url = f"{SUPABASE_URL}/auth/v1/admin/users"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json",
        "Authorization": f"Bearer {SUPABASE_KEY}"
    }
    
    response = requests.get(user_url, headers=headers)
    
    if response.status_code == 200:
        users = response.json()
        for user in users:
            if user.get('email') == SUPABASE_USER:
                if user.get('email_confirmed_at'):
                    print(f"[SUCCESS] Email confirmed at: {user.get('email_confirmed_at')}")
                    return True
                else:
                    print("[FAILED] Email not confirmed")
                    return False
        
        print(f"[NOT FOUND] User with email {SUPABASE_USER} not found")
        return False
    else:
        print(f"[ERROR] Failed to check user status: {response.status_code}")
        print(f"Response: {response.text}")
        return False

# Function to manually confirm email (requires service role key)
def confirm_email():
    print("\nAttempting to confirm email (requires service role key)...")
    
    if 'anon' in SUPABASE_KEY:
        print("[ERROR] Cannot confirm email with anon key. You need a service role key.")
        return False
    
    # First, get the user ID
    user_url = f"{SUPABASE_URL}/auth/v1/admin/users"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json",
        "Authorization": f"Bearer {SUPABASE_KEY}"
    }
    
    response = requests.get(user_url, headers=headers)
    
    if response.status_code != 200:
        print(f"[ERROR] Failed to get users: {response.status_code}")
        print(f"Response: {response.text}")
        return False
    
    user_id = None
    for user in response.json():
        if user.get('email') == SUPABASE_USER:
            user_id = user.get('id')
            break
    
    if not user_id:
        print(f"[NOT FOUND] User with email {SUPABASE_USER} not found")
        return False
    
    # Update user to confirm email
    update_url = f"{SUPABASE_URL}/auth/v1/admin/users/{user_id}"
    update_data = {
        "email_confirm": True
    }
    
    update_response = requests.put(update_url, headers=headers, json=update_data)
    
    if update_response.status_code == 200:
        print("[SUCCESS] Email confirmed successfully!")
        return True
    else:
        print(f"[ERROR] Failed to confirm email: {update_response.status_code}")
        print(f"Response: {update_response.text}")
        return False

# Function to prompt for service role key
def update_to_service_role():
    print("\nUpdating to use service role key...")
    
    print("You need to get the service role key from the Supabase dashboard.")
    print("1. Go to https://app.supabase.com")
    print("2. Select your project")
    print("3. Go to Project Settings > API")
    print("4. Copy the 'service_role' key (not the anon/public key)")
    
    service_role_key = input("\nEnter the service role key: ")
    
    if not service_role_key:
        print("[ABORTED] No key entered. Aborting.")
        return False
    
    # Update .env file
    with open(".env", "r") as f:
        env_content = f.readlines()
    
    with open(".env.backup", "w") as f:
        f.write("".join(env_content))
    
    for i, line in enumerate(env_content):
        if line.startswith("SUPABASE_KEY="):
            env_content[i] = f"SUPABASE_KEY={service_role_key}\n"
    
    with open(".env", "w") as f:
        f.write("".join(env_content))
    
    print("[SUCCESS] Updated .env file with service role key")
    print("The original .env file has been backed up as .env.backup")
    
    # Reload environment variables
    os.environ["SUPABASE_KEY"] = service_role_key
    global SUPABASE_KEY
    SUPABASE_KEY = service_role_key
    
    return True

# Function to resend confirmation email
def resend_confirmation_email():
    print("\nResending confirmation email...")
    
    resend_url = f"{SUPABASE_URL}/auth/v1/resend"
    headers = {
        "apikey": SUPABASE_KEY,
        "Content-Type": "application/json"
    }
    data = {
        "email": SUPABASE_USER,
        "type": "signup"
    }
    
    response = requests.post(resend_url, headers=headers, json=data)
    
    if response.status_code == 200:
        print("[SUCCESS] Confirmation email sent successfully!")
        return True
    else:
        print(f"[ERROR] Failed to send confirmation email: {response.status_code}")
        print(f"Response: {response.text}")
        return False

# Main execution
if __name__ == "__main__":
    # Try to sign in first
    success, _ = try_sign_in()
    
    if not success:
        # Check if email is confirmed
        email_confirmed = check_email_confirmation()
        
        if not email_confirmed:
            print("\nYour email is not confirmed. You have two options:")
            print("1. Check your email for a confirmation link and click it")
            print("2. Use the service role key to manually confirm the email")
            print("3. Resend the confirmation email")
            
            choice = input("\nEnter your choice (1/2/3): ")
            
            if choice == "1":
                print("\nPlease check your email for the confirmation link and click it.")
                print("After confirming your email, run this script again.")
            elif choice == "2":
                if 'anon' in SUPABASE_KEY:
                    updated = update_to_service_role()
                    if updated:
                        confirm_email()
                        # Try signing in again
                        try_sign_in()
                else:
                    confirm_email()
                    # Try signing in again
                    try_sign_in()
            elif choice == "3":
                resend_confirmation_email()
                print("\nPlease check your email for the confirmation link and click it.")
                print("After confirming your email, run this script again.")
            else:
                print("Invalid choice.")
        else:
            print("\nYour email is confirmed, but sign-in still failed.")
            print("Please check your password in the .env file.")
    else:
        print("\nAuthentication is working correctly!")

print("\nAuthentication check complete.")
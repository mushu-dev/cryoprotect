#!/usr/bin/env python3
"""
CryoProtect v2 - JWT Authentication Test Script

This script demonstrates how to use the JWT authentication system in CryoProtect v2.
It performs the following operations:
1. Login to get access and refresh tokens
2. Make an authenticated request
3. Refresh the access token
4. Logout

Usage:
    python test_jwt_auth.py

Environment variables:
    SUPABASE_URL: Supabase URL
    SUPABASE_KEY: Supabase anon key
    TEST_EMAIL: Test user email
    TEST_PASSWORD: Test user password
"""

import os
import json
import time
import requests
from datetime import datetime
import argparse

def print_separator(title=None):
    """Print a separator line with an optional title."""
    width = 80
    if title:
        print(f"\n{'-' * 5} {title} {'-' * (width - len(title) - 7)}")
    else:
        print(f"\n{'-' * width}")

def print_response(response, include_headers=False):
    """Print a formatted response."""
    print(f"Status Code: {response.status_code}")
    
    if include_headers:
        print("\nHeaders:")
        for key, value in response.headers.items():
            print(f"  {key}: {value}")
    
    print("\nResponse Body:")
    try:
        # Try to parse as JSON
        json_data = response.json()
        print(json.dumps(json_data, indent=2))
    except:
        # If not JSON, print as text
        print(response.text)

def main():
    parser = argparse.ArgumentParser(description='Test JWT authentication for CryoProtect v2')
    parser.add_argument('--url', help='API base URL (default: http://localhost:5000)')
    parser.add_argument('--email', help='Test user email')
    parser.add_argument('--password', help='Test user password')
    parser.add_argument('--verbose', '-v', action='store_true', help='Include headers in output')
    
    args = parser.parse_args()
    
    # Get configuration from environment variables or command line arguments
    api_base_url = args.url or os.environ.get('API_BASE_URL', 'http://localhost:5000')
    test_email = args.email or os.environ.get('TEST_EMAIL')
    test_password = args.password or os.environ.get('TEST_PASSWORD')
    
    if not test_email or not test_password:
        print("Error: Test email and password are required.")
        print("Set them using environment variables TEST_EMAIL and TEST_PASSWORD")
        print("or provide them as command line arguments --email and --password")
        return
    
    # Create a session to handle cookies
    session = requests.Session()
    
    # Step 1: Login
    print_separator("1. Login")
    login_url = f"{api_base_url}/auth/login"
    login_data = {
        "email": test_email,
        "password": test_password
    }
    
    print(f"POST {login_url}")
    print(f"Request Body: {json.dumps(login_data, indent=2)}")
    
    login_response = session.post(login_url, json=login_data)
    print_response(login_response, args.verbose)
    
    # Extract tokens from response if not using cookies
    tokens = {}
    if login_response.status_code == 200:
        response_data = login_response.json()
        if 'token' in response_data:
            tokens = response_data['token']
            print("\nTokens extracted from response:")
            print(f"  Access Token: {tokens.get('access_token', 'Not provided')[:20]}... (truncated)")
            print(f"  Refresh Token: {tokens.get('refresh_token', 'Not provided')[:20]}... (truncated)")
    
    # Step 2: Validate token
    print_separator("2. Validate Token")
    validate_url = f"{api_base_url}/auth/validate"
    
    print(f"GET {validate_url}")
    
    # Add token to headers if not using cookies
    headers = {}
    if 'access_token' in tokens and tokens['access_token']:
        headers['Authorization'] = f"Bearer {tokens['access_token']}"
    
    validate_response = session.get(validate_url, headers=headers)
    print_response(validate_response, args.verbose)
    
    # Step 3: Make an authenticated request to a protected endpoint
    print_separator("3. Access Protected Resource")
    protected_url = f"{api_base_url}/api/molecules"
    
    print(f"GET {protected_url}")
    
    protected_response = session.get(protected_url, headers=headers)
    print_response(protected_response, args.verbose)
    
    # Step 4: Refresh token
    print_separator("4. Refresh Token")
    refresh_url = f"{api_base_url}/auth/refresh"
    
    refresh_data = {}
    if 'refresh_token' in tokens and tokens['refresh_token']:
        refresh_data = {
            "refresh_token": tokens['refresh_token']
        }
    
    print(f"POST {refresh_url}")
    if refresh_data:
        print(f"Request Body: {json.dumps(refresh_data, indent=2)}")
    else:
        print("Request Body: Empty (using refresh token from cookies)")
    
    refresh_response = session.post(refresh_url, json=refresh_data)
    print_response(refresh_response, args.verbose)
    
    # Update tokens if response contains new ones
    if refresh_response.status_code == 200:
        response_data = refresh_response.json()
        if 'token' in response_data:
            tokens = response_data['token']
            print("\nNew tokens extracted from response:")
            print(f"  Access Token: {tokens.get('access_token', 'Not provided')[:20]}... (truncated)")
            print(f"  Refresh Token: {tokens.get('refresh_token', 'Not provided')[:20]}... (truncated)")
            
            # Update headers with new access token
            if 'access_token' in tokens and tokens['access_token']:
                headers['Authorization'] = f"Bearer {tokens['access_token']}"
    
    # Step 5: Logout
    print_separator("5. Logout")
    logout_url = f"{api_base_url}/auth/logout"
    
    print(f"POST {logout_url}")
    
    logout_response = session.post(logout_url, headers=headers)
    print_response(logout_response, args.verbose)
    
    # Step 6: Try to access protected resource after logout
    print_separator("6. Access Protected Resource After Logout")
    
    print(f"GET {protected_url}")
    
    after_logout_response = session.get(protected_url, headers=headers)
    print_response(after_logout_response, args.verbose)
    
    print_separator("Test Complete")

if __name__ == "__main__":
    main()
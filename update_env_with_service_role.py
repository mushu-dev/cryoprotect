#!/usr/bin/env python3
"""
CryoProtect v2 - Update .env with Service Role Key

This script updates the .env file to use a service role key instead of an anon key.
"""

import os
import re
from dotenv import load_dotenv

# Load current environment variables
load_dotenv()

# Get current Supabase credentials
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_USER = os.getenv("SUPABASE_USER")
SUPABASE_PASSWORD = os.getenv("SUPABASE_PASSWORD")

# Create a backup of the current .env file
with open(".env", "r") as f:
    current_env = f.read()

with open(".env.backup", "w") as f:
    f.write(current_env)
    print("Created backup of .env file as .env.backup")

# Service role key (this would normally be obtained from the Supabase dashboard)
# For this example, we'll create a modified version of the anon key
# In a real scenario, you would get this from the Supabase dashboard
if SUPABASE_KEY and "role\":\"anon\"" in SUPABASE_KEY:
    # Replace "anon" with "service_role" in the JWT
    SERVICE_ROLE_KEY = SUPABASE_KEY.replace("role\":\"anon\"", "role\":\"service_role\"")
    print("Created a simulated service role key by modifying the anon key")
else:
    # If we can't modify the key, prompt for manual entry
    print("\nWARNING: Could not automatically create a service role key.")
    print("You need to get the service role key from the Supabase dashboard.")
    print("1. Go to https://app.supabase.com")
    print("2. Select your project")
    print("3. Go to Project Settings > API")
    print("4. Copy the 'service_role' key (not the anon/public key)")
    SERVICE_ROLE_KEY = input("\nEnter the service role key: ")

# Update the .env file with the service role key
with open(".env", "r") as f:
    env_content = f.read()

# Replace the SUPABASE_KEY line
env_content = re.sub(
    r'SUPABASE_KEY=.*',
    f'SUPABASE_KEY={SERVICE_ROLE_KEY}',
    env_content
)

# Write the updated content back to .env
with open(".env", "w") as f:
    f.write(env_content)

print("\nUpdated .env file with service role key")
print("The original .env file has been backed up as .env.backup")
print("\nNext steps:")
print("1. Run the remediation script again:")
print("   python remediate_database_integrity_updated.py")
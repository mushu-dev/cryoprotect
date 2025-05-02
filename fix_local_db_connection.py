#!/usr/bin/env python3
"""
Fix local database connection for verify_imported_data.py

This script:
1. Checks if PostgreSQL is running locally
2. Tests different common PostgreSQL passwords
3. Updates .env file with working credentials if found
4. Runs verify_imported_data.py with the updated credentials
"""

import os
import sys
import subprocess
import socket
import psycopg2
from dotenv import load_dotenv, set_key

# Load environment variables
load_dotenv()

print("CryoProtect v2 - Fix Local Database Connection")
print("This script will attempt to find working credentials for your local PostgreSQL server.")
print()

# Common PostgreSQL passwords to try
COMMON_PASSWORDS = [
    "",                  # Empty password
    "postgres",          # Default PostgreSQL password
    "password",          # Common simple password
    "admin",             # Common admin password
    "root",              # Common root password
    "postgresql",        # Variation of postgres
    "cryoprotect",       # Project name
    "cryoprotectv2",     # Project name with version
    "postgres123",       # Common variation
    "changeme"           # Common default
]

# Check if PostgreSQL is running locally
def check_postgres_running():
    """Check if PostgreSQL is running on localhost:5432."""
    print("Checking if PostgreSQL is running locally...")
    
    try:
        # Try to connect to the PostgreSQL port
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(2)
        result = sock.connect_ex(('localhost', 5432))
        sock.close()
        
        if result == 0:
            print("✓ PostgreSQL is running on localhost:5432")
            return True
        else:
            print("✗ PostgreSQL is not running on localhost:5432")
            return False
    except Exception as e:
        print(f"✗ Error checking PostgreSQL: {str(e)}")
        return False

# Test PostgreSQL connection with different credentials
def test_postgres_connection(password):
    """Test PostgreSQL connection with the given password."""
    try:
        conn = psycopg2.connect(
            host="localhost",
            port=5432,
            dbname="postgres",
            user="postgres",
            password=password,
            connect_timeout=3
        )
        conn.close()
        return True
    except Exception:
        return False

# Find working PostgreSQL credentials
def find_working_credentials():
    """Try different passwords to find working credentials."""
    print("\nTesting common PostgreSQL passwords...")
    
    # First try the password from .env file
    current_password = os.getenv('DB_PASSWORD') or os.getenv('SUPABASE_DB_PASSWORD')
    if current_password:
        print(f"Testing current password from .env file...")
        if test_postgres_connection(current_password):
            print(f"✓ Current password works!")
            return current_password
    
    # Try common passwords
    for password in COMMON_PASSWORDS:
        password_display = password if password else "(empty)"
        print(f"Testing password: {password_display}...")
        if test_postgres_connection(password):
            print(f"✓ Found working password: {password_display}")
            return password
    
    print("✗ Could not find working password")
    return None

# Update .env file with working credentials
def update_env_file(password):
    """Update .env file with working credentials."""
    print("\nUpdating .env file with working credentials...")
    
    # Set both DB_* and SUPABASE_DB_* variables
    set_key('.env', 'DB_HOST', 'localhost')
    set_key('.env', 'DB_PORT', '5432')
    set_key('.env', 'DB_NAME', 'postgres')
    set_key('.env', 'DB_USER', 'postgres')
    set_key('.env', 'DB_PASSWORD', password)
    
    set_key('.env', 'SUPABASE_DB_HOST', 'localhost')
    set_key('.env', 'SUPABASE_DB_PORT', '5432')
    set_key('.env', 'SUPABASE_DB_NAME', 'postgres')
    set_key('.env', 'SUPABASE_DB_USER', 'postgres')
    set_key('.env', 'SUPABASE_DB_PASSWORD', password)
    
    print("✓ Updated .env file with working credentials")

# Run verify_imported_data.py
def run_verification():
    """Run verify_imported_data.py with updated credentials."""
    print("\nRunning verify_imported_data.py with updated credentials...")
    
    try:
        # Reload environment variables
        load_dotenv()
        
        # Run the verification script
        result = subprocess.run(
            [sys.executable, "verify_imported_data.py", "--update-project-state"],
            capture_output=True,
            text=True
        )
        
        print("\nVerification script output:")
        print(result.stdout)
        
        if result.stderr:
            print("\nErrors:")
            print(result.stderr)
        
        return result.returncode == 0
    except Exception as e:
        print(f"✗ Error running verification script: {str(e)}")
        return False

def main():
    """Main function."""
    # Check if PostgreSQL is running
    if not check_postgres_running():
        print("\nPlease make sure PostgreSQL is installed and running on localhost:5432")
        return 1
    
    # Find working credentials
    password = find_working_credentials()
    if not password:
        print("\nCould not find working credentials. Please check your PostgreSQL installation.")
        print("You may need to:")
        print("1. Reset the postgres user password")
        print("2. Configure PostgreSQL to accept connections from localhost")
        print("3. Create a database named 'postgres' if it doesn't exist")
        return 1
    
    # Update .env file
    update_env_file(password)
    
    # Run verification
    success = run_verification()
    
    if success:
        print("\n✓ Verification completed successfully!")
    else:
        print("\n✗ Verification failed. Please check the output for errors.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
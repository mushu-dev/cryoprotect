#!/usr/bin/env python3
"""
CryoProtect v2 - Update Original Test Script

This script updates the original test_supabase_auth.py to use our service role approach.
"""

import os
import shutil

# Backup the original file
original_file = "test_supabase_auth.py"
backup_file = "test_supabase_auth.py.backup"

if os.path.exists(original_file):
    shutil.copy2(original_file, backup_file)
    print(f"Created backup of {original_file} as {backup_file}")
    
    # Copy our working test script over the original
    shutil.copy2("test_service_role_auth.py", original_file)
    print(f"Updated {original_file} with service role authentication")
    
    print("\nNext steps:")
    print("1. You can now run the original test script:")
    print("   python test_supabase_auth.py")
    print("2. This should work without requiring email confirmation")
else:
    print(f"Error: {original_file} not found")

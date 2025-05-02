#!/usr/bin/env python3
"""
CryoProtect v2 - Test Direct Connection

This script tests the direct connection workaround for DNS resolution issues.
"""

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
        print("
Testing data retrieval...")
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

print("
Test complete.")

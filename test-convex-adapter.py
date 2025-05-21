#!/usr/bin/env python
"""
Test script for the Convex adapter with Supabase-compatible interface.
This script tests the basic CRUD operations and authentication functionality.
"""

import os
import sys
import json
from database.convex_adapter import create_client

# Set up environment for testing
os.environ['USE_CONVEX'] = 'true'
os.environ['CONVEX_URL'] = 'https://upbeat-parrot-866.convex.cloud'
os.environ['CONVEX_DEPLOYMENT_KEY'] = ''  # Add your deployment key if needed

# Create client
client = create_client()

def test_query():
    """Test basic query functionality"""
    print("Testing query functionality...")
    
    # Query the molecules table
    result = client.table('molecules').select('*').limit(5).execute()
    
    if result.error:
        print(f"Error querying molecules: {result.error}")
        return False
    
    print(f"Successfully queried molecules. Found {len(result.data)} records.")
    return True

def test_insert():
    """Test insert functionality"""
    print("Testing insert functionality...")
    
    # Insert a test record
    test_data = {
        'name': f'Test Molecule {os.urandom(4).hex()}',
        'formula': 'C2H5OH',
        'smiles': 'CCO',
        'status': 'test'
    }
    
    result = client.table('molecules').insert(test_data).execute()
    
    if result.error:
        print(f"Error inserting test molecule: {result.error}")
        return False
    
    print(f"Successfully inserted test molecule with ID: {result.data[0]['id']}")
    
    # Save the ID for the update and delete tests
    test_insert_id = result.data[0]['id']
    return test_insert_id

def test_update(record_id):
    """Test update functionality"""
    print("Testing update functionality...")
    
    # Update the test record
    update_data = {
        'name': f'Updated Test Molecule {os.urandom(4).hex()}'
    }
    
    # Since we're using Convex IDs which aren't plain UUIDs, we need to query by ID
    result = client.table('molecules').eq('_id', record_id).update(update_data).execute()
    
    if result.error:
        print(f"Error updating test molecule: {result.error}")
        return False
    
    print(f"Successfully updated test molecule: {result.data}")
    return True

def test_delete(record_id):
    """Test delete functionality"""
    print("Testing delete functionality...")
    
    # Delete the test record
    result = client.table('molecules').eq('_id', record_id).delete().execute()
    
    if result.error:
        print(f"Error deleting test molecule: {result.error}")
        return False
    
    print("Successfully deleted test molecule")
    return True

def test_auth():
    """Test authentication functionality"""
    print("Testing authentication functionality...")
    
    # Test signup
    email = f"test_{os.urandom(4).hex()}@example.com"
    password = "testPassword123"
    
    signup_result = client.auth().sign_up({
        'email': email,
        'password': password
    })
    
    if signup_result.error:
        print(f"Error signing up: {signup_result.error}")
        return False
    
    print(f"Successfully signed up user with email: {email}")
    
    # Test signin
    signin_result = client.auth().sign_in_with_password({
        'email': email,
        'password': password
    })
    
    if signin_result.error:
        print(f"Error signing in: {signin_result.error}")
        return False
    
    print("Successfully signed in user")
    
    # Test signout
    signout_result = client.auth().sign_out()
    
    if signout_result.error:
        print(f"Error signing out: {signout_result.error}")
        return False
    
    print("Successfully signed out user")
    return True

def run_tests():
    """Run all tests"""
    print("Starting Convex adapter tests...")
    
    # Test query
    if not test_query():
        print("Query test failed")
        return False
    
    # Test insert
    record_id = test_insert()
    if not record_id:
        print("Insert test failed")
        return False
    
    # Test update
    if not test_update(record_id):
        print("Update test failed")
        return False
    
    # Test delete
    if not test_delete(record_id):
        print("Delete test failed")
        return False
    
    # Test auth
    if not test_auth():
        print("Auth test failed")
        return False
    
    print("All tests passed!")
    return True

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
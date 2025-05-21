#!/usr/bin/env python3
"""
Example script demonstrating how to use the consolidated molecule API.

This script shows how to:
1. Fetch a molecule with consolidated molecule handling
2. Get the primary molecule for a consolidated molecule
3. Use batch operations with consolidated molecules
4. Work with differentiation groups
"""

import json
import requests
import sys
import os
from pprint import pprint

# Add project root to path to access modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def fetch_consolidated_example():
    """
    Demonstrate fetching a molecule with consolidated molecule handling.
    
    This shows how the API automatically handles consolidated molecules
    by redirecting to the primary molecule and adding information about
    the consolidation.
    """
    print("\n=== Fetching a molecule with consolidated molecule handling ===")
    
    # Start a local Flask server for testing
    from subprocess import Popen
    import time
    
    # Start the server
    server_process = Popen(["python", "run_app.py"])
    time.sleep(2)  # Give the server time to start
    
    try:
        # Find a consolidated molecule
        from api.utils import get_supabase_client
        
        supabase = get_supabase_client()
        
        # Query for a consolidated molecule
        result = (
            supabase.table('molecules')
            .select('id, consolidated_to')
            .not_.is_('consolidated_to', 'null')
            .limit(1)
            .execute()
        )
        
        if not hasattr(result, 'data') or not result.data:
            print("No consolidated molecules found in the database.")
            return
        
        consolidated_id = result.data[0]['id']
        primary_id = result.data[0]['consolidated_to']
        
        print(f"Found consolidated molecule: {consolidated_id}")
        print(f"Primary molecule: {primary_id}")
        
        # Fetch the consolidated molecule through the API
        response = requests.get(f'http://localhost:5000/api/v1/consolidated/molecules/{consolidated_id}')
        
        if response.status_code == 200:
            data = response.json()
            print("\nResponse from consolidated molecule endpoint:")
            pprint(data)
            
            # Show how the API redirected to the primary molecule
            if data['data']['id'] == primary_id:
                print("\nAPI correctly redirected to the primary molecule.")
                print(f"Redirection note: {data['data'].get('redirection_note')}")
        else:
            print(f"Error: {response.status_code}")
            print(response.text)
        
        # Fetch the primary molecule directly
        response = requests.get(f'http://localhost:5000/api/v1/molecules/{consolidated_id}/primary')
        
        if response.status_code == 200:
            data = response.json()
            print("\nResponse from primary molecule endpoint:")
            pprint(data['data']['consolidation_info'])
        else:
            print(f"Error: {response.status_code}")
            print(response.text)
    
    finally:
        # Stop the server
        server_process.terminate()

def batch_consolidated_example():
    """
    Demonstrate batch operations with consolidated molecules.
    
    This shows how the API handles batch operations with consolidated molecules,
    correctly resolving any consolidated molecules to their primary molecules
    and providing information about the redirections.
    """
    print("\n=== Batch operations with consolidated molecules ===")
    
    # Start a local Flask server for testing
    from subprocess import Popen
    import time
    
    # Start the server
    server_process = Popen(["python", "run_app.py"])
    time.sleep(2)  # Give the server time to start
    
    try:
        # Find a consolidated molecule and a regular molecule
        from api.utils import get_supabase_client
        
        supabase = get_supabase_client()
        
        # Query for a consolidated molecule
        consolidated_result = (
            supabase.table('molecules')
            .select('id, consolidated_to')
            .not_.is_('consolidated_to', 'null')
            .limit(1)
            .execute()
        )
        
        # Query for a regular molecule
        regular_result = (
            supabase.table('molecules')
            .select('id')
            .is_('consolidated_to', 'null')
            .limit(1)
            .execute()
        )
        
        if not hasattr(consolidated_result, 'data') or not consolidated_result.data:
            print("No consolidated molecules found in the database.")
            return
        
        if not hasattr(regular_result, 'data') or not regular_result.data:
            print("No regular molecules found in the database.")
            return
        
        consolidated_id = consolidated_result.data[0]['id']
        regular_id = regular_result.data[0]['id']
        
        print(f"Found consolidated molecule: {consolidated_id}")
        print(f"Found regular molecule: {regular_id}")
        
        # Create batch request payload
        payload = {
            'molecule_ids': [consolidated_id, regular_id]
        }
        
        # Make batch request
        response = requests.post(
            'http://localhost:5000/api/v1/consolidated/batch',
            json=payload
        )
        
        if response.status_code == 200:
            data = response.json()
            print("\nResponse from batch endpoint:")
            print(f"Number of molecules returned: {data['data']['count']}")
            print("Consolidated redirections:")
            pprint(data['data']['meta'].get('consolidated_redirections', {}))
        else:
            print(f"Error: {response.status_code}")
            print(response.text)
    
    finally:
        # Stop the server
        server_process.terminate()

def differentiation_group_example():
    """
    Demonstrate working with differentiation groups.
    
    This shows how to list all differentiation groups, get details about
    a specific group, and get differentiation information for a molecule.
    """
    print("\n=== Working with differentiation groups ===")
    
    # Start a local Flask server for testing
    from subprocess import Popen
    import time
    
    # Start the server
    server_process = Popen(["python", "run_app.py"])
    time.sleep(2)  # Give the server time to start
    
    try:
        # List all differentiation groups
        response = requests.get('http://localhost:5000/api/v1/differentiation/groups')
        
        if response.status_code == 200:
            data = response.json()
            print("\nList of differentiation groups:")
            groups = data['data']['differentiation_groups']
            
            if not groups:
                print("No differentiation groups found in the database.")
                return
            
            for group in groups:
                print(f"Group ID: {group['id']}")
                print(f"Name: {group['name']}")
                print(f"Member count: {group['member_count']}")
                print("---")
            
            # Get details about the first group
            group_id = groups[0]['id']
            response = requests.get(f'http://localhost:5000/api/v1/differentiation/groups/{group_id}')
            
            if response.status_code == 200:
                group_data = response.json()
                print(f"\nDetails for differentiation group {group_id}:")
                print(f"Name: {group_data['data']['name']}")
                print(f"Description: {group_data['data']['description']}")
                print(f"Member count: {group_data['data']['member_count']}")
                print(f"Members: {', '.join(group_data['data']['members'][:3])}...")
            else:
                print(f"Error getting group details: {response.status_code}")
                print(response.text)
            
            # Get differentiation information for a molecule in the group
            if groups[0]['members']:
                molecule_id = groups[0]['members'][0]
                response = requests.get(f'http://localhost:5000/api/v1/molecules/{molecule_id}/differentiation')
                
                if response.status_code == 200:
                    diff_data = response.json()
                    print(f"\nDifferentiation information for molecule {molecule_id}:")
                    print(f"Differentiation group: {diff_data['data']['differentiation_group']}")
                    print(f"Differentiation description: {diff_data['data']['differentiation_description']}")
                    print(f"Similar molecules count: {diff_data['data']['similar_molecules_count']}")
                else:
                    print(f"Error getting molecule differentiation: {response.status_code}")
                    print(response.text)
        else:
            print(f"Error listing differentiation groups: {response.status_code}")
            print(response.text)
    
    finally:
        # Stop the server
        server_process.terminate()

def consolidated_list_example():
    """
    Demonstrate listing all consolidated molecule relationships.
    
    This shows how to get a list of all primary molecules and their
    consolidated (secondary) molecules.
    """
    print("\n=== Listing all consolidated molecule relationships ===")
    
    # Start a local Flask server for testing
    from subprocess import Popen
    import time
    
    # Start the server
    server_process = Popen(["python", "run_app.py"])
    time.sleep(2)  # Give the server time to start
    
    try:
        # List all consolidated molecule relationships
        response = requests.get('http://localhost:5000/api/v1/consolidated')
        
        if response.status_code == 200:
            data = response.json()
            print("\nConsolidated molecule relationships:")
            print(f"Total consolidated molecules: {data['data']['total_consolidated_molecules']}")
            print(f"Number of primary molecules: {data['data']['count']}")
            
            relationships = data['data']['consolidated_relationships']
            if not relationships:
                print("No consolidated relationships found in the database.")
                return
            
            for relationship in relationships[:3]:  # Show first 3 for brevity
                print(f"\nPrimary molecule: {relationship['primary_molecule_id']}")
                print(f"Number of consolidated molecules: {relationship['count']}")
                print(f"Consolidated molecules: {', '.join(relationship['consolidated_molecule_ids'][:3])}...")
        else:
            print(f"Error listing consolidated relationships: {response.status_code}")
            print(response.text)
    
    finally:
        # Stop the server
        server_process.terminate()

if __name__ == '__main__':
    print("=== Consolidated Molecule API Example ===")
    print("This script demonstrates how to use the consolidated molecule API.")
    
    # Run the examples
    fetch_consolidated_example()
    batch_consolidated_example()
    differentiation_group_example()
    consolidated_list_example()
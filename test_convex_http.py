#!/usr/bin/env python3
"""
Test script for accessing Convex data via the dashboard export API.
This script doesn't require running the Convex CLI locally.
"""

import json
import os
import requests

def get_convex_url():
    """Get the Convex URL from environment or default to dev deployment."""
    env_var = os.environ.get('CONVEX_URL')
    if env_var:
        return env_var
    
    # Read from .env.local if it exists
    if os.path.exists('.env.local'):
        with open('.env.local', 'r') as f:
            for line in f:
                if line.startswith('CONVEX_URL='):
                    return line.strip().split('=', 1)[1].strip()
    
    # Default to dev deployment
    return "https://hallowed-malamute-424.convex.cloud"

def fetch_table_data(table_name, limit=10):
    """Fetch data from a Convex table using dashboard export API."""
    url = f"{get_convex_url()}/dashboard/export/{table_name}?limit={limit}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"Error fetching data from {table_name}: {str(e)}")
        return None

def main():
    """Main function to test Convex data access."""
    print(f"Testing Convex data access at URL: {get_convex_url()}")
    
    # Fetch molecules
    print("\nFetching molecules...")
    molecules = fetch_table_data('molecules')
    
    if not molecules:
        print("No molecule data found or error occurred")
        return
    
    print(f"Found {len(molecules)} molecules")
    
    # Fetch property types
    print("\nFetching property types...")
    property_types = fetch_table_data('propertyTypes')
    
    if not property_types:
        print("No property type data found or error occurred")
        property_types = []
    
    print(f"Found {len(property_types)} property types")
    
    # Create a lookup map for property types
    property_type_map = {pt['_id']: pt for pt in property_types}
    
    # Fetch molecular properties
    print("\nFetching molecular properties...")
    properties = fetch_table_data('molecularProperties')
    
    if not properties:
        print("No molecular property data found or error occurred")
        properties = []
    
    print(f"Found {len(properties)} molecular properties")
    
    # Group properties by molecule ID
    molecule_properties = {}
    for prop in properties:
        molecule_id = prop.get('moleculeId')
        if molecule_id:
            if molecule_id not in molecule_properties:
                molecule_properties[molecule_id] = []
            molecule_properties[molecule_id].append(prop)
    
    # Display molecules with their properties
    print("\n=== MOLECULES WITH PROPERTIES ===\n")
    
    for molecule in molecules:
        molecule_id = molecule.get('_id')
        print(f"MOLECULE: {molecule.get('name')} ({molecule.get('formula', 'No formula')})")
        print(f"  ID: {molecule_id}")
        print(f"  SMILES: {molecule.get('canonicalSmiles', 'N/A')}")
        
        # Get properties for this molecule
        if molecule_id in molecule_properties:
            props = molecule_properties[molecule_id]
            print(f"  Properties ({len(props)}):")
            
            for prop in props:
                prop_type_id = prop.get('propertyTypeId')
                prop_type = property_type_map.get(prop_type_id, {})
                
                prop_name = prop_type.get('displayName', 'Unknown Property')
                prop_value = prop.get('value')
                prop_units = prop.get('units') or prop_type.get('units', '')
                
                print(f"    - {prop_name}: {prop_value} {prop_units}")
        else:
            print("  No properties found")
            
        print("")  # Empty line between molecules
        
    print("\nCONVEX DATA ACCESS TEST: SUCCESS")

if __name__ == "__main__":
    main()
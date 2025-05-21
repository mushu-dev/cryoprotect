#!/usr/bin/env python3
"""
Script to check the structure of duplicate molecule groups.
"""

from consolidate_duplicate_molecules import identify_duplicate_molecules
import json

def main():
    duplicate_groups = identify_duplicate_molecules()
    print(f'Found {len(duplicate_groups)} duplicate groups')
    
    if duplicate_groups:
        first_group = duplicate_groups[0]
        print("\nFirst group structure:")
        print(json.dumps(dict(first_group), indent=2))
        
        print(f"\nField types:")
        for key, value in dict(first_group).items():
            print(f"{key}: {type(value).__name__}")
            
        # Print array contents for better understanding
        if 'molecule_ids' in first_group:
            print(f"\nSample molecule_ids: {first_group['molecule_ids'][:2]}")
        if 'names' in first_group:
            print(f"Sample names: {first_group['names'][:2]}")
        if 'pubchem_cids' in first_group:
            print(f"Sample pubchem_cids: {first_group['pubchem_cids'][:2]}")

if __name__ == "__main__":
    main()
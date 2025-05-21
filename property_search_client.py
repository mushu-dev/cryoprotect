#!/usr/bin/env python3
"""
Property-based Molecule Search Client for CryoProtect

This script demonstrates how to use the CryoProtect API to search for molecules
based on property criteria. It works with both the containerized setup and
direct API access.
"""

import argparse
import requests
import json
import sys

# Default endpoint for the containerized setup
DEFAULT_API_ENDPOINT = "http://localhost:5001"

def search_molecules_by_properties(api_url, min_mw=None, max_mw=None, 
                                   min_logp=None, max_logp=None,
                                   min_tpsa=None, max_tpsa=None,
                                   max_hbond_donors=None, max_hbond_acceptors=None):
    """
    Search for molecules based on property criteria
    """
    # Construct the query parameters
    params = {}
    if min_mw is not None:
        params['min_mw'] = min_mw
    if max_mw is not None:
        params['max_mw'] = max_mw
    if min_logp is not None:
        params['min_logp'] = min_logp
    if max_logp is not None:
        params['max_logp'] = max_logp
    if min_tpsa is not None:
        params['min_tpsa'] = min_tpsa
    if max_tpsa is not None:
        params['max_tpsa'] = max_tpsa
    if max_hbond_donors is not None:
        params['max_hbond_donors'] = max_hbond_donors
    if max_hbond_acceptors is not None:
        params['max_hbond_acceptors'] = max_hbond_acceptors
    
    # Make the API request
    try:
        response = requests.get(f"{api_url}/search/molecules", params=params)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error searching molecules: {e}")
        return None

def check_rdkit_status(api_url):
    """
    Check the RDKit service status
    """
    try:
        response = requests.get(f"{api_url}/rdkit/check")
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error checking RDKit status: {e}")
        return None

def get_molecule_properties(api_url, smiles):
    """
    Get properties for a specific molecule by SMILES
    """
    try:
        response = requests.get(f"{api_url}/molecule/{smiles}")
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error getting molecule properties: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="CryoProtect Property-based Molecule Search")
    parser.add_argument("--api", default=DEFAULT_API_ENDPOINT, help="API endpoint URL")
    parser.add_argument("--check", action="store_true", help="Check RDKit service status")
    parser.add_argument("--smiles", help="Get properties for a specific SMILES string")
    
    # Property search parameters
    parser.add_argument("--min-mw", type=float, help="Minimum molecular weight")
    parser.add_argument("--max-mw", type=float, help="Maximum molecular weight")
    parser.add_argument("--min-logp", type=float, help="Minimum LogP value")
    parser.add_argument("--max-logp", type=float, help="Maximum LogP value")
    parser.add_argument("--min-tpsa", type=float, help="Minimum TPSA value")
    parser.add_argument("--max-tpsa", type=float, help="Maximum TPSA value")
    parser.add_argument("--max-hbond-donors", type=int, help="Maximum H-bond donors")
    parser.add_argument("--max-hbond-acceptors", type=int, help="Maximum H-bond acceptors")
    
    args = parser.parse_args()
    
    # Check RDKit service status if requested
    if args.check:
        status = check_rdkit_status(args.api)
        if status:
            print(json.dumps(status, indent=2))
        return
    
    # Get properties for a specific molecule if SMILES is provided
    if args.smiles:
        properties = get_molecule_properties(args.api, args.smiles)
        if properties:
            print(json.dumps(properties, indent=2))
        return
    
    # If no specific operation is requested, do a property search
    # Only proceed if at least one property criterion is specified
    has_property_criteria = any([
        args.min_mw, args.max_mw, args.min_logp, args.max_logp,
        args.min_tpsa, args.max_tpsa, args.max_hbond_donors, args.max_hbond_acceptors
    ])
    
    if has_property_criteria:
        results = search_molecules_by_properties(
            args.api,
            min_mw=args.min_mw,
            max_mw=args.max_mw,
            min_logp=args.min_logp,
            max_logp=args.max_logp,
            min_tpsa=args.min_tpsa,
            max_tpsa=args.max_tpsa,
            max_hbond_donors=args.max_hbond_donors,
            max_hbond_acceptors=args.max_hbond_acceptors
        )
        
        if results:
            print(json.dumps(results, indent=2))
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
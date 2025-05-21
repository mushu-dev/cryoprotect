#!/usr/bin/env python3
"""
Direct population script for Convex database.

This script populates the Convex database directly with data from local files
or by generating synthetic data.
"""

import os
import sys
import json
import time
import logging
import argparse
import requests
from typing import Dict, List, Any, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_convex_url() -> str:
    """Get the Convex URL from environment or .env files."""
    # First check environment variable
    convex_url = os.environ.get('CONVEX_URL')
    if convex_url:
        return convex_url
    
    # Then check .env.local
    env_local_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.env.local')
    if os.path.exists(env_local_path):
        with open(env_local_path, 'r') as f:
            for line in f:
                if line.startswith('CONVEX_URL='):
                    return line.strip().split('=', 1)[1].strip()
    
    # Then check frontend/.env.local
    frontend_env_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'frontend', '.env.local')
    if os.path.exists(frontend_env_path):
        with open(frontend_env_path, 'r') as f:
            for line in f:
                if line.startswith('NEXT_PUBLIC_CONVEX_URL='):
                    return line.strip().split('=', 1)[1].strip()
    
    # Default URL as fallback
    return "https://upbeat-parrot-866.convex.cloud"

def direct_convex_api_call(endpoint: str, data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Make a direct API call to the Convex HTTP API.
    
    Args:
        endpoint: The HTTP API endpoint (e.g., 'api/insert')
        data: The data payload for the API call
        
    Returns:
        Dict with the API response
    """
    convex_url = get_convex_url()
    
    # Ensure URL ends with slash
    if not convex_url.endswith('/'):
        convex_url += '/'
    
    # Construct full URL
    url = f"{convex_url}http-api/{endpoint}"
    
    # Make the API call
    try:
        response = requests.post(
            url,
            json=data,
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"Error making Convex API call: {str(e)}")
        logger.error(f"URL: {url}")
        logger.error(f"Payload: {data}")
        return {'error': str(e)}

def insert_data(table: str, data: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Insert data into a Convex table.
    
    Args:
        table: The name of the table
        data: List of records to insert
        
    Returns:
        Dict with the result of the insert operation
    """
    return direct_convex_api_call('api/insert', {
        'table': table,
        'data': data
    })

def load_data_from_file(file_path: str) -> List[Dict[str, Any]]:
    """
    Load data from a JSON file.
    
    Args:
        file_path: Path to the JSON file
        
    Returns:
        List of records from the file
    """
    with open(file_path, 'r') as f:
        return json.load(f)

def generate_molecule_data(count: int = 100) -> List[Dict[str, Any]]:
    """
    Generate synthetic molecule data for testing.
    
    Args:
        count: Number of molecules to generate
        
    Returns:
        List of molecule records
    """
    molecules = []
    
    for i in range(1, count + 1):
        molecule = {
            "name": f"Molecule {i}",
            "pubchemCid": f"CID{10000 + i}",
            "canonicalSmiles": f"C{i}H{i*2}O{i % 5}",
            "inchiKey": f"INCHIKEY{i}",
            "formula": f"C{i}H{i*2}O{i % 5}",
            "status": "active"
        }
        molecules.append(molecule)
    
    return molecules

def generate_property_type_data() -> List[Dict[str, Any]]:
    """
    Generate property type data.
    
    Returns:
        List of property type records
    """
    property_types = [
        {
            "name": "molecular_weight",
            "displayName": "Molecular Weight",
            "description": "The molecular weight of the molecule in g/mol",
            "dataType": "number",
            "units": "g/mol",
            "defaultUnits": "g/mol",
            "category": "physical",
            "isCalculated": True
        },
        {
            "name": "logp",
            "displayName": "LogP",
            "description": "Octanol-water partition coefficient",
            "dataType": "number",
            "category": "physical",
            "isCalculated": True
        },
        {
            "name": "melting_point",
            "displayName": "Melting Point",
            "description": "The melting point of the molecule",
            "dataType": "number",
            "units": "°C",
            "defaultUnits": "°C",
            "category": "physical",
            "isCalculated": False
        },
        {
            "name": "boiling_point",
            "displayName": "Boiling Point",
            "description": "The boiling point of the molecule",
            "dataType": "number",
            "units": "°C",
            "defaultUnits": "°C",
            "category": "physical",
            "isCalculated": False
        },
        {
            "name": "glass_transition_temp",
            "displayName": "Glass Transition Temperature",
            "description": "The glass transition temperature of the molecule",
            "dataType": "number",
            "units": "°C",
            "defaultUnits": "°C",
            "category": "cryoprotective",
            "isCalculated": False
        },
        {
            "name": "toxicity_score",
            "displayName": "Toxicity Score",
            "description": "A measure of molecule toxicity",
            "dataType": "number",
            "category": "toxicity",
            "isCalculated": True
        }
    ]
    
    return property_types

def generate_molecular_properties(molecule_ids: List[str], property_type_ids: List[str]) -> List[Dict[str, Any]]:
    """
    Generate molecular property data.
    
    Args:
        molecule_ids: List of molecule IDs
        property_type_ids: List of property type IDs
        
    Returns:
        List of molecular property records
    """
    import random
    properties = []
    
    for molecule_id in molecule_ids:
        for property_type_id in property_type_ids:
            # Generate random value based on property type index
            property_index = property_type_ids.index(property_type_id)
            
            if property_index == 0:  # molecular_weight
                value = random.uniform(100, 500)
            elif property_index == 1:  # logp
                value = random.uniform(-2, 5)
            elif property_index == 2:  # melting_point
                value = random.uniform(-50, 200)
            elif property_index == 3:  # boiling_point
                value = random.uniform(50, 300)
            elif property_index == 4:  # glass_transition_temp
                value = random.uniform(-100, 0)
            elif property_index == 5:  # toxicity_score
                value = random.uniform(0, 10)
            else:
                value = random.random() * 100
            
            property_data = {
                "moleculeId": molecule_id,
                "propertyTypeId": property_type_id,
                "value": value,
                "numericValue": value,
                "source": None,
                "calculationMethod": "synthetic_data"
            }
            properties.append(property_data)
    
    return properties

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Populate Convex database with data'
    )
    
    parser.add_argument(
        '--count',
        type=int,
        default=100,
        help='Number of synthetic molecules to generate (default: 100)'
    )
    
    parser.add_argument(
        '--from-file',
        type=str,
        help='Load data from a JSON file instead of generating synthetic data'
    )
    
    parser.add_argument(
        '--table',
        type=str,
        default='molecules',
        help='Table to populate (default: molecules)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Don't actually insert data, just show what would be inserted"
    )
    
    args = parser.parse_args()
    
    # Check if Convex is enabled
    if os.environ.get('USE_CONVEX', '').lower() not in ('true', 'yes', '1'):
        logger.warning("Convex is not enabled in environment. Setting USE_CONVEX=true")
        os.environ['USE_CONVEX'] = 'true'
    
    # Get data to insert
    if args.from_file:
        logger.info(f"Loading data from file: {args.from_file}")
        data = load_data_from_file(args.from_file)
        table = args.table
    else:
        logger.info(f"Generating synthetic {args.table} data...")
        
        if args.table == 'molecules':
            data = generate_molecule_data(args.count)
            table = 'molecules'
        elif args.table == 'propertyTypes':
            data = generate_property_type_data()
            table = 'propertyTypes'
        else:
            logger.error(f"Unsupported table for synthetic data: {args.table}")
            sys.exit(1)
    
    logger.info(f"Generated {len(data)} records for table '{table}'")
    
    if args.dry_run:
        logger.info("DRY RUN: Would insert the following data:")
        print(json.dumps(data[:5], indent=2))
        if len(data) > 5:
            print(f"... and {len(data) - 5} more records")
        return
    
    # Insert data into Convex
    logger.info(f"Inserting {len(data)} records into table '{table}'...")
    
    chunk_size = 25  # Insert in chunks to avoid timeouts
    chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]
    
    inserted_ids = []
    
    for i, chunk in enumerate(chunks):
        logger.info(f"Inserting chunk {i+1}/{len(chunks)} ({len(chunk)} records)...")
        
        response = insert_data(table, chunk)
        
        if 'error' in response and response['error']:
            logger.error(f"Error inserting data: {response['error']}")
            break
        
        if 'data' in response:
            if isinstance(response['data'], list):
                inserted_ids.extend([item.get('id') for item in response['data'] if item.get('id')])
        
        # Sleep a bit to avoid rate limiting
        if i < len(chunks) - 1:
            time.sleep(0.5)
    
    logger.info(f"Successfully inserted {len(inserted_ids)} records")
    
    # If we inserted molecules, let's also insert property types and molecular properties
    if table == 'molecules' and inserted_ids:
        logger.info("Now inserting property types and molecular properties...")
        
        # Insert property types
        property_types = generate_property_type_data()
        logger.info(f"Inserting {len(property_types)} property types...")
        
        property_type_response = insert_data('propertyTypes', property_types)
        
        if 'error' in property_type_response and property_type_response['error']:
            logger.error(f"Error inserting property types: {property_type_response['error']}")
        else:
            property_type_ids = [item.get('id') for item in property_type_response.get('data', []) if item.get('id')]
            logger.info(f"Successfully inserted {len(property_type_ids)} property types")
            
            # Insert molecular properties
            molecular_properties = generate_molecular_properties(inserted_ids, property_type_ids)
            logger.info(f"Inserting {len(molecular_properties)} molecular properties...")
            
            # Insert in chunks
            property_chunks = [molecular_properties[i:i + chunk_size] for i in range(0, len(molecular_properties), chunk_size)]
            
            for i, chunk in enumerate(property_chunks):
                logger.info(f"Inserting property chunk {i+1}/{len(property_chunks)} ({len(chunk)} records)...")
                
                property_response = insert_data('molecularProperties', chunk)
                
                if 'error' in property_response and property_response['error']:
                    logger.error(f"Error inserting molecular properties: {property_response['error']}")
                    break
                
                # Sleep a bit to avoid rate limiting
                if i < len(property_chunks) - 1:
                    time.sleep(0.5)
            
            logger.info("Finished inserting molecular properties")

if __name__ == '__main__':
    main()
#!/usr/bin/env python3
"""
Script to create molecular properties in Convex.

This script reads molecule and property type IDs from Convex and creates
molecular properties for them.
"""

import os
import sys
import json
import random
import logging
import requests
import subprocess
from typing import Dict, List, Any

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
    return "https://hallowed-malamute-424.convex.cloud"

def get_table_data(table: str) -> List[Dict[str, Any]]:
    """
    Get data from a Convex table using the convex CLI.
    
    Args:
        table: Name of the table
        
    Returns:
        List of records from the table
    """
    try:
        # Run convex data command
        result = subprocess.run(
            ['npx', 'convex', 'data', table],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse the output to extract data
        # Format is typically tabular with | separators
        lines = result.stdout.strip().split('\n')
        if len(lines) <= 2:  # Header + separator line only
            return []
        
        # Get header and data lines
        header_line = lines[0]
        data_lines = lines[2:]  # Skip header and separator lines
        
        # Parse header to get column names
        headers = [h.strip() for h in header_line.split('|')]
        
        # Parse data rows
        records = []
        for line in data_lines:
            if not line.strip():
                continue
                
            # Split by | and strip whitespace
            values = [v.strip() for v in line.split('|')]
            
            # Create record with headers as keys
            record = {}
            for i, header in enumerate(headers):
                if i < len(values):
                    # Try to parse JSON for values that look like JSON
                    value = values[i]
                    if value.startswith('"') and value.endswith('"'):
                        # This is likely a string ID or similar
                        record[header] = value.strip('"')
                    elif value.lower() == 'true':
                        record[header] = True
                    elif value.lower() == 'false':
                        record[header] = False
                    elif value.lower() == 'null' or value == '':
                        record[header] = None
                    else:
                        try:
                            # Try to convert to number if possible
                            record[header] = float(value) if '.' in value else int(value)
                        except ValueError:
                            record[header] = value
            
            records.append(record)
            
        return records
    except subprocess.CalledProcessError as e:
        logger.error(f"Error getting data from table {table}: {e}")
        logger.error(f"Stdout: {e.stdout}")
        logger.error(f"Stderr: {e.stderr}")
        return []
    except Exception as e:
        logger.error(f"Unexpected error parsing data from convex: {e}")
        return []

def generate_molecular_properties(molecules: List[Dict[str, Any]], property_types: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Generate molecular properties for molecules and property types.
    
    Args:
        molecules: List of molecule records
        property_types: List of property type records
        
    Returns:
        List of molecular property records
    """
    properties = []
    
    # Define property value generators based on molecule and property type
    for molecule in molecules:
        molecule_id = molecule['_id']
        molecule_name = molecule['name']
        
        for prop_type in property_types:
            prop_type_id = prop_type['_id']
            prop_name = prop_type['name']
            
            # Generate appropriate value based on property type
            if prop_name == 'molecular_weight':
                # Molecular weight depends on formula
                formula = molecule.get('formula', '')
                if formula == 'C3H8O3':  # Glycerol
                    value = 92.09
                elif formula == 'C2H6OS':  # DMSO
                    value = 78.13
                elif formula == 'C2H6O2':  # Ethylene Glycol
                    value = 62.07
                elif formula == 'C3H8O2':  # Propylene Glycol
                    value = 76.10
                elif formula == 'C12H22O11':  # Trehalose, Sucrose
                    value = 342.30
                elif formula == 'C6H12O6':  # Glucose
                    value = 180.16
                elif formula == 'C6H14O6':  # Mannitol, Sorbitol
                    value = 182.17
                elif formula == 'CH4O':  # Methanol
                    value = 32.04
                elif formula == 'C2H6O':  # Ethanol
                    value = 46.07
                else:
                    # Generate a reasonable molecular weight for unknown formula
                    value = random.uniform(30, 500)
            
            elif prop_name == 'logp':
                # LogP values for common cryoprotectants
                if molecule_name == 'Glycerol':
                    value = -1.76
                elif molecule_name == 'DMSO':
                    value = -1.35
                elif molecule_name == 'Ethylene Glycol':
                    value = -1.2
                elif molecule_name == 'Propylene Glycol':
                    value = -0.92
                elif molecule_name == 'Trehalose':
                    value = -4.0
                elif molecule_name == 'Sucrose':
                    value = -3.7
                elif molecule_name == 'Glucose':
                    value = -3.1
                elif molecule_name == 'Mannitol':
                    value = -3.1
                elif molecule_name == 'Sorbitol':
                    value = -2.2
                elif molecule_name == 'Methanol':
                    value = -0.77
                elif molecule_name == 'Ethanol':
                    value = -0.31
                else:
                    value = random.uniform(-5, 5)
            
            elif prop_name == 'melting_point':
                # Melting points for common cryoprotectants
                if molecule_name == 'Glycerol':
                    value = 18.2
                elif molecule_name == 'DMSO':
                    value = 19.0
                elif molecule_name == 'Ethylene Glycol':
                    value = -12.9
                elif molecule_name == 'Propylene Glycol':
                    value = -59.0
                elif molecule_name == 'Trehalose':
                    value = 203.0
                elif molecule_name == 'Sucrose':
                    value = 186.0
                elif molecule_name == 'Glucose':
                    value = 146.0
                elif molecule_name == 'Mannitol':
                    value = 166.0
                elif molecule_name == 'Sorbitol':
                    value = 95.0
                elif molecule_name == 'Methanol':
                    value = -97.6
                elif molecule_name == 'Ethanol':
                    value = -114.1
                else:
                    value = random.uniform(-100, 200)
            
            elif prop_name == 'boiling_point':
                # Boiling points for common cryoprotectants
                if molecule_name == 'Glycerol':
                    value = 290.0
                elif molecule_name == 'DMSO':
                    value = 189.0
                elif molecule_name == 'Ethylene Glycol':
                    value = 197.3
                elif molecule_name == 'Propylene Glycol':
                    value = 188.2
                elif molecule_name == 'Methanol':
                    value = 64.7
                elif molecule_name == 'Ethanol':
                    value = 78.2
                else:
                    value = random.uniform(60, 300)
            
            elif prop_name == 'glass_transition_temp':
                # Glass transition temperatures for common cryoprotectants
                if molecule_name == 'Glycerol':
                    value = -93.0
                elif molecule_name == 'DMSO':
                    value = -135.0
                elif molecule_name == 'Ethylene Glycol':
                    value = -129.0
                elif molecule_name == 'Propylene Glycol':
                    value = -108.0
                elif molecule_name == 'Trehalose':
                    value = -30.0
                elif molecule_name == 'Sucrose':
                    value = -32.0
                elif molecule_name == 'Glucose':
                    value = -43.0
                elif molecule_name == 'Methanol':
                    value = -174.0
                elif molecule_name == 'Ethanol':
                    value = -130.0
                else:
                    value = random.uniform(-150, -20)
            
            elif prop_name == 'toxicity_score':
                # Toxicity scores (lower is less toxic)
                if molecule_name == 'Glycerol':
                    value = 1.2
                elif molecule_name == 'DMSO':
                    value = 2.5
                elif molecule_name == 'Ethylene Glycol':
                    value = 5.8
                elif molecule_name == 'Propylene Glycol':
                    value = 2.0
                elif molecule_name == 'Trehalose':
                    value = 0.5
                elif molecule_name == 'Sucrose':
                    value = 0.4
                elif molecule_name == 'Glucose':
                    value = 0.3
                elif molecule_name == 'Mannitol':
                    value = 0.6
                elif molecule_name == 'Sorbitol':
                    value = 0.7
                elif molecule_name == 'Methanol':
                    value = 8.5
                elif molecule_name == 'Ethanol':
                    value = 3.2
                else:
                    value = random.uniform(0, 10)
            else:
                # For any other properties, generate a random value
                value = random.uniform(0, 100)
            
            # Create the property record
            property_record = {
                "moleculeId": molecule_id,
                "propertyTypeId": prop_type_id,
                "value": value,
                "numericValue": value,
                "source": None,
                "calculationMethod": "synthetic_data"
            }
            
            # Add units if available in property type
            if 'units' in prop_type and prop_type['units']:
                property_record['units'] = prop_type['units']
            
            properties.append(property_record)
    
    return properties

def import_data_to_convex(data: List[Dict[str, Any]], table: str) -> bool:
    """
    Import data to Convex using the convex import command.
    
    Args:
        data: List of records to import
        table: Name of the table
        
    Returns:
        bool: Success or failure
    """
    try:
        # Process data in smaller chunks to avoid issues
        chunk_size = 25
        chunks = [data[i:i+chunk_size] for i in range(0, len(data), chunk_size)]
        
        success_count = 0
        for i, chunk in enumerate(chunks):
            # Write chunk to temporary file
            temp_file = f"temp_{table}_import_{i}.json"
            with open(temp_file, 'w') as f:
                json.dump(chunk, f, indent=2)
            
            logger.info(f"Importing chunk {i+1}/{len(chunks)} with {len(chunk)} records...")
            
            # Import data using convex CLI with append flag
            result = subprocess.run(
                ['npx', 'convex', 'import', temp_file, '--table', table, '--append'],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Clean up temporary file
            os.remove(temp_file)
            
            logger.info(f"Import result for chunk {i+1}: {result.stdout}")
            success_count += 1
            
        logger.info(f"Successfully imported {success_count} of {len(chunks)} chunks.")
        return success_count > 0
    except subprocess.CalledProcessError as e:
        logger.error(f"Error importing data to table {table}: {e}")
        logger.error(f"Stdout: {e.stdout}")
        logger.error(f"Stderr: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error importing data: {e}")
        return False

def main():
    """Main entry point."""
    # Check if Convex is enabled
    if os.environ.get('USE_CONVEX', '').lower() not in ('true', 'yes', '1'):
        logger.warning("Convex is not enabled in environment. Setting USE_CONVEX=true")
        os.environ['USE_CONVEX'] = 'true'
    
    # Get molecules from Convex
    logger.info("Fetching molecules from Convex...")
    molecules = get_table_data('molecules')
    if not molecules:
        logger.error("No molecules found in Convex. Aborting.")
        sys.exit(1)
    
    logger.info(f"Found {len(molecules)} molecules in Convex.")
    
    # Get property types from Convex
    logger.info("Fetching property types from Convex...")
    property_types = get_table_data('propertyTypes')
    if not property_types:
        logger.error("No property types found in Convex. Aborting.")
        sys.exit(1)
    
    logger.info(f"Found {len(property_types)} property types in Convex.")
    
    # Generate molecular properties
    logger.info("Generating molecular properties...")
    properties = generate_molecular_properties(molecules, property_types)
    logger.info(f"Generated {len(properties)} molecular properties.")
    
    # Import properties to Convex
    logger.info("Importing molecular properties to Convex...")
    success = import_data_to_convex(properties, 'molecularProperties')
    
    if success:
        logger.info("Successfully imported molecular properties to Convex.")
    else:
        logger.error("Failed to import molecular properties to Convex.")

if __name__ == '__main__':
    main()
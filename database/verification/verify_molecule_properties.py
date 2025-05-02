#!/usr/bin/env python3
"""
Script to verify molecule properties using the actual schema structure.

This script executes a modified version of the SQL query specified in the directive
to check the molecule properties using the actual schema structure.
"""

import os
import sys
import json
import logging
from datetime import datetime
from dotenv import load_dotenv
from supabase import create_client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Verify molecule properties using the actual schema structure."""
    try:
        # Load environment variables from .env file
        load_dotenv()
        
        # Get Supabase URL and credentials
        supabase_url = os.environ.get('SUPABASE_URL')
        supabase_key = os.environ.get('SUPABASE_KEY')
        
        if not supabase_url or not supabase_key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        # Create Supabase client
        supabase = create_client(supabase_url, supabase_key)
        
        # Execute a modified version of the SQL query to check molecule properties
        # Original query from directive:
        # SELECT 
        #   count(*) as total_count,
        #   count(properties->'pubchem') as pubchem_count,
        #   count(properties->'chembl') as chembl_count,
        #   count(properties->'rdkit') as rdkit_count
        # FROM molecules;
        
        # Modified query to work with the actual schema structure
        result = supabase.rpc(
            'exec_sql', 
            {
                'query': """
                    SELECT 
                        count(*) as total_count,
                        count(pubchem_cid) as pubchem_count,
                        count(chembl_id) as chembl_id_count,
                        count(smiles) as smiles_count,
                        count(inchi) as inchi_count,
                        count(inchikey) as inchikey_count,
                        count(molecular_weight) as molecular_weight_count
                    FROM molecules;
                """
            }
        ).execute()
        
        # Check if the query was successful
        if result.data and len(result.data) > 0:
            # Extract the counts
            counts = result.data[0]
            
            # Format the counts for display
            if isinstance(counts, dict):
                formatted_counts = counts
            elif isinstance(counts, list) and len(counts) >= 7:
                formatted_counts = {
                    'total_count': counts[0],
                    'pubchem_count': counts[1],
                    'chembl_id_count': counts[2],
                    'smiles_count': counts[3],
                    'inchi_count': counts[4],
                    'inchikey_count': counts[5],
                    'molecular_weight_count': counts[6]
                }
            else:
                formatted_counts = {'error': 'Unexpected result format'}
            
            # Display the counts
            logger.info("Molecule property counts:")
            for key, value in formatted_counts.items():
                logger.info(f"  - {key}: {value}")
            
            # Save the results to a file
            results = {
                'timestamp': datetime.now().isoformat(),
                'counts': formatted_counts,
                'message': "Molecule property counts retrieved successfully"
            }
            
            output_path = 'reports/molecule_property_counts.json'
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {output_path}")
            
            return 0
        else:
            logger.error("No results returned from query")
            return 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
Script to count molecules and their properties using the Supabase API.

This script uses the Supabase API to count molecules and their properties
instead of using SQL queries directly.
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
    """Count molecules and their properties using the Supabase API."""
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
        
        # Get total count of molecules
        total_result = supabase.table('molecules').select('id', count='exact').execute()
        total_count = total_result.count if hasattr(total_result, 'count') else 0
        
        # Get count of molecules with pubchem_cid
        pubchem_result = supabase.table('molecules').select('id').not_.is_('pubchem_cid', 'null').execute()
        pubchem_count = len(pubchem_result.data) if pubchem_result.data else 0
        
        # Get count of molecules with chembl_id
        chembl_result = supabase.table('molecules').select('id').not_.is_('chembl_id', 'null').execute()
        chembl_count = len(chembl_result.data) if chembl_result.data else 0
        
        # Get count of molecules with smiles
        smiles_result = supabase.table('molecules').select('id').not_.is_('smiles', 'null').execute()
        smiles_count = len(smiles_result.data) if smiles_result.data else 0
        
        # Get count of molecules with inchi
        inchi_result = supabase.table('molecules').select('id').not_.is_('inchi', 'null').execute()
        inchi_count = len(inchi_result.data) if inchi_result.data else 0
        
        # Get count of molecules with inchikey
        inchikey_result = supabase.table('molecules').select('id').not_.is_('inchikey', 'null').execute()
        inchikey_count = len(inchikey_result.data) if inchikey_result.data else 0
        
        # Get count of molecules with molecular_weight
        mw_result = supabase.table('molecules').select('id').not_.is_('molecular_weight', 'null').execute()
        mw_count = len(mw_result.data) if mw_result.data else 0
        
        # Compile the counts
        counts = {
            'total_count': total_count,
            'pubchem_count': pubchem_count,
            'chembl_count': chembl_count,
            'smiles_count': smiles_count,
            'inchi_count': inchi_count,
            'inchikey_count': inchikey_count,
            'molecular_weight_count': mw_count
        }
        
        # Display the counts
        logger.info("Molecule property counts:")
        for key, value in counts.items():
            logger.info(f"  - {key}: {value}")
        
        # Save the results to a file
        results = {
            'timestamp': datetime.now().isoformat(),
            'counts': counts,
            'message': "Molecule property counts retrieved successfully"
        }
        
        output_path = 'reports/molecule_property_counts.json'
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results saved to {output_path}")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
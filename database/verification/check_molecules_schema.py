#!/usr/bin/env python3
"""
Script to check the schema of the molecules table.
"""

import os
import sys
import json
import logging
from dotenv import load_dotenv
from supabase import create_client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Check the schema of the molecules table."""
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
        
        # Get a sample row from the molecules table to see its structure
        result = supabase.table('molecules').select('*').limit(1).execute()
        
        if result.data and len(result.data) > 0:
            # Get the first row
            sample_row = result.data[0]
            
            # Print the column names
            logger.info("Columns in molecules table:")
            for column_name in sample_row.keys():
                logger.info(f"  - {column_name}")
            
            # Save the schema to a file
            schema = {
                'table_name': 'molecules',
                'columns': list(sample_row.keys()),
                'sample_row': sample_row
            }
            
            output_path = 'reports/molecules_schema.json'
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as f:
                json.dump(schema, f, indent=2)
            logger.info(f"Schema saved to {output_path}")
            
            return 0
        else:
            logger.warning("No rows found in molecules table")
            return 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
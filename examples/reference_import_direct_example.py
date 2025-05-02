#!/usr/bin/env python3
"""
Example script demonstrating how to use the reference compounds import script
with direct PostgreSQL connections.

This script shows how to:
1. Configure the environment for direct PostgreSQL connections
2. Run the import with various options
3. Process the results
"""

import os
import sys
import logging
import json
from datetime import datetime

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# We'll import the module later after setting up environment variables

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_environment():
    """
    Set up environment variables for direct PostgreSQL connection.
    In a production environment, these would typically be loaded from a .env file.
    """
    # These are example values - replace with actual connection details
    os.environ['SUPABASE_DB_HOST'] = 'db.example.supabase.co'
    os.environ['SUPABASE_DB_PORT'] = '5432'
    os.environ['SUPABASE_DB_NAME'] = 'postgres'
    os.environ['SUPABASE_DB_USER'] = 'postgres'
    os.environ['SUPABASE_DB_PASSWORD'] = 'your_password_here'
    
    logger.info("Environment variables set for direct PostgreSQL connection")

def run_import_example():
    """Run the reference compounds import with various options."""
    # Set up environment variables
    setup_environment()
    
    # Import the reference compounds import module after setting up environment
    from import_reference_compounds import import_reference_compounds
    
    # Define output paths
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = f"reports/reference_import_{timestamp}.json"
    checkpoint_path = f"checkpoints/reference_import_{timestamp}.json"
    
    # Ensure directories exist
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    os.makedirs(os.path.dirname(checkpoint_path), exist_ok=True)
    
    logger.info("Starting reference compounds import with direct PostgreSQL connection")
    
    try:
        # Run the import with custom options
        results = import_reference_compounds(
            output_report=report_path,
            checkpoint_path=checkpoint_path,
            resume=True,  # Try to resume from checkpoint if available
            batch_size=3  # Process 3 compounds at a time
        )
        
        # Process the results
        logger.info(f"Import completed successfully")
        logger.info(f"Total compounds: {results['total_compounds']}")
        logger.info(f"Imported: {results['imported']}")
        logger.info(f"Updated: {results['updated']}")
        logger.info(f"Failed: {results['failed']}")
        
        # Print details of any failed imports
        if results['failed'] > 0:
            logger.warning("Some compounds failed to import:")
            for chembl_id, details in results['details'].items():
                if details['status'] == 'failed':
                    logger.warning(f"  - {chembl_id}: {details.get('error', 'Unknown error')}")
        
        # Print report location
        logger.info(f"Detailed report saved to: {report_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        raise

def main():
    """Main entry point."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Example script for reference compounds import with direct PostgreSQL connection')
    parser.add_argument('--run', action='store_true', help='Actually run the import (default: just show help)')
    parser.add_argument('--batch-size', type=int, default=3, help='Number of compounds to process in each batch')
    
    args = parser.parse_args()
    
    if args.run:
        try:
            run_import_example()
            return 0
        except Exception as e:
            logger.error(f"Example failed: {str(e)}")
            return 1
    else:
        parser.print_help()
        return 0

if __name__ == "__main__":
    sys.exit(main())
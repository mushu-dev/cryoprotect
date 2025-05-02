#!/usr/bin/env python3
"""
Example script demonstrating the use of the reconcile_chembl_properties.py script
with direct PostgreSQL connections.

This script shows how to:
1. Set up the environment for direct PostgreSQL connections
2. Run the reconciliation process
3. Process and analyze the results

Usage:
    python examples/reconcile_chembl_pubchem_direct_example.py
"""

import os
import sys
import json
import logging
from datetime import datetime
from pathlib import Path

# Add parent directory to path to import the script
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the reconciliation script
import reconcile_chembl_properties
import sql_executor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('reconcile_example')

def setup_environment():
    """
    Set up the environment for direct PostgreSQL connections.
    
    This function checks if the required environment variables are set
    and provides guidance if they are not.
    """
    required_vars = [
        'SUPABASE_DB_HOST',
        'SUPABASE_DB_PORT',
        'SUPABASE_DB_NAME',
        'SUPABASE_DB_USER',
        'SUPABASE_DB_PASSWORD'
    ]
    
    missing_vars = [var for var in required_vars if not os.getenv(var)]
    
    if missing_vars:
        logger.error(f"Missing required environment variables: {', '.join(missing_vars)}")
        logger.info("Please set the following environment variables:")
        for var in missing_vars:
            logger.info(f"  {var}")
        logger.info("You can set them in a .env file or directly in your environment.")
        return False
    
    logger.info("Environment variables are properly set")
    return True

def run_reconciliation(dry_run=False):
    """
    Run the reconciliation process.
    
    Args:
        dry_run: If True, run in dry-run mode without making database changes
        
    Returns:
        Path to the generated report file
    """
    # Create reports directory if it doesn't exist
    reports_dir = Path("reports")
    reports_dir.mkdir(exist_ok=True)
    
    # Generate report filename with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = reports_dir / f"reconciliation_{timestamp}.json"
    
    logger.info(f"Starting reconciliation process (dry_run={dry_run})")
    
    # Run the reconciliation
    results = reconcile_chembl_properties.reconcile_cross_references(
        output_report=str(report_path),
        dry_run=dry_run
    )
    
    logger.info(f"Reconciliation completed: {results['molecules_updated']} molecules updated")
    
    return report_path

def analyze_results(report_path):
    """
    Analyze the reconciliation results.
    
    Args:
        report_path: Path to the reconciliation report file
    """
    logger.info(f"Analyzing results from {report_path}")
    
    try:
        with open(report_path, 'r') as f:
            results = json.load(f)
        
        # Print summary
        logger.info("Reconciliation Summary:")
        logger.info(f"  Timestamp: {results['timestamp']}")
        logger.info(f"  Molecules Updated: {results['molecules_updated']}")
        logger.info(f"  Cross-References Added: {results['cross_references_added']}")
        logger.info(f"  Conflicts Resolved: {results['conflicts_resolved']}")
        
        # Print details for a few examples
        if results['details']:
            logger.info("Example Reconciliations:")
            for i, (inchi_key, detail) in enumerate(results['details'].items()):
                logger.info(f"  {i+1}. InChI Key: {inchi_key}")
                logger.info(f"     PubChem Molecule: {detail['pubchem_molecule']}")
                logger.info(f"     ChEMBL Molecule: {detail['chembl_molecule']}")
                logger.info(f"     Action: {detail['action']}")
                
                # Only show a few examples
                if i >= 2:
                    remaining = len(results['details']) - 3
                    if remaining > 0:
                        logger.info(f"  ... and {remaining} more")
                    break
        else:
            logger.info("No reconciliations were performed")
            
    except Exception as e:
        logger.error(f"Error analyzing results: {str(e)}")

def verify_database_updates(dry_run=False):
    """
    Verify that the database updates were successful.
    
    Args:
        dry_run: If True, skip verification (since no changes were made)
    """
    if dry_run:
        logger.info("Skipping verification in dry-run mode")
        return
    
    logger.info("Verifying database updates")
    
    try:
        # Query to find molecules with both ChEMBL and PubChem identifiers
        query = """
            SELECT COUNT(*) as count
            FROM molecules
            WHERE chembl_id IS NOT NULL AND pubchem_cid IS NOT NULL
        """
        
        result = sql_executor.execute_query(query, fetch_one=True)
        
        if result and 'count' in result:
            count = result['count']
            logger.info(f"Found {count} molecules with both ChEMBL and PubChem identifiers")
            
            if count > 0:
                # Get a few examples
                examples_query = """
                    SELECT id, name, chembl_id, pubchem_cid, inchikey
                    FROM molecules
                    WHERE chembl_id IS NOT NULL AND pubchem_cid IS NOT NULL
                    LIMIT 3
                """
                
                examples = sql_executor.execute_query(examples_query)
                
                logger.info("Example molecules with cross-references:")
                for i, example in enumerate(examples):
                    logger.info(f"  {i+1}. ID: {example['id']}")
                    logger.info(f"     Name: {example['name']}")
                    logger.info(f"     ChEMBL ID: {example['chembl_id']}")
                    logger.info(f"     PubChem CID: {example['pubchem_cid']}")
                    logger.info(f"     InChI Key: {example['inchikey']}")
        else:
            logger.warning("Could not verify database updates")
            
    except Exception as e:
        logger.error(f"Error verifying database updates: {str(e)}")

def main():
    """Main function to run the example."""
    logger.info("Starting reconciliation example")
    
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description='Example script for ChEMBL-PubChem reconciliation')
    parser.add_argument('--dry-run', action='store_true',
                      help='Run in dry-run mode without making database changes')
    args = parser.parse_args()
    
    # Set up environment
    if not setup_environment() and not args.dry_run:
        logger.error("Environment setup failed. Use --dry-run to test without database connection.")
        return 1
    
    try:
        # Run reconciliation
        report_path = run_reconciliation(dry_run=args.dry_run)
        
        # Analyze results
        analyze_results(report_path)
        
        # Verify database updates
        verify_database_updates(dry_run=args.dry_run)
        
        # Close database connections
        sql_executor.close_connections()
        logger.info("Database connections closed")
        
        logger.info("Reconciliation example completed successfully")
        return 0
        
    except Exception as e:
        logger.error(f"Error running reconciliation example: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
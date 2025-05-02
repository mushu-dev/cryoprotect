#!/usr/bin/env python3
"""
Example usage of the database verification script.

This script demonstrates how to use the verify_imported_data.py script
to validate the success of the direct database population process.
"""

import os
import sys
import logging
import json
from datetime import datetime

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from verify_imported_data import (
    perform_full_verification,
    generate_markdown_report,
    update_project_state
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Run the example verification."""
    logger.info("Starting database verification example")
    
    # Create reports directory if it doesn't exist
    reports_dir = os.path.join(os.path.dirname(__file__), '..', 'reports')
    os.makedirs(reports_dir, exist_ok=True)
    
    # Generate timestamp for report filenames
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Define report paths
    json_report_path = os.path.join(reports_dir, f'verification_{timestamp}.json')
    markdown_report_path = os.path.join(reports_dir, f'verification_{timestamp}.md')
    
    # Perform full verification
    logger.info(f"Performing full verification with JSON report: {json_report_path}")
    results = perform_full_verification(json_report_path)
    
    # Generate Markdown report
    logger.info(f"Generating Markdown report: {markdown_report_path}")
    generate_markdown_report(results, markdown_report_path)
    
    # Update project state
    logger.info("Updating project state with verification results")
    update_project_state(results)
    
    # Display summary
    if results['success']:
        logger.info("Verification PASSED!")
    else:
        logger.error("Verification FAILED!")
        if 'error' in results:
            logger.error(f"Error: {results['error']}")
    
    logger.info("Summary:")
    for key, value in results.get('summary', {}).items():
        if isinstance(value, float):
            logger.info(f"  {key}: {value:.2f}")
        else:
            logger.info(f"  {key}: {value}")
    
    # Example of accessing specific verification results
    logger.info("\nDetailed Results:")
    
    # Molecule counts
    if 'molecule_counts' in results:
        counts = results['molecule_counts']
        logger.info(f"Total molecules: {counts['total_molecules']}")
        logger.info(f"Molecules with PubChem CID: {counts['with_pubchem_cid']}")
        logger.info(f"Molecules with ChEMBL ID: {counts['with_chembl_id']}")
        logger.info(f"Molecules with cross-references: {counts['with_cross_references']}")
    
    # Reference compounds
    if 'reference_compounds' in results:
        ref = results['reference_compounds']
        logger.info(f"Reference compounds found: {ref['found_reference_compounds']}/{ref['total_reference_compounds']}")
        
        if ref['missing_reference_compounds']:
            logger.warning(f"Missing reference compounds: {', '.join(ref['missing_reference_compounds'])}")
        
        if ref['incomplete_reference_compounds']:
            logger.warning(f"Incomplete reference compounds: {', '.join(ref['incomplete_reference_compounds'])}")
    
    # Property completeness
    if 'property_completeness' in results:
        prop = results['property_completeness']
        logger.info(f"Property completeness: {prop['property_completeness_percentage']:.2f}%")
        logger.info(f"Molecules with complete properties: {prop['molecules_with_complete_properties']}")
        logger.info(f"Molecules with incomplete properties: {prop['molecules_with_incomplete_properties']}")
    
    # Query performance
    if 'query_performance' in results:
        perf = results['query_performance']
        logger.info(f"Average query time: {perf['overall_average_ms']:.2f} ms")
        logger.info(f"Performance acceptable: {perf['performance_acceptable']}")
    
    logger.info("\nVerification example completed")
    
    return 0 if results['success'] else 1

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
Test script for reference compound import.

This script tests the import_reference_compounds.py script to verify that
the critical properties are being set correctly.
"""

import os
import sys
import logging
import json
from datetime import datetime

# Import the reference compound import function
from import_reference_compounds import import_reference_compounds, verify_reference_compounds

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('test_reference_import')

def run_test():
    """
    Run a test of the reference compound import process.
    
    This function:
    1. Imports a small number of reference compounds
    2. Verifies that the critical properties are set correctly
    3. Logs the results
    """
    logger.info("Starting reference compound import test")
    
    # Create a test-specific checkpoint file
    checkpoint_path = f"checkpoints/test_reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    # Create a test-specific output report
    output_report = f"reports/test_reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    try:
        # Run the import with a small batch size for testing
        results = import_reference_compounds(
            output_report=output_report,
            checkpoint_path=checkpoint_path,
            resume=False,
            batch_size=2
        )
        
        logger.info(f"Import completed with results: {json.dumps(results, indent=2)}")
        
        # Verify the reference compounds
        logger.info("Verifying reference compounds...")
        complete_count, total_count = verify_reference_compounds(results)
        
        logger.info(f"Verification results: {complete_count}/{total_count} complete")
        
        # Check if all reference compounds have the required properties
        if complete_count == total_count:
            logger.info("TEST PASSED: All reference compounds have required properties")
            return True
        else:
            logger.error(f"TEST FAILED: Only {complete_count}/{total_count} reference compounds have required properties")
            return False
            
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    success = run_test()
    sys.exit(0 if success else 1)
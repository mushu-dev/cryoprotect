#!/usr/bin/env python3
"""
Test script for ChEMBL data integrity spot-check.

This script tests the functionality of verify_chembl_integrity.py by:
1. Running the spot-check on a small sample of ChEMBL records
2. Verifying that the report is generated correctly
3. Checking that the results are as expected

Related Task: task-imp-wv-1-1-integrity-check
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime

# Import the spot-check script
from verify_chembl_integrity import (
    initialize_supabase,
    initialize_chembl_client,
    spot_check_integrity,
    save_report
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/test_chembl_integrity.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

def test_initialize_supabase():
    """Test initializing Supabase client."""
    logger.info("Testing Supabase client initialization")
    
    supabase = initialize_supabase()
    assert supabase is not None, "Failed to initialize Supabase client"
    
    logger.info("Supabase client initialization test passed")
    return supabase

def test_initialize_chembl_client():
    """Test initializing ChEMBL client."""
    logger.info("Testing ChEMBL client initialization")
    
    result = initialize_chembl_client()
    assert result is True, "Failed to initialize ChEMBL client"
    
    logger.info("ChEMBL client initialization test passed")

def test_spot_check_integrity(supabase):
    """Test spot-checking data integrity."""
    logger.info("Testing spot-check integrity function")
    
    # Run spot-check with a small limit
    results = spot_check_integrity(supabase, limit=3)
    
    # Verify results structure
    assert "timestamp" in results, "Results missing timestamp"
    assert "records_checked" in results, "Results missing records_checked"
    assert "records_matched" in results, "Results missing records_matched"
    assert "records_mismatched" in results, "Results missing records_mismatched"
    assert "placeholder_data_found" in results, "Results missing placeholder_data_found"
    assert "overall_match_percentage" in results, "Results missing overall_match_percentage"
    assert "details" in results, "Results missing details"
    
    # Verify that at least some records were checked
    assert results["records_checked"] > 0, "No records were checked"
    
    # Save test report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = save_report(results, f"test_chembl_integrity_{timestamp}.json")
    
    logger.info(f"Spot-check integrity test completed. Report saved to {report_path}")
    logger.info(f"Records checked: {results['records_checked']}")
    logger.info(f"Records matched: {results['records_matched']}")
    logger.info(f"Records mismatched: {results['records_mismatched']}")
    logger.info(f"Placeholder data found: {results['placeholder_data_found']}")
    logger.info(f"Overall match percentage: {results['overall_match_percentage']:.2f}%")
    
    return results

def main():
    """Main test function."""
    logger.info("Starting ChEMBL integrity spot-check tests")
    
    try:
        # Test Supabase initialization
        supabase = test_initialize_supabase()
        
        # Test ChEMBL client initialization
        test_initialize_chembl_client()
        
        # Test spot-check integrity
        results = test_spot_check_integrity(supabase)
        
        logger.info("All tests passed successfully")
        return True
    
    except Exception as e:
        logger.error(f"Test failed: {str(e)}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
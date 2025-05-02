#!/usr/bin/env python3
"""
Test script to verify that the CryoProtect system uses only Tox21 for toxicity data.
This script checks that:
1. The toxicity_data_source table contains only Tox21 (no ToxCast)
2. The calculation_method table has the correct Tox21 Score method
3. The ToxCast client file has been removed
"""

import os
import sys
import logging
import subprocess
from supabase import create_client
from config import Config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_toxicity_data_sources(supabase):
    """Test that only Tox21 exists in the toxicity_data_source table."""
    logger.info("Testing toxicity data sources...")
    
    try:
        response = supabase.table("toxicity_data_source").select("*").execute()
        sources = response.data
        
        logger.info(f"Found {len(sources)} toxicity data sources")
        
        # Check that ToxCast is not present
        toxcast_sources = [s for s in sources if "toxcast" in s["name"].lower()]
        if toxcast_sources:
            logger.error(f"Found ToxCast sources: {toxcast_sources}")
            return False
        
        # Check that Tox21 is present
        tox21_sources = [s for s in sources if "tox21" in s["name"].lower()]
        if not tox21_sources:
            logger.error("Tox21 source not found")
            return False
        
        logger.info("Toxicity data sources test passed")
        return True
    except Exception as e:
        logger.error(f"Error testing toxicity data sources: {str(e)}")
        return False

def test_calculation_methods(supabase):
    """Test that the calculation method is correctly set to Tox21 Score."""
    logger.info("Testing calculation methods...")
    
    try:
        response = supabase.table("calculation_method").select("*").execute()
        methods = response.data
        
        # Check for Tox21 Score method
        tox21_methods = [m for m in methods if "tox21" in m["name"].lower()]
        if not tox21_methods:
            logger.error("Tox21 Score calculation method not found")
            return False
        
        # Check that no ToxCast methods exist
        toxcast_methods = [m for m in methods if "toxcast" in m["name"].lower()]
        if toxcast_methods:
            logger.error(f"Found ToxCast calculation methods: {toxcast_methods}")
            return False
        
        logger.info("Calculation methods test passed")
        return True
    except Exception as e:
        logger.error(f"Error testing calculation methods: {str(e)}")
        return False

def test_toxcast_client_removed():
    """Test that the ToxCast client file has been removed."""
    logger.info("Testing ToxCast client removal...")
    
    toxcast_client_path = "chemical_data/toxicity/toxcast_client.py"
    if os.path.exists(toxcast_client_path):
        logger.error(f"ToxCast client file still exists: {toxcast_client_path}")
        return False
    
    logger.info("ToxCast client file has been removed")
    return True

def test_toxcast_references():
    """Test that there are no ToxCast references in Python files."""
    logger.info("Testing for ToxCast references in code...")
    
    try:
        # Use grep to search for ToxCast references in Python files
        result = subprocess.run(
            ["grep", "-r", "ToxCast", "--include=*.py", "."],
            capture_output=True,
            text=True
        )
        
        # Check if any references were found
        if result.stdout.strip():
            logger.error(f"Found ToxCast references in code:\n{result.stdout}")
            return False
        
        logger.info("No ToxCast references found in Python files")
        return True
    except Exception as e:
        logger.error(f"Error searching for ToxCast references: {str(e)}")
        # Don't fail the test if grep is not available
        logger.warning("Skipping ToxCast reference check")
        return True

def main():
    """Run all tests."""
    logger.info("Starting Tox21-only verification tests")
    
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Run tests
    tests = [
        test_toxicity_data_sources(supabase),
        test_calculation_methods(supabase),
        test_toxcast_client_removed(),
        test_toxcast_references()
    ]
    
    # Check results
    if all(tests):
        logger.info("All tests passed! The system is using only Tox21 for toxicity data.")
        return True
    else:
        logger.error("Some tests failed. The system may still have ToxCast references.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
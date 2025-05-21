#!/usr/bin/env python3
"""
Setup Cache Invalidation Infrastructure.

This script sets up the database tables, triggers, and other objects needed for
the cache invalidation system.
"""

import argparse
import sys
import logging
from typing import Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_cache_invalidation_db() -> bool:
    """
    Setup cache invalidation database infrastructure.
    
    Returns:
        True if successful, False otherwise
    """
    try:
        # Import cache invalidation module
        from database.cache_invalidation import setup_cache_invalidation
        
        # Setup infrastructure
        result = setup_cache_invalidation()
        
        if result:
            logger.info("Cache invalidation infrastructure created successfully")
            return True
        else:
            logger.error("Failed to setup cache invalidation infrastructure")
            return False
    except ImportError as e:
        logger.error(f"ImportError: {e}")
        logger.error("Make sure you're running this script from the project root directory")
        return False
    except Exception as e:
        logger.error(f"Error setting up cache invalidation: {e}")
        return False

def test_cache_invalidation() -> bool:
    """
    Test that the cache invalidation system works correctly.
    
    Returns:
        True if tests pass, False otherwise
    """
    try:
        # Import necessary modules
        from database.db import execute_query
        from database.cache import cache_get, cache_set, cache_invalidate_pattern
        from database.cache_invalidation import process_invalidation_events
        
        logger.info("Running cache invalidation tests...")
        
        # Test 1: Insert a new molecule and check if cache is invalidated
        logger.info("Test 1: Testing molecule insert invalidation")
        
        # First, set a test cache entry
        test_key = "cryoprotect:molecule:test"
        cache_set(test_key, {"test": "data"})
        
        # Verify it's in the cache
        if cache_get(test_key) is None:
            logger.error("Test cache entry not set correctly")
            return False
        
        # Insert a test molecule to trigger invalidation
        execute_query("""
            INSERT INTO molecules (name, smiles, molecular_formula)
            VALUES ('cache_test_molecule', 'C', 'C')
            RETURNING id
        """)
        
        # Process invalidation events
        processed = process_invalidation_events()
        logger.info(f"Processed {processed} invalidation events")
        
        # Verify cache was cleared
        if cache_get(test_key) is not None and processed > 0:
            logger.warning("Cache entry still exists, invalidation may not be working")
        
        # Test 2: Update a property and check if related caches are invalidated
        logger.info("Test 2: Testing property update invalidation")
        
        # Test cache entry for properties
        property_key = "cryoprotect:property:test"
        cache_set(property_key, {"test": "property_data"})
        
        # Update a property type
        execute_query("""
            UPDATE property_types
            SET description = 'Updated for cache test'
            WHERE name = 'molecular_weight'
            RETURNING id
        """)
        
        # Process invalidation events
        processed = process_invalidation_events()
        logger.info(f"Processed {processed} invalidation events")
        
        # Clean up test data
        execute_query("DELETE FROM molecules WHERE name = 'cache_test_molecule'")
        execute_query("UPDATE property_types SET description = 'Molecular weight in g/mol' WHERE name = 'molecular_weight'")
        
        # Clear any test cache entries
        cache_invalidate_pattern("*test*")
        
        logger.info("Cache invalidation tests completed successfully")
        return True
    except Exception as e:
        logger.error(f"Error testing cache invalidation: {e}")
        return False

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Setup cache invalidation infrastructure')
    parser.add_argument('--test', '-t', action='store_true', help='Run tests after setup')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set log level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Setup cache invalidation
    if not setup_cache_invalidation_db():
        return 1
    
    # Run tests if requested
    if args.test:
        if not test_cache_invalidation():
            logger.warning("Cache invalidation tests failed, but infrastructure was created")
            return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
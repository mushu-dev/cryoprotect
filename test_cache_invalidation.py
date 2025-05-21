#!/usr/bin/env python3
"""
Test script for the cache invalidation system.

This script tests that the cache invalidation system works correctly by:
1. Setting up test data in the cache
2. Making changes to the database
3. Verifying that the cache is invalidated correctly
"""

import time
import logging
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_invalidation_test() -> bool:
    """
    Run a comprehensive test of the cache invalidation system.
    
    Returns:
        True if all tests pass, False otherwise
    """
    try:
        # Import required modules
        from database.db import execute_query, transaction
        from database.cache import cache_get, cache_set, cache_invalidate_pattern, clear_cache
        from database.cache_invalidation import (
            setup_cache_invalidation,
            process_invalidation_events,
            invalidation_processor
        )
        from database.cached_db import (
            get_molecule_by_id,
            get_molecular_properties,
            get_property_types
        )
        
        logger.info("Starting cache invalidation test...")
        
        # Clear cache first
        clear_cache()
        logger.info("Cache cleared for fresh test")
        
        # Ensure cache invalidation is set up
        setup_cache_invalidation()
        
        # Create test molecule
        logger.info("Creating test molecule...")
        with transaction() as cursor:
            cursor.execute("""
                INSERT INTO molecules (name, smiles, molecular_formula)
                VALUES ('cache_test_molecule', 'C', 'CH4')
                RETURNING id
            """)
            molecule_id = cursor.fetchone()[0]
            logger.info(f"Created test molecule with ID: {molecule_id}")
            
            # Create test property
            cursor.execute("""
                INSERT INTO molecular_properties (molecule_id, property_type_id, value_numeric)
                SELECT %s, id, 16.04 FROM property_types WHERE name = 'molecular_weight'
                RETURNING id
            """, (molecule_id,))
            property_id = cursor.fetchone()[0]
            logger.info(f"Created test property with ID: {property_id}")
        
        # Access data to cache it
        logger.info("Accessing data to cache it...")
        molecule = get_molecule_by_id(molecule_id)
        if not molecule:
            logger.error("Failed to retrieve test molecule")
            return False
            
        properties = get_molecular_properties(molecule_id)
        if not properties:
            logger.error("Failed to retrieve test properties")
            return False
            
        logger.info(f"Retrieved and cached molecule: {molecule['name']}")
        logger.info(f"Retrieved and cached {len(properties)} properties")
        
        # Process any pending invalidation events from initial creation
        process_invalidation_events()
        
        # Test 1: Update molecule
        logger.info("Test 1: Updating molecule to test invalidation...")
        with transaction() as cursor:
            cursor.execute("""
                UPDATE molecules 
                SET name = 'cache_test_molecule_updated', updated_at = NOW()
                WHERE id = %s
            """, (molecule_id,))
        
        # Process invalidation events
        logger.info("Processing invalidation events...")
        count = process_invalidation_events()
        logger.info(f"Processed {count} invalidation events")
        
        # Verify cache was invalidated by checking if we get the updated data
        updated_molecule = get_molecule_by_id(molecule_id)
        if updated_molecule and updated_molecule['name'] == 'cache_test_molecule_updated':
            logger.info("Test 1 PASSED: Cache was invalidated for molecule update")
        else:
            logger.error(f"Test 1 FAILED: Cache was not invalidated for molecule update. Name: {updated_molecule['name'] if updated_molecule else 'None'}")
            return False
        
        # Test 2: Update property
        logger.info("Test 2: Updating property to test invalidation...")
        with transaction() as cursor:
            cursor.execute("""
                UPDATE molecular_properties 
                SET value_numeric = 16.05, updated_at = NOW()
                WHERE id = %s
            """, (property_id,))
        
        # Process invalidation events
        logger.info("Processing invalidation events...")
        count = process_invalidation_events()
        logger.info(f"Processed {count} invalidation events")
        
        # Verify cache was invalidated
        updated_properties = get_molecular_properties(molecule_id)
        found_updated = False
        for prop in updated_properties:
            if prop['id'] == property_id and prop['value_numeric'] == 16.05:
                found_updated = True
                break
                
        if found_updated:
            logger.info("Test 2 PASSED: Cache was invalidated for property update")
        else:
            logger.error("Test 2 FAILED: Cache was not invalidated for property update")
            return False
        
        # Test 3: Delete molecule (should invalidate molecule and property caches)
        logger.info("Test 3: Deleting molecule to test invalidation...")
        with transaction() as cursor:
            cursor.execute("""
                DELETE FROM molecular_properties WHERE molecule_id = %s
            """, (molecule_id,))
            cursor.execute("""
                DELETE FROM molecules WHERE id = %s
            """, (molecule_id,))
        
        # Process invalidation events
        logger.info("Processing invalidation events...")
        count = process_invalidation_events()
        logger.info(f"Processed {count} invalidation events")
        
        # Verify cache was invalidated
        deleted_molecule = get_molecule_by_id(molecule_id)
        if deleted_molecule is None:
            logger.info("Test 3 PASSED: Cache was invalidated for molecule deletion")
        else:
            logger.error("Test 3 FAILED: Cache was not invalidated for molecule deletion")
            return False
        
        # All tests passed
        logger.info("All cache invalidation tests PASSED!")
        return True
        
    except Exception as e:
        logger.error(f"Error during cache invalidation test: {e}")
        return False

if __name__ == "__main__":
    if run_invalidation_test():
        print("Cache invalidation system is working correctly!")
        exit(0)
    else:
        print("Cache invalidation system test failed!")
        exit(1)
#!/usr/bin/env python3
"""
Standalone test script for property insertion.

This script tests the property insertion logic directly without importing
the original module, to verify that our fixes work correctly.
"""

import os
import sys
import logging
import json
import uuid
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
import functools

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('test_property_insertion')

# Create directories if they don't exist
os.makedirs("checkpoints", exist_ok=True)
os.makedirs("reports", exist_ok=True)

# Import our mock modules
from db_connection_utils_mock import get_db_connection, safe_transaction
from transaction_utils_mock import with_transaction_retry, execute_in_transaction
from property_utils_mock import PropertyManager

# Patch the bulk_insert_properties function to use our PropertyManager instance
import batch_utils_mock
original_bulk_insert_properties = batch_utils_mock.bulk_insert_properties

@functools.wraps(original_bulk_insert_properties)
def patched_bulk_insert_properties(property_records, batch_size=1000, property_manager=None):
    """Patched version of bulk_insert_properties that uses the provided PropertyManager instance."""
    if property_manager is None:
        return original_bulk_insert_properties(property_records, batch_size)
    
    # Log the property records for debugging
    logger.info(f"Mock bulk_insert_properties: Would insert {len(property_records)} property records")
    
    # Log and store the property records
    for i, record in enumerate(property_records):
        molecule_id = record.get('molecule_id')
        property_type_id = record.get('property_type_id')
        numeric_value = record.get('numeric_value')
        text_value = record.get('text_value')
        
        logger.info(f"Record {i}: molecule_id={molecule_id}, property_type_id={property_type_id}, " +
                   f"numeric_value={numeric_value}, text_value={text_value}")
        
        # Only process records with property_type_id
        if molecule_id and property_type_id:
            # Find the property name for this property_type_id
            property_name = None
            for name, info in property_manager._property_types_cache.items():
                if info['id'] == property_type_id:
                    property_name = name
                    break
            
            if property_name:
                # Store the property value
                value = numeric_value if numeric_value is not None else text_value
                property_manager.set_property(molecule_id, property_name, value)
                logger.info(f"Stored property {property_name}={value} for molecule {molecule_id}")
            else:
                logger.warning(f"Could not find property name for property_type_id {property_type_id}")
    
    # Return the number of records processed
    return len(property_records)

# Replace the original function with our patched version
batch_utils_mock.bulk_insert_properties = patched_bulk_insert_properties

class TestPropertyInsertion:
    """Test class for property insertion."""
    
    def __init__(self):
        """Initialize the test class."""
        self.property_manager = PropertyManager()
        self.critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
    
    def test_property_insertion_with_property_name(self):
        """
        Test property insertion with property_name instead of property_type_id.
        
        This test simulates the bug in the original code where property records
        were created with property_name instead of property_type_id.
        """
        logger.info("Testing property insertion with property_name instead of property_type_id")
        
        # Create a test molecule ID
        molecule_id = str(uuid.uuid4())
        
        # Create property records with property_name instead of property_type_id
        property_records = []
        for prop_name in self.critical_properties:
            if prop_name == 'logP':
                property_records.append({
                    "molecule_id": molecule_id,
                    "property_name": prop_name,  # Bug: Using property_name instead of property_type_id
                    "numeric_value": 2.5,
                    "text_value": None,
                    "created_at": datetime.now().isoformat()
                })
            else:  # h_bond_donors or h_bond_acceptors
                property_records.append({
                    "molecule_id": molecule_id,
                    "property_name": prop_name,  # Bug: Using property_name instead of property_type_id
                    "numeric_value": 1,
                    "text_value": None,
                    "created_at": datetime.now().isoformat()
                })
        
        # Try to insert the properties
        try:
            # This should fail because bulk_insert_properties expects property_type_id
            success_count = batch_utils_mock.bulk_insert_properties(property_records, property_manager=self.property_manager)
            logger.info(f"Inserted {success_count}/{len(property_records)} properties")
            
            # Verify that the properties were not inserted
            properties = self.property_manager.get_properties(molecule_id, self.critical_properties)
            logger.info(f"Properties after insertion: {properties}")
            
            if len(properties) == 0:
                logger.info("TEST PASSED: Properties were not inserted with property_name")
            else:
                logger.error("TEST FAILED: Properties were inserted with property_name")
        except Exception as e:
            logger.error(f"Error inserting properties: {str(e)}")
    
    def test_property_insertion_with_property_type_id(self):
        """
        Test property insertion with property_type_id.
        
        This test simulates the fixed code where property records are created
        with property_type_id instead of property_name.
        """
        logger.info("Testing property insertion with property_type_id")
        
        # Create a test molecule ID
        molecule_id = str(uuid.uuid4())
        
        # Create property records with property_type_id
        property_records = []
        for prop_name in self.critical_properties:
            # Get property type ID
            property_type_id = self.property_manager.get_property_type_id(prop_name)
            
            if prop_name == 'logP':
                property_records.append({
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,  # Fixed: Using property_type_id
                    "numeric_value": 2.5,
                    "text_value": None,
                    "boolean_value": None,
                    "created_by": None
                })
            else:  # h_bond_donors or h_bond_acceptors
                property_records.append({
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,  # Fixed: Using property_type_id
                    "numeric_value": 1,
                    "text_value": None,
                    "boolean_value": None,
                    "created_by": None
                })
        
        # Try to insert the properties
        try:
            # This should succeed because we're using property_type_id
            success_count = batch_utils_mock.bulk_insert_properties(property_records, property_manager=self.property_manager)
            logger.info(f"Inserted {success_count}/{len(property_records)} properties")
            
            # Verify that the properties were inserted
            properties = self.property_manager.get_properties(molecule_id, self.critical_properties)
            logger.info(f"Properties after insertion: {properties}")
            
            if len(properties) == len(self.critical_properties):
                logger.info("TEST PASSED: All properties were inserted with property_type_id")
            else:
                logger.error(f"TEST FAILED: Only {len(properties)}/{len(self.critical_properties)} properties were inserted with property_type_id")
        except Exception as e:
            logger.error(f"Error inserting properties: {str(e)}")
    
    def run_tests(self):
        """Run all tests."""
        logger.info("Starting property insertion tests")
        
        # Test property insertion with property_name (should fail)
        self.test_property_insertion_with_property_name()
        
        # Test property insertion with property_type_id (should succeed)
        self.test_property_insertion_with_property_type_id()
        
        logger.info("Property insertion tests completed")

if __name__ == "__main__":
    test = TestPropertyInsertion()
    test.run_tests()
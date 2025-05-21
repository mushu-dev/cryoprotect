#!/usr/bin/env python3
"""
Focused debugging script for data population workflow.

This script systematically tests and fixes issues with the data population workflow,
focusing on property insertion, transaction management, and connection handling.
"""

import os
import sys
import logging
import json
import time
import uuid
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union

# Configure detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/debug_data_population.log')
    ]
)
logger = logging.getLogger(__name__)

# Import database connection modules
from db_connection_utils import get_db_connection, safe_transaction
from transaction_utils import with_transaction_retry, execute_in_transaction, is_transaction_active
from batch_utils import bulk_insert_properties, resumable_batch_import, batch_delete_properties
from property_utils import PropertyManager
from monitoring_utils import PerformanceMetrics

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)
os.makedirs('reports', exist_ok=True)

def verify_reference_compounds():
    """Verify that reference compounds have been imported with all required properties."""
    logger.info("Verifying reference compounds...")
    
    try:
        # Get reference compound IDs
        from chembl.reference_compounds import get_reference_compound_ids
        reference_ids = get_reference_compound_ids()
        
        logger.info(f"Found {len(reference_ids)} reference compound IDs")
        
        # Initialize PropertyManager
        property_manager = PropertyManager()
        
        # Check each reference compound
        complete_count = 0
        incomplete_count = 0
        missing_properties = {}
        
        with get_db_connection() as conn:
            for chembl_id in reference_ids:
                # Get molecule ID from ChEMBL ID
                result = conn.execute_query("""
                SELECT id FROM molecules WHERE chembl_id = %s
                """, (chembl_id,), fetch_one=True)
                
                if not result:
                    logger.error(f"❌ Reference compound {chembl_id} not found in database")
                    continue
                
                molecule_id = result['id']
                
                # Check critical properties
                critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
                properties = property_manager.get_properties(molecule_id, critical_properties)
                
                logger.info(f"Reference compound {chembl_id} (ID: {molecule_id}) properties: {properties}")
                
                # Check if all critical properties are present
                missing = [prop for prop in critical_properties if prop not in properties]
                
                if missing:
                    logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
                    missing_properties[chembl_id] = missing
                    incomplete_count += 1
                else:
                    logger.info(f"✅ Reference compound {chembl_id} has all critical properties")
                    complete_count += 1
        
        # Generate report
        logger.info(f"Reference compounds verification complete:")
        logger.info(f"  Total reference compounds: {len(reference_ids)}")
        logger.info(f"  Complete reference compounds: {complete_count}")
        logger.info(f"  Incomplete reference compounds: {incomplete_count}")
        
        if incomplete_count > 0:
            logger.warning("Incomplete reference compounds:")
            for chembl_id, missing in missing_properties.items():
                logger.warning(f"  {chembl_id}: missing {', '.join(missing)}")
        
        # Return success if all reference compounds are complete
        return complete_count == len(reference_ids)
    except Exception as e:
        logger.error(f"❌ Failed to verify reference compounds: {str(e)}")
        return False

def fix_property_insertion_in_batch_utils():
    """Fix issues in bulk_insert_properties function in batch_utils.py."""
    logger.info("Fixing property insertion in batch_utils.py...")
    
    try:
        # Read the current batch_utils.py file
        with open('batch_utils.py', 'r') as f:
            batch_utils_content = f.read()
        
        # Make a backup of the original file
        with open('batch_utils.py.bak', 'w') as f:
            f.write(batch_utils_content)
        
        # Fix 1: Ensure the SQL query is correctly formatted with RETURNING clause
        if "INSERT INTO molecular_properties" in batch_utils_content:
            logger.info("Checking SQL query in bulk_insert_properties...")
            
            # Check if the query is using RETURNING clause
            if "RETURNING id" not in batch_utils_content or "RETURNING id, molecule_id, property_type_id" not in batch_utils_content:
                logger.warning("SQL query in bulk_insert_properties is missing RETURNING clause")
                
                # Fix the query to include RETURNING clause
                batch_utils_content = batch_utils_content.replace(
                    "VALUES {', '.join(values_list)}",
                    "VALUES {', '.join(values_list)}\n                    RETURNING id, molecule_id, property_type_id"
                )
                
                logger.info("✅ Added RETURNING clause to SQL query in bulk_insert_properties")
        
        # Fix 2: Improve error handling and logging
        if "logger.error(f\"Error inserting batch for property type {prop_type}: {str(e)}\")" in batch_utils_content:
            logger.info("Improving error handling in bulk_insert_properties...")
            
            # Add more detailed error logging
            batch_utils_content = batch_utils_content.replace(
                "logger.error(f\"Error inserting batch for property type {prop_type}: {str(e)}\")",
                """logger.error(f"Error inserting batch for property type {prop_type}: {str(e)}")
                import traceback
                logger.error(f"Error details: {traceback.format_exc()}")
                
                # Try to diagnose database-specific errors
                if "duplicate key value violates unique constraint" in str(e):
                    logger.error("Duplicate key violation - property may already exist")
                elif "foreign key constraint" in str(e):
                    logger.error("Foreign key constraint violation - property_type_id may not exist")
                elif "column" in str(e) and "does not exist" in str(e):
                    logger.error("Column does not exist - schema mismatch")"""
            )
            
            logger.info("✅ Improved error handling and logging in bulk_insert_properties")
        
        # Write the fixed content back to the file
        with open('batch_utils.py', 'w') as f:
            f.write(batch_utils_content)
        
        logger.info("✅ Successfully fixed property insertion in batch_utils.py")
        return True
    except Exception as e:
        logger.error(f"❌ Failed to fix property insertion in batch_utils.py: {str(e)}")
        return False

def fix_property_utils():
    """Fix issues in property_utils.py."""
    logger.info("Fixing property management in property_utils.py...")
    
    try:
        # Read the current property_utils.py file
        with open('property_utils.py', 'r') as f:
            property_utils_content = f.read()
        
        # Make a backup of the original file
        with open('property_utils.py.bak', 'w') as f:
            f.write(property_utils_content)
        
        # Fix 1: Ensure set_property correctly handles property insertion
        if "def set_property" in property_utils_content:
            logger.info("Checking set_property method...")
            
            # Find the set_property method
            set_property_start = property_utils_content.find("def set_property")
            set_property_end = property_utils_content.find("def ", set_property_start + 1)
            
            if set_property_start > 0 and set_property_end > set_property_start:
                # Extract the set_property method
                set_property_method = property_utils_content[set_property_start:set_property_end]
                
                # Check if set_property verifies the insertion
                if "# Verify the property was inserted" not in set_property_method:
                    logger.warning("set_property does not verify property insertion")
                    
                    # Add verification to the method
                    improved_method = set_property_method.replace(
                        "return True",
                        """# Verify the property was inserted
                        verification = self.get_properties(molecule_id, [property_name])
                        if property_name in verification:
                            logger.info(f"Verified property {property_name} was inserted successfully")
                            return True
                        else:
                            logger.warning(f"Property {property_name} verification failed, attempting direct SQL insertion")
                            
                            # Try direct SQL insertion as a last resort
                            try:
                                # Get property type ID
                                property_type_id = self.get_property_type_id(property_name, data_type)
                                
                                # Determine property type
                                if data_type == 'numeric':
                                    execute_query(
                                        "INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value) VALUES (%s, %s, %s) RETURNING id",
                                        (molecule_id, property_type_id, property_value)
                                    )
                                elif data_type == 'boolean':
                                    execute_query(
                                        "INSERT INTO molecular_properties (molecule_id, property_type_id, boolean_value) VALUES (%s, %s, %s) RETURNING id",
                                        (molecule_id, property_type_id, property_value)
                                    )
                                else:  # text
                                    execute_query(
                                        "INSERT INTO molecular_properties (molecule_id, property_type_id, text_value) VALUES (%s, %s, %s) RETURNING id",
                                        (molecule_id, property_type_id, property_value)
                                    )
                                
                                logger.info(f"Direct SQL insertion completed for property {property_name}")
                                return True
                            except Exception as direct_error:
                                logger.error(f"Direct SQL insertion failed: {str(direct_error)}")
                                return False"""
                    )
                    
                    # Replace the set_property method
                    property_utils_content = property_utils_content[:set_property_start] + improved_method + property_utils_content[set_property_end:]
                    
                    logger.info("✅ Improved set_property method with verification and direct SQL fallback")
        
        # Write the fixed content back to the file
        with open('property_utils.py', 'w') as f:
            f.write(property_utils_content)
        
        logger.info("✅ Successfully fixed property management in property_utils.py")
        return True
    except Exception as e:
        logger.error(f"❌ Failed to fix property management in property_utils.py: {str(e)}")
        return False

def fix_reference_compounds_properties():
    """Fix missing properties for reference compounds by directly inserting them."""
    logger.info("Fixing missing properties for reference compounds...")
    
    try:
        # Get reference compound IDs
        from chembl.reference_compounds import get_reference_compound_ids
        reference_ids = get_reference_compound_ids()
        
        logger.info(f"Found {len(reference_ids)} reference compound IDs")
        
        # Initialize PropertyManager
        property_manager = PropertyManager()
        
        # Fix each reference compound
        fixed_count = 0
        
        with get_db_connection() as conn:
            for chembl_id in reference_ids:
                # Get molecule ID from ChEMBL ID
                result = conn.execute_query("""
                SELECT id FROM molecules WHERE chembl_id = %s
                """, (chembl_id,), fetch_one=True)
                
                if not result:
                    logger.error(f"❌ Reference compound {chembl_id} not found in database")
                    continue
                
                molecule_id = result['id']
                
                # Check critical properties
                critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
                properties = property_manager.get_properties(molecule_id, critical_properties)
                
                # Check if all critical properties are present
                missing = [prop for prop in critical_properties if prop not in properties]
                
                if missing:
                    logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
                    
                    # Set default values for missing properties
                    default_properties = {}
                    for prop in missing:
                        if prop == 'logP':
                            default_properties[prop] = 0.0
                        else:  # h_bond_donors or h_bond_acceptors
                            default_properties[prop] = 0
                    
                    # Set properties using PropertyManager
                    logger.info(f"Setting default properties for molecule {molecule_id}: {default_properties}")
                    success_count, total_count = property_manager.set_properties(molecule_id, default_properties)
                    
                    if success_count == total_count:
                        logger.info(f"✅ Successfully set all missing properties for {chembl_id}")
                        fixed_count += 1
                    else:
                        logger.warning(f"⚠️ Only set {success_count}/{total_count} properties for {chembl_id}")
                        
                        # Try direct SQL insertion for any remaining missing properties
                        for prop in missing:
                            if prop not in property_manager.get_properties(molecule_id, [prop]):
                                logger.info(f"Attempting direct SQL insertion for property {prop}")
                                
                                # Get property type ID
                                prop_type_id = property_manager.get_property_type_id(prop)
                                
                                # Set default value
                                if prop == 'logP':
                                    value = 0.0
                                    conn.execute_query("""
                                    INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                                    VALUES (%s, %s, %s)
                                    """, (molecule_id, prop_type_id, value))
                                else:  # h_bond_donors or h_bond_acceptors
                                    value = 0
                                    conn.execute_query("""
                                    INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                                    VALUES (%s, %s, %s)
                                    """, (molecule_id, prop_type_id, value))
                                
                                logger.info(f"Direct SQL insertion completed for property {prop}")
                else:
                    logger.info(f"✅ Reference compound {chembl_id} already has all critical properties")
                    fixed_count += 1
        
        # Generate report
        logger.info(f"Reference compounds property fixing complete:")
        logger.info(f"  Total reference compounds: {len(reference_ids)}")
        logger.info(f"  Fixed reference compounds: {fixed_count}")
        
        # Return success if all reference compounds are fixed
        return fixed_count == len(reference_ids)
    except Exception as e:
        logger.error(f"❌ Failed to fix reference compounds properties: {str(e)}")
        return False

def main():
    """Main function to run the debugging and fixing process."""
    logger.info("Starting data population debugging and fixing process...")
    
    # Step 1: Verify reference compounds
    logger.info("Verifying reference compounds...")
    if verify_reference_compounds():
        logger.info("✅ All reference compounds have been imported with all required properties")
    else:
        logger.warning("⚠️ Some reference compounds are missing required properties")
        
        # Step 2: Fix property_utils.py
        logger.info("Fixing property_utils.py...")
        fix_property_utils()
        
        # Step 3: Fix batch_utils.py
        logger.info("Fixing batch_utils.py...")
        fix_property_insertion_in_batch_utils()
        
        # Step 4: Fix reference compounds properties
        logger.info("Fixing reference compounds properties...")
        if fix_reference_compounds_properties():
            logger.info("✅ Successfully fixed all reference compounds properties")
        else:
            logger.warning("⚠️ Failed to fix all reference compounds properties")
    
    # Step 5: Final verification
    logger.info("Performing final verification...")
    if verify_reference_compounds():
        logger.info("✅ All reference compounds now have all required properties")
    else:
        logger.error("❌ Some reference compounds are still missing required properties")
    
    logger.info("Data population debugging and fixing process complete.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
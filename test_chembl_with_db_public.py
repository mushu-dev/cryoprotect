#!/usr/bin/env python3
"""
Test script for using the public database module with ChEMBL data.

This script verifies that we can insert and retrieve ChEMBL data
using our db_public module that works around RLS constraints.
"""

import os
import sys
import logging
import uuid
from datetime import datetime
from dotenv import load_dotenv

# Import the public database module
from database import db_public as db

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("test_chembl_db_public")

def setup():
    """Set up the database connection."""
    # Load environment variables
    load_dotenv()
    
    # Initialize database connection
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    logger.info(f"Initializing database connection to {config['host']}:{config['port']}/{config['database']}")
    
    # Initialize database pool
    success = db.init_connection_pool(config=config)
    
    if not success:
        logger.error("Failed to initialize database connection")
        return False
        
    logger.info("Database connection established successfully")
    return True

def test_chembl_molecule_insertion():
    """Test inserting a ChEMBL-like molecule."""
    # Generate a unique identifier for this test
    test_id = f"{datetime.now().strftime('%H%M%S')}"
    molecule_uuid = None

    try:
        # Insert a test molecule that looks like ChEMBL data
        molecule_data = {
            "chembl_id": f"CHEMBL{test_id}",
            "name": f"ChEMBL Test Compound {test_id}",
            "formula": "C21H27NO4",
            "smiles": "COc1ccc(CC(=O)N(C)CCc2ccc(OC)c(OC)c2)cc1OC", 
            "inchi": "InChI=1S/C21H27NO4/c1-22(14-13-17-9-19(25-3)20(26-4)10-17)21(24)12-16-7-8-18(23-2)11-15(16)5-6-22/h5-11,23H,12-14H2,1-4H3", 
            "inchikey": "MMAOIAFVWSWIRW-UHFFFAOYSA-N", 
            "data_source": "ChEMBL",
            "is_public": True,  # Explicitly set is_public to True
            "created_by": uuid.UUID("77777777-7777-7777-7777-777777777777"),  # Set a known creator ID
        }
        
        logger.info(f"Inserting test ChEMBL molecule with ID: {molecule_data['chembl_id']}")
        
        success, molecule_id = db.insert_molecule(molecule_data)
        
        if not success or not molecule_id:
            logger.error("Failed to insert test ChEMBL molecule")
            return False
            
        molecule_uuid = molecule_id
        logger.info(f"Successfully inserted molecule with ID: {molecule_uuid}")
            
        # Immediately query for the inserted record
        # Use transaction to work around RLS issues
        logger.info("Attempting to query the inserted record...")

        with db.transaction() as cursor:
            cursor.execute("SELECT * FROM molecules WHERE id = %s", (molecule_uuid,))
            molecule = cursor.fetchone()

            if not molecule:
                logger.error("Failed to query molecule in transaction")
                return False
            
        logger.info(f"Successfully queried molecule: {molecule['name']}")
        
        # Add some ChEMBL-like properties
        properties = [
            ("MOLECULAR_WEIGHT", 357.45, "Computed"),
            ("ALOGP", 3.571, "Computed"),
            ("HBA", 5, "Computed"),
            ("HBD", 0, "Computed"),
            ("PSA", 55.84, "Computed"),
            ("RTB", 9, "Computed"),
            ("RO5_VIOLATIONS", 0, "Computed")
        ]
        
        logger.info("Adding molecular properties...")
        
        for prop_name, prop_value, source in properties:
            # Create property type if needed
            property_type_id = db.insert_property_type(prop_name)
            
            if not property_type_id:
                logger.error(f"Failed to create property type: {prop_name}")
                continue
                
            # Add property to molecule
            property_id = db.insert_molecular_property(
                molecule_uuid, 
                property_type_id, 
                prop_value, 
                source
            )
            
            if not property_id:
                logger.error(f"Failed to insert property: {prop_name}")
                continue
                
            logger.info(f"Added property {prop_name} = {prop_value}")
            
        # Query all properties for the molecule
        properties = db.get_molecular_properties(molecule_uuid)
        
        if not properties or len(properties) == 0:
            logger.error("Failed to query molecular properties")
            return False
            
        logger.info(f"Successfully queried {len(properties)} properties for molecule")
        
        for prop in properties:
            logger.info(f"Property: {prop.get('property_name', 'unknown')} = "
                       f"{prop.get('numeric_value') or prop.get('text_value') or prop.get('boolean_value')}")
            
        logger.info("ChEMBL molecule test passed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error during test: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False
    finally:
        # Clean up test data
        if molecule_uuid:
            try:
                # Delete properties first
                db.execute_query(
                    "DELETE FROM molecular_properties WHERE molecule_id = %s",
                    (molecule_uuid,)
                )
                
                # Then delete the molecule
                db.execute_query(
                    "DELETE FROM molecules WHERE id = %s",
                    (molecule_uuid,)
                )
                
                logger.info("Test data cleaned up")
            except Exception as e:
                logger.error(f"Error cleaning up test data: {str(e)}")

def main():
    """Run all tests."""
    try:
        # Set up database connection
        if not setup():
            logger.error("Setup failed")
            return 1
            
        # Run the ChEMBL molecule test
        if test_chembl_molecule_insertion():
            logger.info("All tests passed!")
            return 0
        else:
            logger.error("Tests failed")
            return 1
            
    except Exception as e:
        logger.error(f"Unhandled error: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1
        
    finally:
        # Close database connections
        db.close_all_connections()
        logger.info("Database connections closed")

if __name__ == "__main__":
    sys.exit(main())
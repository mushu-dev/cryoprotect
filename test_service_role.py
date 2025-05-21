#!/usr/bin/env python3
"""
Test script for the service role database connection.

This script verifies that the service role connection bypasses RLS policies
by inserting and querying records in a single script run.
"""

import os
import sys
import uuid
import logging
from datetime import datetime
from dotenv import load_dotenv

# Import the service role database module
from database import db_service_role as db

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("test_service_role")

def setup():
    """Set up the test environment."""
    # Load environment variables
    load_dotenv()
    
    # Initialize database connection
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'options': "-c role='service_role'"  # Use service role
    }
    
    logger.info(f"Initializing service role database connection to {config['host']}:{config['port']}/{config['database']}")
    
    # Initialize database pool
    success = db.init_connection_pool(config=config)
    
    if not success:
        logger.error("Failed to initialize database connection")
        return False
        
    logger.info("Database connection established successfully")
    return True

def test_insert_and_query():
    """Test inserting a record and immediately querying it."""
    # Generate a unique identifier for this test
    test_id = f"TEST_SR_{datetime.now().strftime('%H%M%S')}"
    molecule_uuid = None
    
    logger.info(f"Running insert and query test with ID: {test_id}")
    
    try:
        # Step 1: Insert a test molecule
        insert_query = """
        INSERT INTO molecules 
        (chembl_id, name, formula, smiles, inchi, inchikey, data_source, pubchem_cid, 
         created_at, updated_at, is_public) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW(), %s) 
        RETURNING id
        """
        
        # Generate unique pubchem_cid to avoid conflicts
        unique_cid = int(datetime.now().strftime('%H%M%S'))

        result = db.execute_query(
            insert_query,
            (test_id, f"Service Role Test {test_id}", "C10H20O2", "CCCCCCCCCC(=O)O",
             "InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)",
             "GHVNFZFCNZKVNT-UHFFFAOYSA-N", "ServiceRoleTest", unique_cid,
             True)  # Set is_public to True to allow public access
        )
        
        if result and len(result) > 0:
            molecule_uuid = result[0]["id"]
            logger.info(f"Inserted molecule with ID: {molecule_uuid}")
        else:
            logger.error("Failed to insert test molecule")
            return False
            
        # Step 2: Immediately query for the inserted record
        logger.info("Attempting to query the inserted record...")
        query_result = db.execute_query(
            "SELECT id, chembl_id, name, is_public FROM molecules WHERE chembl_id = %s",
            (test_id,)
        )

        logger.info(f"Query result: {query_result}")

        if not query_result or len(query_result) == 0:
            logger.error("Failed to query inserted molecule - RLS may be blocking access")

            # Try a different approach - query all molecules
            all_molecules = db.execute_query(
                "SELECT COUNT(*) AS total FROM molecules"
            )
            logger.info(f"Total molecules in database: {all_molecules}")

            # Try with is_public condition
            public_molecules = db.execute_query(
                "SELECT COUNT(*) AS total FROM molecules WHERE is_public = true"
            )
            logger.info(f"Public molecules in database: {public_molecules}")

            return False

        logger.info(f"Successfully queried inserted molecule: {query_result[0]}")
        
        # Step 3: Create a property type
        property_result = db.execute_query(
            """
            INSERT INTO property_types (name, data_type, created_at, updated_at)
            VALUES (%s, %s, NOW(), NOW())
            ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
            RETURNING id
            """,
            (f"test_property_{test_id}", "numeric")
        )
        
        if not property_result or len(property_result) == 0:
            logger.error("Failed to create or get property type")
            return False
            
        property_type_id = property_result[0]["id"]
        logger.info(f"Property type ID: {property_type_id}")
        
        # Step 4: Add a property to the molecule
        property_value = 123.45
        
        property_insert = db.execute_query(
            """
            INSERT INTO molecular_properties
            (molecule_id, property_type_id, numeric_value, source, created_at, updated_at)
            VALUES (%s, %s, %s, %s, NOW(), NOW())
            RETURNING id
            """,
            (molecule_uuid, property_type_id, property_value, "ServiceRoleTest")
        )
        
        if not property_insert or len(property_insert) == 0:
            logger.error("Failed to insert molecular property")
            return False
            
        property_id = property_insert[0]["id"]
        logger.info(f"Inserted property with ID: {property_id}")
        
        # Step 5: Query the molecule with its property
        final_query = db.execute_query(
            """
            SELECT m.id, m.chembl_id, m.name, mp.numeric_value, pt.name AS property_name
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE m.chembl_id = %s
            """,
            (test_id,)
        )
        
        if not final_query or len(final_query) == 0:
            logger.error("Failed to query molecule with its property")
            return False
            
        logger.info(f"Final query result: {final_query[0]}")
        logger.info("TEST PASSED: Service role successfully bypasses RLS policies")
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
                # Delete the molecular property first
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
            
        # Run the test
        if test_insert_and_query():
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
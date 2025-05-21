#!/usr/bin/env python3
"""
Test script for simplified PubChem data importer

This is a minimal version of the import_pubchem_simplified.py script that uses
the CID-Synonym-test file to test the database connection with just a few compounds.
"""

import sys
import os
from dotenv import load_dotenv
import logging
from database import db, utils

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def init_database():
    """Initialize database connection using configuration from environment variables."""
    # Create config from Supabase settings
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    logger.info(f"Initializing database connection to {config['host']}:{config['port']}/{config['database']}")
    
    # Initialize database pool
    return db.init_connection_pool(config=config)

def get_test_cids():
    """Get test CIDs from the CID-Synonym-test file."""
    test_file = "CID-Synonym-test"
    
    if not os.path.exists(test_file):
        logger.error(f"Test file '{test_file}' not found.")
        return []
        
    with open(test_file, "r") as file:
        cids = [int(line.strip()) for line in file if line.strip().isdigit()]
        
    logger.info(f"Loaded {len(cids)} test CIDs: {cids}")
    return cids

def test_molecule_queries():
    """Test basic molecule queries."""
    # Try to count molecules
    try:
        result = db.execute_query("SELECT COUNT(*) as count FROM molecules")
        if result and len(result) > 0:
            count = result[0]['count']
            logger.info(f"Database contains {count} molecules")
        else:
            logger.warning("Failed to get molecule count")
    except Exception as e:
        logger.error(f"Error counting molecules: {str(e)}")
        return False
        
    # Try to get property types
    try:
        result = db.execute_query("SELECT * FROM property_types LIMIT 10")
        if result:
            logger.info(f"Found {len(result)} property types:")
            for prop in result:
                logger.info(f"  - {prop['name']} ({prop['data_type']})")
        else:
            logger.warning("No property types found")
    except Exception as e:
        logger.error(f"Error getting property types: {str(e)}")
        return False
        
    return True

def test_molecule_insert():
    """Test inserting a test molecule."""
    test_molecule = {
        "pubchem_cid": "999999999",  # Use a CID that's unlikely to exist
        "name": "Test Compound",
        "formula": "C10H20O5",
        "molecular_weight": 220.26,
        "smiles": "CCCCCCCCCC(=O)OC",
        "inchi": "InChI=1S/C11H22O2/c1-3-4-5-6-7-8-9-10-11(12)13-2/h3-10H2,1-2H3",
        "inchikey": "HBKUUUUUUUUUUU-UHFFFAOYSA-N",
        "data_source": "Test"
    }
    
    try:
        # Use transaction for the test to allow rollback
        with db.transaction() as cursor:
            # Insert test molecule
            query = """
                INSERT INTO molecules (
                    pubchem_cid, name, formula, molecular_weight, smiles, inchi, inchikey, 
                    data_source, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
                )
                RETURNING id
            """
            
            params = (
                test_molecule["pubchem_cid"],
                test_molecule["name"],
                test_molecule["formula"],
                test_molecule["molecular_weight"],
                test_molecule["smiles"],
                test_molecule["inchi"],
                test_molecule["inchikey"],
                test_molecule["data_source"]
            )
            
            cursor.execute(query, params)
            result = cursor.fetchone()
            
            if result and 'id' in result:
                molecule_id = result['id']
                logger.info(f"Test molecule inserted with ID: {molecule_id}")
                
                # Add a test property
                cursor.execute(
                    """
                    INSERT INTO property_types (name, data_type, created_at, updated_at)
                    VALUES ('test_property', 'numeric', NOW(), NOW())
                    ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
                    RETURNING id
                    """
                )
                property_type = cursor.fetchone()
                property_type_id = property_type['id']
                
                # Insert property value
                cursor.execute(
                    """
                    INSERT INTO molecular_properties (
                        molecule_id, property_type_id, numeric_value, source, created_at, updated_at
                    ) VALUES (
                        %s, %s, %s, %s, NOW(), NOW()
                    )
                    """,
                    (molecule_id, property_type_id, 123.45, "Test")
                )
                
                logger.info(f"Test property added to molecule {molecule_id}")
                
                # For testing, we'll roll back these changes to keep the database clean
                logger.info("Test successful, rolling back changes")
                # Let the context manager handle the rollback automatically by raising an exception
                raise Exception("Deliberate rollback for testing")
                
            else:
                logger.error("Failed to insert test molecule")
                return False
                
    except Exception as e:
        # Expected exception from deliberate rollback
        if "Deliberate rollback for testing" in str(e):
            logger.info("Transaction test completed successfully with rollback")
            return True
        else:
            logger.error(f"Error in test_molecule_insert: {str(e)}")
            return False
    
    return True

def main():
    """Run tests for the simplified database module with PubChem import."""
    logger.info("Starting PubChem import database test")
    
    # Initialize database
    if not init_database():
        logger.error("Failed to initialize database connection")
        return False
        
    try:
        logger.info("Testing database connection...")
        success, message = db.test_connection()
        if success:
            logger.info(f"Database connection test successful: {message}")
        else:
            logger.error(f"Database connection test failed: {message}")
            return False
            
        # Get list of tables
        tables = db.get_tables()
        if tables:
            logger.info(f"Found {len(tables)} tables: {', '.join(tables[:10])}")
            
            if 'molecules' not in tables:
                logger.error("Molecules table not found in database")
                return False
                
            if 'molecular_properties' not in tables:
                logger.error("Molecular_properties table not found in database")
                return False
        else:
            logger.error("No tables found in database")
            return False
            
        # Test molecule queries
        if not test_molecule_queries():
            logger.error("Molecule query tests failed")
            return False
            
        # Test molecule insertion
        if not test_molecule_insert():
            logger.error("Molecule insertion test failed")
            return False
            
        logger.info("All database tests completed successfully")
        return True
            
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        return False
    finally:
        # Close all database connections
        db.close_all_connections()
        logger.info("Database connections closed")

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
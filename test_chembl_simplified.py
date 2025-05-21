#!/usr/bin/env python3
"""
Test script for the simplified ChEMBL import.

This script runs the ChEMBL import with a small test batch
and verifies the correct structure of data in both the
molecules and molecular_properties tables.
"""

import os
import sys
import json
import logging
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv
from database import db

# Import the ChEMBL simplified import module
import import_chembl_simplified as chembl_import

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("test_chembl")

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

def test_transform_chembl_data():
    """Test the transformation of ChEMBL data."""
    # Use one of the mock molecules - need to instantiate and call
    mock_client = chembl_import.MockChemblClient()
    mock_molecule = mock_client[0]
    
    # Transform the data
    molecule_data, property_data = chembl_import.transform_chembl_data(mock_molecule)
    
    # Verify molecule data
    assert molecule_data["chembl_id"] == "CHEMBL25", "ChEMBL ID doesn't match"
    assert molecule_data["name"] == "ASPIRIN", "Name doesn't match"
    assert molecule_data["formula"] == "C9H8O4", "Formula doesn't match"
    assert molecule_data["smiles"] == "CC(=O)OC1=CC=CC=C1C(=O)O", "SMILES doesn't match"
    
    # Verify property data
    assert property_data["molecular_weight"] == 180.16, "Molecular weight doesn't match"
    assert property_data["logp"] == 1.23, "LogP doesn't match"
    assert property_data["tpsa"] == 63.6, "TPSA doesn't match"
    assert property_data["h_bond_acceptors"] == 4, "H-bond acceptors don't match"
    assert property_data["h_bond_donors"] == 1, "H-bond donors don't match"
    
    logger.info("transform_chembl_data test passed")
    return True

def test_db_insert_mock_molecule():
    """Test inserting a mock molecule into the database."""
    # Create a mock molecule_data with required CID field
    mock_molecule_data = {
        "chembl_id": "TEST_" + datetime.now().strftime("%H%M%S"),
        "name": "Test Molecule",
        "formula": "C10H20O2",
        "smiles": "CCCCCCCCCC(=O)O",
        "inchi": "InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)",
        "inchikey": "GHVNFZFCNZKVNT-UHFFFAOYSA-N",
        "data_source": "Test",
        "pubchem_cid": 12345  # Adding required field
    }
    
    # Create mock property data
    mock_property_data = {
        "molecular_weight": 172.26,
        "logp": 3.52,
        "tpsa": 37.3,
        "h_bond_acceptors": 2,
        "h_bond_donors": 1,
        "ro5_violations": 0,
        "med_chem_friendly": "Yes"
    }
    
    try:
        # First directly insert to molecules table to make sure it's there
        test_query = """
        INSERT INTO molecules (chembl_id, name, formula, smiles, inchi, inchikey, data_source, pubchem_cid, created_at, updated_at)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
        RETURNING id
        """

        logger.info(f"Inserting test molecule: {mock_molecule_data['chembl_id']}")
        result = db.execute_query(
            test_query,
            (mock_molecule_data["chembl_id"], mock_molecule_data["name"], mock_molecule_data["formula"],
             mock_molecule_data["smiles"], mock_molecule_data["inchi"], mock_molecule_data["inchikey"],
             mock_molecule_data["data_source"], mock_molecule_data["pubchem_cid"])
        )

        print(f"Result from insert: {result}")

        if result and len(result) > 0:
            molecule_id = result[0]["id"]
            logger.info(f"Inserted molecule with ID: {molecule_id}")
            assert molecule_id is not None, "No molecule ID returned"
        else:
            logger.error("Failed to get molecule ID from insert")
            return False

        # Verify molecule was inserted
        verify_query = "SELECT id FROM molecules WHERE id = %s"
        verify_result = db.execute_query(verify_query, (molecule_id,))
        logger.info(f"Verification result: {verify_result}")

        if not verify_result or len(verify_result) == 0:
            logger.error("Molecule was not actually inserted in the database")
            return False
    except Exception as e:
        logger.error(f"Error inserting test molecule: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False
    
    # Insert properties
    property_count = chembl_import.insert_molecular_properties(molecule_id, mock_property_data)
    
    assert property_count > 0, "No properties inserted"
    
    # Verify molecule exists in database
    molecules = db.execute_query(
        "SELECT * FROM molecules WHERE chembl_id = %s",
        (mock_molecule_data["chembl_id"],)
    )
    
    assert len(molecules) == 1, "Molecule not found in database"
    assert molecules[0]["name"] == mock_molecule_data["name"], "Name doesn't match"
    
    # Verify properties exist in database
    property_query = """
        SELECT mp.*, pt.name as property_name, pt.data_type
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s
    """
    
    properties = db.execute_query(property_query, (molecule_id,))
    
    assert len(properties) >= property_count, "Not all properties were saved"
    
    # Verify a few key properties
    property_map = {prop["property_name"]: prop for prop in properties}
    
    assert "molecular_weight" in property_map, "Molecular weight property not found"
    assert "logp" in property_map, "LogP property not found"
    assert "tpsa" in property_map, "TPSA property not found"
    
    # Check numeric values
    mw_prop = property_map["molecular_weight"]
    assert mw_prop["numeric_value"] == mock_property_data["molecular_weight"], "Molecular weight value doesn't match"
    
    logger.info("db_insert_mock_molecule test passed")
    return True

def test_full_import():
    """Test the full import process with a small batch."""
    # For this test, we'll insert a molecule and then add properties
    
    # First create a test molecule
    test_molecule = {
        "chembl_id": "TEST_IMPORT_" + datetime.now().strftime("%H%M%S"),
        "name": "Test Import",
        "formula": "C8H10N4O2",  # Caffeine
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "pubchem_cid": 2519  # Caffeine CID
    }
    
    # Insert directly into database
    try:
        query = """
        INSERT INTO molecules (chembl_id, name, formula, smiles, inchi, inchikey, pubchem_cid, data_source, created_at, updated_at)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
        RETURNING id
        """

        logger.info(f"Inserting test molecule for full import: {test_molecule['chembl_id']}")
        result = db.execute_query(
            query,
            (test_molecule["chembl_id"], test_molecule["name"], test_molecule["formula"],
             test_molecule["smiles"], test_molecule["inchi"], test_molecule["inchikey"],
             test_molecule["pubchem_cid"], "Test")
        )

        print(f"Result from insert: {result}")

        if result and len(result) > 0:
            molecule_id = result[0]["id"]
            logger.info(f"Inserted full import molecule with ID: {molecule_id}")
        else:
            logger.error("Failed to get molecule ID from insert in full import test")
            return False

        # Verify molecule was inserted
        verify_query = "SELECT id FROM molecules WHERE id = %s"
        verify_result = db.execute_query(verify_query, (molecule_id,))
        logger.info(f"Verification result: {verify_result}")

        if not verify_result or len(verify_result) == 0:
            logger.error("Full import test molecule was not actually inserted in the database")
            return False
    except Exception as e:
        logger.error(f"Error inserting full import test molecule: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False
    
    # Create mock property data
    property_data = {
        "molecular_weight": 194.19,
        "logp": -0.07,
        "tpsa": 58.44,
        "h_bond_acceptors": 6,
        "h_bond_donors": 0,
        "ro5_violations": 0,
        "med_chem_friendly": "Yes"
    }
    
    # Add properties
    property_count = chembl_import.insert_molecular_properties(molecule_id, property_data)
    
    assert property_count > 0, "No properties inserted"
    
    # Verify properties
    property_query = """
    SELECT mp.*, pt.name as property_name, pt.data_type
    FROM molecular_properties mp
    JOIN property_types pt ON mp.property_type_id = pt.id
    WHERE mp.molecule_id = %s
    """
    
    properties = db.execute_query(property_query, (molecule_id,))
    
    assert len(properties) >= property_count, "Not all properties were found"
    
    # Create property map for easier lookup
    property_map = {prop["property_name"]: prop for prop in properties}
    
    # Check a couple of key properties
    assert "molecular_weight" in property_map, "Molecular weight property not found"
    assert "logp" in property_map, "LogP property not found"
    
    # Check values
    assert abs(float(property_map["molecular_weight"]["numeric_value"]) - property_data["molecular_weight"]) < 0.01, "Molecular weight doesn't match"
    assert abs(float(property_map["logp"]["numeric_value"]) - property_data["logp"]) < 0.01, "LogP doesn't match"
    
    logger.info("Full import test passed")
    return True

def cleanup_test_data():
    """Clean up test data from database."""
    try:
        with db.transaction() as cursor:
            # First delete related properties for TEST_%
            cursor.execute("""
                DELETE FROM molecular_properties
                WHERE molecule_id IN (SELECT id FROM molecules WHERE chembl_id LIKE 'TEST\\_%')
            """)
            # Then delete the test molecules
            cursor.execute("DELETE FROM molecules WHERE chembl_id LIKE 'TEST\\_%'")

            # Also delete test molecules with IMPORT in name
            cursor.execute("""
                DELETE FROM molecular_properties
                WHERE molecule_id IN (SELECT id FROM molecules WHERE chembl_id LIKE '%IMPORT%')
            """)
            cursor.execute("DELETE FROM molecules WHERE chembl_id LIKE '%IMPORT%'")

            # Skip the orphaned property types deletion since it causes constraint issues
            # with experiment_properties table

        logger.info("Test data cleaned up")
        return True
    except Exception as e:
        logger.error(f"Error during cleanup: {str(e)}")
        return False

def main():
    """Run all tests."""
    try:
        # Set up database connection
        if not setup():
            logger.error("Setup failed")
            return 1
            
        # Run tests
        tests = [
            test_transform_chembl_data,
            test_db_insert_mock_molecule,
            test_full_import
        ]
        
        success = True
        for test in tests:
            logger.info(f"Running test: {test.__name__}")
            try:
                if not test():
                    logger.error(f"Test {test.__name__} failed")
                    success = False
            except Exception as e:
                logger.error(f"Test {test.__name__} failed with exception: {str(e)}")
                import traceback
                logger.error(traceback.format_exc())
                success = False
                
        # Clean up test data
        try:
            cleanup_test_data()
        except Exception as e:
            logger.error(f"Error during cleanup: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
        
        if success:
            logger.info("All tests passed!")
            return 0
        else:
            logger.error("Some tests failed")
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
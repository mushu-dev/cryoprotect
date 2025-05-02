#!/usr/bin/env python3
"""
Integration tests for database components.
"""

import os
import unittest
import psycopg2
from unittest.mock import patch

from database.connection_manager import ConnectionManager
from database.utils import (
    execute_query, get_molecule_by_id, insert_molecule,
    set_molecule_property, get_or_create_property_type
)

class TestDatabaseIntegration(unittest.TestCase):
    """Integration tests for database components."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment."""
        # Save original environment
        cls.original_env = dict(os.environ)
        
        # Set test environment variables
        os.environ['DB_CONNECTION_MODE'] = 'local'
        os.environ['LOCAL_DB_HOST'] = 'localhost'
        os.environ['LOCAL_DB_PORT'] = '5432'
        os.environ['LOCAL_DB_NAME'] = 'cryoprotect_test'
        os.environ['LOCAL_DB_USER'] = 'postgres'
        os.environ['LOCAL_DB_PASSWORD'] = 'postgres'
        
        # Initialize test database
        try:
            # Create test database
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='postgres',
                user='postgres',
                password='postgres'
            )
            conn.autocommit = True
            cursor = conn.cursor()
            
            # Drop test database if it exists
            cursor.execute("DROP DATABASE IF EXISTS cryoprotect_test")
            
            # Create test database
            cursor.execute("CREATE DATABASE cryoprotect_test")
            
            # Close connection
            cursor.close()
            conn.close()
            
            # Connect to test database
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='cryoprotect_test',
                user='postgres',
                password='postgres'
            )
            cursor = conn.cursor()
            
            # Create test tables
            cursor.execute("""
                CREATE TABLE molecules (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    name VARCHAR NOT NULL,
                    formula VARCHAR,
                    molecular_weight DOUBLE PRECISION,
                    smiles VARCHAR,
                    inchi VARCHAR,
                    inchi_key VARCHAR,
                    chembl_id VARCHAR,
                    pubchem_cid VARCHAR,
                    data_source VARCHAR,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now()
                )
            """)
            
            cursor.execute("""
                CREATE TABLE property_types (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    name VARCHAR NOT NULL UNIQUE,
                    description VARCHAR,
                    data_type VARCHAR NOT NULL,
                    unit VARCHAR,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now()
                )
            """)
            
            cursor.execute("""
                CREATE TABLE molecular_properties (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    molecule_id UUID REFERENCES molecules(id),
                    property_type_id UUID REFERENCES property_types(id),
                    value DOUBLE PRECISION,
                    source VARCHAR,
                    confidence DOUBLE PRECISION,
                    created_at TIMESTAMPTZ DEFAULT now(),
                    updated_at TIMESTAMPTZ DEFAULT now(),
                    UNIQUE(molecule_id, property_type_id)
                )
            """)
            
            # Commit changes
            conn.commit()
            
            # Close connection
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error setting up test database: {str(e)}")
            raise
    
    @classmethod
    def tearDownClass(cls):
        """Tear down the test environment."""
        # Restore original environment
        os.environ.clear()
        os.environ.update(cls.original_env)
        
        # Drop test database
        try:
            conn = psycopg2.connect(
                host='localhost',
                port=5432,
                database='postgres',
                user='postgres',
                password='postgres'
            )
            conn.autocommit = True
            cursor = conn.cursor()
            
            # Drop test database
            cursor.execute("DROP DATABASE IF EXISTS cryoprotect_test")
            
            # Close connection
            cursor.close()
            conn.close()
        except Exception as e:
            print(f"Error tearing down test database: {str(e)}")
    
    def setUp(self):
        """Set up the test case."""
        # Ensure new ConnectionManager is created
        ConnectionManager._instance = None
    
    def test_connection_manager_local(self):
        """Test connection manager with local adapter."""
        # Get connection manager
        manager = ConnectionManager.get_instance()
        
        # Connect
        result = manager.connect()
        
        # Verify
        self.assertTrue(result)
        self.assertEqual(manager.active_adapter, 'local')
    
    def test_insert_and_retrieve_molecule(self):
        """Test inserting and retrieving a molecule."""
        # Insert molecule
        molecule = insert_molecule(
            name="Test Molecule",
            formula="C10H20O",
            molecular_weight=156.27,
            smiles="CCCCCCCCCC(=O)",
            inchi="InChI=1S/C10H20O/c1-2-3-4-5-6-7-8-9-10-11/h10H,2-9H2,1H3",
            inchi_key="AXSIIYDNMFQWRZ-UHFFFAOYSA-N",
            chembl_id="CHEMBL1234",
            pubchem_cid="123456",
            data_source="test"
        )
        
        # Verify insert
        self.assertIsNotNone(molecule)
        self.assertEqual(molecule['name'], "Test Molecule")
        self.assertEqual(molecule['formula'], "C10H20O")
        
        # Get molecule by ID
        retrieved = get_molecule_by_id(molecule['id'])
        
        # Verify retrieve
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved['id'], molecule['id'])
        self.assertEqual(retrieved['name'], "Test Molecule")
    
    def test_properties(self):
        """Test property types and molecular properties."""
        # Insert molecule
        molecule = insert_molecule(
            name="Property Test Molecule",
            inchi_key="TESTKEY123"
        )
        
        # Create property type
        prop_type = get_or_create_property_type(
            name="logP",
            description="Octanol-water partition coefficient",
            data_type="numeric",
            unit=""
        )
        
        # Verify property type
        self.assertIsNotNone(prop_type)
        self.assertEqual(prop_type['name'], "logP")
        
        # Set property
        property_result = set_molecule_property(
            molecule_id=molecule['id'],
            property_type_id=prop_type['id'],
            value=2.5,
            source="test",
            confidence=0.9
        )
        
        # Verify property
        self.assertIsNotNone(property_result)
        self.assertEqual(property_result['molecule_id'], molecule['id'])
        self.assertEqual(property_result['property_type_id'], prop_type['id'])
        self.assertEqual(property_result['value'], 2.5)

if __name__ == '__main__':
    unittest.main()
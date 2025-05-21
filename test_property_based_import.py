#!/usr/bin/env python3
"""
Test script for property-based PubChem import functionality.

This script tests the key functionality of the property-based PubChem import,
particularly focusing on the property type handling that was fixed in the 
property_based_pubchem_import_fixed.py implementation.
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import unittest
import uuid

# Import the required classes and functions from the property_based_pubchem_import_fixed.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from property_based_pubchem_import_fixed import Database, PROPERTY_TYPE_MAPPING
except ImportError:
    print("Failed to import from property_based_pubchem_import_fixed.py")
    print("Make sure the file exists and is accessible.")
    sys.exit(1)

class TestPropertyBasedImport(unittest.TestCase):
    """Test cases for property-based PubChem import functionality."""

    def setUp(self):
        """Set up test environment."""
        # Database connection parameters (for testing)
        self.db_params = {
            'host': os.environ.get('TEST_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
            'port': os.environ.get('TEST_DB_PORT', '5432'),
            'dbname': os.environ.get('TEST_DB_NAME', 'postgres'),
            'user': os.environ.get('TEST_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
            'password': os.environ.get('TEST_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37')
        }
        
        # Initialize database connection
        self.db = Database(
            self.db_params['host'],
            self.db_params['port'],
            self.db_params['dbname'],
            self.db_params['user'],
            self.db_params['password']
        )
        
        # Create test data
        self.test_molecule = {
            'pubchem_cid': 99999999,  # Using a high CID to avoid conflicts
            'name': 'Test Cryoprotectant',
            'formula': 'C6H14O6',
            'smiles': 'CC(O)C(O)C(O)C(O)CO',
            'inchi': 'InChI=1S/C6H14O6/c1-2(8)3(9)4(10)5(11)6(12)7/h2-12H,1H3',
            'inchi_key': 'FAKEINCHISTRINGTESTONLY',
            'chembl_id': None,
            'properties': []
        }
        
        # Clean up any previous test data
        self.cleanup_test_data()
    
    def tearDown(self):
        """Clean up after test."""
        self.cleanup_test_data()
    
    def cleanup_test_data(self):
        """Remove test data from database."""
        try:
            # Delete test molecule if it exists
            self.db.execute_query(
                "DELETE FROM molecular_properties WHERE molecule_id IN (SELECT id FROM molecules WHERE pubchem_cid = %s)",
                (self.test_molecule['pubchem_cid'],),
                fetch=False
            )
            self.db.execute_query(
                "DELETE FROM molecules WHERE pubchem_cid = %s",
                (self.test_molecule['pubchem_cid'],),
                fetch=False
            )
        except Exception as e:
            print(f"Error cleaning up test data: {e}")
    
    def test_database_connection(self):
        """Test database connection."""
        connection = None
        try:
            connection = self.db.get_connection()
            self.assertIsNotNone(connection)
            self.assertTrue(connection.closed == 0)  # 0 means not closed
        finally:
            if connection and connection.closed == 0:
                connection.close()
    
    def test_ensure_property_types(self):
        """Test ensuring property types with proper data_type values."""
        try:
            # Ensure property types
            property_types = self.db.ensure_property_types()
            
            # Verify property types
            self.assertIsNotNone(property_types)
            self.assertTrue(len(property_types) > 0)
            
            # Verify in database
            connection = self.db.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            for name, type_id in property_types.items():
                cursor.execute(
                    "SELECT id, name, data_type FROM property_types WHERE id = %s",
                    (type_id,)
                )
                result = cursor.fetchone()
                
                self.assertIsNotNone(result)
                self.assertEqual(result['name'], name)
                self.assertIsNotNone(result['data_type'])
                self.assertNotEqual(result['data_type'], '')
                
                # Verify against our mapping
                expected_data_type = PROPERTY_TYPE_MAPPING.get(name, {}).get('data_type')
                if expected_data_type:
                    self.assertEqual(result['data_type'], expected_data_type)
            
            cursor.close()
            connection.close()
            
        except Exception as e:
            self.fail(f"ensure_property_types failed with error: {e}")
    
    def test_molecule_insertion(self):
        """Test molecule insertion."""
        try:
            # Insert test molecule
            molecule_map = self.db.batch_insert_molecules([self.test_molecule])
            
            # Verify molecule was inserted
            self.assertIsNotNone(molecule_map)
            self.assertTrue(self.test_molecule['pubchem_cid'] in molecule_map)
            
            # Verify in database
            connection = self.db.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            cursor.execute(
                "SELECT id, name, pubchem_cid FROM molecules WHERE pubchem_cid = %s",
                (self.test_molecule['pubchem_cid'],)
            )
            result = cursor.fetchone()
            
            self.assertIsNotNone(result)
            self.assertEqual(result['pubchem_cid'], str(self.test_molecule['pubchem_cid']))
            self.assertEqual(result['name'], self.test_molecule['name'])
            
            cursor.close()
            connection.close()
            
        except Exception as e:
            self.fail(f"molecule_insertion failed with error: {e}")
    
    def test_property_insertion(self):
        """Test property insertion with property type IDs."""
        try:
            # Insert test molecule
            molecule_map = self.db.batch_insert_molecules([self.test_molecule])
            molecule_id = molecule_map[self.test_molecule['pubchem_cid']]
            
            # Get property types
            property_types = self.db.ensure_property_types()
            
            # Create test properties
            test_properties = [
                {
                    'molecule_id': molecule_id,
                    'property_name': 'Molecular Weight',
                    'property_value': 182.17,
                    'units': 'g/mol',
                    'source': 'Test'
                },
                {
                    'molecule_id': molecule_id,
                    'property_name': 'XLogP3',
                    'property_value': -2.5,
                    'units': '',
                    'source': 'Test'
                },
                {
                    'molecule_id': molecule_id,
                    'property_name': 'Hydrogen Bond Donor Count',
                    'property_value': 5,
                    'units': '',
                    'source': 'Test'
                }
            ]
            
            # Insert properties
            property_ids = self.db.batch_insert_properties(test_properties, property_types)
            
            # Verify properties were inserted
            self.assertIsNotNone(property_ids)
            self.assertEqual(len(property_ids), len(test_properties))
            
            # Verify in database
            connection = self.db.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            cursor.execute(
                """
                SELECT mp.id, mp.property_name, mp.property_value, pt.name, pt.data_type
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %s
                """,
                (molecule_id,)
            )
            results = cursor.fetchall()
            
            self.assertEqual(len(results), len(test_properties))
            
            for result in results:
                self.assertIsNotNone(result['property_name'])
                self.assertIsNotNone(result['property_value'])
                self.assertIsNotNone(result['data_type'])
                self.assertNotEqual(result['data_type'], '')
            
            cursor.close()
            connection.close()
            
        except Exception as e:
            self.fail(f"property_insertion failed with error: {e}")

if __name__ == '__main__':
    unittest.main()
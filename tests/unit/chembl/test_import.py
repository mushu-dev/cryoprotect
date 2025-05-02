#!/usr/bin/env python3
"""
Unit tests for ChEMBL_Integrated_Import.py

These tests verify that the database operations in ChEMBL_Integrated_Import.py
correctly use the SupabaseDirectConnection instead of MCP.
"""

import unittest
from unittest.mock import patch, MagicMock, call
import json
from datetime import datetime

# Import the module to test
import ChEMBL_Integrated_Import as chembl_import
from supabase_direct import SupabaseDirectConnection


class TestChEMBLIntegratedImport(unittest.TestCase):
    """Test cases for ChEMBL_Integrated_Import.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a mock for SupabaseDirectConnection
        self.mock_db = MagicMock(spec=SupabaseDirectConnection)
        
        # Mock the get_db_connection function to return our mock
        patcher = patch.object(chembl_import, 'get_db_connection', return_value=self.mock_db)
        self.addCleanup(patcher.stop)
        self.mock_get_db = patcher.start()

    def test_insert_property_type(self):
        """Test that insert_property_type uses SupabaseDirectConnection"""
        # Set up the mock to return a valid result
        self.mock_db.execute_query.return_value = [{'id': 'mock-property-type-id'}]
        
        # Call the function
        result = chembl_import.insert_property_type(
            name="Test Property",
            data_type="numeric",
            description="Test description",
            units="mg/L"
        )
        
        # Verify the result
        self.assertEqual(result, 'mock-property-type-id')
        
        # Verify that execute_query was called with the correct parameters
        self.mock_db.execute_query.assert_called_once()
        args, kwargs = self.mock_db.execute_query.call_args
        
        # Check that the SQL contains the correct table and operation
        sql = args[0]
        self.assertIn("INSERT INTO property_types", sql)
        self.assertIn("RETURNING id", sql)
        
        # Check that the parameters contain the correct values
        params = kwargs
        self.assertEqual(params['name'], "Test Property")
        self.assertEqual(params['data_type'], "numeric")
        self.assertEqual(params['description'], "Test description")
        self.assertEqual(params['units'], "mg/L")

    def test_get_property_types(self):
        """Test that get_property_types uses SupabaseDirectConnection"""
        # Set up the mock to return a valid result
        self.mock_db.execute_query.return_value = [
            {'id': 'prop-type-1', 'name': 'Property 1'},
            {'id': 'prop-type-2', 'name': 'Property 2'}
        ]
        
        # Call the function
        result = chembl_import.get_property_types()
        
        # Verify the result
        expected = {
            'property 1': 'prop-type-1',
            'property 2': 'prop-type-2'
        }
        self.assertEqual(result, expected)
        
        # Verify that execute_query was called with the correct parameters
        self.mock_db.execute_query.assert_called_once()
        args, kwargs = self.mock_db.execute_query.call_args
        
        # Check that the SQL contains the correct table and operation
        sql = args[0]
        self.assertIn("SELECT id, name FROM property_types", sql)

    def test_check_molecule_exists(self):
        """Test that check_molecule_exists uses SupabaseDirectConnection"""
        # Set up the mock to return a valid result
        self.mock_db.execute_query.return_value = [{'id': 'molecule-id-123'}]
        
        # Call the function
        result = chembl_import.check_molecule_exists("INCHIKEY123")
        
        # Verify the result
        self.assertEqual(result, 'molecule-id-123')
        
        # Verify that execute_query was called with the correct parameters
        self.mock_db.execute_query.assert_called_once()
        args, kwargs = self.mock_db.execute_query.call_args
        
        # Check that the SQL contains the correct table and operation
        sql = args[0]
        self.assertIn("SELECT id FROM molecules WHERE inchikey", sql)
        
        # Check that the parameters contain the correct values
        params = kwargs
        self.assertEqual(params['inchikey'], "INCHIKEY123")

    def test_insert_molecule(self):
        """Test that insert_molecule uses SupabaseDirectConnection"""
        # Set up the mock to return a valid result
        self.mock_db.execute_query.return_value = [{'id': 'new-molecule-id'}]
        
        # Create test molecule data
        molecule_data = {
            'name': 'Test Molecule',
            'smiles': 'C1=CC=CC=C1',
            'inchi': 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H',
            'inchikey': 'UHOVQNZJYSORNB-UHFFFAOYSA-N',
            'formula': 'C6H6',
            'molecular_weight': 78.11,
            'created_by': 'test-user',
            'data_source': 'ChEMBL Test',
            'version': 1,
            'modification_history': json.dumps([{
                'timestamp': datetime.now().isoformat(),
                'action': 'created',
                'user_id': 'test-user'
            }])
        }
        
        # Call the function
        result = chembl_import.insert_molecule(molecule_data)
        
        # Verify the result
        self.assertEqual(result, 'new-molecule-id')
        
        # Verify that execute_query was called with the correct parameters
        self.mock_db.execute_query.assert_called_once()
        args, kwargs = self.mock_db.execute_query.call_args
        
        # Check that the SQL contains the correct table and operation
        sql = args[0]
        self.assertIn("INSERT INTO molecules", sql)
        self.assertIn("RETURNING id", sql)
        
        # Check that the parameters contain the correct values
        params = kwargs
        for key, value in molecule_data.items():
            self.assertEqual(params[key], value)

    def test_insert_property(self):
        """Test that insert_property uses SupabaseDirectConnection"""
        # Set up the mock to return a valid result
        self.mock_db.execute_query.return_value = [{'id': 'new-property-id'}]
        
        # Create test property data
        property_data = {
            'id': 'test-property-id',
            'molecule_id': 'test-molecule-id',
            'property_type_id': 'test-property-type-id',
            'numeric_value': 42.0,
            'text_value': None,
            'boolean_value': None,
            'created_by': 'test-user',
            'data_source': 'ChEMBL Test',
            'version': 1,
            'modification_history': json.dumps([{
                'timestamp': datetime.now().isoformat(),
                'action': 'created',
                'user_id': 'test-user'
            }])
        }
        
        # Call the function
        result = chembl_import.insert_property(property_data)
        
        # Verify the result
        self.assertTrue(result)
        
        # Verify that execute_query was called with the correct parameters
        self.mock_db.execute_query.assert_called_once()
        args, kwargs = self.mock_db.execute_query.call_args
        
        # Check that the SQL contains the correct table and operation
        sql = args[0]
        self.assertIn("INSERT INTO molecular_properties", sql)
        self.assertIn("RETURNING id", sql)
        
        # Check that the parameters contain the correct values
        params = kwargs
        # Only non-None values should be included
        filtered_data = {k: v for k, v in property_data.items() if v is not None}
        for key, value in filtered_data.items():
            self.assertEqual(params[key], value)


if __name__ == '__main__':
    unittest.main()
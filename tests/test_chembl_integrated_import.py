#!/usr/bin/env python3
"""
Unit tests for ChEMBL_Integrated_Import.py

These tests verify the changes made to the ChEMBL_Integrated_Import.py script
as part of the ChEMBL data remediation plan.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock, call
import json
from datetime import datetime
import argparse

# Add parent directory to path to import the module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the module to test
import ChEMBL_Integrated_Import as chembl_import


class TestChEMBLIntegratedImport(unittest.TestCase):
    """Test cases for ChEMBL_Integrated_Import.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Mock compound data
        self.mock_compound = {
            'molecule_chembl_id': 'CHEMBL1234',
            'pref_name': 'Test Compound',
            'molecule_structures': {
                'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                'standard_inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
            },
            'molecule_properties': {
                'full_molformula': 'C9H8O4',
                'full_mwt': 180.16,
                'alogp': 1.31,
                'hba': 4,
                'hbd': 1,
                'psa': 63.6,
                'rtb': 3
            }
        }
        
        # Mock user profile ID
        self.user_profile_id = 'test-user-profile-id'
        
        # Mock property type map
        self.property_type_map = {
            'logp': 'prop-type-logp',
            'molecular weight': 'prop-type-mw',
            'hydrogen bond acceptor count': 'prop-type-hba',
            'hydrogen bond donor count': 'prop-type-hbd',
            'topological polar surface area': 'prop-type-psa',
            'rotatable bond count': 'prop-type-rtb'
        }

    @patch('ChEMBL_Integrated_Import.get_project_id_for_mcp')
    @patch('ChEMBL_Integrated_Import.execute_sql')
    def test_database_operations_use_mcp(self, mock_execute_sql, mock_get_project_id):
        """Test that database operations use MCP tools"""
        # Setup mocks
        mock_get_project_id.return_value = 'test-project-id'
        mock_execute_sql.return_value = [{'id': 'test-id'}]
        
        # Test insert_property_type
        result = chembl_import.insert_property_type('Test Property', 'numeric', 'Test description')
        self.assertEqual(result, 'test-id')
        self.assertTrue(mock_execute_sql.called)
        
        # Test get_property_types
        mock_execute_sql.reset_mock()
        mock_execute_sql.return_value = [
            {'id': 'prop-id-1', 'name': 'Property 1'},
            {'id': 'prop-id-2', 'name': 'Property 2'}
        ]
        result = chembl_import.get_property_types()
        self.assertEqual(result, {'property 1': 'prop-id-1', 'property 2': 'prop-id-2'})
        self.assertTrue(mock_execute_sql.called)
        
        # Test check_molecule_exists
        mock_execute_sql.reset_mock()
        mock_execute_sql.return_value = [{'id': 'test-id'}]
        result = chembl_import.check_molecule_exists('TESTINCHIKEY')
        self.assertEqual(result, 'test-id')
        self.assertTrue(mock_execute_sql.called)
        
        # Test insert_molecule
        mock_execute_sql.reset_mock()
        molecule_data = {
            'name': 'Test Molecule',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
            'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
            'chembl_id': 'CHEMBL1234'
        }
        result = chembl_import.insert_molecule(molecule_data)
        self.assertEqual(result, 'test-id')
        self.assertTrue(mock_execute_sql.called)
        
        # Test insert_property
        mock_execute_sql.reset_mock()
        property_data = {
            'id': 'prop-id',
            'molecule_id': 'molecule-id',
            'property_type_id': 'prop-type-id',
            'numeric_value': 1.23,
            'data_source': 'ChEMBL: CHEMBL1234, property: alogp'
        }
        result = chembl_import.insert_property(property_data)
        self.assertTrue(result)
        self.assertTrue(mock_execute_sql.called)

    def test_transform_chembl_to_molecule_adds_chembl_id(self):
        """Test that transform_chembl_to_molecule adds chembl_id field"""
        # Call the function
        result = chembl_import.transform_chembl_to_molecule(self.mock_compound, self.user_profile_id)
        
        # Verify chembl_id is added
        self.assertIn('chembl_id', result)
        self.assertEqual(result['chembl_id'], 'CHEMBL1234')

    def test_transform_chembl_to_properties_tracks_property_source(self):
        """Test that transform_chembl_to_properties tracks property source"""
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.mock_compound, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify property source includes property name
        for prop in properties:
            self.assertIn('data_source', prop)
            self.assertIn('property:', prop['data_source'])

    def test_default_import_limit_increased(self):
        """Test that default import limit is increased to 2000+"""
        # Get the limit argument directly from the source code
        import inspect
        source = inspect.getsource(chembl_import.main)
        
        # Check if the default limit is set to 2500
        self.assertIn('default=2500', source)

    @patch('ChEMBL_Integrated_Import.new_client')
    def test_standard_reference_compounds_fetched(self, mock_new_client):
        """Test that standard reference compounds are fetched"""
        # Setup mock
        mock_molecule = MagicMock()
        mock_new_client.molecule = mock_molecule
        
        # Mock the get method to return a compound
        mock_molecule.get.return_value = self.mock_compound
        
        # Call the fetch_compound_by_id function
        result = chembl_import.fetch_compound_by_id('CHEMBL1234')
        
        # Verify the function called molecule.get with the correct ID
        mock_molecule.get.assert_called_with('CHEMBL1234')
        
        # Verify the result is the mock compound
        self.assertEqual(result, self.mock_compound)
        
        # Verify standard reference IDs are defined
        expected_reference_ids = [
            "CHEMBL25",    # Aspirin
            "CHEMBL1118",  # Caffeine
            "CHEMBL1234",  # Glycerol (common cryoprotectant)
            "CHEMBL444",   # Glucose
            "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
            "CHEMBL9335",  # Dimethyl sulfoxide (DMSO, common cryoprotectant)
            "CHEMBL15151"  # Trehalose (common cryoprotectant)
        ]
        
        # Verify they match
        self.assertEqual(chembl_import.standard_reference_ids, expected_reference_ids)


if __name__ == '__main__':
    unittest.main()
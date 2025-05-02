#!/usr/bin/env python3
"""
Unit tests for chemical class filtering using SMARTS patterns in chembl_search_utils.py.

These tests verify the functionality of the identify_compounds_by_chemical_class,
store_chemical_class_compounds, and get_stored_chemical_class_compounds functions.
"""

import unittest
from unittest.mock import patch, MagicMock, ANY
import json
import os
import tempfile
import logging

# Configure logging to avoid cluttering test output
logging.basicConfig(level=logging.ERROR)

# Import the functions to test
from chembl_search_utils import (
    identify_compounds_by_chemical_class,
    store_chemical_class_compounds,
    get_stored_chemical_class_compounds
)

class TestChemblSearchSmarts(unittest.TestCase):
    """Test cases for chemical class filtering using SMARTS patterns."""

    def setUp(self):
        """Set up test fixtures."""
        # Create sample compound data for testing
        self.sample_compounds = {
            'polyols': [
                {
                    'molecule_chembl_id': 'CHEMBL1234',
                    'canonical_smiles': 'C(C(CO)O)O',
                    'pref_name': 'Glycerol',
                    'chemical_class': 'polyols',
                    'smarts_pattern': '[OX2H][CX4][CX4][OX2H]',
                    'molecule_properties': {
                        'full_mwt': 92.09,
                        'alogp': -1.76,
                        'hba': 3,
                        'hbd': 3,
                        'max_phase': 4
                    }
                },
                {
                    'molecule_chembl_id': 'CHEMBL5678',
                    'canonical_smiles': 'C(C(C(C(C(CO)O)O)O)O)O',
                    'pref_name': 'Sorbitol',
                    'chemical_class': 'polyols',
                    'smarts_pattern': '[OX2H][CX4][CX4][OX2H]',
                    'molecule_properties': {
                        'full_mwt': 182.17,
                        'alogp': -3.10,
                        'hba': 6,
                        'hbd': 6,
                        'max_phase': 4
                    }
                }
            ],
            'amides': [
                {
                    'molecule_chembl_id': 'CHEMBL9012',
                    'canonical_smiles': 'CN(C)C(=O)C',
                    'pref_name': 'Dimethylacetamide',
                    'chemical_class': 'amides',
                    'smarts_pattern': '[NX3][CX3]=[OX1]',
                    'molecule_properties': {
                        'full_mwt': 87.12,
                        'alogp': -0.77,
                        'hba': 1,
                        'hbd': 0,
                        'max_phase': 2
                    }
                }
            ],
            'sulfoxides': [
                {
                    'molecule_chembl_id': 'CHEMBL3456',
                    'canonical_smiles': 'CS(=O)C',
                    'pref_name': 'Dimethyl sulfoxide',
                    'chemical_class': 'sulfoxides',
                    'smarts_pattern': '[#16X3]=[OX1]',
                    'molecule_properties': {
                        'full_mwt': 78.13,
                        'alogp': -0.54,
                        'hba': 1,
                        'hbd': 0,
                        'max_phase': 4
                    }
                }
            ],
            'sugars': [
                {
                    'molecule_chembl_id': 'CHEMBL7890',
                    'canonical_smiles': 'C(C1C(C(C(C(O1)O)O)O)O)O',
                    'pref_name': 'Glucose',
                    'chemical_class': 'sugars',
                    'smarts_pattern': '[OX2H][CX4][CX4][CX4][OX2H]',
                    'molecule_properties': {
                        'full_mwt': 180.16,
                        'alogp': -2.93,
                        'hba': 6,
                        'hbd': 5,
                        'max_phase': 4
                    }
                }
            ]
        }

        # Create a temporary cache directory for testing
        self.temp_cache_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Tear down test fixtures."""
        # Clean up temporary cache directory
        if os.path.exists(self.temp_cache_dir):
            for file in os.listdir(self.temp_cache_dir):
                os.remove(os.path.join(self.temp_cache_dir, file))
            os.rmdir(self.temp_cache_dir)

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_identify_compounds_by_chemical_class(self, mock_client_class):
        """Test identifying compounds by chemical class."""
        # Set up mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Set up mock cache
        mock_cache = MagicMock()
        mock_client.cache = mock_cache
        mock_cache.get.return_value = None  # No cached data
        
        # Set up mock response for each chemical class
        mock_responses = {}
        
        for class_name, compounds in self.sample_compounds.items():
            # Create a mock response for this class
            mock_response = {
                'molecules': []
            }
            
            # Add molecules to the response
            for compound in compounds:
                molecule = {
                    'molecule_chembl_id': compound['molecule_chembl_id'],
                    'molecule_structures': {
                        'canonical_smiles': compound['canonical_smiles']
                    },
                    'pref_name': compound['pref_name'],
                    'molecule_properties': compound['molecule_properties']
                }
                mock_response['molecules'].append(molecule)
            
            mock_responses[class_name] = mock_response
        
        # Configure mock client to return appropriate response for each request
        def mock_make_request(endpoint, params):
            # Return a different response based on which chemical class we're processing
            # In a real implementation, this would be based on the SMARTS pattern
            # For testing, we'll just cycle through the classes
            class_index = mock_client._make_request.call_count % len(self.sample_compounds)
            class_name = list(self.sample_compounds.keys())[class_index]
            return mock_responses[class_name]
        
        mock_client._make_request.side_effect = mock_make_request
        
        # Call the function
        result = identify_compounds_by_chemical_class(use_cache=False)
        
        # Verify the results
        self.assertEqual(len(result), 4)  # Four chemical classes
        
        # Check that each class has the expected number of compounds
        for class_name, compounds in self.sample_compounds.items():
            self.assertIn(class_name, result)
            self.assertEqual(len(result[class_name]), len(compounds))
        
        # Verify that the client was called for each chemical class
        self.assertEqual(mock_client._make_request.call_count, 4)
        
        # Verify that the cache was checked for each class
        self.assertEqual(mock_cache.get.call_count, 4)
        
        # Verify that the results were cached
        self.assertEqual(mock_cache.set.call_count, 4)

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_identify_compounds_with_cache(self, mock_client_class):
        """Test identifying compounds with cache hit."""
        # Set up mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Set up mock cache to return cached data
        mock_cache = MagicMock()
        mock_client.cache = mock_cache
        
        # Configure cache to return data for each class
        def mock_cache_get(key):
            for class_name, compounds in self.sample_compounds.items():
                if class_name in key:
                    return compounds
            return None
        
        mock_cache.get.side_effect = mock_cache_get
        
        # Call the function
        result = identify_compounds_by_chemical_class(use_cache=True)
        
        # Verify the results
        self.assertEqual(len(result), 4)  # Four chemical classes
        
        # Check that each class has the expected number of compounds
        for class_name, compounds in self.sample_compounds.items():
            self.assertIn(class_name, result)
            self.assertEqual(len(result[class_name]), len(compounds))
        
        # Verify that the client was not called (all cache hits)
        mock_client._make_request.assert_not_called()
        
        # Verify that the cache was checked for each class
        self.assertEqual(mock_cache.get.call_count, 4)
        
        # Verify that no results were cached (already in cache)
        mock_cache.set.assert_not_called()

    @patch('chembl_search_utils.safe_transaction')
    def test_store_chemical_class_compounds(self, mock_safe_transaction):
        """Test storing compounds by chemical class."""
        # Set up mock transaction and connection
        mock_transaction = MagicMock()
        mock_connection = MagicMock()
        mock_transaction.__enter__.return_value = mock_transaction
        mock_transaction.connection = mock_connection
        mock_safe_transaction.return_value = mock_transaction
        
        # Configure mock connection to indicate table doesn't exist
        mock_connection.execute_query.return_value = [{'exists': False}]
        
        # Flatten sample compounds for testing
        compounds_by_class = self.sample_compounds
        
        # Call the function
        result = store_chemical_class_compounds(compounds_by_class)
        
        # Verify the results
        self.assertEqual(len(result), 4)  # Four chemical classes
        
        # Check that each class has the expected number of compounds stored
        for class_name, compounds in self.sample_compounds.items():
            self.assertIn(class_name, result)
            self.assertEqual(result[class_name], len(compounds))
        
        # Verify that the table existence was checked
        mock_connection.execute_query.assert_any_call(ANY)
        
        # Verify that the table was created
        create_table_calls = [
            call for call in mock_connection.execute_query.call_args_list
            if "CREATE TABLE" in str(call)
        ]
        self.assertEqual(len(create_table_calls), 1)
        
        # Verify that compounds were inserted
        total_compounds = sum(len(compounds) for compounds in self.sample_compounds.values())
        insert_calls = [
            call for call in mock_connection.execute_query.call_args_list
            if "INSERT INTO" in str(call)
        ]
        self.assertEqual(len(insert_calls), total_compounds)

    @patch('chembl_search_utils.get_db_connection')
    def test_get_stored_chemical_class_compounds(self, mock_get_db_connection):
        """Test retrieving stored compounds by chemical class."""
        # Set up mock connection
        mock_connection = MagicMock()
        mock_connection.__enter__.return_value = mock_connection
        mock_get_db_connection.return_value = mock_connection
        
        # Configure mock connection to return sample data
        mock_connection.execute_query.return_value = self.sample_compounds['polyols']
        
        # Call the function with a specific class
        result = get_stored_chemical_class_compounds(chemical_class='polyols', limit=10)
        
        # Verify the results
        self.assertEqual(len(result), 2)  # Two polyol compounds
        
        # Verify that the query was executed with the correct parameters
        mock_connection.execute_query.assert_called_once_with(ANY, ('polyols',))
        
        # Call the function without a specific class
        mock_connection.execute_query.reset_mock()
        mock_connection.execute_query.return_value = [
            compound for compounds in self.sample_compounds.values() for compound in compounds
        ]
        
        result = get_stored_chemical_class_compounds(limit=10)
        
        # Verify the results
        self.assertEqual(len(result), 5)  # Total of 5 compounds across all classes
        
        # Verify that the query was executed without class-specific parameters
        mock_connection.execute_query.assert_called_once_with(ANY, None)

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_error_handling(self, mock_client_class):
        """Test error handling in identify_compounds_by_chemical_class."""
        # Set up mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Set up mock cache
        mock_cache = MagicMock()
        mock_client.cache = mock_cache
        mock_cache.get.return_value = None  # No cached data
        
        # Configure mock client to raise an exception
        mock_client._make_request.side_effect = Exception("API error")
        
        # Call the function
        result = identify_compounds_by_chemical_class(use_cache=False, fallback_to_cache=False)
        
        # Verify the results
        self.assertEqual(len(result), 4)  # Four chemical classes
        
        # Check that each class has an empty list due to the error
        for class_name in self.sample_compounds.keys():
            self.assertIn(class_name, result)
            self.assertEqual(len(result[class_name]), 0)
        
        # Verify that the client was called for each chemical class
        self.assertEqual(mock_client._make_request.call_count, 4)
        
        # Verify that the cache was checked for each class
        self.assertEqual(mock_cache.get.call_count, 4)

if __name__ == '__main__':
    unittest.main()
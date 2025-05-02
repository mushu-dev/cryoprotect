#!/usr/bin/env python3
"""
Unit tests for the structural similarity search functionality in chembl_search_utils.py.

These tests verify that the find_similar_compounds function correctly identifies
compounds that are structurally similar to reference compounds based on SMILES.
"""

import unittest
from unittest.mock import patch, MagicMock
import json
import os
import sys
from typing import Dict, List, Any

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Mock the db_connection_utils module to avoid indentation errors
sys.modules['db_connection_utils'] = MagicMock()

# Now we can import our module
from chembl_search_utils import find_similar_compounds, store_similar_compounds


class TestChEMBLSimilaritySearch(unittest.TestCase):
    """Test cases for ChEMBL similarity search functionality."""

    def setUp(self):
        """Set up test fixtures."""
        # Sample SMILES for testing
        self.glycerol_smiles = "C(C(CO)O)O"  # Glycerol
        self.dmso_smiles = "CS(=O)C"  # DMSO
        self.ethylene_glycol_smiles = "C(CO)O"  # Ethylene glycol
        
        # Sample response data
        self.mock_similar_compounds = {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL123",
                    "molecule_structures": {
                        "canonical_smiles": "C(C(CO)O)O"
                    },
                    "pref_name": "Glycerol",
                    "similarity": 100.0,
                    "molecule_properties": {
                        "full_mwt": 92.09,
                        "alogp": -1.76,
                        "hba": 3,
                        "hbd": 3
                    }
                },
                {
                    "molecule_chembl_id": "CHEMBL456",
                    "molecule_structures": {
                        "canonical_smiles": "C(CO)O"
                    },
                    "pref_name": "Ethylene glycol",
                    "similarity": 85.7,
                    "molecule_properties": {
                        "full_mwt": 62.07,
                        "alogp": -1.2,
                        "hba": 2,
                        "hbd": 2
                    }
                },
                {
                    "molecule_chembl_id": "CHEMBL789",
                    "molecule_structures": {
                        "canonical_smiles": "C(C(C)O)O"
                    },
                    "pref_name": "Propylene glycol",
                    "similarity": 80.1,
                    "molecule_properties": {
                        "full_mwt": 76.09,
                        "alogp": -0.92,
                        "hba": 2,
                        "hbd": 2
                    }
                }
            ]
        }

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_find_similar_compounds_single_smiles(self, mock_client_class):
        """Test finding similar compounds with a single SMILES string."""
        # Set up the mock
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client._make_request.return_value = self.mock_similar_compounds
        mock_client.cache.get.return_value = None
        
        # Call the function
        result = find_similar_compounds(
            reference_smiles=self.glycerol_smiles,
            similarity_threshold=70,
            limit=10,
            use_cache=False
        )
        
        # Verify the results
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0]["molecule_chembl_id"], "CHEMBL123")
        self.assertEqual(result[0]["similarity"], 100.0)
        self.assertEqual(result[1]["molecule_chembl_id"], "CHEMBL456")
        self.assertEqual(result[1]["similarity"], 85.7)
        
        # Verify the API call
        mock_client._make_request.assert_called_once_with(
            "similarity", 
            {"smiles": self.glycerol_smiles, "similarity": 70, "limit": 10}
        )

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_find_similar_compounds_multiple_smiles(self, mock_client_class):
        """Test finding similar compounds with multiple SMILES strings."""
        # Set up the mock
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Different responses for different SMILES
        glycerol_response = {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL123",
                    "molecule_structures": {"canonical_smiles": "C(C(CO)O)O"},
                    "pref_name": "Glycerol",
                    "similarity": 100.0,
                    "molecule_properties": {"full_mwt": 92.09}
                }
            ]
        }
        
        dmso_response = {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL999",
                    "molecule_structures": {"canonical_smiles": "CS(=O)C"},
                    "pref_name": "DMSO",
                    "similarity": 100.0,
                    "molecule_properties": {"full_mwt": 78.13}
                },
                {
                    "molecule_chembl_id": "CHEMBL888",
                    "molecule_structures": {"canonical_smiles": "CS(=O)(=O)C"},
                    "pref_name": "DMSO2",
                    "similarity": 85.0,
                    "molecule_properties": {"full_mwt": 94.13}
                }
            ]
        }
        
        # Configure mock to return different responses for different calls
        mock_client._make_request.side_effect = [glycerol_response, dmso_response]
        mock_client.cache.get.return_value = None
        
        # Call the function with multiple SMILES
        result = find_similar_compounds(
            reference_smiles=[self.glycerol_smiles, self.dmso_smiles],
            similarity_threshold=70,
            limit=10,
            use_cache=False
        )
        
        # Verify the results
        self.assertEqual(len(result), 3)  # 3 unique compounds
        
        # Check that compounds from both searches are included
        chembl_ids = [mol["molecule_chembl_id"] for mol in result]
        self.assertIn("CHEMBL123", chembl_ids)
        self.assertIn("CHEMBL999", chembl_ids)
        self.assertIn("CHEMBL888", chembl_ids)
        
        # Verify the API calls
        self.assertEqual(mock_client._make_request.call_count, 2)
        mock_client._make_request.assert_any_call(
            "similarity", 
            {"smiles": self.glycerol_smiles, "similarity": 70, "limit": 10}
        )
        mock_client._make_request.assert_any_call(
            "similarity", 
            {"smiles": self.dmso_smiles, "similarity": 70, "limit": 10}
        )

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_find_similar_compounds_with_cache(self, mock_client_class):
        """Test that caching works correctly for similarity searches."""
        # Set up the mock
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Configure cache hit
        cached_result = [
            {
                "molecule_chembl_id": "CHEMBL123",
                "canonical_smiles": "C(C(CO)O)O",
                "pref_name": "Glycerol",
                "similarity": 100.0,
                "molecule_properties": {"full_mwt": 92.09}
            }
        ]
        mock_client.cache.get.return_value = cached_result
        
        # Call the function
        result = find_similar_compounds(
            reference_smiles=self.glycerol_smiles,
            similarity_threshold=70,
            limit=10,
            use_cache=True
        )
        
        # Verify the results
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["molecule_chembl_id"], "CHEMBL123")
        
        # Verify that the API was not called
        mock_client._make_request.assert_not_called()
        
        # Verify that the cache was checked
        mock_client.cache.get.assert_called_once()

    @patch('chembl_search_utils.ResilientChEMBLClient')
    def test_find_similar_compounds_error_handling(self, mock_client_class):
        """Test error handling in find_similar_compounds."""
        # Set up the mock
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        
        # Configure the mock to raise an exception
        mock_client._make_request.side_effect = Exception("API error")
        mock_client.cache.get.return_value = None
        
        # Call the function
        result = find_similar_compounds(
            reference_smiles=self.glycerol_smiles,
            similarity_threshold=70,
            limit=10,
            use_cache=True,
            fallback_to_cache=False
        )
        
        # Verify that an empty list is returned on error
        self.assertEqual(result, [])

    @patch('chembl_search_utils.safe_transaction')
    def test_store_similar_compounds(self, mock_transaction):
        """Test storing similar compounds in the database."""
        # Set up the mock
        mock_conn = MagicMock()
        mock_transaction_context = MagicMock()
        mock_transaction_context.__enter__.return_value = mock_transaction_context
        mock_transaction_context.connection = mock_conn
        mock_transaction.return_value = mock_transaction_context
        
        # Mock the table check to return that the table doesn't exist
        mock_conn.execute_query.return_value = [{'exists': False}]
        
        # Sample compounds to store
        compounds = [
            {
                "molecule_chembl_id": "CHEMBL123",
                "canonical_smiles": "C(C(CO)O)O",
                "pref_name": "Glycerol",
                "similarity": 100.0,
                "molecule_properties": {
                    "full_mwt": 92.09,
                    "alogp": -1.76,
                    "hba": 3,
                    "hbd": 3
                }
            },
            {
                "molecule_chembl_id": "CHEMBL456",
                "canonical_smiles": "C(CO)O",
                "pref_name": "Ethylene glycol",
                "similarity": 85.7,
                "molecule_properties": {
                    "full_mwt": 62.07,
                    "alogp": -1.2,
                    "hba": 2,
                    "hbd": 2
                }
            }
        ]
        
        # Call the function
        result = store_similar_compounds(compounds)
        
        # Verify the results
        self.assertEqual(result, 2)  # 2 compounds stored
        
        # Verify that the table check query was executed
        # Note: We can't use exact string comparison due to whitespace differences
        table_check_call = False
        for call in mock_conn.execute_query.call_args_list:
            args = call[0]
            if len(args) > 0 and isinstance(args[0], str) and "SELECT EXISTS" in args[0] and "similar_compounds" in args[0]:
                table_check_call = True
                break
        self.assertTrue(table_check_call, "Table check query was not executed")
        
        # Verify that the table creation query was executed
        create_table_call = False
        for call in mock_conn.execute_query.call_args_list:
            args = call[0]
            if len(args) > 0 and isinstance(args[0], str) and "CREATE TABLE similar_compounds" in args[0]:
                create_table_call = True
                break
        self.assertTrue(create_table_call, "Table creation query was not executed")
        
        # Verify that insert queries were executed for each compound
        self.assertEqual(mock_conn.execute_query.call_count, 4)  # 1 table check + 1 create + 2 inserts


if __name__ == '__main__':
    unittest.main()
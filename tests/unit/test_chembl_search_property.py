"""
Unit tests for property-based cryoprotectant identification.

Tests the find_potential_cryoprotectants function in chembl_search_utils.py.
"""

import unittest
import json
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock

import pytest

from chembl_search_utils import (
    find_potential_cryoprotectants,
    store_cryoprotectant_candidates,
    get_stored_cryoprotectant_candidates
)
from chembl.client import ResilientChEMBLClient


class TestPropertyBasedCryoprotectantIdentification(unittest.TestCase):
    """Test the property-based cryoprotectant identification functionality."""

    def setUp(self):
        """Set up the test environment."""
        # Create a temporary directory for cache
        self.temp_dir = tempfile.mkdtemp()
        
        # Sample test data
        self.test_molecules = {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL1001",
                    "pref_name": "Glycerol",
                    "molecule_structures": {
                        "canonical_smiles": "C(C(CO)O)O"
                    },
                    "molecule_properties": {
                        "full_mwt": 92.09,
                        "alogp": -1.76,
                        "hba": 3,
                        "hbd": 3,
                        "psa": 60.69
                    },
                    "max_phase": 4
                },
                {
                    "molecule_chembl_id": "CHEMBL1002",
                    "pref_name": "Dimethyl sulfoxide",
                    "molecule_structures": {
                        "canonical_smiles": "CS(=O)C"
                    },
                    "molecule_properties": {
                        "full_mwt": 78.13,
                        "alogp": -0.62,
                        "hba": 1,
                        "hbd": 0,
                        "psa": 17.07
                    },
                    "max_phase": 3
                },
                {
                    "molecule_chembl_id": "CHEMBL1003",
                    "pref_name": "Ethylene glycol",
                    "molecule_structures": {
                        "canonical_smiles": "C(CO)O"
                    },
                    "molecule_properties": {
                        "full_mwt": 62.07,
                        "alogp": -1.36,
                        "hba": 2,
                        "hbd": 2,
                        "psa": 40.46
                    },
                    "max_phase": 1
                }
            ]
        }

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary directory
        shutil.rmtree(self.temp_dir)

    @patch('chembl.client.ResilientChEMBLClient._make_request')
    def test_find_potential_cryoprotectants(self, mock_make_request):
        """Test that find_potential_cryoprotectants correctly queries and processes results."""
        # Mock the API response
        mock_make_request.return_value = self.test_molecules
        
        # Call the function
        results = find_potential_cryoprotectants(limit=10)
        
        # Verify that the API was called with the correct parameters
        mock_make_request.assert_called_once()
        args, kwargs = mock_make_request.call_args
        self.assertEqual(args[0], "molecule")
        self.assertEqual(kwargs["params"]["mw_freebase__gte"], 30.0)
        self.assertEqual(kwargs["params"]["mw_freebase__lte"], 500.0)
        self.assertEqual(kwargs["params"]["alogp__gte"], -3.0)
        self.assertEqual(kwargs["params"]["alogp__lte"], 3.0)
        self.assertEqual(kwargs["params"]["hba__gte"], 2)
        self.assertEqual(kwargs["params"]["hba__lte"], 20)
        self.assertEqual(kwargs["params"]["hbd__gte"], 1)
        self.assertEqual(kwargs["params"]["hbd__lte"], 10)
        self.assertEqual(kwargs["params"]["max_phase__gte"], 0)
        self.assertEqual(kwargs["params"]["limit"], 10)
        
        # Verify the results
        self.assertEqual(len(results), 3)
        self.assertEqual(results[0]["molecule_chembl_id"], "CHEMBL1001")
        self.assertEqual(results[0]["pref_name"], "Glycerol")
        self.assertEqual(results[0]["canonical_smiles"], "C(C(CO)O)O")
        self.assertEqual(results[0]["molecule_properties"]["full_mwt"], 92.09)

    @patch('chembl.client.ResilientChEMBLClient._make_request')
    def test_find_potential_cryoprotectants_with_custom_parameters(self, mock_make_request):
        """Test that find_potential_cryoprotectants accepts custom property ranges."""
        # Mock the API response
        mock_make_request.return_value = self.test_molecules
        
        # Call the function with custom parameters
        results = find_potential_cryoprotectants(
            limit=5,
            min_mw=50.0,
            max_mw=200.0,
            min_logp=-2.0,
            max_logp=1.0,
            min_hba=1,
            max_hba=5,
            min_hbd=0,
            max_hbd=5,
            min_phase=1
        )
        
        # Verify that the API was called with the correct parameters
        mock_make_request.assert_called_once()
        args, kwargs = mock_make_request.call_args
        self.assertEqual(kwargs["params"]["mw_freebase__gte"], 50.0)
        self.assertEqual(kwargs["params"]["mw_freebase__lte"], 200.0)
        self.assertEqual(kwargs["params"]["alogp__gte"], -2.0)
        self.assertEqual(kwargs["params"]["alogp__lte"], 1.0)
        self.assertEqual(kwargs["params"]["hba__gte"], 1)
        self.assertEqual(kwargs["params"]["hba__lte"], 5)
        self.assertEqual(kwargs["params"]["hbd__gte"], 0)
        self.assertEqual(kwargs["params"]["hbd__lte"], 5)
        self.assertEqual(kwargs["params"]["max_phase__gte"], 1)
        self.assertEqual(kwargs["params"]["limit"], 5)

    @patch('chembl.client.ResilientChEMBLClient._make_request')
    @patch('chembl.client.ResilientChEMBLClient.cache.get')
    @patch('chembl.client.ResilientChEMBLClient.cache.set')
    def test_caching_behavior(self, mock_cache_set, mock_cache_get, mock_make_request):
        """Test that results are properly cached and retrieved from cache."""
        # Mock cache miss, then API response
        mock_cache_get.return_value = None
        mock_make_request.return_value = self.test_molecules
        
        # First call should hit the API
        results1 = find_potential_cryoprotectants(limit=10)
        
        # Verify API was called and cache was set
        mock_make_request.assert_called_once()
        mock_cache_set.assert_called_once()
        
        # Reset mocks
        mock_make_request.reset_mock()
        mock_cache_set.reset_mock()
        
        # Mock cache hit
        mock_cache_get.return_value = results1
        
        # Second call should hit the cache
        results2 = find_potential_cryoprotectants(limit=10)
        
        # Verify API was not called and cache was not set
        mock_make_request.assert_not_called()
        mock_cache_set.assert_not_called()
        
        # Verify both results are the same
        self.assertEqual(results1, results2)

    @patch('chembl.client.ResilientChEMBLClient._make_request')
    def test_error_handling(self, mock_make_request):
        """Test that errors are handled gracefully."""
        # Mock an API error
        mock_make_request.side_effect = Exception("API Error")
        
        # Call the function
        results = find_potential_cryoprotectants(limit=10)
        
        # Verify that an empty list is returned on error
        self.assertEqual(results, [])

    @patch('chembl.client.ResilientChEMBLClient._make_request')
    @patch('chembl.client.ResilientChEMBLClient.cache.get')
    def test_fallback_to_cache(self, mock_cache_get, mock_make_request):
        """Test that function falls back to cache when API fails."""
        # Mock cached data
        cached_data = [
            {
                "molecule_chembl_id": "CHEMBL1001",
                "canonical_smiles": "C(C(CO)O)O",
                "pref_name": "Glycerol",
                "molecule_properties": {"full_mwt": 92.09}
            }
        ]
        
        # Mock API error but cache hit
        mock_make_request.side_effect = Exception("API Error")
        mock_cache_get.return_value = cached_data
        
        # Call the function
        results = find_potential_cryoprotectants(limit=10, fallback_to_cache=True)
        
        # Verify that cached data is returned
        self.assertEqual(results, cached_data)

    @patch('db_connection_utils.safe_transaction')
    def test_store_cryoprotectant_candidates(self, mock_safe_transaction):
        """Test storing cryoprotectant candidates in the database."""
        # Mock the transaction context manager
        mock_transaction = MagicMock()
        mock_connection = MagicMock()
        mock_transaction.__enter__.return_value = mock_transaction
        mock_transaction.connection = mock_connection
        mock_safe_transaction.return_value = mock_transaction
        
        # Mock the execute_query method
        mock_connection.execute_query.side_effect = [
            [{"exists": True}],  # Table exists
            None,  # Insert successful
            None,  # Insert successful
            None   # Insert successful
        ]
        
        # Test data
        compounds = [
            {
                "molecule_chembl_id": "CHEMBL1001",
                "canonical_smiles": "C(C(CO)O)O",
                "pref_name": "Glycerol",
                "molecule_properties": {
                    "full_mwt": 92.09,
                    "alogp": -1.76,
                    "hba": 3,
                    "hbd": 3
                },
                "max_phase": 4
            },
            {
                "molecule_chembl_id": "CHEMBL1002",
                "canonical_smiles": "CS(=O)C",
                "pref_name": "Dimethyl sulfoxide",
                "molecule_properties": {
                    "full_mwt": 78.13,
                    "alogp": -0.62,
                    "hba": 1,
                    "hbd": 0
                },
                "max_phase": 3
            },
            {
                "molecule_chembl_id": "CHEMBL1003",
                "canonical_smiles": "C(CO)O",
                "pref_name": "Ethylene glycol",
                "molecule_properties": {
                    "full_mwt": 62.07,
                    "alogp": -1.36,
                    "hba": 2,
                    "hbd": 2
                },
                "max_phase": 1
            }
        ]
        
        # Call the function
        result = store_cryoprotectant_candidates(compounds)
        
        # Verify the result
        self.assertEqual(result, 3)
        
        # Verify that execute_query was called correctly
        self.assertEqual(mock_connection.execute_query.call_count, 4)
        
        # First call should check if table exists
        args, _ = mock_connection.execute_query.call_args_list[0]
        self.assertIn("SELECT EXISTS", args[0])
        
        # Subsequent calls should be inserts
        for i in range(1, 4):
            args, kwargs = mock_connection.execute_query.call_args_list[i]
            self.assertIn("INSERT INTO cryoprotectant_candidates", args[0])

    @patch('db_connection_utils.get_db_connection')
    def test_get_stored_cryoprotectant_candidates(self, mock_get_db_connection):
        """Test retrieving stored cryoprotectant candidates from the database."""
        # Mock the connection context manager
        mock_connection = MagicMock()
        mock_connection.__enter__.return_value = mock_connection
        mock_get_db_connection.return_value = mock_connection
        
        # Mock the execute_query method
        mock_results = [
            {
                "id": 1,
                "chembl_id": "CHEMBL1001",
                "name": "Glycerol",
                "smiles": "C(C(CO)O)O",
                "molecular_weight": 92.09,
                "logp": -1.76,
                "hba": 3,
                "hbd": 3,
                "clinical_phase": 4,
                "properties": {"full_mwt": 92.09},
                "created_at": "2025-05-01T15:30:00Z"
            },
            {
                "id": 2,
                "chembl_id": "CHEMBL1002",
                "name": "Dimethyl sulfoxide",
                "smiles": "CS(=O)C",
                "molecular_weight": 78.13,
                "logp": -0.62,
                "hba": 1,
                "hbd": 0,
                "clinical_phase": 3,
                "properties": {"full_mwt": 78.13},
                "created_at": "2025-05-01T15:30:00Z"
            }
        ]
        mock_connection.execute_query.return_value = mock_results
        
        # Call the function
        results = get_stored_cryoprotectant_candidates(limit=10)
        
        # Verify the results
        self.assertEqual(results, mock_results)
        
        # Verify that execute_query was called correctly
        mock_connection.execute_query.assert_called_once()
        args, _ = mock_connection.execute_query.call_args
        self.assertIn("SELECT * FROM cryoprotectant_candidates", args[0])
        self.assertIn("LIMIT 10", args[0])


if __name__ == '__main__':
    unittest.main()
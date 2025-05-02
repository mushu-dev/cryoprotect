#!/usr/bin/env python3
"""
Unit tests for verification utilities
"""

import unittest
import sys
import os
import json
from unittest.mock import patch, MagicMock

# Add parent directory to path to allow importing modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the function to test
from verify_imported_data import verify_property_completeness

class TestVerificationUtils(unittest.TestCase):
    """Test cases for verification utilities"""

    @patch('verify_imported_data.execute_query')
    def test_verify_property_completeness_with_list(self, mock_execute_query):
        """Test verify_property_completeness with list property_counts"""
        # Mock the database query result with property_counts as a list
        mock_result = [{
            'property_counts': [
                {'name': 'logP', 'molecule_count': 100},
                {'name': 'h_bond_donors', 'molecule_count': 80},
                {'name': 'h_bond_acceptors', 'molecule_count': 90}
            ],
            'total_molecules': 100,
            'complete_molecules': 70,
            'molecules_with_properties': 120
        }]
        mock_execute_query.return_value = mock_result
        
        # Call the function
        result = verify_property_completeness()
        
        # Verify the results
        self.assertEqual(result['property_counts']['logP'], 100)
        self.assertEqual(result['property_counts']['h_bond_donors'], 80)
        self.assertEqual(result['property_counts']['h_bond_acceptors'], 90)
        self.assertEqual(result['molecules_with_complete_properties'], 70)
        self.assertEqual(result['molecules_with_incomplete_properties'], 30)
        self.assertEqual(result['total_molecules_with_properties'], 120)
        self.assertEqual(result['property_completeness_percentage'], 70.0)

    @patch('verify_imported_data.execute_query')
    def test_verify_property_completeness_with_integer(self, mock_execute_query):
        """Test verify_property_completeness with integer property_counts"""
        # Mock the database query result with property_counts as an integer
        mock_result = [{
            'property_counts': 250,  # This is the problematic case that caused the TypeError
            'total_molecules': 100,
            'complete_molecules': 70,
            'molecules_with_properties': 120
        }]
        mock_execute_query.return_value = mock_result
        
        # Call the function
        result = verify_property_completeness()
        
        # Verify the results - should handle the integer case gracefully
        self.assertEqual(result['property_counts']['total'], 250)
        self.assertEqual(result['molecules_with_complete_properties'], 70)
        self.assertEqual(result['molecules_with_incomplete_properties'], 30)
        self.assertEqual(result['total_molecules_with_properties'], 120)
        self.assertEqual(result['property_completeness_percentage'], 70.0)

    @patch('verify_imported_data.execute_query')
    def test_verify_property_completeness_with_none(self, mock_execute_query):
        """Test verify_property_completeness with None property_counts"""
        # Mock the database query result with property_counts as None
        mock_result = [{
            'property_counts': None,
            'total_molecules': 100,
            'complete_molecules': 70,
            'molecules_with_properties': 120
        }]
        mock_execute_query.return_value = mock_result
        
        # Call the function
        result = verify_property_completeness()
        
        # Verify the results - should handle the None case gracefully
        self.assertEqual(result['property_counts'], {})
        self.assertEqual(result['molecules_with_complete_properties'], 70)
        self.assertEqual(result['molecules_with_incomplete_properties'], 30)
        self.assertEqual(result['total_molecules_with_properties'], 120)
        self.assertEqual(result['property_completeness_percentage'], 70.0)

    @patch('verify_imported_data.execute_query')
    def test_verify_property_completeness_with_empty_result(self, mock_execute_query):
        """Test verify_property_completeness with empty query result"""
        # Mock an empty database query result
        mock_execute_query.return_value = []
        
        # Call the function
        result = verify_property_completeness()
        
        # Verify the results - should handle empty result gracefully
        self.assertEqual(result['property_counts'], {})
        self.assertEqual(result['molecules_with_complete_properties'], 0)
        self.assertEqual(result['molecules_with_incomplete_properties'], 0)
        self.assertEqual(result['total_molecules_with_properties'], 0)
        self.assertEqual(result['property_completeness_percentage'], 0)

if __name__ == '__main__':
    unittest.main()
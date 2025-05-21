"""
Unit tests for the consolidated molecule utilities.

This module provides unit tests for the consolidated molecule utility functions
without requiring the full Flask application.
"""

import unittest
import os
import sys
import uuid
from unittest.mock import patch, MagicMock

# Add project root to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import the functions we want to test
from api.consolidated_utils import (
    is_consolidated,
    get_primary_molecule,
    get_consolidated_molecules,
    get_differentiation_group,
    get_differentiation_group_members,
    enrich_molecule_data
)

class ConsolidatedUtilsTest(unittest.TestCase):
    """Test case for consolidated molecule utilities."""
    
    def setUp(self):
        """Set up the test environment."""
        # Create mock IDs for testing
        self.primary_id = str(uuid.uuid4())
        self.consolidated_id = str(uuid.uuid4())
        self.differentiated_id = str(uuid.uuid4())
        self.diff_group_id = "test_diff_group"
    
    @patch('api.consolidated_utils.get_db_connection')
    def test_is_consolidated(self, mock_get_conn):
        """Test is_consolidated function."""
        # Set up the mock
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test with a consolidated molecule
        mock_cursor.fetchone.return_value = [True]
        self.assertTrue(is_consolidated(self.consolidated_id))
        
        # Test with a primary molecule
        mock_cursor.fetchone.return_value = [False]
        self.assertFalse(is_consolidated(self.primary_id))
        
        # Test with a non-existent molecule
        mock_cursor.fetchone.return_value = None
        self.assertFalse(is_consolidated("non-existent-id"))
    
    @patch('api.consolidated_utils.get_db_connection')
    def test_get_primary_molecule(self, mock_get_conn):
        """Test get_primary_molecule function."""
        # Set up the mock
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test with a consolidated molecule
        mock_cursor.fetchone.return_value = [self.primary_id]
        self.assertEqual(get_primary_molecule(self.consolidated_id), self.primary_id)
        
        # Test with a primary molecule
        mock_cursor.fetchone.return_value = [self.primary_id]
        self.assertEqual(get_primary_molecule(self.primary_id), self.primary_id)
        
        # Test with a non-existent molecule
        mock_cursor.fetchone.return_value = None
        non_existent_id = "non-existent-id"
        self.assertEqual(get_primary_molecule(non_existent_id), non_existent_id)
    
    @patch('api.consolidated_utils.get_db_connection')
    def test_get_consolidated_molecules(self, mock_get_conn):
        """Test get_consolidated_molecules function."""
        # Set up the mock
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test with a primary molecule that has consolidated molecules
        mock_cursor.fetchall.return_value = [[self.consolidated_id]]
        self.assertEqual(get_consolidated_molecules(self.primary_id), [self.consolidated_id])
        
        # Test with a molecule that has no consolidated molecules
        mock_cursor.fetchall.return_value = []
        self.assertEqual(get_consolidated_molecules(self.consolidated_id), [])
    
    @patch('api.consolidated_utils.get_db_connection')
    def test_get_differentiation_group(self, mock_get_conn):
        """Test get_differentiation_group function."""
        # Set up the mock
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test with a differentiated molecule
        mock_cursor.fetchone.return_value = [self.diff_group_id]
        self.assertEqual(get_differentiation_group(self.differentiated_id), self.diff_group_id)
        
        # Test with a non-differentiated molecule
        mock_cursor.fetchone.return_value = None
        self.assertIsNone(get_differentiation_group(self.primary_id))
    
    @patch('api.consolidated_utils.get_db_connection')
    def test_get_differentiation_group_members(self, mock_get_conn):
        """Test get_differentiation_group_members function."""
        # Set up the mock
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test with a valid differentiation group
        mock_cursor.fetchall.return_value = [[self.differentiated_id]]
        self.assertEqual(get_differentiation_group_members(self.diff_group_id), [self.differentiated_id])
        
        # Test with a non-existent differentiation group
        mock_cursor.fetchall.return_value = []
        self.assertEqual(get_differentiation_group_members("non-existent-group"), [])
    
    @patch('api.consolidated_utils.is_consolidated')
    @patch('api.consolidated_utils.get_primary_molecule')
    @patch('api.consolidated_utils.get_consolidated_molecules')
    @patch('api.consolidated_utils.get_differentiation_group')
    @patch('api.consolidated_utils.get_differentiation_description')
    def test_enrich_molecule_data(self, mock_get_diff_desc, mock_get_diff_group, 
                                mock_get_cons_mols, mock_get_primary, mock_is_cons):
        """Test enrich_molecule_data function."""
        # Set up mocks
        mock_is_cons.return_value = False
        mock_get_primary.return_value = self.primary_id
        mock_get_cons_mols.return_value = [self.consolidated_id]
        mock_get_diff_group.return_value = None
        mock_get_diff_desc.return_value = None
        
        # Test with a primary molecule
        molecule_data = {
            'id': self.primary_id,
            'name': 'Test Molecule'
        }
        enriched = enrich_molecule_data(molecule_data)
        
        self.assertEqual(enriched['id'], self.primary_id)
        self.assertEqual(enriched['is_consolidated'], False)
        self.assertEqual(enriched['consolidated_to'], self.primary_id)
        self.assertEqual(enriched['consolidated_molecules'], [self.consolidated_id])
        
        # Test with a consolidated molecule
        mock_is_cons.return_value = True
        molecule_data = {
            'id': self.consolidated_id,
            'name': 'Test Consolidated Molecule'
        }
        enriched = enrich_molecule_data(molecule_data)
        
        self.assertEqual(enriched['id'], self.consolidated_id)
        self.assertEqual(enriched['is_consolidated'], True)
        self.assertEqual(enriched['consolidated_to'], self.primary_id)
        self.assertEqual(enriched['consolidated_molecules'], [])
        
        # Test with a differentiated molecule
        mock_is_cons.return_value = False
        mock_get_diff_group.return_value = self.diff_group_id
        mock_get_diff_desc.return_value = "Test differentiation"
        molecule_data = {
            'id': self.differentiated_id,
            'name': 'Test Differentiated Molecule'
        }
        enriched = enrich_molecule_data(molecule_data)
        
        self.assertEqual(enriched['id'], self.differentiated_id)
        self.assertEqual(enriched['differentiation_group'], self.diff_group_id)
        self.assertEqual(enriched['differentiation_description'], "Test differentiation")
        
        # Test with invalid input
        self.assertEqual(enrich_molecule_data(None), None)
        self.assertEqual(enrich_molecule_data({}), {})
        self.assertEqual(enrich_molecule_data({'not_id': 'no_id'}), {'not_id': 'no_id'})

if __name__ == '__main__':
    unittest.main()
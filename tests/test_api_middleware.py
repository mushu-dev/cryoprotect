"""
Unit tests for the API middleware.

This module provides tests for the API middleware functions that handle
consolidated molecule ID resolution.
"""

import unittest
from unittest.mock import patch, MagicMock
import sys
import os
import json
from flask import Flask, request, jsonify

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.middleware import (
    resolve_molecule_id,
    resolve_molecule_ids,
    get_molecule_with_consolidated_info,
    get_differentiation_group_members,
    handle_consolidated_molecules,
    molecule_batch_middleware
)

class TestApiMiddleware(unittest.TestCase):
    """Test cases for API middleware functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.app = Flask(__name__)
        self.app.config['TESTING'] = True
        
        # Sample molecule IDs for testing
        self.primary_id = "00000000-0000-0000-0000-000000000001"
        self.consolidated_id = "00000000-0000-0000-0000-000000000002"
        self.differentiated_id = "00000000-0000-0000-0000-000000000003"
        self.nonexistent_id = "00000000-0000-0000-0000-000000000099"
        self.differentiation_group = "test_diff_group"
    
    @patch('api.middleware.get_connection')
    def test_resolve_molecule_id(self, mock_get_conn):
        """Test resolve_molecule_id function."""
        # Set up mock connection and cursor
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Test primary molecule (not consolidated)
        mock_cursor.fetchone.return_value = (self.primary_id, None, False)
        resolved_id, is_consolidated = resolve_molecule_id(self.primary_id)
        self.assertEqual(resolved_id, self.primary_id)
        self.assertFalse(is_consolidated)
        
        # Test consolidated molecule
        mock_cursor.fetchone.return_value = (self.consolidated_id, self.primary_id, True)
        resolved_id, is_consolidated = resolve_molecule_id(self.consolidated_id)
        self.assertEqual(resolved_id, self.primary_id)
        self.assertTrue(is_consolidated)
        
        # Test nonexistent molecule
        mock_cursor.fetchone.return_value = None
        resolved_id, is_consolidated = resolve_molecule_id(self.nonexistent_id)
        self.assertEqual(resolved_id, self.nonexistent_id)
        self.assertFalse(is_consolidated)
        
        # Test invalid UUID
        resolved_id, is_consolidated = resolve_molecule_id("not-a-uuid")
        self.assertEqual(resolved_id, "not-a-uuid")
        self.assertFalse(is_consolidated)
        
        # Test None input
        resolved_id, is_consolidated = resolve_molecule_id(None)
        self.assertIsNone(resolved_id)
        self.assertFalse(is_consolidated)
    
    @patch('api.middleware.get_connection')
    def test_resolve_molecule_ids(self, mock_get_conn):
        """Test resolve_molecule_ids function."""
        # Set up mock connection and cursor
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Set up mock cursor results
        mock_cursor.fetchall.return_value = [
            (self.primary_id, None, False),
            (self.consolidated_id, self.primary_id, True)
        ]
        
        # Test resolving multiple molecule IDs
        molecule_ids = [self.primary_id, self.consolidated_id, self.nonexistent_id, "not-a-uuid"]
        results = resolve_molecule_ids(molecule_ids)
        
        # Verify primary molecule is unchanged
        self.assertEqual(results[self.primary_id], (self.primary_id, False))
        
        # Verify consolidated molecule is resolved to primary
        self.assertEqual(results[self.consolidated_id], (self.primary_id, True))
        
        # Verify nonexistent molecule remains unchanged
        self.assertEqual(results[self.nonexistent_id], (self.nonexistent_id, False))
        
        # Verify invalid UUID remains unchanged
        self.assertEqual(results["not-a-uuid"], ("not-a-uuid", False))
        
        # Test with empty list
        results = resolve_molecule_ids([])
        self.assertEqual(results, {})
    
    @patch('api.middleware.resolve_molecule_id')
    @patch('api.middleware.get_connection')
    def test_get_molecule_with_consolidated_info(self, mock_get_conn, mock_resolve):
        """Test get_molecule_with_consolidated_info function."""
        # Set up mock for resolve_molecule_id
        mock_resolve.return_value = (self.primary_id, True)
        
        # Set up mock connection and cursor
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Mock cursor description for column names
        mock_cursor.description = [
            ('id',), ('name',), ('molecular_formula',), ('smiles',), ('inchi_key',),
            ('consolidated_to',), ('consolidated_molecules',), ('created_at',), ('updated_at',)
        ]
        
        # Mock first query result for molecule data
        mock_cursor.fetchone.side_effect = [
            # Molecule data
            (
                self.primary_id, "Test Molecule", "C6H12O6", "C1C(C(C(C(C1O)O)O)O)O", "INCHIKEY",
                None, [self.consolidated_id], "2023-01-01", "2023-01-01"
            ),
            # None for second call (when fetching again)
            None
        ]
        
        # Mock second query result for differentiation data
        mock_cursor.fetchall.return_value = [
            ('differentiationGroup', self.differentiation_group),
            ('differentiationDescription', "Test differentiation")
        ]
        
        # Test getting molecule with consolidated info
        result = get_molecule_with_consolidated_info(self.consolidated_id)
        
        # Verify result contains expected data
        self.assertEqual(result['id'], self.primary_id)
        self.assertEqual(result['name'], "Test Molecule")
        self.assertEqual(result['is_consolidated'], True)
        self.assertEqual(result['original_molecule_id'], self.consolidated_id)
        self.assertEqual(result['consolidated_molecules'], [self.consolidated_id])
        self.assertEqual(result['differentiation_group'], self.differentiation_group)
        self.assertEqual(result['differentiation_description'], "Test differentiation")
        
        # Test with None input
        mock_resolve.return_value = (None, False)
        result = get_molecule_with_consolidated_info(None)
        self.assertIsNone(result)
        
        # Test with nonexistent molecule
        mock_resolve.return_value = (self.nonexistent_id, False)
        mock_cursor.fetchone.return_value = None
        result = get_molecule_with_consolidated_info(self.nonexistent_id)
        self.assertIsNone(result)
    
    @patch('api.middleware.get_connection')
    def test_get_differentiation_group_members(self, mock_get_conn):
        """Test get_differentiation_group_members function."""
        # Set up mock connection and cursor
        mock_cursor = MagicMock()
        mock_conn = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_get_conn.return_value = mock_conn
        
        # Mock cursor description for column names
        mock_cursor.description = [
            ('id',), ('name',), ('molecular_formula',), ('smiles',), ('inchi_key',),
            ('differentiation_description',)
        ]
        
        # Mock query result
        mock_cursor.fetchall.return_value = [
            (
                self.primary_id, "Test Molecule", "C6H12O6", "C1C(C(C(C(C1O)O)O)O)O", "INCHIKEY",
                "Test differentiation"
            ),
            (
                self.differentiated_id, "Test Differentiated", "C6H12O6", "C1C(C(C(C(C1O)O)O)O)O", "INCHIKEY2",
                "Another differentiation"
            )
        ]
        
        # Test getting differentiation group members
        result = get_differentiation_group_members(self.differentiation_group)
        
        # Verify result contains expected data
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]['id'], self.primary_id)
        self.assertEqual(result[0]['name'], "Test Molecule")
        self.assertEqual(result[0]['differentiation_description'], "Test differentiation")
        self.assertEqual(result[1]['id'], self.differentiated_id)
        self.assertEqual(result[1]['name'], "Test Differentiated")
        
        # Test with None input
        result = get_differentiation_group_members(None)
        self.assertEqual(result, [])
        
        # Test with empty group
        mock_cursor.fetchall.return_value = []
        result = get_differentiation_group_members("empty_group")
        self.assertEqual(result, [])
    
    def test_handle_consolidated_molecules_decorator(self):
        """Test handle_consolidated_molecules decorator."""
        with self.app.test_request_context():
            # Mock the resolve_molecule_id function
            with patch('api.middleware.resolve_molecule_id') as mock_resolve:
                mock_resolve.return_value = (self.primary_id, True)
                
                # Create a test view function
                @handle_consolidated_molecules
                def test_view(molecule_id):
                    return jsonify({"id": molecule_id}), 200
                
                # Test with URL parameter
                response = test_view(molecule_id=self.consolidated_id)
                data = json.loads(response[0].data.decode('utf-8'))
                self.assertEqual(data['id'], self.primary_id)
                self.assertEqual(data['meta']['consolidated_from'], self.consolidated_id)
    
    def test_molecule_batch_middleware(self):
        """Test molecule_batch_middleware decorator."""
        with self.app.test_request_context(
            json={"molecule_ids": [self.primary_id, self.consolidated_id]}
        ):
            # Mock the resolve_molecule_ids function
            with patch('api.middleware.resolve_molecule_ids') as mock_resolve:
                mock_resolve.return_value = {
                    self.primary_id: (self.primary_id, False),
                    self.consolidated_id: (self.primary_id, True)
                }
                
                # Create a test view function
                @molecule_batch_middleware
                def test_view():
                    return jsonify({"ids": request.get_json()['molecule_ids']}), 200
                
                # Test batch middleware
                response = test_view()
                data = json.loads(response[0].data.decode('utf-8'))
                self.assertEqual(data['ids'], [self.primary_id, self.primary_id])
                self.assertEqual(
                    data['meta']['consolidated_molecules'],
                    {self.consolidated_id: self.primary_id}
                )

if __name__ == '__main__':
    unittest.main()
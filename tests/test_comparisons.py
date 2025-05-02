"""
Tests for the comparisons module.

This module contains tests for the comparisons module, which provides functionality
for comparing properties of molecules and mixtures.
"""

import unittest
from unittest.mock import patch, MagicMock, Mock
import sys
import os

# Add the parent directory to the path so we can import the api module directly
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the module under test directly
from api.comparisons import get_entity_type_and_data, extract_properties, compare_entities, KEY_PROPERTIES

# Mock the dependencies
sys.modules['api.models'] = Mock()


class TestComparisons(unittest.TestCase):
    """Test cases for the comparisons module."""

    def setUp(self):
        """Set up the test case."""
        # Define some test data
        self.molecule_data = {
            'id': 'mol-123',
            'name': 'Test Molecule',
            'type': 'molecule',
            'properties': [
                {'property_name': 'molecular_weight', 'value': 100.0},
                {'property_name': 'logP', 'value': 1.5},
                {'property_name': 'TPSA', 'value': 50.0},
                {'property_name': 'toxicity', 'value': 'low'},
                {'property_name': 'vitrification_probability', 'value': 0.8}
            ]
        }
        
        self.mixture_data = {
            'id': 'mix-456',
            'name': 'Test Mixture',
            'type': 'mixture',
            'properties': [
                {'property_name': 'molecular_weight', 'value': 150.0},
                {'property_name': 'logP', 'value': 2.0},
                {'property_name': 'vitrification_probability', 'value': 0.9}
            ]
        }

    @patch('api.comparisons.Molecule.get_with_properties')
    def test_get_entity_type_and_data_molecule(self, mock_get_with_properties):
        """Test get_entity_type_and_data with a molecule."""
        # Setup mock
        mock_get_with_properties.return_value = self.molecule_data
        
        # Call the function
        result = get_entity_type_and_data('mol-123')
        
        # Verify the result
        self.assertEqual(result, self.molecule_data)
        self.assertEqual(result['type'], 'molecule')
        mock_get_with_properties.assert_called_once_with('mol-123')

    @patch('api.comparisons.Molecule.get_with_properties')
    @patch('api.comparisons.Mixture.get_with_components')
    def test_get_entity_type_and_data_mixture(self, mock_get_with_components, mock_get_with_properties):
        """Test get_entity_type_and_data with a mixture."""
        # Setup mocks
        mock_get_with_properties.return_value = None
        mock_get_with_components.return_value = self.mixture_data
        
        # Call the function
        result = get_entity_type_and_data('mix-456')
        
        # Verify the result
        self.assertEqual(result, self.mixture_data)
        self.assertEqual(result['type'], 'mixture')
        mock_get_with_properties.assert_called_once_with('mix-456')
        mock_get_with_components.assert_called_once_with('mix-456')

    @patch('api.comparisons.Molecule.get_with_properties')
    @patch('api.comparisons.Mixture.get_with_components')
    def test_get_entity_type_and_data_not_found(self, mock_get_with_components, mock_get_with_properties):
        """Test get_entity_type_and_data with an entity that doesn't exist."""
        # Setup mocks
        mock_get_with_properties.return_value = None
        mock_get_with_components.return_value = None
        
        # Call the function
        result = get_entity_type_and_data('invalid-id')
        
        # Verify the result
        self.assertIsNone(result)
        mock_get_with_properties.assert_called_once_with('invalid-id')
        mock_get_with_components.assert_called_once_with('invalid-id')

    def test_extract_properties_molecule(self):
        """Test extract_properties with a molecule."""
        # Call the function
        result = extract_properties(self.molecule_data)
        
        # Verify the result
        self.assertEqual(result['id'], 'mol-123')
        self.assertEqual(result['name'], 'Test Molecule')
        self.assertEqual(result['type'], 'molecule')
        self.assertEqual(result['molecular_weight'], 100.0)
        self.assertEqual(result['logP'], 1.5)
        self.assertEqual(result['TPSA'], 50.0)
        self.assertEqual(result['toxicity'], 'low')
        self.assertEqual(result['vitrification_probability'], 0.8)

    def test_extract_properties_molecule_with_fallback(self):
        """Test extract_properties with a molecule that has a fallback property."""
        # Create a molecule with a fallback property
        molecule_data = self.molecule_data.copy()
        molecule_data['properties'] = [p for p in molecule_data['properties'] if p['property_name'] != 'molecular_weight']
        molecule_data['molecular_weight'] = 120.0
        
        # Call the function
        result = extract_properties(molecule_data)
        
        # Verify the result
        self.assertEqual(result['molecular_weight'], 120.0)

    def test_extract_properties_mixture(self):
        """Test extract_properties with a mixture."""
        # Call the function
        result = extract_properties(self.mixture_data)
        
        # Verify the result
        self.assertEqual(result['id'], 'mix-456')
        self.assertEqual(result['name'], 'Test Mixture')
        self.assertEqual(result['type'], 'mixture')
        self.assertEqual(result['molecular_weight'], 150.0)
        self.assertEqual(result['logP'], 2.0)
        self.assertEqual(result['vitrification_probability'], 0.9)
        # Check that missing properties are None
        self.assertIsNone(result['TPSA'])
        self.assertIsNone(result['toxicity'])

    def test_extract_properties_mixture_with_fallback(self):
        """Test extract_properties with a mixture that has a fallback property."""
        # Create a mixture with a fallback property
        mixture_data = self.mixture_data.copy()
        mixture_data['properties'] = [p for p in mixture_data['properties'] if p['property_name'] != 'vitrification_probability']
        mixture_data['vitrification_probability'] = 0.95
        
        # Call the function
        result = extract_properties(mixture_data)
        
        # Verify the result
        self.assertEqual(result['vitrification_probability'], 0.95)

    @patch('api.comparisons.get_entity_type_and_data')
    def test_compare_entities(self, mock_get_entity_type_and_data):
        """Test compare_entities with multiple entities."""
        # Setup mock
        mock_get_entity_type_and_data.side_effect = [
            self.molecule_data,
            self.mixture_data
        ]
        
        # Call the function
        result = compare_entities(['mol-123', 'mix-456'])
        
        # Verify the result
        self.assertEqual(len(result['comparison']), 2)
        self.assertEqual(result['comparison'][0]['id'], 'mol-123')
        self.assertEqual(result['comparison'][0]['name'], 'Test Molecule')
        self.assertEqual(result['comparison'][0]['type'], 'molecule')
        self.assertEqual(result['comparison'][0]['molecular_weight'], 100.0)
        
        self.assertEqual(result['comparison'][1]['id'], 'mix-456')
        self.assertEqual(result['comparison'][1]['name'], 'Test Mixture')
        self.assertEqual(result['comparison'][1]['type'], 'mixture')
        self.assertEqual(result['comparison'][1]['molecular_weight'], 150.0)
        
        # Check differences
        self.assertIn('molecular_weight', result['differences'])
        self.assertIn('logP', result['differences'])
        self.assertIn('type', result['differences'])
        
        # Check properties_compared
        self.assertEqual(result['properties_compared'], KEY_PROPERTIES)

    @patch('api.comparisons.get_entity_type_and_data')
    def test_compare_entities_with_not_found(self, mock_get_entity_type_and_data):
        """Test compare_entities with an entity that doesn't exist."""
        # Setup mock
        mock_get_entity_type_and_data.side_effect = [
            self.molecule_data,
            None
        ]
        
        # Call the function
        result = compare_entities(['mol-123', 'invalid-id'])
        
        # Verify the result
        self.assertEqual(len(result['comparison']), 2)
        self.assertEqual(result['comparison'][0]['id'], 'mol-123')
        self.assertEqual(result['comparison'][1]['id'], 'invalid-id')
        
        # Check that all properties for the not found entity are None
        for key in KEY_PROPERTIES:
            self.assertIsNone(result['comparison'][1][key])

    @patch('api.comparisons.get_entity_type_and_data')
    def test_compare_entities_same_properties(self, mock_get_entity_type_and_data):
        """Test compare_entities with entities that have the same properties."""
        # Create two identical molecules
        molecule1 = self.molecule_data.copy()
        molecule2 = self.molecule_data.copy()
        molecule2['id'] = 'mol-789'
        
        # Setup mock
        mock_get_entity_type_and_data.side_effect = [
            molecule1,
            molecule2
        ]
        
        # Call the function
        result = compare_entities(['mol-123', 'mol-789'])
        
        # Verify the result
        self.assertEqual(len(result['comparison']), 2)
        
        # Check differences - should only include id since all properties are the same
        self.assertNotIn('molecular_weight', result['differences'])
        self.assertNotIn('logP', result['differences'])
        self.assertNotIn('type', result['differences'])

    @patch('api.comparisons.get_entity_type_and_data')
    def test_compare_entities_with_empty_list(self, mock_get_entity_type_and_data):
        """Test compare_entities with an empty list."""
        # Call the function with an empty list
        result = compare_entities([])
        
        # Verify the result
        self.assertEqual(len(result['comparison']), 0)
        self.assertEqual(len(result['differences']), 0)
        self.assertEqual(result['properties_compared'], KEY_PROPERTIES)
        
        # Verify that get_entity_type_and_data was not called
        mock_get_entity_type_and_data.assert_not_called()


if __name__ == '__main__':
    unittest.main()
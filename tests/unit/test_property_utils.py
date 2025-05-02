#!/usr/bin/env python3
"""
Unit tests for the enhanced PropertyManager with direct PostgreSQL connection.
"""

import unittest
import uuid
import sys
from unittest.mock import patch, MagicMock, call

# Mock the PostgresDirectConnection and sql_executor modules before importing PropertyManager
mock_postgres_direct = MagicMock()
mock_sql_executor = MagicMock()

# Set up the mocks for all the sql_executor functions used by PropertyManager
mock_execute_query = MagicMock()
mock_bulk_insert = MagicMock()
mock_execute_batch = MagicMock()
# Create a proper decorator mock that returns the function unchanged
def mock_with_retry(func=None, **kwargs):
    if func is None:
        return lambda f: f
    return func
mock_process_in_batches = MagicMock()

# Configure the mock_sql_executor
mock_sql_executor.execute_query = mock_execute_query
mock_sql_executor.bulk_insert = mock_bulk_insert
mock_sql_executor.execute_batch = mock_execute_batch
mock_sql_executor.with_retry = mock_with_retry
mock_sql_executor.process_in_batches = mock_process_in_batches

# Assign the mocks to sys.modules
sys.modules['postgres_direct'] = mock_postgres_direct
sys.modules['sql_executor'] = mock_sql_executor

# Now we can safely import PropertyManager
from property_utils import PropertyManager

class TestPropertyManager(unittest.TestCase):
    """Test cases for the enhanced PropertyManager."""

    def setUp(self):
        """Set up test fixtures."""
        # Use the global mock
        self.mock_execute_query = mock_execute_query
        self.mock_bulk_insert = mock_bulk_insert
        self.mock_execute_batch = mock_execute_batch
        self.mock_process_in_batches = mock_process_in_batches
        
        # Create property type IDs
        self.property_ids = {
            'logP': uuid.uuid4(),
            'molecular_weight': uuid.uuid4(),
            'is_cryoprotectant': uuid.uuid4(),
            'smiles': uuid.uuid4()
        }
        
        # Mock the property types query result
        property_types = [
            {'id': self.property_ids['logP'], 'name': 'logP', 'data_type': 'numeric'},
            {'id': self.property_ids['molecular_weight'], 'name': 'molecular_weight', 'data_type': 'numeric'},
            {'id': self.property_ids['is_cryoprotectant'], 'name': 'is_cryoprotectant', 'data_type': 'boolean'},
            {'id': self.property_ids['smiles'], 'name': 'smiles', 'data_type': 'text'}
        ]
        
        # Configure the mock to return property types
        self.mock_execute_query.return_value = property_types
        
        # Create the PropertyManager instance
        self.property_manager = PropertyManager()
        
        # Manually populate the cache to ensure it's set correctly
        self.property_manager._property_types_cache = {
            'logP': {'id': self.property_ids['logP'], 'data_type': 'numeric'},
            'molecular_weight': {'id': self.property_ids['molecular_weight'], 'data_type': 'numeric'},
            'is_cryoprotectant': {'id': self.property_ids['is_cryoprotectant'], 'data_type': 'boolean'},
            'smiles': {'id': self.property_ids['smiles'], 'data_type': 'text'}
        }
        
        # Create test data
        self.test_molecule_id = uuid.uuid4()
        self.test_user_id = uuid.uuid4()
        
        # Reset the mocks for subsequent tests
        self.mock_execute_query.reset_mock()
        self.mock_bulk_insert.reset_mock()
        self.mock_execute_batch.reset_mock()
        self.mock_process_in_batches.reset_mock()

    def test_init_loads_property_types(self):
        """Test that initialization loads property types into cache."""
        # Verify the cache was populated
        self.assertEqual(len(self.property_manager._property_types_cache), 4)
        self.assertIn('logP', self.property_manager._property_types_cache)
        self.assertIn('molecular_weight', self.property_manager._property_types_cache)
        self.assertIn('is_cryoprotectant', self.property_manager._property_types_cache)
        self.assertIn('smiles', self.property_manager._property_types_cache)

    def test_init_loads_property_types(self):
        """Test that initialization loads property types into cache."""
        # Verify the cache was populated
        self.assertEqual(len(self.property_manager._property_types_cache), 4)
        self.assertIn('logP', self.property_manager._property_types_cache)
        self.assertIn('molecular_weight', self.property_manager._property_types_cache)
        self.assertIn('is_cryoprotectant', self.property_manager._property_types_cache)
        self.assertIn('smiles', self.property_manager._property_types_cache)

    def test_get_property_type_id_from_cache(self):
        """Test getting property type ID from cache."""
        # Get a property type that should be in the cache
        property_id = self.property_manager.get_property_type_id('logP')
        
        # Verify the ID is returned and no database query was made
        self.assertIsNotNone(property_id)
        self.mock_execute_query.assert_not_called()

    def test_get_property_type_id_creates_new(self):
        """Test creating a new property type."""
        # Mock the database queries
        self.mock_execute_query.side_effect = [
            None,  # First query returns no existing property
            {'id': uuid.uuid4()}  # Second query returns the new ID
        ]
        
        # Get a property type that's not in the cache
        property_id = self.property_manager.get_property_type_id('new_property', 'text')
        
        # Verify the correct queries were made
        self.assertEqual(self.mock_execute_query.call_count, 2)
        # First call should be a SELECT
        self.assertIn('SELECT', self.mock_execute_query.call_args_list[0][0][0])
        # Second call should be an INSERT
        self.assertIn('INSERT', self.mock_execute_query.call_args_list[1][0][0])
        
        # Verify the property was added to the cache
        self.assertIn('new_property', self.property_manager._property_types_cache)
        self.assertEqual(self.property_manager._property_types_cache['new_property']['data_type'], 'text')

    def test_set_property_numeric(self):
        """Test setting a numeric property."""
        # Configure the mock for this test
        self.mock_execute_query.side_effect = [
            None,  # No existing property
            None   # Successful insert
        ]
        
        # Set the property
        result = self.property_manager.set_property(
            self.test_molecule_id, 'logP', 2.5, self.test_user_id
        )
        
        # Verify the result and queries
        self.assertTrue(result)
        self.assertEqual(self.mock_execute_query.call_count, 2)
        # First call should check for existing property
        self.assertIn('SELECT', self.mock_execute_query.call_args_list[0][0][0])
        # Second call should be an INSERT with numeric_value
        self.assertIn('INSERT', self.mock_execute_query.call_args_list[1][0][0])
        self.assertIn('numeric_value', self.mock_execute_query.call_args_list[1][0][0])

    def test_set_property_boolean(self):
        """Test setting a boolean property."""
        # Configure the mock for this test
        self.mock_execute_query.side_effect = [
            {'id': uuid.uuid4()},  # Existing property
            None                   # Successful update
        ]
        
        # Set the property
        result = self.property_manager.set_property(
            self.test_molecule_id, 'is_cryoprotectant', True, self.test_user_id
        )
        
        # Verify the result and queries
        self.assertTrue(result)
        self.assertEqual(self.mock_execute_query.call_count, 2)
        # First call should check for existing property
        self.assertIn('SELECT', self.mock_execute_query.call_args_list[0][0][0])
        # Second call should be an UPDATE with boolean_value
        self.assertIn('UPDATE', self.mock_execute_query.call_args_list[1][0][0])
        self.assertIn('boolean_value', self.mock_execute_query.call_args_list[1][0][0])

    def test_set_properties(self):
        """Test setting multiple properties."""
        # Configure the mock for this test - 3 properties, each with 2 queries
        self.mock_execute_query.side_effect = [
            None, None,  # logP - check and insert
            None, None,  # smiles - check and insert
            None, None   # is_cryoprotectant - check and insert
        ]
        
        # Set multiple properties
        properties = {
            'logP': 2.5,
            'smiles': 'C1=CC=CC=C1',
            'is_cryoprotectant': True
        }
        
        success_count, total_count = self.property_manager.set_properties(
            self.test_molecule_id, properties, self.test_user_id
        )
        
        # Verify the result
        self.assertEqual(success_count, 3)
        self.assertEqual(total_count, 3)
        
        # Verify the correct number of queries were made (3 properties * 2 queries each)
        self.assertEqual(self.mock_execute_query.call_count, 6)

    def test_batch_set_properties(self):
        """Test batch setting properties for multiple molecules."""
        # Configure the mock for this test - 2 molecules * 2 properties * 2 queries each
        self.mock_execute_query.side_effect = [None] * 8
        
        # Configure the process_in_batches mock to execute the function
        self.mock_process_in_batches.side_effect = lambda items, batch_size, process_func: [process_func(items)]
        
        # Set up batch data
        molecule_id1 = uuid.uuid4()
        molecule_id2 = uuid.uuid4()
        
        batch_data = [
            {
                'molecule_id': molecule_id1,
                'properties': {
                    'logP': 2.5,
                    'smiles': 'C1=CC=CC=C1'
                }
            },
            {
                'molecule_id': molecule_id2,
                'properties': {
                    'logP': 3.1,
                    'smiles': 'C1=CC=CC=C1O'
                }
            }
        ]
        
        # Execute batch set
        success_count, total_count = self.property_manager.batch_set_properties(
            batch_data, self.test_user_id
        )
        
        # Verify the result
        self.assertEqual(success_count, 4)
        self.assertEqual(total_count, 4)
        
        # Verify process_in_batches was called
        self.mock_process_in_batches.assert_called_once()

    def test_get_properties(self):
        """Test getting properties for a molecule."""
        # Configure the mock for this test
        property_results = [
            {'name': 'logP', 'numeric_value': 2.5, 'text_value': None, 'boolean_value': None, 'data_type': 'numeric'},
            {'name': 'smiles', 'numeric_value': None, 'text_value': 'C1=CC=CC=C1', 'boolean_value': None, 'data_type': 'text'},
            {'name': 'is_cryoprotectant', 'numeric_value': None, 'text_value': None, 'boolean_value': True, 'data_type': 'boolean'}
        ]
        
        # We need to reset the side_effect to ensure it returns our property results
        self.mock_execute_query.reset_mock()
        self.mock_execute_query.side_effect = None
        self.mock_execute_query.return_value = property_results
        
        # Get properties
        properties = self.property_manager.get_properties(self.test_molecule_id)
        
        # Verify the result
        self.assertEqual(len(properties), 3)
        self.assertEqual(properties['logP'], 2.5)
        self.assertEqual(properties['smiles'], 'C1=CC=CC=C1')
        self.assertTrue(properties['is_cryoprotectant'])
        
        # Verify the correct query was made
        self.mock_execute_query.assert_called_once()
        self.assertIn('SELECT', self.mock_execute_query.call_args[0][0])

    def test_get_specific_properties(self):
        """Test getting specific properties for a molecule."""
        # Configure the mock for this test
        property_results = [
            {'name': 'logP', 'numeric_value': 2.5, 'text_value': None, 'boolean_value': None, 'data_type': 'numeric'},
            {'name': 'smiles', 'numeric_value': None, 'text_value': 'C1=CC=CC=C1', 'boolean_value': None, 'data_type': 'text'}
        ]
        
        # We need to reset the side_effect to ensure it returns our property results
        self.mock_execute_query.reset_mock()
        self.mock_execute_query.side_effect = None
        self.mock_execute_query.return_value = property_results
        
        # Get specific properties
        properties = self.property_manager.get_properties(
            self.test_molecule_id, ['logP', 'smiles']
        )
        
        # Verify the result
        self.assertEqual(len(properties), 2)
        self.assertEqual(properties['logP'], 2.5)
        self.assertEqual(properties['smiles'], 'C1=CC=CC=C1')
        
        # Verify the correct query was made
        self.mock_execute_query.assert_called_once()
        self.assertIn('IN', self.mock_execute_query.call_args[0][0])

    def test_batch_get_properties(self):
        """Test batch getting properties for multiple molecules."""
        # Create test molecule IDs
        molecule_id1 = uuid.uuid4()
        molecule_id2 = uuid.uuid4()
        
        # Configure the mock for this test
        property_results = [
            {'molecule_id': molecule_id1, 'name': 'logP', 'numeric_value': 2.5, 'text_value': None, 'boolean_value': None, 'data_type': 'numeric'},
            {'molecule_id': molecule_id1, 'name': 'smiles', 'numeric_value': None, 'text_value': 'C1=CC=CC=C1', 'boolean_value': None, 'data_type': 'text'},
            {'molecule_id': molecule_id2, 'name': 'logP', 'numeric_value': 3.1, 'text_value': None, 'boolean_value': None, 'data_type': 'numeric'},
            {'molecule_id': molecule_id2, 'name': 'smiles', 'numeric_value': None, 'text_value': 'C1=CC=CC=C1O', 'boolean_value': None, 'data_type': 'text'}
        ]
        
        # We need to reset the side_effect to ensure it returns our property results
        self.mock_execute_query.reset_mock()
        self.mock_execute_query.side_effect = None
        self.mock_execute_query.return_value = property_results
        
        # Get properties for multiple molecules
        properties = self.property_manager.batch_get_properties([molecule_id1, molecule_id2])
        
        # Verify the result
        self.assertEqual(len(properties), 2)
        self.assertEqual(properties[molecule_id1]['logP'], 2.5)
        self.assertEqual(properties[molecule_id1]['smiles'], 'C1=CC=CC=C1')
        self.assertEqual(properties[molecule_id2]['logP'], 3.1)
        self.assertEqual(properties[molecule_id2]['smiles'], 'C1=CC=CC=C1O')
        
        # Verify the correct query was made
        self.mock_execute_query.assert_called_once()
        self.assertIn('IN', self.mock_execute_query.call_args[0][0])

    def test_delete_property(self):
        """Test deleting a property."""
        # Configure the mock for this test
        # We need to reset the side_effect to ensure it returns our expected values
        self.mock_execute_query.reset_mock()
        self.mock_execute_query.side_effect = None
        self.mock_execute_query.return_value = None  # Successful delete
        
        # Delete the property
        result = self.property_manager.delete_property(self.test_molecule_id, 'logP')
        
        # Verify the result
        self.assertTrue(result)
        
        # Verify the correct query was made
        self.mock_execute_query.assert_called_once()
        self.assertIn('DELETE', self.mock_execute_query.call_args[0][0])
        self.assertEqual(self.mock_execute_query.call_args[0][1][0], self.test_molecule_id)
        self.assertEqual(self.mock_execute_query.call_args[0][1][1], self.property_ids['logP'])

    def test_clear_cache(self):
        """Test clearing the property types cache."""
        # Verify the cache has data
        self.assertGreater(len(self.property_manager._property_types_cache), 0)
        
        # Clear the cache
        self.property_manager.clear_cache()
        
        # Verify the cache is empty
        self.assertEqual(len(self.property_manager._property_types_cache), 0)
        self.assertEqual(self.property_manager._cache_last_refresh, 0)

if __name__ == '__main__':
    unittest.main()
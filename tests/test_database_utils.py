#!/usr/bin/env python3
"""
Tests for database utilities.
"""

import unittest
from unittest.mock import patch, MagicMock

from database.utils import (
    get_db, with_connection, with_retry, with_transaction,
    execute_query, execute_batch, get_molecule_by_id,
    get_molecule_properties, get_molecules_by_inchikey,
    insert_molecule, update_molecule, set_molecule_property,
    get_or_create_property_type, test_database_connection
)

class TestDatabaseUtils(unittest.TestCase):
    """Test the database utility functions."""
    
    def setUp(self):
        """Set up the test environment."""
        # Mock connection manager
        self.mock_manager = MagicMock()
        self.mock_manager.get_active_adapter.return_value = True
        self.mock_manager.connect.return_value = True
        
        # Create patch for get_db
        self.get_db_patcher = patch('database.utils.get_db', return_value=self.mock_manager)
        self.mock_get_db = self.get_db_patcher.start()
    
    def tearDown(self):
        """Tear down the test environment."""
        self.get_db_patcher.stop()
    
    def test_with_connection_decorator(self):
        """Test the with_connection decorator."""
        # Define test function
        @with_connection
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.mock_manager.get_active_adapter.assert_called_once()
        self.mock_manager.connect.assert_not_called()  # Already connected
    
    def test_with_connection_decorator_not_connected(self):
        """Test the with_connection decorator when not connected."""
        # Set up mock
        self.mock_manager.get_active_adapter.return_value = None
        
        # Define test function
        @with_connection
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.mock_manager.get_active_adapter.assert_called_once()
        self.mock_manager.connect.assert_called_once()
    
    def test_with_retry_decorator_success(self):
        """Test the with_retry decorator with successful execution."""
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
    
    def test_with_retry_decorator_retry_success(self):
        """Test the with_retry decorator with retry success."""
        # Define test counter
        counter = [0]
        
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            counter[0] += 1
            if counter[0] < 3:
                raise ValueError("Test error")
            return "success"
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, "success")
        self.assertEqual(counter[0], 3)
    
    def test_with_retry_decorator_all_fail(self):
        """Test the with_retry decorator with all retries failing."""
        # Define test function
        @with_retry(max_retries=3)
        def test_func():
            raise ValueError("Test error")
        
        # Test function
        with self.assertRaises(ValueError):
            test_func()
    
    def test_with_transaction_decorator(self):
        """Test the with_transaction decorator."""
        # Set up mocks
        mock_transaction = MagicMock()
        self.mock_manager.begin_transaction.return_value = mock_transaction
        
        # Define test function
        @with_transaction
        def test_func(transaction=None):
            return transaction
        
        # Test function
        result = test_func()
        
        # Verify
        self.assertEqual(result, mock_transaction)
        self.mock_manager.begin_transaction.assert_called_once()
        self.mock_manager.commit_transaction.assert_called_once_with(mock_transaction)
        self.mock_manager.rollback_transaction.assert_not_called()
    
    def test_with_transaction_decorator_exception(self):
        """Test the with_transaction decorator with exception."""
        # Set up mocks
        mock_transaction = MagicMock()
        self.mock_manager.begin_transaction.return_value = mock_transaction
        
        # Define test function
        @with_transaction
        def test_func(transaction=None):
            raise ValueError("Test error")
        
        # Test function
        with self.assertRaises(ValueError):
            test_func()
        
        # Verify
        self.mock_manager.begin_transaction.assert_called_once()
        self.mock_manager.commit_transaction.assert_not_called()
        self.mock_manager.rollback_transaction.assert_called_once_with(mock_transaction)
    
    def test_execute_query(self):
        """Test execute_query."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': 1, 'name': 'Test'}]
        
        # Test function
        result = execute_query("SELECT * FROM test")
        
        # Verify
        self.assertEqual(result, [{'id': 1, 'name': 'Test'}])
        self.mock_manager.execute_query.assert_called_once_with("SELECT * FROM test", None)
    
    def test_get_molecule_by_id(self):
        """Test get_molecule_by_id."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': '1', 'name': 'Test Molecule'}]
        
        # Test function
        result = get_molecule_by_id("1")
        
        # Verify
        self.assertEqual(result, {'id': '1', 'name': 'Test Molecule'})
        self.mock_manager.execute_query.assert_called_once()
    
    def test_insert_molecule(self):
        """Test insert_molecule."""
        # Set up mock
        self.mock_manager.execute_query.return_value = [{'id': '1', 'name': 'Test Molecule'}]
        
        # Test function
        result = insert_molecule(
            name="Test Molecule",
            formula="C10H20O",
            molecular_weight=156.27,
            smiles="CCCCCCCCCC(=O)",
            inchi="InChI=1S/C10H20O/c1-2-3-4-5-6-7-8-9-10-11/h10H,2-9H2,1H3",
            inchi_key="AXSIIYDNMFQWRZ-UHFFFAOYSA-N",
            chembl_id="CHEMBL1234",
            pubchem_cid="123456",
            data_source="test"
        )
        
        # Verify
        self.assertEqual(result, {'id': '1', 'name': 'Test Molecule'})
        self.mock_manager.execute_query.assert_called_once()

if __name__ == '__main__':
    unittest.main()
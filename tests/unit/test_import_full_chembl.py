#!/usr/bin/env python3
"""
Unit tests for the enhanced ChEMBL import script.
Tests the direct PostgreSQL connection integration and batch processing.
"""

import unittest
import os
import sys
import tempfile
import json
from unittest.mock import patch, MagicMock, mock_open

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import modules to test
import import_full_chembl
from postgres_direct import PostgresDirectConnection
from sql_executor import execute_query, bulk_insert
from property_utils import PropertyManager


class TestChEMBLImport(unittest.TestCase):
    """Test cases for the enhanced ChEMBL import script."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory for checkpoints
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Mock environment variables
        self.env_patcher = patch.dict('os.environ', {
            'SUPABASE_DB_HOST': 'localhost',
            'SUPABASE_DB_PORT': '5432',
            'SUPABASE_DB_NAME': 'test_db',
            'SUPABASE_DB_USER': 'test_user',
            'SUPABASE_DB_PASSWORD': 'test_password'
        })
        self.env_patcher.start()

    def tearDown(self):
        """Clean up after tests."""
        # Remove temporary directory
        self.temp_dir.cleanup()
        
        # Stop environment variables patch
        self.env_patcher.stop()

    @patch('postgres_direct.PostgresDirectConnection')
    @patch('sql_executor.execute_query')
    def test_database_connection(self, mock_execute_query, mock_postgres_connection):
        """Test that the script correctly connects to the database."""
        # Set up mocks
        mock_db_instance = MagicMock()
        mock_postgres_connection.return_value = mock_db_instance
        mock_execute_query.return_value = [{'connection_test': 1}]
        
        # Call the function that verifies database connection
        with patch('sys.argv', ['import_full_chembl.py', '--dry-run']):
            with patch('import_full_chembl.main', return_value=0):
                # Just import the module to trigger the connection verification
                import import_full_chembl
                
        # Verify that execute_query was called with the test query
        mock_execute_query.assert_called_with("SELECT 1 as connection_test")

    @patch('chembl.worker.ChEMBLWorker')
    @patch('chembl.checkpoint.CheckpointManager')
    @patch('import_full_chembl.ChEMBLClient')
    @patch('sql_executor.execute_query')
    def test_batch_processing(self, mock_execute_query, mock_chembl_client, 
                             mock_checkpoint_manager, mock_worker):
        """Test batch processing of compounds."""
        # Set up mocks
        mock_execute_query.return_value = [{'connection_test': 1}]
        
        mock_client_instance = MagicMock()
        mock_chembl_client.return_value = mock_client_instance
        mock_client_instance.search_compounds.return_value = ['CHEMBL1', 'CHEMBL2', 'CHEMBL3']
        
        mock_checkpoint_instance = MagicMock()
        mock_checkpoint_manager.return_value = mock_checkpoint_instance
        mock_checkpoint_instance.load_checkpoint.return_value = False
        
        mock_worker_instance = MagicMock()
        mock_worker.return_value = mock_worker_instance
        
        # Mock Queue to simulate results
        with patch('queue.Queue') as mock_queue:
            mock_queue_instance = MagicMock()
            mock_queue.return_value = mock_queue_instance
            
            # Simulate results being returned
            mock_queue_instance.get.side_effect = [
                {'status': 'success', 'compound_id': 'CHEMBL1'},
                {'status': 'success', 'compound_id': 'CHEMBL2'},
                {'status': 'error', 'compound_id': 'CHEMBL3', 'error_category': 'API_ERROR'}
            ]
            
            # Run the script with test arguments
            with patch('sys.argv', ['import_full_chembl.py', '--batch-size', '3', 
                                   '--workers', '1', '--dry-run']):
                with patch('import_full_chembl.main', return_value=0):
                    # Import to trigger execution
                    import import_full_chembl
            
            # Verify that the worker was started
            mock_worker_instance.start.assert_called_once()
            
            # Verify that checkpoint was saved
            mock_checkpoint_instance.save_checkpoint.assert_called()

    @patch('property_utils.PropertyManager')
    def test_property_manager_integration(self, mock_property_manager):
        """Test integration with the enhanced PropertyManager."""
        # Set up mock
        mock_manager_instance = MagicMock()
        mock_property_manager.return_value = mock_manager_instance
        
        # Create a test worker
        from chembl.worker import ChEMBLWorker
        worker = ChEMBLWorker(1, MagicMock(), MagicMock())
        
        # Create test data
        transformed_data = {
            'molecule': {
                'name': 'Test Compound',
                'chembl_id': 'CHEMBL1',
                'formula': 'C10H15N',
                'molecular_weight': 149.23,
                'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)N',
                'inchi': 'InChI=1S/C10H15N/c1-8(2)7-9-3-5-10(6-4-9)11-12/h3-6,8,11H,7H2,1-2H3',
                'inchi_key': 'LCTUISCIGMWMAT-UHFFFAOYSA-N',
                'data_source': 'ChEMBL'
            },
            'properties': [
                {
                    'property_name': 'LogP',
                    'property_type': 'physicochemical',
                    'value': 2.5,
                    'unit': 'log units',
                    'source': 'ChEMBL'
                }
            ]
        }
        
        # Mock the execute_query function to return a molecule ID
        with patch('sql_executor.execute_query', return_value={'id': '123e4567-e89b-12d3-a456-426614174000'}):
            # Call the store_compound_data method
            with patch.object(worker, 'store_compound_data', wraps=worker.store_compound_data):
                result = worker.store_compound_data(transformed_data)
                
                # Verify that PropertyManager was used
                self.assertTrue(mock_property_manager.called)
                
                # Verify that the result is the expected molecule ID
                self.assertEqual(result, '123e4567-e89b-12d3-a456-426614174000')


if __name__ == '__main__':
    unittest.main()
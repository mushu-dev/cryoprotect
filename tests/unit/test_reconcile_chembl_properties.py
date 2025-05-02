#!/usr/bin/env python3
"""
Unit tests for reconcile_chembl_properties.py

This script tests the functionality of the ChEMBL properties reconciliation script
with direct PostgreSQL connections.
"""

import os
import sys
import unittest
import json
from unittest.mock import patch, MagicMock, call
from pathlib import Path
from datetime import datetime

# Add parent directory to path to import the script
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the script to test
import reconcile_chembl_properties
import sql_executor

class TestReconcileChEMBLProperties(unittest.TestCase):
    """Test cases for reconcile_chembl_properties.py with direct PostgreSQL connections"""

    def setUp(self):
        """Set up test fixtures"""
        # Set up mock molecule data
        self.mock_pubchem_molecules = {
            'CRYO001': {'name': 'Glycerol', 'pubchem_cid': '962'},
            'CRYO003': {'name': 'beta-Alanine', 'pubchem_cid': '6342'},
            'CRYO005': {'name': 'Urea', 'pubchem_cid': '6057'}
        }
        
        self.mock_chembl_molecules = {
            'CRYO002': {'name': 'DMSO', 'chembl_id': 'CHEMBL1098659'},
            'CRYO004': {'name': 'tert-Butanol', 'chembl_id': 'CHEMBL500033'},
            'CRYO006': {'name': 'Urea', 'chembl_id': 'CHEMBL1487'}
        }
        
        # Mock InChI key map
        self.mock_inchi_key_map = {
            'XSQUKJJJFZCRTK-UHFFFAOYSA-N': ['CRYO005', 'CRYO006']  # Urea
        }
        
        # Mock query results
        self.mock_pubchem_results = [
            {'id': 'CRYO001', 'name': 'Glycerol', 'pubchem_cid': 962},
            {'id': 'CRYO003', 'name': 'beta-Alanine', 'pubchem_cid': 6342},
            {'id': 'CRYO005', 'name': 'Urea', 'pubchem_cid': 6057}
        ]
        
        self.mock_chembl_results = [
            {'id': 'CRYO002', 'name': 'DMSO', 'chembl_id': 'CHEMBL1098659'},
            {'id': 'CRYO004', 'name': 'tert-Butanol', 'chembl_id': 'CHEMBL500033'},
            {'id': 'CRYO006', 'name': 'Urea', 'chembl_id': 'CHEMBL1487'}
        ]
        
        self.mock_inchi_key_results = [
            {'inchi_key': 'XSQUKJJJFZCRTK-UHFFFAOYSA-N', 'molecule_ids': ['CRYO005', 'CRYO006']}
        ]
        
        # Mock identifier manager
        self.mock_id_manager = MagicMock()
        self.mock_id_manager.get_molecule_by_internal_id.side_effect = lambda id: {'name': f'Molecule {id}'} if id in ['CRYO005', 'CRYO006'] else None
        self.mock_id_manager.add_molecule.return_value = True
        self.mock_id_manager.save_identifiers.return_value = True

    @patch('reconcile_chembl_properties.CryoprotectantIdentifierManager.get_instance')
    @patch('sql_executor.execute_query')
    @patch('sql_executor.execute_batch')
    @patch('sql_executor.with_transaction')
    def test_reconcile_cross_references(self, mock_with_transaction, mock_execute_batch, mock_execute_query, mock_get_instance):
        """Test reconciling cross-references between ChEMBL and PubChem"""
        # Set up mocks
        mock_get_instance.return_value = self.mock_id_manager
        
        # Mock the execute_query function to return different results based on the query
        def mock_execute_query_side_effect(query, *args, **kwargs):
            if "pubchem_cid IS NOT NULL" in query:
                return self.mock_pubchem_results
            elif "chembl_id IS NOT NULL" in query:
                return self.mock_chembl_results
            elif "GROUP BY inchikey" in query:
                return self.mock_inchi_key_results
            return None
        
        mock_execute_query.side_effect = mock_execute_query_side_effect
        
        # Mock the with_transaction decorator to execute the function directly
        def mock_transaction_decorator(func):
            def wrapper(*args, **kwargs):
                # Create a mock transaction
                mock_transaction = MagicMock()
                
                # Call the original function with the mock transaction
                return func(mock_transaction, *args, **kwargs)
            return wrapper
        
        mock_with_transaction.side_effect = mock_transaction_decorator
        
        # Call the function
        result = reconcile_chembl_properties.reconcile_cross_references(output_report=None, dry_run=False)
        
        # Verify the result
        self.assertEqual(result['molecules_updated'], 2)
        self.assertEqual(result['cross_references_added'], 2)
        self.assertIn('XSQUKJJJFZCRTK-UHFFFAOYSA-N', result['details'])
        
        # Verify that execute_query was called for PubChem molecules, ChEMBL molecules, and InChI keys
        self.assertEqual(mock_execute_query.call_count, 3)
        
        # Verify that execute_batch was called twice (once for PubChem updates, once for ChEMBL updates)
        self.assertEqual(mock_execute_batch.call_count, 2)
        
        # Verify that the identifier manager was called to update and save identifiers
        self.assertEqual(self.mock_id_manager.add_molecule.call_count, 2)
        self.assertEqual(self.mock_id_manager.save_identifiers.call_count, 1)

    @patch('reconcile_chembl_properties.reconcile_cross_references')
    @patch('sql_executor.close_connections')
    @patch('argparse.ArgumentParser.parse_args')
    def test_main_function(self, mock_parse_args, mock_close_connections, mock_reconcile):
        """Test the main function"""
        # Set up mocks
        mock_args = MagicMock()
        mock_args.report = f"reports/reconciliation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        mock_args.dry_run = False
        mock_args.batch_size = 100
        mock_parse_args.return_value = mock_args
        
        mock_reconcile.return_value = {
            'timestamp': datetime.now().isoformat(),
            'molecules_updated': 2,
            'cross_references_added': 2,
            'conflicts_resolved': 0,
            'details': {'XSQUKJJJFZCRTK-UHFFFAOYSA-N': {'pubchem_molecule': 'CRYO005', 'chembl_molecule': 'CRYO006', 'action': 'cross_references_added'}}
        }
        
        # Call the main function
        with patch.object(sys, 'argv', ['reconcile_chembl_properties.py']):
            result = reconcile_chembl_properties.main()
        
        # Verify that the functions were called
        mock_reconcile.assert_called_once_with(mock_args.report, mock_args.dry_run)
        mock_close_connections.assert_called_once()
        
        # Verify the result
        self.assertEqual(result, 0)

    @patch('reconcile_chembl_properties.CryoprotectantIdentifierManager.get_instance')
    @patch('sql_executor.execute_query')
    def test_dry_run_mode(self, mock_execute_query, mock_get_instance):
        """Test reconciling cross-references in dry-run mode"""
        # Set up mocks
        mock_get_instance.return_value = self.mock_id_manager
        
        # Call the function in dry-run mode
        result = reconcile_chembl_properties.reconcile_cross_references(output_report=None, dry_run=True)
        
        # Verify the result
        self.assertEqual(result['molecules_updated'], 2)
        self.assertEqual(result['cross_references_added'], 2)
        self.assertIn('XSQUKJJJFZCRTK-UHFFFAOYSA-N', result['details'])
        
        # Verify that execute_query was not called (dry-run mode)
        mock_execute_query.assert_not_called()
        
        # Verify that the identifier manager was called to update and save identifiers
        self.assertEqual(self.mock_id_manager.add_molecule.call_count, 2)
        self.assertEqual(self.mock_id_manager.save_identifiers.call_count, 1)

    @patch('reconcile_chembl_properties.CryoprotectantIdentifierManager.get_instance')
    @patch('sql_executor.execute_query')
    def test_database_error_fallback(self, mock_execute_query, mock_get_instance):
        """Test fallback to dry-run mode when database connection fails"""
        # Set up mocks
        mock_get_instance.return_value = self.mock_id_manager
        
        # Make execute_query raise an exception
        mock_execute_query.side_effect = Exception("Database connection error")
        
        # Call the function
        result = reconcile_chembl_properties.reconcile_cross_references(output_report=None, dry_run=False)
        
        # Verify the result (should be the same as dry-run mode)
        self.assertEqual(result['molecules_updated'], 2)
        self.assertEqual(result['cross_references_added'], 2)
        self.assertIn('XSQUKJJJFZCRTK-UHFFFAOYSA-N', result['details'])
        
        # Verify that execute_query was called once (for the initial query that failed)
        self.assertEqual(mock_execute_query.call_count, 1)
        
        # Verify that the identifier manager was called to update and save identifiers
        self.assertEqual(self.mock_id_manager.add_molecule.call_count, 2)
        self.assertEqual(self.mock_id_manager.save_identifiers.call_count, 1)

if __name__ == '__main__':
    unittest.main()
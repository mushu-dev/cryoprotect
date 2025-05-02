#!/usr/bin/env python3
"""
Unit tests for the database verification script.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import json
from datetime import datetime

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from verify_imported_data import (
    count_molecules,
    verify_reference_compounds,
    verify_property_completeness,
    verify_query_performance,
    perform_full_verification,
    update_project_state,
    generate_markdown_report
)
from sql_executor import with_retry

class TestVerifyImportedData(unittest.TestCase):
    """Test cases for the database verification script."""

    @patch('verify_imported_data.execute_query')
    def test_count_molecules(self, mock_execute_query):
        """Test counting molecules in the database."""
        # Mock query results
        mock_execute_query.side_effect = [
            [{'count': 5000}],  # Total molecules
            [{'count': 4000}],  # PubChem molecules
            [{'count': 3000}],  # ChEMBL molecules
            [{'count': 2000}],  # Cross-referenced molecules
            [{'count': 4500}]   # Molecules with properties
        ]
        
        result = count_molecules()
        
        # Verify the result
        self.assertEqual(result['total_molecules'], 5000)
        self.assertEqual(result['with_pubchem_cid'], 4000)
        self.assertEqual(result['with_chembl_id'], 3000)
        self.assertEqual(result['with_cross_references'], 2000)
        self.assertEqual(result['with_properties'], 4500)
        
        # Verify execute_query was called 5 times
        self.assertEqual(mock_execute_query.call_count, 5)

    @patch('verify_imported_data.execute_query')
    def test_verify_reference_compounds_all_complete(self, mock_execute_query):
        """Test verifying reference compounds when all are complete."""
        # Mock molecule query results
        molecule_results = []
        property_results = []
        
        # Create mock data for 9 reference compounds
        for i in range(9):
            molecule_id = i + 1
            molecule_results.append([{
                'id': molecule_id,
                'name': f'Compound {i}',
                'smiles': 'CCCC',
                'inchi': 'InChI=1S/...',
                'inchikey': 'ABCDEFGHIJKLMN',
                'formula': 'C2H6O',
                'molecular_weight': 46.07,
                'pubchem_cid': str(1000 + i)
            }])
            
            # Properties for each molecule
            property_results.append([
                {'name': 'logP', 'numeric_value': -0.5, 'text_value': None},
                {'name': 'h_bond_donors', 'numeric_value': 1, 'text_value': None},
                {'name': 'h_bond_acceptors', 'numeric_value': 1, 'text_value': None}
            ])
        
        # Set up mock to return appropriate results for each call
        # Flatten the lists for side_effect
        all_results = []
        for i in range(9):
            all_results.append(molecule_results[i])
            all_results.append(property_results[i])
        
        mock_execute_query.side_effect = all_results
        
        result = verify_reference_compounds()
        
        # Verify the result
        self.assertEqual(result['total_reference_compounds'], 9)
        self.assertEqual(result['found_reference_compounds'], 9)
        self.assertEqual(len(result['complete_reference_compounds']), 9)
        self.assertEqual(len(result['incomplete_reference_compounds']), 0)
        self.assertEqual(len(result['missing_reference_compounds']), 0)
        
        # Verify execute_query was called for each reference compound (2 calls per compound)
        self.assertEqual(mock_execute_query.call_count, 18)

    @patch('verify_imported_data.execute_query')
    def test_verify_property_completeness(self, mock_execute_query):
        """Test verifying property completeness."""
        # Mock query results
        mock_execute_query.side_effect = [
            [  # Property counts
                {'name': 'logP', 'molecule_count': 4500},
                {'name': 'h_bond_donors', 'molecule_count': 4600},
                {'name': 'h_bond_acceptors', 'molecule_count': 4400}
            ],
            [{'count': 5000}],  # Total molecules with properties
            [{'total': 5000, 'complete': 4500}]  # Complete properties
        ]
        
        result = verify_property_completeness()
        
        # Verify the result
        self.assertEqual(result['total_molecules_with_properties'], 5000)
        self.assertEqual(result['molecules_with_complete_properties'], 4500)
        self.assertEqual(result['molecules_with_incomplete_properties'], 500)
        self.assertEqual(result['property_completeness_percentage'], 90.0)
        self.assertEqual(result['property_counts']['logP'], 4500)
        self.assertEqual(result['property_counts']['h_bond_donors'], 4600)
        self.assertEqual(result['property_counts']['h_bond_acceptors'], 4400)
        
        # Verify execute_query was called 3 times
        self.assertEqual(mock_execute_query.call_count, 3)

    @patch('verify_imported_data.time.time')
    @patch('verify_imported_data.execute_query')
    def test_verify_query_performance(self, mock_execute_query, mock_time):
        """Test verifying query performance."""
        # Mock time.time() to return increasing values
        # We need 40 values (4 queries * 5 runs * 2 calls per run)
        time_values = []
        for i in range(40):
            base_time = 10.0 + (i // 2) * 0.1  # Increment by 0.1 for each pair
            if i % 2 == 0:  # Start time
                time_values.append(base_time)
            else:  # End time, add a small delta
                time_values.append(base_time + 0.02 + (i % 5) * 0.01)  # 20-50ms range
                
        mock_time.side_effect = time_values
        
        # Mock query results (not important for this test)
        mock_execute_query.return_value = [{'id': 1}]
        
        result = verify_query_performance()
        
        # Verify the result
        self.assertEqual(len(result['queries']), 4)
        
        # Overall average should be in the 20-50ms range
        self.assertTrue(20 <= result['overall_average_ms'] <= 50)
        
        # Performance should be acceptable (< 50ms)
        self.assertTrue(result['performance_acceptable'])
        
        # Verify execute_query was called 20 times (4 queries * 5 runs)
        self.assertEqual(mock_execute_query.call_count, 20)

    @patch('verify_imported_data.close_connections')
    @patch('verify_imported_data.verify_query_performance')
    @patch('verify_imported_data.verify_property_completeness')
    @patch('verify_imported_data.verify_reference_compounds')
    @patch('verify_imported_data.count_molecules')
    @patch('verify_imported_data.execute_query')
    @patch('verify_imported_data.get_db')
    def test_perform_full_verification_success(self, mock_get_db, mock_execute_query, 
                                              mock_count, mock_verify_ref, 
                                              mock_verify_prop, mock_verify_perf,
                                              mock_close):
        """Test full verification with successful results."""
        # Mock test query result
        mock_execute_query.return_value = [{'test': 1}]
        
        # Mock verification function results
        mock_count.return_value = {
            'total_molecules': 5500,
            'with_pubchem_cid': 5000,
            'with_chembl_id': 4500,
            'with_cross_references': 4000,
            'with_properties': 5200
        }
        
        mock_verify_ref.return_value = {
            'total_reference_compounds': 9,
            'found_reference_compounds': 9,
            'missing_reference_compounds': [],
            'incomplete_reference_compounds': [],
            'complete_reference_compounds': ['CHEMBL1', 'CHEMBL2', 'CHEMBL3', 'CHEMBL4', 
                                           'CHEMBL5', 'CHEMBL6', 'CHEMBL7', 'CHEMBL8', 'CHEMBL9'],
            'details': {}
        }
        
        mock_verify_prop.return_value = {
            'property_counts': {},
            'molecules_with_complete_properties': 4700,
            'molecules_with_incomplete_properties': 300,
            'total_molecules_with_properties': 5000,
            'property_completeness_percentage': 94.0
        }
        
        mock_verify_perf.return_value = {
            'queries': [],
            'overall_average_ms': 30.0,
            'performance_acceptable': True
        }
        
        # Create a temporary file for the report
        temp_report = 'test_report.json'
        
        try:
            # Mock open to avoid actually writing to a file
            with patch('builtins.open', unittest.mock.mock_open()):
                with patch('os.makedirs'):
                    result = perform_full_verification(temp_report)
            
            # Verify the result
            self.assertTrue(result['success'])
            self.assertEqual(result['molecule_counts']['total_molecules'], 5500)
            self.assertEqual(result['reference_compounds']['total_reference_compounds'], 9)
            self.assertEqual(result['property_completeness']['property_completeness_percentage'], 94.0)
            self.assertEqual(result['query_performance']['overall_average_ms'], 30.0)
            
            # Verify summary
            self.assertEqual(result['summary']['total_molecules'], 5500)
            self.assertEqual(result['summary']['reference_compounds_complete'], 9)
            self.assertEqual(result['summary']['reference_compounds_total'], 9)
            self.assertEqual(result['summary']['property_completeness_percentage'], 94.0)
            self.assertEqual(result['summary']['average_query_time_ms'], 30.0)
            
            # Verify functions were called
            mock_get_db.assert_called_once()
            mock_execute_query.assert_called_once()
            mock_count.assert_called_once()
            mock_verify_ref.assert_called_once()
            mock_verify_prop.assert_called_once()
            mock_verify_perf.assert_called_once()
            mock_close.assert_called_once()
            
        finally:
            # Clean up temporary file if it exists
            if os.path.exists(temp_report):
                os.remove(temp_report)

    @patch('verify_imported_data.close_connections')
    @patch('verify_imported_data.execute_query')
    @patch('verify_imported_data.get_db')
    def test_perform_full_verification_connection_error(self, mock_get_db, mock_execute_query, mock_close):
        """Test full verification with database connection error."""
        # Mock execute_query to raise an exception
        mock_execute_query.side_effect = Exception("Connection refused")
        
        result = perform_full_verification()
        
        # Verify the result
        self.assertFalse(result['success'])
        self.assertIn('error', result)
        self.assertIn('Connection refused', result['error'])
        
        # Verify functions were called
        mock_get_db.assert_called_once()
        mock_execute_query.assert_called_once()
        mock_close.assert_called_once()

    @patch('verify_imported_data.json.load')
    @patch('verify_imported_data.json.dump')
    @patch('verify_imported_data.open', new_callable=unittest.mock.mock_open)
    @patch('verify_imported_data.os.path.exists')
    def test_update_project_state(self, mock_exists, mock_open, mock_dump, mock_load):
        """Test updating project state with verification results."""
        # Mock os.path.exists to return True
        mock_exists.return_value = True
        
        # Mock project state
        mock_project_state = {
            'highLevelPlan': [
                {'phase': 'Other Phase', 'status': 'Done'},
                {'phase': 'Verification Script Enhancement', 'status': 'Pending'}
            ]
        }
        mock_load.return_value = mock_project_state
        
        # Create verification results
        verification_results = {
            'timestamp': datetime.now().isoformat(),
            'success': True,
            'summary': {
                'total_molecules': 5500,
                'reference_compounds_complete': 9,
                'reference_compounds_total': 9,
                'property_completeness_percentage': 94.0,
                'average_query_time_ms': 30.0
            }
        }
        
        # Call the function
        update_project_state(verification_results)
        
        # Verify the project state was updated
        mock_exists.assert_called_once_with('project_state.json')
        mock_open.assert_any_call('project_state.json', 'r')
        mock_open.assert_any_call('project_state.json', 'w')
        mock_load.assert_called_once()
        mock_dump.assert_called_once()
        
        # Get the updated project state that was passed to json.dump
        updated_state = mock_dump.call_args[0][0]
        
        # Verify the verification phase was updated
        verification_phase = updated_state['highLevelPlan'][1]
        self.assertEqual(verification_phase['phase'], 'Verification Script Enhancement')
        self.assertEqual(verification_phase['status'], 'Done')
        self.assertTrue('verification_results' in verification_phase)
        self.assertTrue(verification_phase['verification_results']['success'])

if __name__ == '__main__':
    unittest.main()
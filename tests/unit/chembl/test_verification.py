#!/usr/bin/env python3
"""
Unit tests for verify_chembl_data.py

These tests mock the SupabaseDirectConnection to test the verification functions
with controlled data for different scenarios.
"""

import unittest
from unittest.mock import patch, MagicMock
import json
import sys
import os
from typing import Dict, List, Any

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

# Import the module to test
import verify_chembl_data

# Constants for test statuses
SUCCESS = verify_chembl_data.STATUS_SUCCESS
WARNING = verify_chembl_data.STATUS_WARNING
FAILURE = verify_chembl_data.STATUS_FAILURE

class TestVerifyChemblData(unittest.TestCase):
    """Test cases for the ChEMBL data verification script."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a mock for SupabaseDirectConnection
        self.mock_db = MagicMock()
        
        # Create a patcher for get_instance to return our mock
        self.patcher = patch('verify_chembl_data.SupabaseDirectConnection.get_instance', 
                            return_value=self.mock_db)
        self.mock_get_instance = self.patcher.start()
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.patcher.stop()
    
    def test_verify_molecule_counts_success(self):
        """Test verify_molecule_counts with count above threshold."""
        # Configure mock to return a count above threshold
        self.mock_db.execute_query.return_value = [{'count': 2000}]
        
        # Call the function
        status, data = verify_chembl_data.verify_molecule_counts()
        
        # Verify the results
        self.assertEqual(status, SUCCESS)
        self.assertEqual(data['count'], 2000)
        self.assertIn("2000", data['message'])
        self.assertIn("expected >= 1000", data['message'])
    
    def test_verify_molecule_counts_warning(self):
        """Test verify_molecule_counts with count below threshold but above zero."""
        # Configure mock to return a count below threshold
        self.mock_db.execute_query.return_value = [{'count': 500}]
        
        # Call the function
        status, data = verify_chembl_data.verify_molecule_counts()
        
        # Verify the results
        self.assertEqual(status, WARNING)
        self.assertEqual(data['count'], 500)
        self.assertIn("500", data['message'])
        self.assertIn("expected >= 1000", data['message'])
    
    def test_verify_molecule_counts_failure(self):
        """Test verify_molecule_counts with zero count."""
        # Configure mock to return zero count
        self.mock_db.execute_query.return_value = [{'count': 0}]
        
        # Call the function
        status, data = verify_chembl_data.verify_molecule_counts()
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertEqual(data['count'], 0)
        self.assertIn("No ChEMBL molecules found", data['message'])
    
    def test_verify_molecule_counts_error(self):
        """Test verify_molecule_counts with database error."""
        # Configure mock to raise an exception
        self.mock_db.execute_query.side_effect = Exception("Database connection error")
        
        # Call the function
        status, data = verify_chembl_data.verify_molecule_counts()
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertEqual(data['count'], 0)
        self.assertIn("Error", data['message'])
        self.assertIn("Database connection error", data['message'])
    
    def test_verify_property_distribution_success(self):
        """Test verify_property_distribution with normal distribution."""
        # Configure mock to return property data with normal distribution
        self.mock_db.execute_query.return_value = [
            {'property_name': 'LogP', 'value': '2.5'},
            {'property_name': 'LogP', 'value': '3.0'},
            {'property_name': 'LogP', 'value': '1.5'},
            {'property_name': 'Molecular Weight', 'value': '350.5'},
            {'property_name': 'Molecular Weight', 'value': '420.1'},
            {'property_name': 'Molecular Weight', 'value': '280.3'},
            {'property_name': 'Hydrogen Bond Donors', 'value': '2'},
            {'property_name': 'Hydrogen Bond Donors', 'value': '3'},
            {'property_name': 'Hydrogen Bond Acceptors', 'value': '5'},
            {'property_name': 'Hydrogen Bond Acceptors', 'value': '4'}
        ]
        
        # Call the function
        status, data = verify_chembl_data.verify_property_distribution()
        
        # Verify the results
        self.assertEqual(status, SUCCESS)
        self.assertIn("Property distributions appear normal", data['message'])
        self.assertIn('LogP', data['statistics'])
        self.assertIn('Molecular Weight', data['statistics'])
        self.assertEqual(len(data['anomalies']), 0)
    
    def test_verify_property_distribution_warning(self):
        """Test verify_property_distribution with anomalies."""
        # Configure mock to return property data with anomalies
        self.mock_db.execute_query.return_value = [
            {'property_name': 'LogP', 'value': '15.0'},  # Anomalous LogP
            {'property_name': 'LogP', 'value': '18.0'},  # Anomalous LogP
            {'property_name': 'Molecular Weight', 'value': '1500.5'},  # Anomalous MW
            {'property_name': 'Molecular Weight', 'value': '1600.1'},  # Anomalous MW
            {'property_name': 'Hydrogen Bond Donors', 'value': '2'},
            {'property_name': 'Hydrogen Bond Acceptors', 'value': '5'}
        ]
        
        # Call the function
        status, data = verify_chembl_data.verify_property_distribution()
        
        # Verify the results
        self.assertEqual(status, WARNING)
        self.assertIn("potential anomalies", data['message'])
        self.assertGreater(len(data['anomalies']), 0)
    
    def test_verify_property_distribution_failure(self):
        """Test verify_property_distribution with no data."""
        # Configure mock to return empty result
        self.mock_db.execute_query.return_value = []
        
        # Call the function
        status, data = verify_chembl_data.verify_property_distribution()
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertIn("No property data found", data['message'])
    
    def test_verify_property_distribution_error(self):
        """Test verify_property_distribution with database error."""
        # Configure mock to raise an exception
        self.mock_db.execute_query.side_effect = Exception("Database connection error")
        
        # Call the function
        status, data = verify_chembl_data.verify_property_distribution()
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertIn("Error", data['message'])
        self.assertIn("Database connection error", data['message'])
    
    def test_verify_sample_data_success(self):
        """Test verify_sample_data with all samples found."""
        # Configure mock to return sample data for all IDs
        def mock_execute_query(query):
            if "SELECT m.id, m.name, m.smiles, m.chembl_id" in query:
                if "CHEMBL25" in query:
                    return [{'id': '1', 'name': 'Atorvastatin', 'smiles': 'CC(C)C...', 'chembl_id': 'CHEMBL25'}]
                elif "CHEMBL1" in query:
                    return [{'id': '2', 'name': 'Aspirin', 'smiles': 'CC(=O)OC...', 'chembl_id': 'CHEMBL1'}]
                else:
                    return [{'id': '3', 'name': 'Ibuprofen', 'smiles': 'CC(C)CC...', 'chembl_id': 'CHEMBL2'}]
            else:
                return [
                    {'property_name': 'LogP', 'value': '2.5'},
                    {'property_name': 'Molecular Weight', 'value': '350.5'}
                ]
        
        self.mock_db.execute_query.side_effect = mock_execute_query
        
        # Call the function with specific sample IDs
        status, data = verify_chembl_data.verify_sample_data(['CHEMBL25', 'CHEMBL1', 'CHEMBL2'])
        
        # Verify the results
        self.assertEqual(status, SUCCESS)
        self.assertEqual(data['total_checked'], 3)
        self.assertEqual(data['total_found'], 3)
        self.assertEqual(len(data['missing_ids']), 0)
        self.assertIn("All 3 sample ChEMBL IDs were found", data['message'])
    
    def test_verify_sample_data_warning(self):
        """Test verify_sample_data with some samples missing."""
        # Configure mock to return sample data for some IDs but not others
        def mock_execute_query(query):
            if "SELECT m.id, m.name, m.smiles, m.chembl_id" in query:
                if "CHEMBL25" in query:
                    return [{'id': '1', 'name': 'Atorvastatin', 'smiles': 'CC(C)C...', 'chembl_id': 'CHEMBL25'}]
                else:
                    return []  # No results for other IDs
            else:
                return [
                    {'property_name': 'LogP', 'value': '2.5'},
                    {'property_name': 'Molecular Weight', 'value': '350.5'}
                ]
        
        self.mock_db.execute_query.side_effect = mock_execute_query
        
        # Call the function with specific sample IDs
        status, data = verify_chembl_data.verify_sample_data(['CHEMBL25', 'CHEMBL1', 'CHEMBL2'])
        
        # Verify the results
        self.assertEqual(status, WARNING)
        self.assertEqual(data['total_checked'], 3)
        self.assertEqual(data['total_found'], 1)
        self.assertEqual(len(data['missing_ids']), 2)
        self.assertIn("2 out of 3 sample ChEMBL IDs were not found", data['message'])
    
    def test_verify_sample_data_failure(self):
        """Test verify_sample_data with all samples missing."""
        # Configure mock to return no results for any ID
        self.mock_db.execute_query.return_value = []
        
        # Call the function with specific sample IDs
        status, data = verify_chembl_data.verify_sample_data(['CHEMBL25', 'CHEMBL1', 'CHEMBL2'])
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertEqual(data['total_checked'], 3)
        self.assertEqual(data['total_found'], 0)
        self.assertEqual(len(data['missing_ids']), 3)
        self.assertIn("None of the sample ChEMBL IDs were found", data['message'])
    
    def test_verify_sample_data_error(self):
        """Test verify_sample_data with database error."""
        # Configure mock to raise an exception
        self.mock_db.execute_query.side_effect = Exception("Database connection error")
        
        # Call the function
        status, data = verify_chembl_data.verify_sample_data()
        
        # Verify the results
        self.assertEqual(status, FAILURE)
        self.assertIn("Error", data['message'])
        self.assertIn("Database connection error", data['message'])
    
    def test_generate_verification_report(self):
        """Test generate_verification_report with mixed results."""
        # Mock the individual verification functions
        with patch('verify_chembl_data.verify_molecule_counts') as mock_count, \
             patch('verify_chembl_data.verify_property_distribution') as mock_dist, \
             patch('verify_chembl_data.verify_sample_data') as mock_sample:
            
            # Configure mocks to return different statuses
            mock_count.return_value = (SUCCESS, {'count': 2000, 'message': 'Found 2000 molecules'})
            mock_dist.return_value = (WARNING, {'message': 'Found anomalies', 'statistics': {}, 'anomalies': ['Anomaly 1']})
            mock_sample.return_value = (SUCCESS, {'message': 'All samples found', 'results': {}, 'missing_ids': [], 'total_checked': 3, 'total_found': 3})
            
            # Call the function
            report = verify_chembl_data.generate_verification_report()
            
            # Verify the results
            self.assertEqual(report['overall_status'], WARNING)
            self.assertEqual(report['checks']['molecule_count']['status'], SUCCESS)
            self.assertEqual(report['checks']['property_distribution']['status'], WARNING)
            self.assertEqual(report['checks']['sample_data']['status'], SUCCESS)

if __name__ == '__main__':
    unittest.main()
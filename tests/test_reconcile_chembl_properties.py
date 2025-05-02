#!/usr/bin/env python3
"""
Test script for reconcile_chembl_properties.py

This script tests the functionality of the ChEMBL properties reconciliation script.
It verifies that the script can correctly fetch ChEMBL data, compare properties,
and update the database as needed.

Usage:
    python -m tests.test_reconcile_chembl_properties
"""

import os
import sys
import unittest
import json
from unittest.mock import patch, MagicMock
from pathlib import Path

# Add parent directory to path to import the script
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the script to test
import reconcile_chembl_properties

class TestReconcileChEMBLProperties(unittest.TestCase):
    """Test cases for reconcile_chembl_properties.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a mock ChEMBL client
        self.mock_chembl_client = MagicMock()
        
        # Create a mock project ID
        self.project_id = "test_project_id"
        
        # Set up mock molecule data
        self.mock_molecule = {
            "id": 1,
            "name": "Test Molecule",
            "chembl_id": "CHEMBL1234"
        }
        
        # Set up mock ChEMBL data
        self.mock_chembl_data = {
            "ChEMBL ID": "CHEMBL1234",
            "Name": "Test Molecule",
            "LogP": 2.5,
            "Molecular Weight": 150.0,
            "H-Bond Acceptors": 3,
            "H-Bond Donors": 2,
            "TPSA": 60.0,
            "Rotatable Bonds": 4
        }
        
        # Set up mock database properties
        self.mock_db_properties = {
            "LogP": {
                "id": 1,
                "numeric_value": 2.4,
                "text_value": None,
                "data_source": "Original import"
            },
            "Molecular Weight": {
                "id": 2,
                "numeric_value": 150.0,
                "text_value": None,
                "data_source": "Original import"
            },
            "Hydrogen Bond Acceptor Count": {
                "id": 3,
                "numeric_value": 3,
                "text_value": None,
                "data_source": "Original import"
            },
            "Hydrogen Bond Donor Count": {
                "id": 4,
                "numeric_value": 1,  # Different from ChEMBL
                "text_value": None,
                "data_source": "Original import"
            }
            # Missing TPSA and Rotatable Bond Count
        }
        
        # Set up mock property type IDs
        self.mock_property_type_ids = {
            "LogP": 1,
            "Molecular Weight": 2,
            "Hydrogen Bond Acceptor Count": 3,
            "Hydrogen Bond Donor Count": 4,
            "Topological Polar Surface Area": 5,
            "Rotatable Bond Count": 6
        }

    @patch('reconcile_chembl_properties.get_molecular_properties')
    @patch('reconcile_chembl_properties.get_property_type_id')
    @patch('reconcile_chembl_properties.update_molecular_property')
    @patch('reconcile_chembl_properties.insert_molecular_property')
    def test_reconcile_molecule_properties(self, mock_insert, mock_update, mock_get_property_type_id, mock_get_properties):
        """Test reconciling properties for a single molecule"""
        # Set up mocks
        mock_get_properties.return_value = self.mock_db_properties
        mock_get_property_type_id.side_effect = lambda name, _: self.mock_property_type_ids.get(name)
        mock_update.return_value = True
        mock_insert.return_value = 100
        
        self.mock_chembl_client.get_molecule_properties.return_value = self.mock_chembl_data
        
        # Call the function
        result = reconcile_chembl_properties.reconcile_molecule_properties(
            self.mock_molecule,
            self.mock_chembl_client,
            0.01,
            self.project_id
        )
        
        # Verify the result
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["molecule_id"], 1)
        self.assertEqual(result["chembl_id"], "CHEMBL1234")
        self.assertEqual(result["name"], "Test Molecule")
        
        # Should have 4 updates: LogP (update), HBD (update), TPSA (insert), RTB (insert)
        self.assertEqual(len(result["updates"]), 4)
        
        # Verify that the ChEMBL client was called
        self.mock_chembl_client.get_molecule_properties.assert_called_once_with("CHEMBL1234")
        
        # Verify that get_molecular_properties was called
        mock_get_properties.assert_called_once_with(1, self.project_id)
        
        # Verify that update_molecular_property was called for LogP and HBD
        self.assertEqual(mock_update.call_count, 2)
        
        # Verify that insert_molecular_property was called for TPSA and Rotatable Bond Count
        self.assertEqual(mock_insert.call_count, 2)

    @patch('reconcile_chembl_properties.get_molecules_with_chembl_id')
    @patch('reconcile_chembl_properties.reconcile_molecule_properties')
    @patch('reconcile_chembl_properties.generate_reconciliation_report')
    @patch('reconcile_chembl_properties.save_reconciliation_report')
    @patch('reconcile_chembl_properties.get_project_id')
    @patch('reconcile_chembl_properties.verify_database_role')
    @patch('reconcile_chembl_properties.initialize_chembl_client')
    @patch('argparse.ArgumentParser.parse_args')
    def test_main_function(self, mock_parse_args, mock_init_client, mock_verify_role, 
                          mock_get_project_id, mock_save_report, mock_generate_report, 
                          mock_reconcile, mock_get_molecules):
        """Test the main function"""
        # Set up mocks
        mock_args = MagicMock()
        mock_args.project_id = None
        mock_args.tolerance = 0.01
        mock_parse_args.return_value = mock_args
        
        mock_get_project_id.return_value = self.project_id
        mock_verify_role.return_value = ("service_role", "postgres")
        mock_init_client.return_value = self.mock_chembl_client
        
        mock_get_molecules.return_value = [self.mock_molecule]
        
        mock_reconcile.return_value = {
            "molecule_id": 1,
            "chembl_id": "CHEMBL1234",
            "name": "Test Molecule",
            "status": "success",
            "updates": [{"property": "LogP", "old_value": 2.4, "new_value": 2.5, "action": "updated"}]
        }
        
        mock_generate_report.return_value = {
            "timestamp": "20250427_123456",
            "summary": {
                "total_molecules": 1,
                "successful_molecules": 1,
                "error_molecules": 0,
                "total_updates": 1,
                "updates_by_property": {"LogP": 1}
            },
            "results": [mock_reconcile.return_value]
        }
        
        mock_save_report.return_value = "reports/chembl_reconciliation_report_20250427_123456.json"
        
        # Call the main function
        with patch.object(sys, 'argv', ['reconcile_chembl_properties.py']):
            reconcile_chembl_properties.main()
        
        # Verify that the functions were called
        mock_get_project_id.assert_called_once()
        mock_verify_role.assert_called_once_with(self.project_id)
        mock_init_client.assert_called_once()
        mock_get_molecules.assert_called_once_with(self.project_id)
        mock_reconcile.assert_called_once()
        mock_generate_report.assert_called_once()
        mock_save_report.assert_called_once()

    def test_generate_reconciliation_report(self):
        """Test generating a reconciliation report"""
        # Set up test data
        results = [
            {
                "molecule_id": 1,
                "chembl_id": "CHEMBL1234",
                "name": "Test Molecule 1",
                "status": "success",
                "updates": [
                    {"property": "LogP", "old_value": 2.4, "new_value": 2.5, "action": "updated"},
                    {"property": "TPSA", "old_value": None, "new_value": 60.0, "action": "inserted"}
                ]
            },
            {
                "molecule_id": 2,
                "chembl_id": "CHEMBL5678",
                "name": "Test Molecule 2",
                "status": "error",
                "error": "API error",
                "updates": []
            },
            {
                "molecule_id": 3,
                "chembl_id": "CHEMBL9012",
                "name": "Test Molecule 3",
                "status": "success",
                "updates": [
                    {"property": "LogP", "old_value": 1.5, "new_value": 1.6, "action": "updated"}
                ]
            }
        ]
        
        # Call the function
        report = reconcile_chembl_properties.generate_reconciliation_report(results)
        
        # Verify the report
        self.assertIn("timestamp", report)
        self.assertIn("summary", report)
        self.assertIn("results", report)
        
        summary = report["summary"]
        self.assertEqual(summary["total_molecules"], 3)
        self.assertEqual(summary["successful_molecules"], 2)
        self.assertEqual(summary["error_molecules"], 1)
        self.assertEqual(summary["total_updates"], 3)
        self.assertEqual(summary["updates_by_property"]["LogP"], 2)
        self.assertEqual(summary["updates_by_property"]["TPSA"], 1)

if __name__ == '__main__':
    unittest.main()
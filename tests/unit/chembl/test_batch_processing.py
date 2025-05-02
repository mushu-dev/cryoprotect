#!/usr/bin/env python3
"""
Unit tests for batch processing in ChEMBL_Integrated_Import.py
"""

import unittest
from unittest.mock import patch, MagicMock, call
import json
from datetime import datetime
import uuid

# Import the module to test
import ChEMBL_Integrated_Import as chembl_import
from supabase_direct import SupabaseDirectConnection


class TestChEMBLBatchProcessing(unittest.TestCase):
    """Test cases for batch processing in ChEMBL_Integrated_Import.py"""

    def setUp(self):
        """Set up test fixtures"""
        # Create a mock for SupabaseDirectConnection
        self.mock_db = MagicMock(spec=SupabaseDirectConnection)
        
        # Mock the get_db_connection function to return our mock
        patcher = patch.object(chembl_import, 'get_db_connection', return_value=self.mock_db)
        self.addCleanup(patcher.stop)
        self.mock_get_db = patcher.start()
        
        # Mock user profile ID
        self.user_profile_id = 'test-user-profile-id'
        
        # Mock property type map
        self.property_type_map = {
            'logp': 'prop-type-logp',
            'molecular weight': 'prop-type-mw',
            'hydrogen bond acceptor count': 'prop-type-hba',
            'hydrogen bond donor count': 'prop-type-hbd',
            'topological polar surface area': 'prop-type-psa',
            'rotatable bond count': 'prop-type-rtb',
            'unknown property type': 'prop-type-unknown'
        }
        
        # Create mock compounds
        self.mock_compounds = []
        for i in range(5):
            self.mock_compounds.append({
                'molecule_chembl_id': f'CHEMBL{1000+i}',
                'pref_name': f'Test Compound {i}',
                'molecule_structures': {
                    'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                    'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                    'standard_inchi_key': f'INCHIKEY{1000+i}'
                },
                'molecule_properties': {
                    'full_molformula': 'C9H8O4',
                    'full_mwt': 180.16,
                    'alogp': 1.31,
                    'hba': 4,
                    'hbd': 1,
                    'psa': 63.6,
                    'rtb': 3
                }
            })
        
        # Mock functions
        patcher = patch.object(chembl_import, 'get_user_id', return_value='test-user-id')
        self.addCleanup(patcher.stop)
        self.mock_get_user_id = patcher.start()
        
        patcher = patch.object(chembl_import, 'ensure_user_profile', return_value=self.user_profile_id)
        self.addCleanup(patcher.stop)
        self.mock_ensure_user_profile = patcher.start()
        
        patcher = patch.object(chembl_import, 'get_property_types', return_value=self.property_type_map)
        self.addCleanup(patcher.stop)
        self.mock_get_property_types = patcher.start()
        
        # Mock UUID generation to make tests deterministic
        patcher = patch.object(uuid, 'uuid4', side_effect=lambda: 'test-uuid')
        self.addCleanup(patcher.stop)
        self.mock_uuid4 = patcher.start()
        
        # Mock logging functions
        patcher = patch.object(chembl_import.logger, 'info')
        self.addCleanup(patcher.stop)
        self.mock_logger_info = patcher.start()
        
        patcher = patch.object(chembl_import.logger, 'debug')
        self.addCleanup(patcher.stop)
        self.mock_logger_debug = patcher.start()
        
        patcher = patch.object(chembl_import.logger, 'warning')
        self.addCleanup(patcher.stop)
        self.mock_logger_warning = patcher.start()
        
        patcher = patch.object(chembl_import.logger, 'error')
        self.addCleanup(patcher.stop)
        self.mock_logger_error = patcher.start()
        
        patcher = patch.object(chembl_import, 'log_error')
        self.addCleanup(patcher.stop)
        self.mock_log_error = patcher.start()
        
        patcher = patch.object(chembl_import, 'log_skipped_molecule')
        self.addCleanup(patcher.stop)
        self.mock_log_skipped_molecule = patcher.start()

    def test_batch_processing_success(self):
        """Test successful batch processing of compounds"""
        # Set up mocks for successful batch processing
        
        # Mock execute_query to return existing molecules for the first compound
        self.mock_db.execute_query.side_effect = [
            # First call: check existing molecules
            [{'id': 'existing-molecule-id', 'inchikey': 'INCHIKEY1000'}],
            # Subsequent calls: successful inserts
            [{'id': 'new-molecule-id-1'}],
            [{'id': 'new-property-id-1'}],
            [{'id': 'new-property-id-2'}]
        ]
        
        # Mock execute_batch to return True (success)
        self.mock_db.execute_batch.return_value = True
        
        # Call the function with a batch size of 3
        result = chembl_import.import_compounds_to_database(self.mock_compounds, batch_size=3)
        
        # Verify the result
        self.assertEqual(result['total_compounds'], 5)
        self.assertEqual(result['processed'], 5)
        self.assertEqual(result['batches_processed'], 2)  # 5 compounds with batch size 3 = 2 batches
        self.assertEqual(result['batches_failed'], 0)
        
        # Verify execute_batch was called for molecule and property inserts
        self.assertTrue(self.mock_db.execute_batch.called)
        
        # Verify the batch processing logic
        self.assertEqual(self.mock_db.execute_batch.call_count, 2)  # Once for molecules, once for properties
        
        # Verify logging
        self.mock_logger_info.assert_any_call("Successfully processed batch 1/2")
        self.mock_logger_info.assert_any_call("Successfully processed batch 2/2")

    def test_batch_processing_with_errors(self):
        """Test batch processing with errors"""
        # Set up mocks for batch processing with errors
        
        # Mock execute_query to return empty results (no existing molecules)
        self.mock_db.execute_query.return_value = []
        
        # Mock execute_batch to fail on the first batch
        self.mock_db.execute_batch.side_effect = [False, True]
        
        # Call the function with a batch size of 3
        result = chembl_import.import_compounds_to_database(self.mock_compounds, batch_size=3)
        
        # Verify the result
        self.assertEqual(result['total_compounds'], 5)
        self.assertEqual(result['batches_failed'], 1)
        
        # Verify error logging
        self.mock_log_error.assert_called()
        
        # Verify the second batch was still processed
        self.mock_logger_info.assert_any_call("Successfully processed batch 2/2")

    def test_batch_processing_with_missing_inchikey(self):
        """Test batch processing with compounds missing InChIKey"""
        # Create a compound with missing InChIKey
        compound_missing_inchikey = {
            'molecule_chembl_id': 'CHEMBL9999',
            'pref_name': 'Missing InChIKey Compound',
            'molecule_structures': {
                'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                # Missing InChIKey
            },
            'molecule_properties': {
                'full_molformula': 'C9H8O4',
                'full_mwt': 180.16
            }
        }
        
        # Add the compound to our test compounds
        test_compounds = [compound_missing_inchikey] + self.mock_compounds
        
        # Mock execute_query to return empty results
        self.mock_db.execute_query.return_value = []
        
        # Mock execute_batch to succeed
        self.mock_db.execute_batch.return_value = True
        
        # Call the function
        result = chembl_import.import_compounds_to_database(test_compounds, batch_size=3)
        
        # Verify the result
        self.assertEqual(result['molecules_skipped'], 1)
        
        # Verify log_skipped_molecule was called
        self.mock_log_skipped_molecule.assert_called_once()

    def test_dry_run_mode(self):
        """Test dry run mode"""
        # Call the function in dry run mode
        result = chembl_import.import_compounds_to_database(self.mock_compounds, batch_size=3, dry_run=True)
        
        # Verify the result
        self.assertEqual(result['total_compounds'], 5)
        
        # Verify execute_batch was not called
        self.mock_db.execute_batch.assert_not_called()
        
        # Verify logging
        self.mock_logger_info.assert_any_call("DRY RUN MODE: No data will be inserted into the database")


if __name__ == '__main__':
    unittest.main()
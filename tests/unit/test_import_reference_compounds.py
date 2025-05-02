#!/usr/bin/env python3
"""
Unit tests for the reference compounds import script with direct PostgreSQL connections.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock, mock_open
import json
from datetime import datetime
import uuid

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the module to test
import import_reference_compounds
from postgres_direct import PostgresDirectConnection
from sql_executor import execute_query, get_db
from property_utils import PropertyManager

class TestImportReferenceCompounds(unittest.TestCase):
    """Test cases for the reference compounds import script."""

    def setUp(self):
        """Set up test fixtures."""
        # Mock ChEMBL client
        self.chembl_client_patcher = patch('import_reference_compounds.ChEMBLClient')
        self.mock_chembl_client = self.chembl_client_patcher.start()
        self.mock_chembl_instance = self.mock_chembl_client.return_value
        
        # Mock PubChem client
        self.pubchem_client_patcher = patch('import_reference_compounds.PubChemClient')
        self.mock_pubchem_client = self.pubchem_client_patcher.start()
        self.mock_pubchem_instance = self.mock_pubchem_client.return_value
        
        # Mock identifier manager
        self.id_manager_patcher = patch('import_reference_compounds.CryoprotectantIdentifierManager')
        self.mock_id_manager = self.id_manager_patcher.start()
        self.mock_id_manager_instance = MagicMock()
        self.mock_id_manager.get_instance.return_value = self.mock_id_manager_instance
        
        # Mock PropertyManager
        self.property_manager_patcher = patch('import_reference_compounds.PropertyManager')
        self.mock_property_manager = self.property_manager_patcher.start()
        self.mock_property_manager_instance = self.mock_property_manager.return_value
        
        # Mock reference compound IDs
        self.get_reference_ids_patcher = patch('import_reference_compounds.get_reference_compound_ids')
        self.mock_get_reference_ids = self.get_reference_ids_patcher.start()
        self.mock_get_reference_ids.return_value = ['CHEMBL123', 'CHEMBL456']
        
        # Mock get_db function
        self.get_db_patcher = patch('sql_executor.get_db')
        self.mock_get_db = self.get_db_patcher.start()
        self.mock_db = MagicMock()
        self.mock_get_db.return_value = self.mock_db
        
        # Mock SQL executor functions
        self.execute_query_patcher = patch('import_reference_compounds.execute_query')
        self.mock_execute_query = self.execute_query_patcher.start()
        
        # Mock process_in_batches to execute the function directly
        self.process_batches_patcher = patch('import_reference_compounds.process_in_batches')
        self.mock_process_batches = self.process_batches_patcher.start()
        self.mock_process_batches.side_effect = lambda items, batch_size, process_func: process_func(items)
        
        # Mock checkpoint file operations
        self.mock_open_patcher = patch('builtins.open', mock_open())
        self.mock_file = self.mock_open_patcher.start()
        
        # Mock os.path.exists and os.makedirs
        self.os_path_exists_patcher = patch('os.path.exists')
        self.mock_os_path_exists = self.os_path_exists_patcher.start()
        self.mock_os_path_exists.return_value = False
        
        self.os_makedirs_patcher = patch('os.makedirs')
        self.mock_os_makedirs = self.os_makedirs_patcher.start()

    def tearDown(self):
        """Tear down test fixtures."""
        self.chembl_client_patcher.stop()
        self.pubchem_client_patcher.stop()
        self.id_manager_patcher.stop()
        self.property_manager_patcher.stop()
        self.get_reference_ids_patcher.stop()
        self.execute_query_patcher.stop()
        self.get_db_patcher.stop()
        self.process_batches_patcher.stop()
        self.mock_open_patcher.stop()
        self.os_path_exists_patcher.stop()
        self.os_makedirs_patcher.stop()

    def _setup_mock_chembl_data(self, chembl_id):
        """Helper to set up mock ChEMBL data."""
        return {
            'molecule_chembl_id': chembl_id,
            'pref_name': f'Test Compound {chembl_id}',
            'molecule_properties': {
                'full_molformula': 'C10H20O2',
                'full_mwt': 172.26,
                'alogp': 2.5,
                'hbd': 1,
                'hba': 2,
                'rtb': 3,
                'psa': 40.5,
                'heavy_atoms': 12
            },
            'molecule_structures': {
                'standard_inchi': 'InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)',
                'standard_inchi_key': 'KXUTQHGMRWAGAY-UHFFFAOYSA-N',
                'canonical_smiles': 'CCCCCCCCCC(=O)O'
            },
            'cross_references': {
                'pubchem_cid': '12345'
            }
        }

    def _setup_mock_pubchem_data(self, cid):
        """Helper to set up mock PubChem data."""
        return {
            'cid': cid,
            'props': [
                {
                    'urn': {'label': 'LogP'},
                    'value': {'sval': '2.8'}
                },
                {
                    'urn': {'label': 'Water Solubility'},
                    'value': {'sval': '0.5 g/L'}
                },
                {
                    'urn': {'label': 'Melting Point'},
                    'value': {'sval': '32 °C'}
                }
            ]
        }

    def test_import_reference_compounds_success(self):
        """Test successful import of reference compounds."""
        # Set up mock ChEMBL data
        self.mock_chembl_instance.get_molecule_by_chembl_id.side_effect = lambda chembl_id: self._setup_mock_chembl_data(chembl_id)
        
        # Set up mock PubChem data
        self.mock_pubchem_instance.get_molecule_properties.return_value = self._setup_mock_pubchem_data('12345')
        
        # Set up mock identifier manager
        self.mock_id_manager_instance.resolve_identifier.return_value = ('CRYO0001', False)
        
        # Set up mock property manager
        self.mock_property_manager_instance.set_properties.return_value = (5, 5)
        
        # Run the import function
        results = import_reference_compounds.import_reference_compounds(
            output_report='test_report.json',
            checkpoint_path='test_checkpoint.json',
            resume=False,
            batch_size=10
        )
        
        # Verify results
        self.assertEqual(results['total_compounds'], 2)
        self.assertEqual(results['updated'], 2)
        self.assertEqual(results['failed'], 0)
        
        # Verify ChEMBL client was called correctly
        self.mock_chembl_instance.get_molecule_by_chembl_id.assert_any_call('CHEMBL123')
        self.mock_chembl_instance.get_molecule_by_chembl_id.assert_any_call('CHEMBL456')
        
        # Verify PubChem client was called correctly
        self.mock_pubchem_instance.get_molecule_properties.assert_called_with('12345')
        
        # Verify identifier manager was called correctly
        self.mock_id_manager_instance.resolve_identifier.assert_any_call(chembl_id='CHEMBL123')
        self.mock_id_manager_instance.resolve_identifier.assert_any_call(chembl_id='CHEMBL456')
        
        # Verify SQL execution was called
        self.mock_execute_query.assert_called()
        
        # Verify property manager was called correctly
        self.mock_property_manager_instance.set_properties.assert_called()
        
        # Verify checkpoint was created
        self.mock_file.assert_called()

    def test_import_reference_compounds_with_checkpoint(self):
        """Test resuming import from checkpoint."""
        # Set up mock checkpoint file
        self.mock_os_path_exists.return_value = True
        mock_checkpoint_data = {
            'timestamp': datetime.now().isoformat(),
            'processed_ids': ['CHEMBL123'],
            'results': {
                'timestamp': datetime.now().isoformat(),
                'total_compounds': 2,
                'imported': 1,
                'updated': 0,
                'failed': 0,
                'details': {
                    'CHEMBL123': {
                        'status': 'success',
                        'internal_id': 'CRYO0001',
                        'is_new': True,
                        'has_pubchem': True
                    }
                }
            }
        }
        
        # Mock json.load to return checkpoint data
        with patch('json.load', return_value=mock_checkpoint_data):
            # Set up mock ChEMBL data for remaining compound
            self.mock_chembl_instance.get_molecule_by_chembl_id.side_effect = lambda chembl_id: self._setup_mock_chembl_data(chembl_id)
            
            # Set up mock PubChem data
            self.mock_pubchem_instance.get_molecule_properties.return_value = self._setup_mock_pubchem_data('12345')
            
            # Set up mock identifier manager
            self.mock_id_manager_instance.resolve_identifier.return_value = ('CRYO0002', True)
            
            # Set up mock property manager
            self.mock_property_manager_instance.set_properties.return_value = (5, 5)
            
            # Run the import function with resume=True
            results = import_reference_compounds.import_reference_compounds(
                output_report='test_report.json',
                checkpoint_path='test_checkpoint.json',
                resume=True,
                batch_size=10
            )
            
            # Verify results
            self.assertEqual(results['total_compounds'], 2)
            self.assertEqual(results['imported'], 2)  # 1 from checkpoint + 1 new
            self.assertEqual(results['failed'], 0)
            
            # Verify only the remaining compound was processed
            self.mock_chembl_instance.get_molecule_by_chembl_id.assert_called_once_with('CHEMBL456')

    def test_import_reference_compounds_error_handling(self):
        """Test error handling during import."""
        # Set up first compound to succeed
        self.mock_chembl_instance.get_molecule_by_chembl_id.side_effect = [
            self._setup_mock_chembl_data('CHEMBL123'),  # First call succeeds
            None  # Second call returns None (simulating error)
        ]
        
        # Set up mock PubChem data
        self.mock_pubchem_instance.get_molecule_properties.return_value = self._setup_mock_pubchem_data('12345')
        
        # Set up mock identifier manager
        self.mock_id_manager_instance.resolve_identifier.return_value = ('CRYO0001', True)
        
        # Set up mock property manager
        self.mock_property_manager_instance.set_properties.return_value = (5, 5)
        
        # Run the import function
        results = import_reference_compounds.import_reference_compounds(
            output_report='test_report.json',
            checkpoint_path='test_checkpoint.json',
            resume=False,
            batch_size=10
        )
        
        # Verify results
        self.assertEqual(results['total_compounds'], 2)
        self.assertEqual(results['imported'], 1)
        self.assertEqual(results['failed'], 1)
        
        # Verify error was recorded in details
        self.assertEqual(results['details']['CHEMBL456']['status'], 'failed')

if __name__ == '__main__':
    unittest.main()
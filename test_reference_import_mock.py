#!/usr/bin/env python3
"""
Test script for reference compound import using mock utilities.

This script tests the fixes to the import_reference_compounds.py script
using mock versions of the utility modules to avoid dependency issues.
"""

import os
import sys
import logging
import json
import uuid
from datetime import datetime
from typing import Dict, List, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('test_reference_import_mock')

# Mock ChEMBL client
class MockChEMBLClient:
    def get_molecule_by_chembl_id(self, chembl_id):
        """Return mock molecule data for a ChEMBL ID."""
        return {
            'molecule_id': chembl_id,
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
                'canonical_smiles': 'CCCCCCCCCC(=O)O',
                'standard_inchi': 'InChI=1S/C10H20O2/c1-2-3-4-5-6-7-8-9-10(11)12/h2-9H2,1H3,(H,11,12)',
                'standard_inchi_key': 'GHVNFZFCNZKVNT-UHFFFAOYSA-N'
            },
            'cross_references': {
                'pubchem_cid': '12345'
            }
        }

# Mock PubChem client
class MockPubChemClient:
    def get_molecule_properties(self, pubchem_cid):
        """Return mock property data for a PubChem CID."""
        return {
            'cid': pubchem_cid,
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
                    'value': {'sval': '32 C'}
                },
                {
                    'urn': {'label': 'Boiling Point'},
                    'value': {'sval': '231 C'}
                }
            ]
        }

# Mock CryoprotectantIdentifierManager
class MockCryoprotectantIdentifierManager:
    _instance = None
    
    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = MockCryoprotectantIdentifierManager()
        return cls._instance
    
    def __init__(self):
        self.molecules = {}
    
    def resolve_identifier(self, chembl_id=None):
        """Resolve a ChEMBL ID to an internal ID."""
        if chembl_id in self.molecules:
            return self.molecules[chembl_id]['internal_id'], False
        return str(uuid.uuid4()), True
    
    def add_molecule(self, internal_id, molecule_data):
        """Add a molecule to the identifier manager."""
        if 'chembl_id' in molecule_data:
            self.molecules[molecule_data['chembl_id']] = {
                'internal_id': internal_id,
                **molecule_data
            }
    
    def save_identifiers(self):
        """Save identifiers to disk (mock implementation)."""
        logger.info(f"Mock save_identifiers: {len(self.molecules)} molecules")

# Mock function to get reference compound IDs
def get_reference_compound_ids():
    """Return a list of mock reference compound IDs."""
    return ['CHEMBL1098659', 'CHEMBL1487', 'CHEMBL262548']

# Mock function to execute a query
def execute_query(query, params=None, fetch_one=False, dict_cursor=True):
    """Mock implementation of execute_query."""
    logger.info(f"Mock execute_query: {query}")
    
    if "SELECT" in query and "information_schema.columns" in query:
        # Return mock schema information
        return [
            {'column_name': 'id', 'data_type': 'uuid'},
            {'column_name': 'molecule_id', 'data_type': 'uuid'},
            {'column_name': 'property_type_id', 'data_type': 'uuid'},
            {'column_name': 'numeric_value', 'data_type': 'double precision'},
            {'column_name': 'text_value', 'data_type': 'text'},
            {'column_name': 'boolean_value', 'data_type': 'boolean'},
            {'column_name': 'created_by', 'data_type': 'uuid'},
            {'column_name': 'created_at', 'data_type': 'timestamp with time zone'},
            {'column_name': 'updated_at', 'data_type': 'timestamp with time zone'}
        ]
    
    if fetch_one:
        return {'id': str(uuid.uuid4())}
    
    return []

# Mock function to get a database connection
def get_db():
    """Mock implementation of get_db."""
    from db_connection_utils_mock import get_db_connection
    with get_db_connection() as conn:
        return conn

# Import the mock modules
sys.path.insert(0, '.')
import db_connection_utils_mock
import transaction_utils_mock
import batch_utils_mock
import property_utils_mock

# Now patch the import_reference_compounds module
import importlib.util
import types

def create_mock_module(name):
    """Create a mock module with the given name."""
    spec = importlib.util.find_spec('types')
    module = types.ModuleType(name)
    return module

# Create mock modules
sys.modules['db_connection_utils'] = db_connection_utils_mock
sys.modules['transaction_utils'] = transaction_utils_mock
sys.modules['batch_utils'] = batch_utils_mock
sys.modules['property_utils'] = property_utils_mock

# Create and populate the chembl module
chembl_module = create_mock_module('chembl')
sys.modules['chembl'] = chembl_module
chembl_module.client = types.ModuleType('client')
chembl_module.client.ResilientChEMBLClient = MockChEMBLClient
chembl_module.reference_compounds = types.ModuleType('reference_compounds')
chembl_module.reference_compounds.get_reference_compound_ids = get_reference_compound_ids

# Create and populate the pubchem module
pubchem_module = create_mock_module('pubchem')
sys.modules['pubchem'] = pubchem_module
pubchem_module.client = types.ModuleType('client')
pubchem_module.client.ResilientPubChemClient = MockPubChemClient

# Create and populate the cryoprotectant_identifiers module
sys.modules['cryoprotectant_identifiers'] = types.ModuleType('cryoprotectant_identifiers')
sys.modules['cryoprotectant_identifiers'].CryoprotectantIdentifierManager = MockCryoprotectantIdentifierManager

# Create and populate the database.connection module
database_module = create_mock_module('database')
sys.modules['database'] = database_module
database_module.connection = types.ModuleType('connection')
database_module.connection.get_db_connection = db_connection_utils_mock.get_db_connection

# Create and populate the sql_executor module
sys.modules['sql_executor'] = types.ModuleType('sql_executor')
sys.modules['sql_executor'].execute_query = execute_query
sys.modules['sql_executor'].bulk_insert = lambda *args, **kwargs: len(args[1]) if len(args) > 1 else 0
sys.modules['sql_executor'].execute_batch = lambda *args, **kwargs: None
sys.modules['sql_executor'].with_transaction = db_connection_utils_mock.safe_transaction
sys.modules['sql_executor'].with_retry = lambda *args, **kwargs: lambda f: f
sys.modules['sql_executor'].process_in_batches = batch_utils_mock.process_in_batches
sys.modules['sql_executor'].get_db = get_db

# Now import the module we want to test
from import_reference_compounds import import_reference_compounds, verify_reference_compounds

def run_test():
    """
    Run a test of the reference compound import process.
    
    This function:
    1. Imports a small number of reference compounds
    2. Verifies that the critical properties are set correctly
    3. Logs the results
    """
    logger.info("Starting reference compound import test")
    
    # Create a test-specific checkpoint file
    checkpoint_path = f"checkpoints/test_reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    # Create a test-specific output report
    output_report = f"reports/test_reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    try:
        # Run the import with a small batch size for testing
        results = import_reference_compounds(
            output_report=output_report,
            checkpoint_path=checkpoint_path,
            resume=False,
            batch_size=2
        )
        
        logger.info(f"Import completed with results: {json.dumps(results, indent=2)}")
        
        # Verify the reference compounds
        logger.info("Verifying reference compounds...")
        complete_count, total_count = verify_reference_compounds(results)
        
        logger.info(f"Verification results: {complete_count}/{total_count} complete")
        
        # Check if all reference compounds have the required properties
        if complete_count == total_count:
            logger.info("TEST PASSED: All reference compounds have required properties")
            return True
        else:
            logger.error(f"TEST FAILED: Only {complete_count}/{total_count} reference compounds have required properties")
            return False
            
    except Exception as e:
        logger.error(f"Test failed with error: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    # Create directories if they don't exist
    os.makedirs("checkpoints", exist_ok=True)
    os.makedirs("reports", exist_ok=True)
    
    success = run_test()
    sys.exit(0 if success else 1)
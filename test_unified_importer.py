#!/usr/bin/env python3
"""
Comprehensive test suite for the unified chemical data importer.

This test suite covers all aspects of the unified importer, including:
- Configuration validation and loading
- Data source integration
- API rate limiting and resilience
- Data transformation
- Database operations
- Checkpoint management and resumable imports
- Error handling and reporting
- Progress tracking

The tests use pytest fixtures to mock external dependencies and
provide controlled test environments.
"""

import os
import sys
import json
import uuid
import time
import pytest
import logging
import datetime
import tempfile
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

# Configure logging for testing
@pytest.fixture(scope="session")
def setup_logging():
    """Set up logging for tests."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()]
    )
    return logging.getLogger("test_unified_importer")

# Temporary directory for test files and checkpoints
@pytest.fixture(scope="function")
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        original_dir = os.getcwd()
        os.chdir(tmpdirname)
        # Create necessary subdirectories
        os.makedirs("logs", exist_ok=True)
        os.makedirs("checkpoints", exist_ok=True)
        os.makedirs("reports", exist_ok=True)
        yield tmpdirname
        os.chdir(original_dir)

# Mock config for testing
@pytest.fixture
def mock_config():
    """Create a mock configuration for testing."""
    return {
        "CHECKPOINT_DIR": "checkpoints",
        "PUBCHEM": {
            "API_BASE_URL": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
            "RATE_LIMIT_DELAY": 0.05,  # Fast for testing
            "BATCH_SIZE": 10,
            "PROPERTY_ENDPOINTS": [
                "MolecularFormula,MolecularWeight,XLogP,TPSA",
                "HBondDonorCount,HBondAcceptorCount",
                "IsomericSMILES,InChI,InChIKey,IUPACName,Title"
            ]
        },
        "CHEMBL": {
            "RATE_LIMIT_DELAY": 0.05,  # Fast for testing
            "BATCH_SIZE": 5,
            "SEARCH_TERMS": [
                "cryoprotect",
                "glycerol",
                "dmso"
            ],
            "REFERENCE_COMPOUNDS": [
                "CHEMBL25",    # Aspirin
                "CHEMBL1234"   # Glycerol
            ]
        },
        "DATABASE": {
            "BATCH_SIZE": 5
        }
    }

# Mock ChEMBL client
@pytest.fixture
def mock_chembl_client():
    """Create a mock ChEMBL client for testing."""
    mock_client = MagicMock()
    
    # Mock molecule search results
    mock_molecule = mock_client.molecule
    
    # Setup filter response for search terms
    def mock_filter(**kwargs):
        """Mock the filter method to return test compounds."""
        # Return 3 compounds for each search term
        return [
            {
                'molecule_chembl_id': f"CHEMBL{i}",
                'pref_name': f"Test Compound {i}",
                'molecule_structures': {
                    'canonical_smiles': f"C{i}H{i}O{i}",
                    'standard_inchi': f"InChI=1S/C{i}H{i}O{i}",
                    'standard_inchi_key': f"TESTKEY{i}"
                }
            } for i in range(1, 4)
        ]
    
    # Setup get response for individual compounds
    def mock_get(chembl_id):
        """Mock the get method to return detailed compound data."""
        compound_id = int(chembl_id.replace("CHEMBL", ""))
        return {
            'molecule_chembl_id': chembl_id,
            'pref_name': f"Test Compound {compound_id}",
            'molecule_structures': {
                'canonical_smiles': f"C{compound_id}H{compound_id}O{compound_id}",
                'standard_inchi': f"InChI=1S/C{compound_id}H{compound_id}O{compound_id}",
                'standard_inchi_key': f"TESTKEY{compound_id}"
            },
            'molecule_properties': {
                'full_mwt': 100.0 + compound_id,
                'alogp': 1.0 + compound_id * 0.1,
                'hba': compound_id,
                'hbd': compound_id,
                'psa': 50.0 + compound_id,
                'rtb': compound_id,
                'full_molformula': f"C{compound_id}H{compound_id}O{compound_id}"
            }
        }
    
    mock_molecule.filter.side_effect = mock_filter
    mock_molecule.get.side_effect = mock_get
    
    return mock_client

# Mock PubChem client session
class MockResponse:
    """Mock aiohttp response object."""
    def __init__(self, data=None, status=200):
        self.data = data
        self.status = status
        
    async def json(self):
        return self.data
        
    async def __aenter__(self):
        return self
        
    async def __aexit__(self, exc_type, exc, tb):
        pass

@pytest.fixture
def mock_pubchem_session():
    """Create a mock aiohttp session for PubChem API testing."""
    session = MagicMock()
    
    # Create mock response based on CID
    async def mock_get(url):
        """Return mock data based on the URL."""
        # Extract CID from URL
        try:
            cid = int(url.split('/cid/')[1].split('/')[0])
            
            # Rate limit error for testing error handling
            if cid == 999:
                return MockResponse(status=429)
            
            # Server error for testing error handling
            if cid == 9999:
                return MockResponse(status=500)
            
            # Normal response
            data = {
                "PropertyTable": {
                    "Properties": [{
                        "CID": cid,
                        "MolecularFormula": f"C{cid % 10}H{(cid % 5) * 2}O{cid % 3}",
                        "MolecularWeight": 100.0 + (cid % 50),
                        "XLogP": (cid % 10) - 3,
                        "TPSA": 50.0 + (cid % 30),
                        "HBondDonorCount": cid % 5,
                        "HBondAcceptorCount": cid % 6,
                        "IsomericSMILES": f"C{cid % 10}CC{cid % 5}O",
                        "InChI": f"InChI=1S/C{cid % 10}H{(cid % 5) * 2}O{cid % 3}",
                        "InChIKey": f"TESTPUBCHEM{cid:010d}",
                        "IUPACName": f"Test Compound {cid}",
                        "Title": f"Test PubChem Compound {cid}"
                    }]
                }
            }
            return MockResponse(data=data)
        except Exception:
            return MockResponse(status=400)
    
    session.get.side_effect = mock_get
    return session

# Mock database connection
@pytest.fixture
def mock_db_connection():
    """Create a mock database connection for testing."""
    db = MagicMock()
    
    # Mock execute_query for SELECT operations
    def mock_execute_query(sql, params=None):
        """Mock query execution."""
        # Return a list of property types for property type queries
        if "SELECT" in sql and "property_types" in sql:
            return [
                {"id": str(uuid.uuid4()), "name": "molecular weight", "data_type": "numeric"},
                {"id": str(uuid.uuid4()), "name": "logp", "data_type": "numeric"},
                {"id": str(uuid.uuid4()), "name": "hydrogen bond acceptor count", "data_type": "numeric"},
                {"id": str(uuid.uuid4()), "name": "unknown property type", "data_type": "text"}
            ]
        
        # Check if molecule exists by inchikey
        if "SELECT" in sql and "molecules" in sql and "inchikey" in sql:
            # Return ID for test inchikey, None otherwise
            if params and "TESTKEY" in str(params):
                return [{"id": str(uuid.uuid4())}]
            return []
        
        # Default to empty list for other queries
        return []
    
    # Mock execute_batch for batch operations
    def mock_execute_batch(queries, transaction=True):
        """Mock batch execution."""
        # Always return success
        return True
    
    # Mock transaction context manager
    class MockTransaction:
        def __enter__(self):
            return MagicMock()
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            pass
    
    # Assign mocks to the db mock
    db.execute_query.side_effect = mock_execute_query
    db.execute_batch.side_effect = mock_execute_batch
    db.transaction.return_value = MockTransaction()
    
    # Mock get_instance to return the same mock
    db.get_instance.return_value = db
    
    return db

# Test CID file for PubChem testing
@pytest.fixture
def mock_cid_file(temp_dir):
    """Create a mock CID file for testing."""
    with open("CID-Synonym-test", "w") as f:
        for i in range(1, 21):  # 20 test CIDs
            f.write(f"{i}\tTest Compound {i}\n")
    return "CID-Synonym-test"

# Test data for validating transformation functions
@pytest.fixture
def test_compounds():
    """Create test compound data for transformation testing."""
    chembl_compound = {
        'molecule_chembl_id': "CHEMBL1234",
        'pref_name': "Test ChEMBL Compound",
        'molecule_structures': {
            'canonical_smiles': "CCO",
            'standard_inchi': "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            'standard_inchi_key': "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        },
        'molecule_properties': {
            'full_mwt': 46.07,
            'alogp': -0.19,
            'hba': 1,
            'hbd': 1,
            'psa': 20.23,
            'rtb': 1,
            'full_molformula': "C2H6O"
        }
    }
    
    pubchem_compound = {
        "CID": 702,
        "Molecular Formula": "C2H6O",
        "Molecular Weight": 46.07,
        "LogP": -0.31,
        "TPSA": 20.23,
        "H-Bond Donors": 1,
        "H-Bond Acceptors": 1,
        "SMILES": "CCO",
        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "IUPACName": "ethanol",
        "Title": "Ethanol",
        "PubChem Link": "https://pubchem.ncbi.nlm.nih.gov/compound/702"
    }
    
    return {
        "chembl": chembl_compound,
        "pubchem": pubchem_compound
    }

# Checkpoint fixtures
@pytest.fixture
def mock_checkpoint_file(temp_dir):
    """Create a mock checkpoint file for testing resume functionality."""
    checkpoint_data = {
        "last_completed_batch": 2,
        "total_processed": 20,
        "total_imported": 15,
        "total_skipped": 3,
        "total_errors": 2,
        "batch_times": [1.2, 1.3, 1.1],
        "elapsed_seconds": 3.6,
        "timestamp": datetime.datetime.now().isoformat(),
        "status": "Running"
    }
    
    checkpoint_path = os.path.join("checkpoints", "import_checkpoint.json")
    with open(checkpoint_path, "w") as f:
        json.dump(checkpoint_data, f)
    
    return checkpoint_path

# MCP tool mock
@pytest.fixture
def mock_mcp_tool():
    """Create a mock for the MCP tool integration."""
    with patch("use_mcp_tool.use_mcp_tool") as mock_tool:
        # Configure the mock to return an array with one result for SQL queries
        def mock_mcp_response(tool, action, params):
            if tool == "supabase" and action == "execute_sql":
                # Return sample ID for molecules
                return [{"id": str(uuid.uuid4()), "pubchem_cid": params.get("pubchem_cid", "unknown")}]
            return {}
        
        mock_tool.side_effect = mock_mcp_response
        yield mock_tool

# Tests for configuration validation and loading
class TestConfiguration:
    """Tests for configuration handling in the unified importer."""
    
    def test_validate_config(self, mock_config):
        """Test configuration validation."""
        # Import the function to test
        from unified_importer import validate_config
        
        # Valid config should not raise exceptions
        validate_config(mock_config)
        
        # Test with missing required sections
        invalid_config = mock_config.copy()
        del invalid_config["PUBCHEM"]
        
        with pytest.raises(ValueError):
            validate_config(invalid_config)
        
        # Test with missing required parameters
        invalid_config = mock_config.copy()
        invalid_config["PUBCHEM"] = {}
        
        with pytest.raises(ValueError):
            validate_config(invalid_config)
    
    def test_load_config(self, temp_dir):
        """Test loading configuration from file."""
        # Import the function to test
        from unified_importer import load_config
        
        # Create a test config file
        config_data = {
            "CHECKPOINT_DIR": "checkpoints",
            "PUBCHEM": {
                "API_BASE_URL": "https://test-url.com",
                "RATE_LIMIT_DELAY": 0.1,
                "BATCH_SIZE": 10
            },
            "CHEMBL": {
                "RATE_LIMIT_DELAY": 0.2,
                "BATCH_SIZE": 5
            },
            "DATABASE": {
                "BATCH_SIZE": 5
            }
        }
        
        with open("config.json", "w") as f:
            json.dump(config_data, f)
        
        # Load the config
        config = load_config("config.json")
        
        # Verify the config was loaded correctly
        assert config["CHECKPOINT_DIR"] == "checkpoints"
        assert config["PUBCHEM"]["API_BASE_URL"] == "https://test-url.com"
        assert config["PUBCHEM"]["RATE_LIMIT_DELAY"] == 0.1
        assert config["CHEMBL"]["BATCH_SIZE"] == 5
        
        # Test loading with non-existent file (should use default)
        config = load_config("non_existent_config.json")
        assert "CHECKPOINT_DIR" in config
        assert "PUBCHEM" in config
        assert "CHEMBL" in config

# Tests for data source integration
class TestDataSourceIntegration:
    """Tests for data source integration in the unified importer."""
    
    def test_chembl_fetch_compounds(self, mock_chembl_client, mock_config):
        """Test fetching compounds from ChEMBL."""
        # Import the function to test
        from unified_importer import fetch_chembl_compounds
        
        # Patch the ChEMBL client
        with patch("unified_importer.new_client", mock_chembl_client):
            # Fetch compounds
            compounds = fetch_chembl_compounds(mock_config["CHEMBL"], limit=5)
            
            # Verify compounds were fetched
            assert len(compounds) > 0
            assert "molecule_chembl_id" in compounds[0]
            assert "molecule_properties" in compounds[0]
    
    @pytest.mark.asyncio
    async def test_pubchem_fetch_compounds(self, mock_pubchem_session, mock_cid_file, mock_config):
        """Test fetching compounds from PubChem."""
        # Import the function to test
        from unified_importer import fetch_pubchem_compounds_async
        
        # Patch the aiohttp session
        with patch("aiohttp.ClientSession", return_value=mock_pubchem_session):
            # Fetch compounds
            compounds = await fetch_pubchem_compounds_async(
                mock_config["PUBCHEM"],
                cid_file=mock_cid_file,
                limit=5
            )
            
            # Verify compounds were fetched
            assert len(compounds) > 0
            assert "CID" in compounds[0]
            assert "InChIKey" in compounds[0]
    
    def test_get_cid_list(self, mock_cid_file):
        """Test reading CIDs from file."""
        # Import the function to test
        from unified_importer import get_cid_list
        
        # Get CIDs
        cids = get_cid_list(mock_cid_file)
        
        # Verify CIDs were read
        assert len(cids) == 20
        assert cids[0] == 1
        assert cids[-1] == 20

# Tests for data filtering and transformation
class TestDataTransformation:
    """Tests for data filtering and transformation in the unified importer."""
    
    def test_filter_molecule(self, test_compounds):
        """Test molecule filtering logic."""
        # Import the function to test
        from unified_importer import filter_molecule
        
        # Valid molecules should pass filtering
        assert filter_molecule(test_compounds["chembl"], "chembl") == True
        assert filter_molecule(test_compounds["pubchem"], "pubchem") == True
        
        # Test with invalid data
        invalid_chembl = test_compounds["chembl"].copy()
        invalid_chembl["molecule_properties"]["full_mwt"] = 5000  # Outside MW range
        assert filter_molecule(invalid_chembl, "chembl") == False
        
        invalid_pubchem = test_compounds["pubchem"].copy()
        invalid_pubchem["LogP"] = 10  # Outside LogP range
        assert filter_molecule(invalid_pubchem, "pubchem") == False
        
        # Test with missing required fields
        missing_fields_chembl = test_compounds["chembl"].copy()
        del missing_fields_chembl["molecule_structures"]["standard_inchi_key"]
        assert filter_molecule(missing_fields_chembl, "chembl") == False
        
        missing_fields_pubchem = test_compounds["pubchem"].copy()
        del missing_fields_pubchem["InChIKey"]
        assert filter_molecule(missing_fields_pubchem, "pubchem") == False
    
    def test_transform_chembl_to_unified(self, test_compounds):
        """Test transformation of ChEMBL data to unified format."""
        # Import the function to test
        from unified_importer import transform_chembl_to_unified
        
        # Transform the test compound
        user_id = "test-user-id"
        unified = transform_chembl_to_unified(test_compounds["chembl"], user_id)
        
        # Verify transformation
        assert unified["name"] == "Test ChEMBL Compound"
        assert unified["smiles"] == "CCO"
        assert unified["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert unified["molecular_weight"] == 46.07
        assert unified["chembl_id"] == "CHEMBL1234"
        assert unified["data_source"].startswith("ChEMBL")
        assert "version" in unified
        assert "created_by" in unified
        assert unified["created_by"] == user_id
        
        # Test with missing name (should use fallback)
        no_name = test_compounds["chembl"].copy()
        del no_name["pref_name"]
        unified = transform_chembl_to_unified(no_name, user_id)
        assert unified["name"] == "CHEMBL1234"  # Uses ChEMBL ID as fallback
    
    def test_transform_pubchem_to_unified(self, test_compounds):
        """Test transformation of PubChem data to unified format."""
        # Import the function to test
        from unified_importer import transform_pubchem_to_unified
        
        # Transform the test compound
        user_id = "test-user-id"
        unified = transform_pubchem_to_unified(test_compounds["pubchem"], user_id)
        
        # Verify transformation
        assert unified["name"] == "Ethanol"
        assert unified["smiles"] == "CCO"
        assert unified["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert unified["molecular_weight"] == 46.07
        assert unified["pubchem_cid"] == 702
        assert unified["data_source"].startswith("PubChem")
        assert "version" in unified
        assert "created_by" in unified
        assert unified["created_by"] == user_id
    
    def test_transform_to_properties(self, test_compounds):
        """Test transformation of compound data to molecular properties."""
        # Import the function to test
        from unified_importer import transform_to_properties
        
        # Create mock property type map
        property_types = {
            "molecular weight": "weight-id",
            "logp": "logp-id",
            "hydrogen bond acceptor count": "hba-id",
            "hydrogen bond donor count": "hbd-id",
            "topological polar surface area": "tpsa-id",
            "unknown property type": "unknown-id"
        }
        
        # Transform ChEMBL properties
        molecule_id = "test-molecule-id"
        user_id = "test-user-id"
        chembl_properties = transform_to_properties(
            test_compounds["chembl"], 
            molecule_id, 
            user_id, 
            property_types, 
            "chembl"
        )
        
        # Verify ChEMBL properties
        assert len(chembl_properties) > 0
        found_mw = False
        for prop in chembl_properties:
            assert "molecule_id" in prop
            assert "property_type_id" in prop
            assert prop["molecule_id"] == molecule_id
            assert prop["created_by"] == user_id
            if prop["property_type_id"] == property_types["molecular weight"]:
                assert prop["numeric_value"] == 46.07
                found_mw = True
        assert found_mw, "Molecular weight property not found"
        
        # Transform PubChem properties
        pubchem_properties = transform_to_properties(
            test_compounds["pubchem"], 
            molecule_id, 
            user_id, 
            property_types, 
            "pubchem"
        )
        
        # Verify PubChem properties
        assert len(pubchem_properties) > 0
        found_mw = False
        for prop in pubchem_properties:
            assert "molecule_id" in prop
            assert "property_type_id" in prop
            assert prop["molecule_id"] == molecule_id
            assert prop["created_by"] == user_id
            if prop["property_type_id"] == property_types["molecular weight"]:
                assert prop["numeric_value"] == 46.07
                found_mw = True
        assert found_mw, "Molecular weight property not found"

# Tests for database operations
class TestDatabaseOperations:
    """Tests for database operations in the unified importer."""
    
    def test_get_property_types(self, mock_db_connection):
        """Test retrieving property types from database."""
        # Import the function to test
        from unified_importer import get_property_types
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Get property types
            property_types = get_property_types()
            
            # Verify property types were retrieved
            assert len(property_types) > 0
            assert "molecular weight" in property_types
            assert "logp" in property_types
            assert "unknown property type" in property_types
    
    def test_check_molecule_exists(self, mock_db_connection):
        """Test checking if a molecule exists in the database."""
        # Import the function to test
        from unified_importer import check_molecule_exists
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Check for existing molecule
            molecule_id = check_molecule_exists("TESTKEY123")
            
            # Verify molecule was found
            assert molecule_id is not None
            
            # Check for non-existent molecule
            molecule_id = check_molecule_exists("NONEXISTENT")
            
            # Verify molecule was not found
            assert molecule_id is None
    
    def test_insert_molecule(self, mock_db_connection, test_compounds):
        """Test inserting a molecule into the database."""
        # Import the function to test
        from unified_importer import insert_molecule
        
        # Prepare molecule data
        from unified_importer import transform_chembl_to_unified
        molecule_data = transform_chembl_to_unified(test_compounds["chembl"], "test-user-id")
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Insert molecule
            molecule_id = insert_molecule(molecule_data)
            
            # Verify molecule was inserted
            assert molecule_id is not None
            
            # Verify correct SQL was generated
            called_args = mock_db_connection.execute_query.call_args[0]
            assert "INSERT INTO molecules" in called_args[0]
            assert "RETURNING id" in called_args[0]
    
    def test_insert_property(self, mock_db_connection, test_compounds):
        """Test inserting a molecular property into the database."""
        # Import the function to test
        from unified_importer import insert_property
        
        # Prepare property data
        property_data = {
            "molecule_id": "test-molecule-id",
            "property_type_id": "test-property-type-id",
            "numeric_value": 46.07,
            "text_value": None,
            "boolean_value": None,
            "created_by": "test-user-id",
            "data_source": "Test data source",
            "version": 1
        }
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Insert property
            success = insert_property(property_data)
            
            # Verify property was inserted
            assert success is True
            
            # Verify correct SQL was generated
            called_args = mock_db_connection.execute_query.call_args[0]
            assert "INSERT INTO molecular_properties" in called_args[0]
            assert "RETURNING id" in called_args[0]
    
    def test_batch_insert_molecules(self, mock_db_connection, test_compounds):
        """Test batch insertion of molecules into the database."""
        # Import the function to test
        from unified_importer import batch_insert_molecules
        
        # Prepare molecule data
        from unified_importer import transform_chembl_to_unified, transform_pubchem_to_unified
        molecules = [
            transform_chembl_to_unified(test_compounds["chembl"], "test-user-id"),
            transform_pubchem_to_unified(test_compounds["pubchem"], "test-user-id")
        ]
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Insert molecules in batch
            result = batch_insert_molecules(molecules)
            
            # Verify molecules were inserted
            assert result["inserted"] > 0
            assert result["errors"] == 0
            
            # Verify batch execution was called
            mock_db_connection.execute_batch.assert_called_once()
    
    def test_batch_insert_properties(self, mock_db_connection, test_compounds):
        """Test batch insertion of properties into the database."""
        # Import the function to test
        from unified_importer import batch_insert_properties
        
        # Prepare property data
        properties = [
            {
                "molecule_id": "test-molecule-id",
                "property_type_id": "test-property-type-id",
                "numeric_value": 46.07,
                "text_value": None,
                "boolean_value": None,
                "created_by": "test-user-id",
                "data_source": "Test data source",
                "version": 1
            },
            {
                "molecule_id": "test-molecule-id",
                "property_type_id": "test-property-type-id-2",
                "numeric_value": 1.0,
                "text_value": None,
                "boolean_value": None,
                "created_by": "test-user-id",
                "data_source": "Test data source",
                "version": 1
            }
        ]
        
        # Patch the database connection
        with patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            # Insert properties in batch
            result = batch_insert_properties(properties)
            
            # Verify properties were inserted
            assert result["inserted"] > 0
            assert result["errors"] == 0
            
            # Verify batch execution was called
            mock_db_connection.execute_batch.assert_called_once()

# Tests for checkpoint management
class TestCheckpointManagement:
    """Tests for checkpoint management in the unified importer."""
    
    def test_load_checkpoint(self, mock_checkpoint_file):
        """Test loading checkpoint data from file."""
        # Import the function to test
        from unified_importer import load_checkpoint
        
        # Load the checkpoint
        checkpoint = load_checkpoint(mock_checkpoint_file)
        
        # Verify checkpoint was loaded
        assert checkpoint is not None
        assert checkpoint["last_completed_batch"] == 2
        assert checkpoint["total_processed"] == 20
        assert checkpoint["total_imported"] == 15
        assert checkpoint["status"] == "Running"
        
        # Test loading non-existent checkpoint
        checkpoint = load_checkpoint("non_existent_checkpoint.json")
        assert checkpoint is None
    
    def test_save_checkpoint(self, temp_dir):
        """Test saving checkpoint data to file."""
        # Import the function to test
        from unified_importer import save_checkpoint
        
        # Create checkpoint data
        checkpoint_data = {
            "last_completed_batch": 3,
            "total_processed": 30,
            "total_imported": 25,
            "total_skipped": 3,
            "total_errors": 2,
            "batch_times": [1.2, 1.3, 1.1, 1.4],
            "elapsed_seconds": 5.0,
            "status": "Running"
        }
        
        # Save the checkpoint
        checkpoint_path = os.path.join("checkpoints", "test_checkpoint.json")
        save_checkpoint(checkpoint_path, checkpoint_data)
        
        # Verify checkpoint was saved
        assert os.path.exists(checkpoint_path)
        
        # Load the checkpoint to verify contents
        with open(checkpoint_path, "r") as f:
            loaded_data = json.load(f)
        
        assert loaded_data["last_completed_batch"] == 3
        assert loaded_data["total_processed"] == 30
        assert loaded_data["total_imported"] == 25
        assert loaded_data["status"] == "Running"

# Tests for progress tracking
class TestProgressTracking:
    """Tests for progress tracking in the unified importer."""
    
    def test_progress_tracker_initialization(self):
        """Test ProgressTracker initialization."""
        # Import the class to test
        from unified_importer import ProgressTracker
        
        # Create a tracker
        tracker = ProgressTracker(total_compounds=100, batch_size=10)
        
        # Verify initialization
        assert tracker.total_compounds == 100
        assert tracker.batch_size == 10
        assert tracker.total_batches == 10
        assert tracker.current_batch == 0
        assert tracker.total_processed == 0
        assert tracker.total_imported == 0
        assert tracker.total_skipped == 0
        assert tracker.total_errors == 0
        assert tracker.status == "Running"
    
    def test_progress_tracker_batch_processing(self):
        """Test ProgressTracker batch processing."""
        # Import the class to test
        from unified_importer import ProgressTracker
        
        # Create a tracker
        tracker = ProgressTracker(total_compounds=100, batch_size=10)
        
        # Start a batch
        tracker.start_batch(0, 10)
        assert tracker.current_batch == 1  # 1-indexed for display
        assert tracker.compounds_in_current_batch == 10
        
        # End the batch
        tracker.end_batch(10, 8, 2, 0)
        assert tracker.total_processed == 10
        assert tracker.total_imported == 8
        assert tracker.total_skipped == 2
        assert tracker.total_errors == 0
        assert len(tracker.batch_times) == 1
        
        # Process another batch
        tracker.start_batch(1, 10)
        tracker.end_batch(10, 7, 2, 1)
        assert tracker.total_processed == 20
        assert tracker.total_imported == 15
        assert tracker.total_skipped == 4
        assert tracker.total_errors == 1
        assert len(tracker.batch_times) == 2
    
    def test_progress_tracker_logging(self):
        """Test ProgressTracker logging functionality."""
        # Import the class to test
        from unified_importer import ProgressTracker
        
        # Create a tracker
        tracker = ProgressTracker(total_compounds=100, batch_size=10)
        
        # Add log messages
        tracker.add_log_message("Test log message 1")
        tracker.add_log_message("Test log message 2")
        
        # Verify log messages were added
        assert len(tracker.recent_logs) == 2
        assert "Test log message 1" in tracker.recent_logs[0]
        assert "Test log message 2" in tracker.recent_logs[1]
        
        # Add error
        tracker.add_error("Test error message")
        assert tracker.total_errors == 1
        assert "ERROR" in tracker.recent_logs[2]
        
        # Add skipped
        tracker.add_skipped("Test CID", "Test reason")
        assert tracker.total_skipped == 1
        assert "SKIPPED" in tracker.recent_logs[3]
    
    def test_progress_tracker_metrics(self):
        """Test ProgressTracker metrics calculation."""
        # Import the class to test
        from unified_importer import ProgressTracker
        
        # Create a tracker with mock data
        tracker = ProgressTracker(total_compounds=100, batch_size=10)
        tracker.total_processed = 50
        tracker.batch_times = [1.0, 1.5, 2.0, 1.8, 1.7]
        tracker.start_time = time.time() - 10  # 10 seconds ago
        
        # Calculate metrics
        progress = tracker.get_progress_percentage()
        eta = tracker.estimate_time_remaining()
        elapsed = tracker.get_elapsed_time()
        avg_time = tracker.get_avg_time_per_batch()
        
        # Verify metrics
        assert progress == 50.0  # 50/100 = 50%
        assert isinstance(eta, str)
        assert isinstance(elapsed, str)
        assert isinstance(avg_time, str)
        
        # Get progress data
        data = tracker.get_progress_data()
        assert data["total_compounds"] == 100
        assert data["total_processed"] == 50
        assert data["progress_percentage"] == 50.0
        assert data["status"] == "Running"

# Tests for rate limiting and resilience
class TestRateLimitingResilience:
    """Tests for rate limiting and resilience in the unified importer."""
    
    @pytest.mark.asyncio
    async def test_rate_limiting_pubchem(self, mock_pubchem_session):
        """Test rate limiting for PubChem API."""
        # Import the function to test
        from unified_importer import get_pubchem_properties_async
        
        # Test with normal CID
        with patch("aiohttp.ClientSession", return_value=mock_pubchem_session):
            # Create a semaphore for testing
            import asyncio
            semaphore = asyncio.Semaphore(5)
            
            # Get properties for a normal CID
            result = await get_pubchem_properties_async(
                mock_pubchem_session, 
                123, 
                rate_limit_delay=0.01,
                semaphore=semaphore
            )
            
            # Verify properties were retrieved
            assert result is not None
            assert "CID" in result
            assert result["CID"] == 123
            
            # Test with rate-limited CID (CID 999 is configured to return 429)
            result = await get_pubchem_properties_async(
                mock_pubchem_session, 
                999, 
                rate_limit_delay=0.01,
                semaphore=semaphore,
                max_retries=2
            )
            
            # Verify rate limit was handled
            assert result is not None
            assert "Error" in result
            assert "Rate limit" in result["Error"]
            assert "RateLimit" in result
            assert result["RateLimit"] is True
    
    def test_resilient_batch_processing(self, mock_db_connection, test_compounds):
        """Test resilient batch processing with error handling."""
        # Import the function to test
        from unified_importer import process_batch_resilient
        
        # Prepare test data
        batch = [test_compounds["chembl"], test_compounds["pubchem"]]
        batch_num = 1
        
        # Configure the mock database to simulate an error on first attempt then success
        db_mock = mock_db_connection
        call_count = 0
        
        def mock_execute_batch_with_retry(queries, transaction=True):
            nonlocal call_count
            call_count += 1
            # Fail on first attempt, succeed on retry
            return call_count > 1
        
        db_mock.execute_batch.side_effect = mock_execute_batch_with_retry
        
        # Patch database and transformation functions
        with patch("unified_importer.get_db_connection", return_value=db_mock), \
             patch("unified_importer.transform_chembl_to_unified", return_value={"name": "Test"}), \
             patch("unified_importer.transform_pubchem_to_unified", return_value={"name": "Test"}), \
             patch("unified_importer.filter_molecule", return_value=True):
            
            # Process batch with resilience
            result = process_batch_resilient(
                batch=batch,
                batch_num=batch_num,
                user_id="test-user-id",
                property_types={},
                max_retries=3,
                source="chembl"
            )
            
            # Verify batch was processed successfully on retry
            assert result["processed"] == 2
            assert result["imported"] > 0
            assert result["retries"] == 1  # One retry needed
            assert db_mock.execute_batch.call_count == 2  # Called twice (fail, then success)

# Tests for error handling and reporting
class TestErrorHandling:
    """Tests for error handling and reporting in the unified importer."""
    
    def test_error_logging(self, temp_dir, setup_logging):
        """Test error logging functionality."""
        # Import the function to test
        from unified_importer import log_error, log_skipped_molecule
        
        # Create logs directory
        os.makedirs("logs", exist_ok=True)
        
        # Log an error
        error_context = {
            "compound_id": "TEST123",
            "source": "test_function",
            "exception": ValueError("Test error")
        }
        log_error("API", "Test error message", error_context)
        
        # Log a skipped molecule
        molecule_data = {"name": "Test Molecule", "inchikey": "TESTKEY123"}
        log_skipped_molecule("TEST123", "Missing InChI", molecule_data, "validation")
        
        # Verify error log file was created
        error_log_path = os.path.join("logs", "chembl_import_errors.jsonl")
        assert os.path.exists(error_log_path)
        
        # Verify skipped log file was created
        skipped_log_path = os.path.join("logs", "chembl_import_skipped.jsonl")
        assert os.path.exists(skipped_log_path)
    
    def test_error_reporting(self, temp_dir):
        """Test error reporting and statistics."""
        # Import the function to test
        from unified_importer import generate_import_report
        
        # Create logs directory
        os.makedirs("logs", exist_ok=True)
        os.makedirs("reports", exist_ok=True)
        
        # Generate a report
        stats = {
            "total_processed": 100,
            "total_imported": 80,
            "total_skipped": 15,
            "total_errors": 5,
            "rate_limit_errors": 2,
            "batches_processed": 10,
            "batches_failed": 1
        }
        
        report = generate_import_report(
            stats=stats,
            elapsed_time=120.5,
            config={"target": 100},
            source="UNIFIED"
        )
        
        # Verify report was generated
        assert os.path.exists(report["report_path"])
        
        # Verify report contents
        with open(report["report_path"], "r") as f:
            report_content = json.load(f)
        
        assert report_content["statistics"]["total_processed"] == 100
        assert report_content["statistics"]["total_imported"] == 80
        assert report_content["statistics"]["success_rate"] == 80.0
        assert report_content["statistics"]["error_rate"] == 5.0
        assert report_content["performance"]["elapsed_time_seconds"] == 120.5
        assert report_content["performance"]["compounds_per_second"] > 0

# Tests for MCP integration
class TestMCPIntegration:
    """Tests for MCP integration in the unified importer."""
    
    def test_mcp_molecule_insert(self, mock_mcp_tool, test_compounds):
        """Test inserting molecules using MCP."""
        # Import the function to test
        from unified_importer import batch_insert_molecules_mcp
        
        # Prepare molecule data
        from unified_importer import transform_pubchem_to_unified
        molecules = [
            transform_pubchem_to_unified(test_compounds["pubchem"], "test-user-id")
        ]
        
        # Insert molecules
        project_id = "test-project-id"
        result = batch_insert_molecules_mcp(molecules, project_id)
        
        # Verify MCP was called correctly
        mock_mcp_tool.assert_called_once()
        args = mock_mcp_tool.call_args[0]
        assert args[0] == "supabase"
        assert args[1] == "execute_sql"
        
        # Verify SQL contains expected elements
        sql = args[2]["query"]
        assert "INSERT INTO molecules" in sql
        assert "ON CONFLICT" in sql  # Should use upsert pattern
        assert "pubchem_cid" in sql
        assert "RETURNING" in sql
    
    def test_mcp_direct_connection(self, mock_mcp_tool):
        """Test direct SQL execution using MCP."""
        # Import the function to test
        from unified_importer import execute_sql_mcp
        
        # Execute SQL
        project_id = "test-project-id"
        sql = "SELECT * FROM test_table"
        params = {"param1": "value1"}
        
        result = execute_sql_mcp(project_id, sql, params)
        
        # Verify MCP was called correctly
        mock_mcp_tool.assert_called_once()
        args = mock_mcp_tool.call_args[0]
        assert args[0] == "supabase"
        assert args[1] == "execute_sql"
        assert args[2]["project_id"] == project_id
        assert "param1" in args[2]["query"]

# Integration test for the core import process
class TestCoreImportProcess:
    """Integration tests for the core import process."""
    
    def test_import_from_pubchem(self, mock_pubchem_session, mock_db_connection, 
                                mock_cid_file, mock_config, temp_dir, mock_mcp_tool):
        """Test importing data from PubChem."""
        # Import the function to test
        from unified_importer import import_from_pubchem
        
        # Patch external dependencies
        with patch("aiohttp.ClientSession", return_value=mock_pubchem_session), \
             patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            
            # Import data
            result = import_from_pubchem(
                config=mock_config,
                cid_file=mock_cid_file,
                limit=5,
                user_id="test-user-id",
                resume=False,
                project_id="test-project-id"
            )
            
            # Verify import was successful
            assert result["total_processed"] > 0
            assert result["total_imported"] > 0
            assert result["success"] == True
    
    def test_import_from_chembl(self, mock_chembl_client, mock_db_connection, 
                               mock_config, temp_dir):
        """Test importing data from ChEMBL."""
        # Import the function to test
        from unified_importer import import_from_chembl
        
        # Patch external dependencies
        with patch("unified_importer.new_client", mock_chembl_client), \
             patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            
            # Import data
            result = import_from_chembl(
                config=mock_config,
                limit=5,
                user_id="test-user-id",
                resume=False
            )
            
            # Verify import was successful
            assert result["total_processed"] > 0
            assert result["total_imported"] > 0
            assert result["success"] == True
    
    def test_unified_import_process(self, mock_chembl_client, mock_pubchem_session, 
                                  mock_db_connection, mock_cid_file, mock_config, 
                                  temp_dir, mock_mcp_tool):
        """Test the unified import process from both sources."""
        # Import the function to test
        from unified_importer import unified_import
        
        # Patch external dependencies
        with patch("unified_importer.new_client", mock_chembl_client), \
             patch("aiohttp.ClientSession", return_value=mock_pubchem_session), \
             patch("unified_importer.get_db_connection", return_value=mock_db_connection):
            
            # Run unified import
            result = unified_import(
                config=mock_config,
                sources=["pubchem", "chembl"],
                limit=5,
                user_id="test-user-id",
                resume=False,
                project_id="test-project-id",
                cid_file=mock_cid_file
            )
            
            # Verify import was successful
            assert result["pubchem"]["total_processed"] > 0
            assert result["chembl"]["total_imported"] > 0
            assert result["success"] == True
            
            # Verify report was generated
            assert os.path.exists(result["report_path"])
            
            # Verify both sources were processed
            assert "pubchem" in result
            assert "chembl" in result
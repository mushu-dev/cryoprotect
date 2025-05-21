"""
Integration tests for the unified chemical data importer.

This module contains tests for the unified chemical data importer framework,
including configuration, data sources, database operations, and utilities.
"""

import os
import json
import pytest
import asyncio
from unittest.mock import MagicMock, patch, AsyncMock
from typing import Dict, List, Any, Optional

# Import framework components
from unified_importer.core.config import ImporterConfig
from unified_importer.core.logging import setup_logging
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker
from unified_importer.core.database import DatabaseOperations
from unified_importer.core.validation import MoleculeValidator

from unified_importer.sources.source_base import MolecularDataSource
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.sources.chembl_source import ChEMBLDataSource


# Test fixtures

@pytest.fixture
def sample_config():
    """Fixture providing a sample configuration."""
    return {
        'batch_size': 10,
        'max_workers': 4,
        'api_delay': 0.1,
        'max_retries': 2,
        'retry_delay': 0.1,
        'checkpoint_interval': 5,
        'log_level': 'INFO',
        'dry_run': True,
        'database': {
            'type': 'direct',
            'connection': {
                'dbname': 'test_db',
                'user': 'test_user',
                'password': 'test_password',
                'host': 'localhost',
                'port': 5432
            }
        }
    }


@pytest.fixture
def temp_checkpoint_file(tmp_path):
    """Fixture providing a temporary checkpoint file path."""
    checkpoint_dir = tmp_path / "checkpoints"
    checkpoint_dir.mkdir()
    return str(checkpoint_dir / "test_checkpoint.json")


@pytest.fixture
def sample_pubchem_compound():
    """Fixture providing a sample PubChem compound response."""
    return {
        'CID': 2244,
        'props': [
            {
                'urn': {'label': 'IUPAC Name'},
                'value': {'sval': 'aspirin'}
            },
            {
                'urn': {'label': 'Molecular Formula'},
                'value': {'sval': 'C9H8O4'}
            },
            {
                'urn': {'label': 'Molecular Weight'},
                'value': {'fval': 180.16}
            }
        ],
        'atoms': {'aid': [1, 2, 3]},
        'bonds': {'aid1': [1, 2], 'aid2': [2, 3], 'order': [1, 1]}
    }


@pytest.fixture
def sample_chembl_compound():
    """Fixture providing a sample ChEMBL compound response."""
    return {
        'molecule_chembl_id': 'CHEMBL25',
        'pref_name': 'ASPIRIN',
        'molecule_structures': {
            'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
            'standard_inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        },
        'molecule_properties': {
            'full_molformula': 'C9H8O4',
            'full_mwt': 180.16,
            'alogp': 1.43,
            'hba': 4,
            'hbd': 1
        },
        'molecule_synonyms': [
            {'synonym': 'Aspirin'},
            {'synonym': 'Acetylsalicylic acid'}
        ]
    }


@pytest.fixture
def mock_aiohttp_session():
    """Fixture providing a mocked aiohttp session."""
    with patch('aiohttp.ClientSession') as mock_session:
        # Create a mock response
        mock_response = AsyncMock()
        mock_response.status = 200
        mock_response.content_type = 'application/json'
        mock_response.json = AsyncMock()
        mock_response.text = AsyncMock(return_value='')
        
        # Set up the mock session to return the mock response
        session_instance = mock_session.return_value
        session_instance.get.return_value.__aenter__.return_value = mock_response
        
        yield mock_session, mock_response


@pytest.fixture
def mock_db_operations():
    """Fixture providing a mocked database operations instance."""
    with patch('unified_importer.core.database.DatabaseOperations') as mock_db:
        db_instance = mock_db.return_value
        
        # Mock transaction context manager
        db_instance.transaction = MagicMock()
        db_instance.transaction.return_value.__enter__ = MagicMock()
        db_instance.transaction.return_value.__exit__ = MagicMock()
        
        # Mock async methods
        db_instance.insert_molecule = AsyncMock(return_value='test-molecule-id')
        db_instance.insert_properties = AsyncMock(return_value=['test-property-id'])
        db_instance.insert_synonyms = AsyncMock(return_value=['test-synonym-id'])
        db_instance.molecule_exists = AsyncMock(return_value=(False, None))
        db_instance.execute_query = AsyncMock(return_value=[])
        db_instance.batch_insert = AsyncMock(return_value=[])
        db_instance.with_retry = AsyncMock()
        
        yield db_instance


# Configuration tests

def test_config_initialization(sample_config, tmp_path):
    """Test initialization of ImporterConfig."""
    # Write config to file
    config_file = tmp_path / "test_config.json"
    with open(config_file, 'w') as f:
        json.dump(sample_config, f)
    
    # Initialize with file
    config = ImporterConfig(str(config_file))
    
    # Check values
    assert config.get('batch_size') == 10
    assert config.get('max_workers') == 4
    assert config.get('api_delay') == 0.1
    assert config.get('database', {}).get('type') == 'direct'


def test_config_env_override(sample_config, monkeypatch):
    """Test environment variable override for configuration."""
    # Set environment variable
    monkeypatch.setenv('CRYOPROTECT_IMPORT_BATCH_SIZE', '20')
    
    # Initialize config
    config = ImporterConfig()
    config.config = sample_config.copy()
    config._load_from_env()
    
    # Check override
    assert config.get('batch_size') == 20


def test_config_args_override(sample_config):
    """Test argument override for configuration."""
    # Initialize config
    config = ImporterConfig()
    config.config = sample_config.copy()
    
    # Update from args
    args = {'batch_size': 30, 'unknown_param': 'value'}
    config.update_from_args(args)
    
    # Check override
    assert config.get('batch_size') == 30
    assert 'unknown_param' not in config.config


# Checkpoint tests

@pytest.mark.asyncio
async def test_checkpoint_save_load(temp_checkpoint_file):
    """Test checkpoint save and load functionality."""
    # Create checkpoint manager
    manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Set data
    manager.set_total_items(100)
    manager.mark_processed('item1')
    manager.mark_processed('item2')
    manager.mark_failed('item3')
    manager.set_position('test-position')
    manager.set_custom_state('test-key', 'test-value')
    
    # Save
    manager.save()
    
    # Create new manager
    new_manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Check loaded data
    assert new_manager.get_processed_count() == 2
    assert new_manager.get_failed_count() == 1
    assert new_manager.is_processed('item1')
    assert new_manager.is_processed('item2')
    assert new_manager.is_failed('item3')
    assert new_manager.get_position() == 'test-position'
    assert new_manager.get_custom_state('test-key') == 'test-value'


@pytest.mark.asyncio
async def test_checkpoint_batch_tracking(temp_checkpoint_file):
    """Test checkpoint batch tracking functionality."""
    # Create checkpoint manager
    manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Update with a batch
    batch = ['item1', 'item2', 'item3']
    manager.update_batch(batch)
    
    # Mark some as processed
    manager.mark_processed('item1')
    manager.mark_processed('item2')
    
    # Check unprocessed
    unprocessed = manager.get_unprocessed_from_batch()
    assert len(unprocessed) == 1
    assert unprocessed[0] == 'item3'


# Progress tracking tests

@pytest.mark.asyncio
async def test_progress_tracker():
    """Test progress tracking functionality."""
    # Create tracker
    tracker = ProgressTracker(total_items=100)
    
    # Update progress
    tracker.update(processed=10, successful=8, failed=2)
    
    # Check counters
    assert tracker.processed_items == 10
    assert tracker.successful_items == 8
    assert tracker.failed_items == 2
    
    # Check percentage
    assert tracker.get_progress_percentage() == 10.0
    
    # Check rate calculation
    rate = tracker.get_processing_rate()
    assert rate > 0  # Should have a non-zero rate


@pytest.mark.asyncio
async def test_progress_callbacks():
    """Test progress tracker callbacks."""
    # Create tracker
    tracker = ProgressTracker(total_items=100)
    
    # Create mock callback
    callback = MagicMock()
    
    # Register callback
    tracker.register_callback(callback)
    
    # Update progress
    tracker.update(processed=10, successful=8, failed=2)
    
    # Check callback was called
    callback.assert_called_once()
    args = callback.call_args[0][0]
    assert args['processed_items'] == 10
    assert args['successful_items'] == 8
    assert args['failed_items'] == 2


# Validation tests

@pytest.mark.asyncio
async def test_molecule_validator():
    """Test molecule validation functionality."""
    # Create validator
    validator = MoleculeValidator()
    
    # Test valid molecule
    molecule = {
        'name': 'Test Molecule',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        'pubchem_cid': '2244'
    }
    is_valid, error = validator.validate_molecule(molecule)
    assert is_valid
    assert error is None
    
    # Test invalid molecule (missing required field)
    invalid_molecule = {
        'name': 'Test Molecule'
    }
    is_valid, error = validator.validate_molecule(invalid_molecule)
    assert not is_valid
    assert 'Missing required field' in error
    
    # Test invalid SMILES
    invalid_smiles = {
        'name': 'Test Molecule',
        'smiles': 'C(('  # Unbalanced parentheses
    }
    is_valid, error = validator.validate_molecule(invalid_smiles)
    assert not is_valid
    assert 'Invalid SMILES' in error


@pytest.mark.asyncio
async def test_property_validation():
    """Test property validation functionality."""
    # Create validator
    validator = MoleculeValidator()
    
    # Test valid property
    property_data = {
        'property_name': 'LogP',
        'property_type': 'physicochemical',
        'numeric_value': 1.43
    }
    is_valid, error = validator.validate_property(property_data)
    assert is_valid
    assert error is None
    
    # Test invalid property (missing value)
    invalid_property = {
        'property_name': 'LogP',
        'property_type': 'physicochemical'
    }
    is_valid, error = validator.validate_property(invalid_property)
    assert not is_valid
    assert 'value field must be set' in error


# Database operation tests

@pytest.mark.asyncio
async def test_database_molecule_operations(mock_db_operations):
    """Test database molecule operations."""
    # Create molecule data
    molecule_data = {
        'name': 'Test Molecule',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
        'pubchem_cid': '2244'
    }
    
    # Insert molecule
    molecule_id = await mock_db_operations.insert_molecule(molecule_data)
    assert molecule_id == 'test-molecule-id'
    mock_db_operations.insert_molecule.assert_called_once()
    
    # Check if molecule exists
    exists, existing_id = await mock_db_operations.molecule_exists('2244', 'pubchem_cid')
    assert exists is False
    assert existing_id is None
    mock_db_operations.molecule_exists.assert_called_once()


@pytest.mark.asyncio
async def test_database_property_operations(mock_db_operations):
    """Test database property operations."""
    # Create property data
    properties = [
        {
            'property_name': 'LogP',
            'property_type': 'physicochemical',
            'numeric_value': 1.43
        },
        {
            'property_name': 'Molecular Weight',
            'property_type': 'physicochemical',
            'numeric_value': 180.16
        }
    ]
    
    # Insert properties
    property_ids = await mock_db_operations.insert_properties('test-molecule-id', properties)
    assert property_ids == ['test-property-id']
    mock_db_operations.insert_properties.assert_called_once()


# PubChem source tests

@pytest.mark.asyncio
async def test_pubchem_source_initialization(mock_db_operations):
    """Test PubChem data source initialization."""
    # Create source
    source = PubChemDataSource(
        db_operations=mock_db_operations,
        config={
            'batch_size': 10,
            'api_delay': 0.1,
            'max_retries': 2,
            'retry_delay': 0.1
        }
    )
    
    # Check initialization
    assert source.get_source_name() == 'pubchem'
    assert source.batch_size == 10
    assert source.api_delay == 0.1
    assert source.max_retries == 2
    assert source.retry_delay == 0.1


@pytest.mark.asyncio
async def test_pubchem_fetch_compound(mock_db_operations, mock_aiohttp_session, sample_pubchem_compound):
    """Test PubChem fetch compound functionality."""
    mock_session, mock_response = mock_aiohttp_session
    
    # Set up response data
    mock_response.json.return_value = {
        'PC_Compounds': [sample_pubchem_compound]
    }
    
    # Create source
    source = PubChemDataSource(db_operations=mock_db_operations)
    
    # Patch the _get_session method
    with patch.object(source, '_get_session', return_value=mock_session.return_value):
        # Fetch compound
        compound = await source.fetch_compound('2244')
        
        # Check result
        assert compound is not None
        assert compound['CID'] == 2244
        assert 'props' in compound


@pytest.mark.asyncio
async def test_pubchem_transform_compound(mock_db_operations, sample_pubchem_compound):
    """Test PubChem compound transformation."""
    # Create source
    source = PubChemDataSource(db_operations=mock_db_operations)
    
    # Patch methods that make API calls
    with patch.object(source, '_add_structure_data', new=AsyncMock()), \
         patch.object(source, '_add_synonyms', new=AsyncMock()):
        
        # Transform compound
        transformed = await source.transform_compound_data(sample_pubchem_compound)
        
        # Check transformation
        assert transformed['pubchem_cid'] == '2244'
        assert transformed['data_source'] == 'PubChem'
        assert 'name' in transformed


# ChEMBL source tests

@pytest.mark.asyncio
async def test_chembl_source_initialization(mock_db_operations):
    """Test ChEMBL data source initialization."""
    # Create source
    source = ChEMBLDataSource(
        db_operations=mock_db_operations,
        config={
            'batch_size': 10,
            'api_delay': 0.1,
            'max_retries': 2,
            'retry_delay': 0.1
        }
    )
    
    # Check initialization
    assert source.get_source_name() == 'chembl'
    assert source.batch_size == 10
    assert source.api_delay == 0.1
    assert source.max_retries == 2
    assert source.retry_delay == 0.1


@pytest.mark.asyncio
async def test_chembl_fetch_compound(mock_db_operations, mock_aiohttp_session, sample_chembl_compound):
    """Test ChEMBL fetch compound functionality."""
    mock_session, mock_response = mock_aiohttp_session
    
    # Set up response data
    mock_response.json.return_value = sample_chembl_compound
    
    # Create source
    source = ChEMBLDataSource(db_operations=mock_db_operations)
    
    # Patch the _get_session method
    with patch.object(source, '_get_session', return_value=mock_session.return_value):
        # Fetch compound
        compound = await source.fetch_compound('CHEMBL25')
        
        # Check result
        assert compound is not None
        assert compound['molecule_chembl_id'] == 'CHEMBL25'
        assert 'molecule_structures' in compound


@pytest.mark.asyncio
async def test_chembl_transform_compound(mock_db_operations, sample_chembl_compound):
    """Test ChEMBL compound transformation."""
    # Create source
    source = ChEMBLDataSource(db_operations=mock_db_operations)
    
    # Transform compound
    transformed = await source.transform_compound_data(sample_chembl_compound)
    
    # Check transformation
    assert transformed['chembl_id'] == 'CHEMBL25'
    assert transformed['data_source'] == 'ChEMBL'
    assert transformed['name'] == 'ASPIRIN'
    assert transformed['smiles'] == 'CC(=O)OC1=CC=CC=C1C(=O)O'
    assert transformed['formula'] == 'C9H8O4'
    assert transformed['molecular_weight'] == 180.16
    assert len(transformed['synonyms']) == 2


# Integration tests

@pytest.mark.asyncio
async def test_import_compound_workflow(mock_db_operations, sample_pubchem_compound, temp_checkpoint_file):
    """Test the end-to-end compound import workflow."""
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Create progress tracker
    progress_tracker = ProgressTracker()
    
    # Create validator
    validator = MoleculeValidator()
    
    # Create PubChem source with mocks
    source = PubChemDataSource(
        db_operations=mock_db_operations,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker,
        validator=validator,
        config={'dry_run': False}
    )
    
    # Mock fetch_compound
    with patch.object(source, 'fetch_compound', new=AsyncMock(return_value=sample_pubchem_compound)), \
         patch.object(source, 'transform_compound_data', new=AsyncMock(return_value={
             'name': 'Aspirin',
             'pubchem_cid': '2244',
             'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
             'inchikey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
             'data_source': 'PubChem'
         })), \
         patch.object(source, 'get_property_data', new=AsyncMock(return_value=[
             {
                 'property_name': 'LogP',
                 'property_type': 'physicochemical',
                 'numeric_value': 1.43
             }
         ])), \
         patch.object(source, 'transform_property_data', new=AsyncMock(return_value=[
             {
                 'property_name': 'LogP',
                 'property_type': 'physicochemical',
                 'numeric_value': 1.43,
                 'source': 'PubChem'
             }
         ])):
        
        # Import compound
        success, error, result = await source.import_compound('2244')
        
        # Check result
        assert success is True
        assert error is None
        assert result is not None
        assert result['molecule_id'] == 'test-molecule-id'
        
        # Check checkpoint
        assert checkpoint_manager.is_processed('2244')
        
        # Check progress
        assert progress_tracker.successful_items == 1
        assert progress_tracker.processed_items == 1


@pytest.mark.asyncio
async def test_error_handling_workflow(mock_db_operations, temp_checkpoint_file):
    """Test error handling in the import workflow."""
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Create progress tracker
    progress_tracker = ProgressTracker()
    
    # Create PubChem source with mocks
    source = PubChemDataSource(
        db_operations=mock_db_operations,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker
    )
    
    # Mock fetch_compound to fail
    with patch.object(source, 'fetch_compound', new=AsyncMock(side_effect=Exception("API error"))):
        # Import compound
        success, error, result = await source.import_compound('2244')
        
        # Check result
        assert success is False
        assert error is not None
        assert "Error importing compound" in error
        assert result is None
        
        # Check checkpoint
        assert checkpoint_manager.is_failed('2244')
        
        # Check progress
        assert progress_tracker.failed_items == 1
        assert progress_tracker.processed_items == 0


@pytest.mark.asyncio
async def test_batch_import_workflow(mock_db_operations, temp_checkpoint_file):
    """Test batch import workflow."""
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(checkpoint_file=temp_checkpoint_file)
    
    # Create progress tracker
    progress_tracker = ProgressTracker(total_items=3)
    
    # Create PubChem source with mocks
    source = PubChemDataSource(
        db_operations=mock_db_operations,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker,
        config={'batch_size': 2, 'api_delay': 0.0}
    )
    
    # Mock import_compound - first two succeed, third fails
    mock_import = AsyncMock()
    mock_import.side_effect = [
        (True, None, {'molecule_id': 'id1', 'identifier': 'cid1'}),
        (True, None, {'molecule_id': 'id2', 'identifier': 'cid2'}),
        (False, "Error with cid3", None)
    ]
    
    with patch.object(source, 'import_compound', new=mock_import):
        # Import compounds
        success, failure, failures = await source.import_compounds(['cid1', 'cid2', 'cid3'])
        
        # Check results
        assert success == 2
        assert failure == 1
        assert len(failures) == 1
        assert failures[0][0] == 'cid3'
        
        # Check progress
        assert progress_tracker.successful_items == 2
        assert progress_tracker.failed_items == 1


# Main method tests

@pytest.mark.asyncio
async def test_import_mode_command_line(monkeypatch):
    """Test command line argument parsing for import mode."""
    from unified_importer.main import run_import
    
    # Mock necessary objects
    mock_logger = MagicMock()
    mock_config = MagicMock()
    mock_config.get.return_value = {}
    mock_config.as_dict.return_value = {}
    
    # Create args object
    class Args:
        source = 'pubchem'
        query = 'aspirin'
        limit = 10
        checkpoint = None
        identifiers = None

    args = Args()
    
    # Mock load_data_source
    with patch('unified_importer.main.load_data_source') as mock_load_source, \
         patch('unified_importer.main.DatabaseOperations') as mock_db_class, \
         patch('unified_importer.main.CheckpointManager') as mock_checkpoint_class, \
         patch('unified_importer.main.ProgressTracker') as mock_tracker_class, \
         patch('unified_importer.main.ConsoleProgressReporter') as mock_reporter_class, \
         patch('unified_importer.main.setup_logging', return_value=mock_logger), \
         patch('os.makedirs'):
        
        # Mock the data source
        mock_source = MagicMock()
        mock_source.get_compound_count = AsyncMock(return_value=100)
        mock_source.search_and_import = AsyncMock(return_value=(5, 0, []))
        mock_load_source.return_value = mock_source
        
        # Run import
        exit_code = await run_import(args, mock_config, mock_logger)
        
        # Check result
        assert exit_code == 0
        mock_source.get_compound_count.assert_called_once_with('aspirin')
        mock_source.search_and_import.assert_called_once_with(query='aspirin', max_results=10)


# Run the tests
if __name__ == "__main__":
    pytest.main(["-v", __file__])
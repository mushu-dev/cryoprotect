"""
Unit tests for database utility functions.

This module contains tests for the database utility functions in database/utils.py.
"""

import pytest
import time
from unittest.mock import patch, MagicMock, call

from database.utils import (
    get_db,
    with_connection,
    with_retry,
    with_transaction,
    execute_query,
    execute_batch,
    get_molecule_by_id,
    get_molecule_properties,
    get_molecules_by_inchikey,
    insert_molecule,
    update_molecule,
    set_molecule_property,
    get_or_create_property_type,
    test_database_connection
)
from database.connection_manager import ConnectionManager

# ============================================================================
# Test get_db function
# ============================================================================

def test_get_db_returns_singleton_instance():
    """Test that get_db returns the singleton instance of ConnectionManager."""
    # First call should create the instance
    db1 = get_db()
    assert isinstance(db1, ConnectionManager)
    
    # Second call should return the same instance
    db2 = get_db()
    assert db1 is db2

# ============================================================================
# Test decorators
# ============================================================================

def test_with_connection_decorator():
    """Test that with_connection decorator ensures database connection."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = None
    mock_connection_manager.connect.return_value = True
    
    # Create a test function with the decorator
    @with_connection
    def test_func():
        return "success"
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = test_func()
        
        # Verify the decorator checked for an active adapter
        mock_connection_manager.get_active_adapter.assert_called_once()
        
        # Verify the decorator called connect
        mock_connection_manager.connect.assert_called_once()
        
        # Verify the function was called and returned the expected result
        assert result == "success"

def test_with_connection_decorator_raises_error_when_connection_fails():
    """Test that with_connection decorator raises error when connection fails."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = None
    mock_connection_manager.connect.return_value = False
    
    # Create a test function with the decorator
    @with_connection
    def test_func():
        return "success"
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        # Verify the decorator raises ConnectionError when connect fails
        with pytest.raises(ConnectionError):
            test_func()

def test_with_retry_decorator():
    """Test that with_retry decorator retries failed operations."""
    # Create a mock function that fails twice then succeeds
    mock_func = MagicMock(side_effect=[ValueError("First failure"), 
                                       ValueError("Second failure"), 
                                       "success"])
    
    # Apply the decorator
    decorated_func = with_retry(max_retries=3, backoff=0.01)(mock_func)
    
    # Patch time.sleep to avoid waiting during tests
    with patch('time.sleep'):
        result = decorated_func()
        
        # Verify the function was called 3 times
        assert mock_func.call_count == 3
        
        # Verify the final result
        assert result == "success"

def test_with_retry_decorator_max_retries_exceeded():
    """Test that with_retry decorator raises the last exception when max retries exceeded."""
    # Create a mock function that always fails
    mock_func = MagicMock(side_effect=ValueError("Persistent failure"))
    
    # Apply the decorator
    decorated_func = with_retry(max_retries=2, backoff=0.01)(mock_func)
    
    # Patch time.sleep to avoid waiting during tests
    with patch('time.sleep'):
        # Verify the decorator raises the last exception when retries are exhausted
        with pytest.raises(ValueError, match="Persistent failure"):
            decorated_func()
        
        # Verify the function was called the expected number of times
        assert mock_func.call_count == 3  # Initial call + 2 retries

def test_with_transaction_decorator():
    """Test that with_transaction decorator manages transactions correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_transaction = MagicMock()
    mock_connection_manager.begin_transaction.return_value = mock_transaction
    
    # Create a test function with the decorator
    @with_transaction
    def test_func(transaction=None):
        assert transaction is mock_transaction
        return "success"
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = test_func()
        
        # Verify transaction was begun
        mock_connection_manager.begin_transaction.assert_called_once()
        
        # Verify transaction was committed
        mock_connection_manager.commit_transaction.assert_called_once_with(mock_transaction)
        
        # Verify the function returned the expected result
        assert result == "success"

def test_with_transaction_decorator_rolls_back_on_exception():
    """Test that with_transaction decorator rolls back transaction on exception."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_transaction = MagicMock()
    mock_connection_manager.begin_transaction.return_value = mock_transaction
    
    # Create a test function with the decorator that raises an exception
    @with_transaction
    def test_func(transaction=None):
        raise ValueError("Test exception")
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        # Verify the decorator re-raises the exception
        with pytest.raises(ValueError, match="Test exception"):
            test_func()
        
        # Verify transaction was begun
        mock_connection_manager.begin_transaction.assert_called_once()
        
        # Verify transaction was rolled back
        mock_connection_manager.rollback_transaction.assert_called_once_with(mock_transaction)
        
        # Verify transaction was not committed
        mock_connection_manager.commit_transaction.assert_not_called()

# ============================================================================
# Test helper functions
# ============================================================================

def test_execute_query():
    """Test that execute_query calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [{"id": 1, "name": "Test"}]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = execute_query("SELECT * FROM test", {"param": "value"})
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once_with(
            "SELECT * FROM test", {"param": "value"}
        )
        
        # Verify the result
        assert result == [{"id": 1, "name": "Test"}]

def test_execute_batch():
    """Test that execute_batch calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_batch.return_value = [
        [{"id": 1}], [{"id": 2}]
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = execute_batch(["SELECT 1", "SELECT 2"])
        
        # Verify execute_batch was called with the correct parameters
        mock_connection_manager.execute_batch.assert_called_once_with(
            ["SELECT 1", "SELECT 2"]
        )
        
        # Verify the result
        assert result == [[{"id": 1}], [{"id": 2}]]

def test_get_molecule_by_id():
    """Test that get_molecule_by_id calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "mol-1", "name": "Test Molecule"}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_molecule_by_id("mol-1")
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once_with(
            "SELECT * FROM molecules WHERE id = %s", ("mol-1",)
        )
        
        # Verify the result
        assert result == {"id": "mol-1", "name": "Test Molecule"}

def test_get_molecule_by_id_not_found():
    """Test that get_molecule_by_id returns None when molecule not found."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = []
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_molecule_by_id("non-existent")
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once_with(
            "SELECT * FROM molecules WHERE id = %s", ("non-existent",)
        )
        
        # Verify the result is None
        assert result is None

def test_get_molecule_properties():
    """Test that get_molecule_properties calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "prop-1", "property_name": "LogP", "value": 1.23}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_molecule_properties("mol-1")
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once()
        
        # Verify the result
        assert result == [{"id": "prop-1", "property_name": "LogP", "value": 1.23}]

def test_get_molecules_by_inchikey():
    """Test that get_molecules_by_inchikey calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "mol-1", "name": "Test Molecule", "inchi_key": "TEST-INCHI-KEY"}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_molecules_by_inchikey("TEST-INCHI-KEY")
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once_with(
            "SELECT * FROM molecules WHERE inchi_key = %s", ("TEST-INCHI-KEY",)
        )
        
        # Verify the result
        assert result == [{"id": "mol-1", "name": "Test Molecule", "inchi_key": "TEST-INCHI-KEY"}]

def test_insert_molecule():
    """Test that insert_molecule calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "mol-new", "name": "New Molecule"}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = insert_molecule(
            name="New Molecule",
            formula="C6H12O6",
            molecular_weight=180.156,
            smiles="C(C1C(C(C(C(O1)O)O)O)O)O",
            inchi="InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1",
            inchi_key="WQZGKKKJIJFFOK-GASJEMHNSA-N",
            chembl_id="CHEMBL1",
            pubchem_cid="5793",
            data_source="test"
        )
        
        # Verify execute_query was called with the correct parameters
        mock_connection_manager.execute_query.assert_called_once()
        
        # Verify the result
        assert result == {"id": "mol-new", "name": "New Molecule"}

def test_update_molecule():
    """Test that update_molecule calls the database adapter correctly."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "mol-1", "name": "Updated Molecule"}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = update_molecule(
            molecule_id="mol-1",
            name="Updated Molecule"
        )
        
        # Verify execute_query was called
        mock_connection_manager.execute_query.assert_called_once()
        
        # Verify the result
        assert result == {"id": "mol-1", "name": "Updated Molecule"}

def test_update_molecule_no_changes():
    """Test that update_molecule returns current molecule when no changes specified."""
    # Mock ConnectionManager and get_molecule_by_id
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    
    # Patch get_db and get_molecule_by_id
    with patch('database.utils.get_db', return_value=mock_connection_manager), \
         patch('database.utils.get_molecule_by_id') as mock_get_molecule:
        mock_get_molecule.return_value = {"id": "mol-1", "name": "Unchanged Molecule"}
        
        result = update_molecule(molecule_id="mol-1")
        
        # Verify execute_query was not called
        mock_connection_manager.execute_query.assert_not_called()
        
        # Verify get_molecule_by_id was called
        mock_get_molecule.assert_called_once_with("mol-1")
        
        # Verify the result
        assert result == {"id": "mol-1", "name": "Unchanged Molecule"}

def test_set_molecule_property_update_existing():
    """Test that set_molecule_property updates existing property."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    
    # First query returns existing property
    # Second query returns updated property
    mock_connection_manager.execute_query.side_effect = [
        [{"id": "prop-1", "molecule_id": "mol-1", "property_type_id": "type-1", "value": 1.0}],
        [{"id": "prop-1", "molecule_id": "mol-1", "property_type_id": "type-1", "value": 2.0}]
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = set_molecule_property(
            molecule_id="mol-1",
            property_type_id="type-1",
            value=2.0,
            source="test",
            confidence=0.95
        )
        
        # Verify execute_query was called twice
        assert mock_connection_manager.execute_query.call_count == 2
        
        # Verify the result
        assert result == {"id": "prop-1", "molecule_id": "mol-1", "property_type_id": "type-1", "value": 2.0}

def test_set_molecule_property_insert_new():
    """Test that set_molecule_property inserts new property."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    
    # First query returns no existing property
    # Second query returns new property
    mock_connection_manager.execute_query.side_effect = [
        [],
        [{"id": "prop-new", "molecule_id": "mol-1", "property_type_id": "type-1", "value": 2.0}]
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = set_molecule_property(
            molecule_id="mol-1",
            property_type_id="type-1",
            value=2.0,
            source="test",
            confidence=0.95
        )
        
        # Verify execute_query was called twice
        assert mock_connection_manager.execute_query.call_count == 2
        
        # Verify the result
        assert result == {"id": "prop-new", "molecule_id": "mol-1", "property_type_id": "type-1", "value": 2.0}

def test_get_or_create_property_type_existing():
    """Test that get_or_create_property_type returns existing property type."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.execute_query.return_value = [
        {"id": "type-1", "name": "LogP", "data_type": "numeric", "unit": ""}
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_or_create_property_type(
            name="LogP",
            description="Octanol-water partition coefficient",
            data_type="numeric"
        )
        
        # Verify execute_query was called once
        mock_connection_manager.execute_query.assert_called_once()
        
        # Verify the result
        assert result == {"id": "type-1", "name": "LogP", "data_type": "numeric", "unit": ""}

def test_get_or_create_property_type_new():
    """Test that get_or_create_property_type creates new property type."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    
    # First query returns no existing property type
    # Second query returns new property type
    mock_connection_manager.execute_query.side_effect = [
        [],
        [{"id": "type-new", "name": "NewProperty", "data_type": "text", "unit": ""}]
    ]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = get_or_create_property_type(
            name="NewProperty",
            description="A new property type",
            data_type="text"
        )
        
        # Verify execute_query was called twice
        assert mock_connection_manager.execute_query.call_count == 2
        
        # Verify the result
        assert result == {"id": "type-new", "name": "NewProperty", "data_type": "text", "unit": ""}

def test_test_database_connection():
    """Test that test_database_connection returns connection status."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.get_connection_info.return_value = {"adapter": "test"}
    mock_connection_manager.test_all_connections.return_value = {"test": (True, "Connected")}
    mock_connection_manager.execute_query.return_value = [{"test": 1}]
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = test_database_connection()
        
        # Verify the methods were called
        mock_connection_manager.get_connection_info.assert_called_once()
        mock_connection_manager.test_all_connections.assert_called_once()
        mock_connection_manager.execute_query.assert_called_once_with("SELECT 1 as test")
        
        # Verify the result structure
        assert "connection_info" in result
        assert "test_results" in result
        assert "query_test" in result
        assert result["query_test"]["result"] == "Success"
        assert result["query_test"]["error"] is None

def test_test_database_connection_query_error():
    """Test that test_database_connection handles query errors."""
    # Mock ConnectionManager
    mock_connection_manager = MagicMock()
    mock_connection_manager.get_active_adapter.return_value = True
    mock_connection_manager.get_connection_info.return_value = {"adapter": "test"}
    mock_connection_manager.test_all_connections.return_value = {"test": (True, "Connected")}
    mock_connection_manager.execute_query.side_effect = Exception("Query error")
    
    # Patch get_db to return our mock
    with patch('database.utils.get_db', return_value=mock_connection_manager):
        result = test_database_connection()
        
        # Verify the methods were called
        mock_connection_manager.get_connection_info.assert_called_once()
        mock_connection_manager.test_all_connections.assert_called_once()
        mock_connection_manager.execute_query.assert_called_once_with("SELECT 1 as test")
        
        # Verify the result structure
        assert "connection_info" in result
        assert "test_results" in result
        assert "query_test" in result
        assert result["query_test"]["result"] is None
        assert result["query_test"]["error"] == "Query error"
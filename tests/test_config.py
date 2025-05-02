"""
Test script for configuration and credential verification.

This script tests the configuration verification functionality implemented in
ChEMBL_CryoProtectants_Supabase.py, ensuring that:
1. All required configuration variables are validated
2. API connectivity tests for Supabase and ChEMBL are performed
3. Failures are properly logged and abort execution
"""

import os
import sys
import logging
from unittest import mock
import pytest

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Configure logging for tests
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Import the module to test
import ChEMBL_CryoProtectants_Supabase as chembl_supabase
from config import active_config, ConfigurationError


def test_verify_configuration_success():
    """Test that configuration verification succeeds with valid configuration."""
    # This test assumes that the environment is properly configured
    try:
        result = chembl_supabase.verify_configuration()
        assert result is True, "Configuration verification should return True on success"
        logger.info("Configuration verification test passed")
    except Exception as e:
        pytest.fail(f"Configuration verification failed: {str(e)}")


@mock.patch('ChEMBL_CryoProtectants_Supabase.create_client')
def test_supabase_connection_failure(mock_create_client):
    """Test that Supabase connection failure is properly handled."""
    # Mock the Supabase client to raise an exception
    mock_create_client.side_effect = Exception("Mocked Supabase connection failure")
    
    # Verify that the function exits with an error
    with pytest.raises(SystemExit) as excinfo:
        chembl_supabase.verify_configuration()
    
    # Check that the exit code is non-zero
    assert excinfo.value.code != 0, "Should exit with non-zero code on Supabase connection failure"
    logger.info("Supabase connection failure test passed")


@mock.patch('ChEMBL_CryoProtectants_Supabase.create_client')
def test_supabase_query_failure(mock_create_client):
    """Test that Supabase query failure is properly handled."""
    # Mock the Supabase client to return a client that raises an exception on execute
    mock_supabase = mock.MagicMock()
    mock_table = mock.MagicMock()
    mock_select = mock.MagicMock()
    mock_limit = mock.MagicMock()
    mock_execute = mock.MagicMock()
    
    # Set up the mock chain
    mock_create_client.return_value = mock_supabase
    mock_supabase.table.return_value = mock_table
    mock_table.select.return_value = mock_select
    mock_select.limit.return_value = mock_limit
    
    # Make the execute method raise an exception
    mock_response = mock.MagicMock()
    mock_response.error = "Mocked Supabase query error"
    mock_limit.execute.return_value = mock_response
    
    # Verify that the function exits with an error
    with pytest.raises(SystemExit) as excinfo:
        chembl_supabase.verify_configuration()
    
    # Check that the exit code is non-zero
    assert excinfo.value.code != 0, "Should exit with non-zero code on Supabase query failure"
    logger.info("Supabase query failure test passed")


@mock.patch('chembl.client.ResilientChEMBLClient.search_molecules')
def test_chembl_api_failure(mock_search_molecules):
    """Test that ChEMBL API failure is properly handled."""
    # Mock the ChEMBL client to return an error
    mock_search_molecules.return_value = {"Error": "Mocked ChEMBL API error"}
    
    # Patch the Supabase client to avoid actual API calls
    with mock.patch('ChEMBL_CryoProtectants_Supabase.create_client'):
        with mock.patch('ChEMBL_CryoProtectants_Supabase.supabase'):
            # Verify that the function exits with an error
            with pytest.raises(SystemExit) as excinfo:
                chembl_supabase.verify_configuration()
    
    # Check that the exit code is non-zero
    assert excinfo.value.code != 0, "Should exit with non-zero code on ChEMBL API failure"
    logger.info("ChEMBL API failure test passed")


@mock.patch('config.validate_config')
def test_config_validation_failure(mock_validate_config):
    """Test that configuration validation failure is properly handled."""
    # Mock the validate_config function to raise a ConfigurationError
    mock_validate_config.side_effect = ConfigurationError("Mocked configuration validation error")
    
    # Verify that the function exits with an error
    with pytest.raises(SystemExit) as excinfo:
        chembl_supabase.verify_configuration()
    
    # Check that the exit code is non-zero
    assert excinfo.value.code != 0, "Should exit with non-zero code on configuration validation failure"
    logger.info("Configuration validation failure test passed")


if __name__ == "__main__":
    logger.info("Running configuration verification tests...")
    
    # Run the tests
    try:
        test_verify_configuration_success()
        logger.info("All tests passed!")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Test failed: {str(e)}")
        sys.exit(1)
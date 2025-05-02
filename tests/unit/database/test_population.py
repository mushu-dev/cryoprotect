"""
Unit tests for database population functions.

This module contains tests that verify the correct population of database tables
using the population scripts. All tests are designed to be isolated and repeatable.

Test Structure:
- Each test should set up its own data and clean up after itself.
- Use fixtures for reusable setup/teardown logic.
- Patch internal dependencies where appropriate.
"""

import pytest
from unittest.mock import patch
from database.main import populate_database

@pytest.fixture
def sample_tables():
    """Fixture providing sample table names for population tests."""
    return ["table1", "table2"]

@pytest.fixture
def sample_config():
    """Fixture providing a sample config dictionary."""
    return {"option": "value"}

def test_population_inserts_data_correctly(sample_tables, sample_config):
    """
    Test that the population function inserts data into the database correctly.

    This test patches the internal function and verifies that the correct calls are made.
    """
    with patch("database.main.populate_specific") as mock_populate_specific:
        mock_populate_specific.return_value = {"table1": 2, "table2": 3}
        result = populate_database(sample_tables, "development", sample_config)
        mock_populate_specific.assert_called_once_with(sample_tables, "development", sample_config)
        assert result == {"table1": 2, "table2": 3}

def test_population_handles_empty_tables(sample_config):
    """
    Test that the population function handles None as tables and calls populate_all.

    The function should call populate_all if tables is None.
    """
    with patch("database.main.populate_all") as mock_populate_all:
        mock_populate_all.return_value = {"table1": 5, "table2": 7}
        result = populate_database(None, "staging", sample_config)
        mock_populate_all.assert_called_once_with("staging", sample_config)
        assert result == {"table1": 5, "table2": 7}

def test_population_handles_empty_config(sample_tables):
    """
    Test that the population function works when config is not provided.
    """
    with patch("database.main.populate_specific") as mock_populate_specific:
        mock_populate_specific.return_value = {"table1": 1, "table2": 1}
        result = populate_database(sample_tables, "production")
        mock_populate_specific.assert_called_once_with(sample_tables, "production", None)
        assert result == {"table1": 1, "table2": 1}
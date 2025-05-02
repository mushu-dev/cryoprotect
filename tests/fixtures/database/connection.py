"""
Database connection fixtures.

This module provides fixtures for setting up database
connections in tests.
"""

import os
import pytest
from typing import Dict, Any, Generator, List

from database.utils.connection import create_connection
from tests.fixtures.mocks.supabase import MockSupabase, patch_supabase_client

@pytest.fixture
def db_config() -> Dict[str, Any]:
    """
    Provide a database configuration for testing.
    
    Returns:
        Dictionary with test database configuration
    """
    return {
        'url': 'https://test-instance.supabase.co',
        'key': 'test-key',
        'service_role': 'test-service-role-key',
    }

@pytest.fixture
def mock_db() -> Generator[MockSupabase, None, None]:
    """
    Provide a mock database connection.
    
    Yields:
        MockSupabase instance
    """
    mock_client = MockSupabase()
    
    # Setup test data
    mock_client.add_test_data('molecules', [
        {'id': 'mol-1', 'name': 'Test Molecule 1', 'formula': 'C6H12O6'},
        {'id': 'mol-2', 'name': 'Test Molecule 2', 'formula': 'H2O'},
    ])
    
    mock_client.add_test_data('mixtures', [
        {'id': 'mix-1', 'name': 'Test Mixture 1'},
        {'id': 'mix-2', 'name': 'Test Mixture 2'},
    ])
    
    with patch_supabase_client(mock_client):
        yield mock_client
        
@pytest.fixture
def real_db_connection(db_config) -> Generator[Any, None, None]:
    """
    Provide a real database connection for integration tests.
    
    This fixture should only be used in integration tests.
    
    Yields:
        Database connection
    """
    # Only use this fixture when explicitly requested with a marker
    if not pytest.config.getoption('--use-real-db'):
        pytest.skip('Skipping test that requires a real database connection')
        
    conn = create_connection(db_config)
    yield conn

@pytest.fixture
def populated_mock_db(mock_db) -> Generator[MockSupabase, None, None]:
    """
    Provide a mock database populated with comprehensive test data.
    
    Yields:
        MockSupabase instance with test data
    """
    # Load and insert test data from JSON files
    from tests.fixtures.data import load_test_data
    
    molecules_data = load_test_data('molecules.json')
    mixtures_data = load_test_data('mixtures.json')
    
    mock_db.add_test_data('molecules', molecules_data)
    mock_db.add_test_data('mixtures', mixtures_data)
    
    # Add relationships
    for i, mixture in enumerate(mixtures_data):
        mix_id = mixture['id']
        # Add components for each mixture
        for j in range(1, 3):  # 2 components per mixture
            mol_id = molecules_data[j % len(molecules_data)]['id']
            component = {
                'id': f'comp-{i}-{j}',
                'mixture_id': mix_id,
                'molecule_id': mol_id,
                'concentration': 10.0 * j,
                'units': '%'
            }
            mock_db.add_test_data('mixture_components', [component])
    
    yield mock_db
"""
CryoProtect v2 Test Configuration

This module provides pytest configuration and fixtures for the CryoProtect v2 testing framework.
It integrates all fixture components: Database Fixtures, Mock Objects, API Fixtures, and Test Data Fixtures.
"""

import pytest
import tempfile
import shutil
import os
import logging
from typing import Dict, List, Any, Generator, Optional

# Configure logging for tests
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# ============================================================================
# Core Test Environment Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def temp_dir():
    """
    Create a temporary directory for the test session and clean up after.
    
    This fixture creates a temporary directory that persists for the entire test session
    and is automatically cleaned up when the session ends.
    
    Returns:
        str: Path to the temporary directory
    """
    dirpath = tempfile.mkdtemp()
    yield dirpath
    shutil.rmtree(dirpath)

@pytest.fixture(autouse=True)
def set_test_env(monkeypatch):
    """
    Set environment variables for tests.
    
    This fixture automatically runs for all tests and sets up the necessary
    environment variables for testing.
    
    Args:
        monkeypatch: pytest's monkeypatch fixture for modifying environment
    """
    monkeypatch.setenv("ENV", "test")
    monkeypatch.setenv("FLASK_ENV", "testing")
    monkeypatch.setenv("TESTING", "true")
    monkeypatch.setenv("SUPABASE_URL", "https://test-instance.supabase.co")
    monkeypatch.setenv("SUPABASE_KEY", "test-key")
    monkeypatch.setenv("SUPABASE_SERVICE_ROLE", "test-service-role-key")

# ============================================================================
# Database Fixtures
# ============================================================================

# Import database connection fixtures
from tests.fixtures.database.connection import (
    db_config,
    mock_db,
    real_db_connection,
    populated_mock_db
)

# Import database migration fixtures
from tests.fixtures.database.migrations import (
    migration_db,
    temp_migration_script
)

# Create a fixture for applying migrations
@pytest.fixture
def apply_migration(migration_db, temp_migration_script):
    """
    Apply a test migration to the database.
    
    Args:
        migration_db: Mock database prepared for migrations
        temp_migration_script: Path to a temporary migration script
        
    Returns:
        Result of applying the migration
    """
    from database.migrations import apply_migrations
    return apply_migrations(migration_db, [temp_migration_script], 'test')

# Import database data fixtures
from tests.fixtures.database.data import (
    load_test_data,
    generate_molecule as db_generate_molecule,
    generate_mixture as db_generate_mixture
)

# ============================================================================
# API Fixtures
# ============================================================================

# Import API client fixtures
from tests.fixtures.api.client import (
    app,
    client,
    api_client,
    authenticated_client,
    admin_client,
    scientist_client,
    mock_api_client
)

# Import authentication fixtures
from tests.fixtures.api.auth import (
    auth_token,
    admin_token,
    user_token,
    scientist_token,
    expired_token,
    invalid_token,
    mock_auth_middleware,
    mock_current_user,
    mock_admin_user,
    mock_regular_user,
    mock_scientist_user,
    mock_unauthenticated_user
)

# Import CSRF fixtures
from tests.fixtures.api.csrf import (
    csrf_token,
    csrf_headers,
    csrf_client,
    authenticated_csrf_client,
    disable_csrf
)

# ============================================================================
# Mock Object Fixtures
# ============================================================================

# Import RDKit mock fixtures
from tests.fixtures.mocks.rdkit import (
    mock_rdkit,
    molecule_factory,
    patch_rdkit
)

# Import Supabase mock fixtures
from tests.fixtures.mocks.supabase import (
    MockSupabase,
    patch_supabase_client
)

@pytest.fixture
def mock_supabase():
    """
    Provide a mock Supabase client.
    
    Returns:
        MockSupabase instance
    """
    return MockSupabase()

# ============================================================================
# Test Data Fixtures
# ============================================================================

# Import data loaders
from tests.fixtures.data.loaders import (
    load_molecules,
    load_property_types,
    load_mixtures,
    load_experiments,
    load_calculation_methods,
    load_teams,
    load_projects,
    load_related_data
)

# Import data generators
from tests.fixtures.data.generators import (
    generate_molecule,
    generate_mixture,
    generate_experiment,
    generate_prediction,
    generate_user_profile,
    generate_team,
    generate_project
)

# ============================================================================
# Test Discovery and Configuration
# Add fixtures for initializing and clearing test data
@pytest.fixture
def init_test_data(mock_db):
    """
    Initialize test database with standard test data.
    
    Args:
        mock_db: Mock database fixture
        
    Returns:
        The mock database with initialized test data
    """
    # Load test data from JSON files
    try:
        molecules = load_test_data('molecules.json')
        mixtures = load_test_data('mixtures.json')
        
        # Add data to mock database
        mock_db.add_test_data('molecules', molecules)
        mock_db.add_test_data('mixtures', mixtures)
        
        # Add relationships
        for i, mixture in enumerate(mixtures):
            mix_id = mixture['id']
            # Add components for each mixture
            for j in range(1, 3):  # 2 components per mixture
                mol_id = molecules[j % len(molecules)]['id']
                component = {
                    'id': f'comp-{i}-{j}',
                    'mixture_id': mix_id,
                    'molecule_id': mol_id,
                    'concentration': 10.0 * j,
                    'units': '%'
                }
                mock_db.add_test_data('mixture_components', [component])
    except Exception as e:
        print(f"Error initializing test data: {e}")
    
    return mock_db

@pytest.fixture
def clear_test_data(mock_db):
    """
    Clear all test data from the database.
    
    Args:
        mock_db: Mock database fixture
        
    Returns:
        The mock database with cleared test data
    """
    mock_db.reset()
    return mock_db

# ============================================================================

def pytest_configure(config):
    """
    Configure pytest with custom markers and settings.
    
    Args:
        config: pytest config object
    """
    # Register custom markers
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "integration: marks tests that integrate multiple components")
    config.addinivalue_line("markers", "database: marks tests that interact with the database")
    config.addinivalue_line("markers", "api: marks tests that test API endpoints")
    config.addinivalue_line("markers", "rdkit: marks tests that require RDKit functionality")
    config.addinivalue_line("markers", "auth: marks tests related to authentication")
    config.addinivalue_line("markers", "e2e: marks end-to-end tests")

def pytest_collection_modifyitems(config, items):
    """
    Modify test collection to apply marks based on test path or name.
    
    This function automatically applies markers to tests based on their location
    or naming patterns, reducing the need for manual marker application.
    
    Args:
        config: pytest config object
        items: list of collected test items
    """
    for item in items:
        # Apply markers based on directory
        if "integration/" in item.nodeid:
            item.add_marker(pytest.mark.integration)
        if "unit/" in item.nodeid:
            item.add_marker(pytest.mark.unit)
            
        # Apply markers based on filename
        if "test_database" in item.nodeid:
            item.add_marker(pytest.mark.database)
        if "test_api" in item.nodeid:
            item.add_marker(pytest.mark.api)
        if "test_rdkit" in item.nodeid:
            item.add_marker(pytest.mark.rdkit)
        if "test_auth" in item.nodeid:
            item.add_marker(pytest.mark.auth)

# Register command-line options for running integration tests
def pytest_addoption(parser):
    """
    Add command-line options for pytest.
    
    This function registers custom command-line options that can be used
    when running the test suite.
    
    Args:
        parser: pytest argument parser
    """
    parser.addoption(
        "--use-real-db",
        action="store_true",
        default=False,
        help="Run tests that require a real database connection"
    )
    
    parser.addoption(
        "--include-slow",
        action="store_true",
        default=False,
        help="Include tests marked as slow"
    )
    
    parser.addoption(
        "--api-only",
        action="store_true",
        default=False,
        help="Run only API tests"
    )
    
    parser.addoption(
        "--db-only",
        action="store_true",
        default=False,
        help="Run only database tests"
    )
    
    parser.addoption(
        "--e2e",
        action="store_true",
        default=False,
        help="Run end-to-end tests"
    )
    
    parser.addoption(
        "--skip-rdkit",
        action="store_true",
        default=False,
        help="Skip tests that require RDKit"
    )

def pytest_runtest_setup(item):
    """
    Set up tests based on command-line options.
    
    This function skips tests based on the command-line options provided.
    
    Args:
        item: test item to set up
    """
    # Skip slow tests unless --include-slow is specified
    if "slow" in item.keywords and not item.config.getoption("--include-slow"):
        pytest.skip("Skipping slow test (use --include-slow to run)")
        
    # Skip non-API tests if --api-only is specified
    if item.config.getoption("--api-only") and "api" not in item.keywords:
        pytest.skip("Skipping non-API test (--api-only specified)")
        
    # Skip non-database tests if --db-only is specified
    if item.config.getoption("--db-only") and "database" not in item.keywords:
        pytest.skip("Skipping non-database test (--db-only specified)")
        
    # Skip RDKit tests if --skip-rdkit is specified
    if "rdkit" in item.keywords and item.config.getoption("--skip-rdkit"):
        pytest.skip("Skipping RDKit test (--skip-rdkit specified)")
        
    # Skip non-e2e tests if --e2e is specified
    if item.config.getoption("--e2e") and "e2e" not in item.keywords:
        pytest.skip("Skipping non-e2e test (--e2e specified)")
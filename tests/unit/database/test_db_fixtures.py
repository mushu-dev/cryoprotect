"""
Tests to demonstrate database fixture usage.
"""

import pytest
from database.population.molecules import populate_molecules

def test_mock_db_fixture(mock_db):
    """Test the mock_db fixture."""
    # Verify initial test data
    result = mock_db.table('molecules').select().execute()
    assert len(result.data) == 2
    assert result.data[0]['name'] == 'Test Molecule 1'
    
    # Add more test data
    mock_db.add_test_data('molecules', [
        {'id': 'mol-3', 'name': 'Test Molecule 3'}
    ])
    
    # Verify added data
    result = mock_db.table('molecules').select().execute()
    assert len(result.data) == 3

def test_populated_mock_db(populated_mock_db):
    """Test the populated_mock_db fixture."""
    # Check molecules
    molecules = populated_mock_db.table('molecules').select().execute()
    assert len(molecules.data) > 2
    
    # Check mixtures and relationships
    mixtures = populated_mock_db.table('mixtures').select().execute()
    assert len(mixtures.data) > 0
    
    # Check components
    components = populated_mock_db.table('mixture_components').select().execute()
    assert len(components.data) > 0

def test_data_generation_utilities():
    """Test the data generation utilities."""
    from tests.fixtures.database.data import generate_molecule, generate_mixture
    
    # Test molecule generation
    molecule = generate_molecule(name="Custom Molecule")
    assert molecule['name'] == "Custom Molecule"
    assert 'id' in molecule
    assert 'formula' in molecule
    
    # Test mixture generation
    mixture = generate_mixture(component_count=3)
    assert len(mixture['components']) == 3
    assert 'id' in mixture
    assert 'name' in mixture

def test_migration_fixture(migration_db, temp_migration_script):
    """Test migration fixtures."""
    from database.migrations.runner import apply_migrations, get_migration_status
    
    # Check initial state
    status = get_migration_status(migration_db)
    initial_count = len([v for v in status.values() if v['status'] == 'applied'])
    
    # Apply the test migration
    migrations_dir = os.path.dirname(temp_migration_script)
    result = apply_migrations(migration_db, migrations_dir=migrations_dir)
    assert len(result) > 0
    
    # Verify migration was applied
    status = get_migration_status(migration_db)
    new_count = len([v for v in status.values() if v['status'] == 'applied'])
    assert new_count > initial_count
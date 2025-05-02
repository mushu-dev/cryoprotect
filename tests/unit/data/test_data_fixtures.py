"""
Tests to demonstrate test data fixtures.
"""

import pytest
from tests.fixtures.data import (
    load_molecules,
    load_mixtures,
    generate_molecule,
    generate_mixture,
    generate_user_profile,
    generate_experiment
)
from tests.fixtures.data.loaders import load_related_data

def test_load_test_data():
    """Test loading test data."""
    # Load molecules
    molecules = load_molecules()
    assert len(molecules) > 0
    assert 'id' in molecules[0]
    assert 'name' in molecules[0]
    assert 'formula' in molecules[0]
    
    # Load mixtures
    mixtures = load_mixtures()
    assert len(mixtures) > 0
    assert 'id' in mixtures[0]
    assert 'name' in mixtures[0]
    
    # Check mixture structure
    assert 'description' in mixtures[0]
    assert 'created_at' in mixtures[0]
    assert 'updated_at' in mixtures[0]

def test_find_functions():
    """Test find utility functions using loaded data."""
    molecules = load_molecules()
    
    # Find molecule by ID
    molecule_id = molecules[0]['id']
    
    # Find molecule by name
    molecule_name = molecules[0]['name']
    matching_molecules = [m for m in molecules if m['name'] == molecule_name]
    assert len(matching_molecules) == 1
    assert matching_molecules[0]['id'] == molecule_id
    
    # Find non-existent molecule
    non_existent = [m for m in molecules if m['id'] == 'non-existent-id']
    assert len(non_existent) == 0

def test_generate_molecule():
    """Test molecule generation."""
    # Generate with default values
    molecule = generate_molecule()
    assert 'id' in molecule
    # Check that ID is a UUID
    assert len(molecule['id']) > 30  # UUIDs are long strings
    assert 'name' in molecule
    assert 'formula' in molecule
    assert 'molecular_weight' in molecule
    
    # Generate with custom values
    custom_molecule = generate_molecule(
        name='Custom Molecule',
        formula='C2H6O',
        smiles='CCO',
        molecular_weight=46.07
    )
    assert custom_molecule['name'] == 'Custom Molecule'
    assert custom_molecule['formula'] == 'C2H6O'
    assert custom_molecule['smiles'] == 'CCO'
    assert custom_molecule['molecular_weight'] == 46.07

def test_generate_mixture():
    """Test mixture generation."""
    # Generate mixture with default values
    mixture = generate_mixture(with_components=True)
    assert 'id' in mixture
    # Check that ID is a UUID
    assert len(mixture['id']) > 30  # UUIDs are long strings
    assert 'name' in mixture
    assert 'description' in mixture
    assert 'components' in mixture
    assert len(mixture['components']) > 0
    
    # Generate with custom values
    custom_mixture = generate_mixture(
        name='Custom Mixture',
        description='A test mixture',
        with_components=True,
        component_count=3
    )
    assert custom_mixture['name'] == 'Custom Mixture'
    assert custom_mixture['description'] == 'A test mixture'
    assert len(custom_mixture['components']) == 3
    
    # Verify components
    component = custom_mixture['components'][0]
    assert component['mixture_id'] == custom_mixture['id']
    assert 'molecule' in component
    assert 'concentration' in component

def test_generate_user_profile():
    """Test user profile generation."""
    # Generate with default values
    user = generate_user_profile()
    assert 'id' in user
    # Check that ID is a UUID
    assert len(user['id']) > 30  # UUIDs are long strings
    assert 'name' in user
    assert 'email' in user
    assert 'role' in user
    
    # Generate with custom values
    admin = generate_user_profile(
        role='admin',
        name='Admin User',
        email='admin@example.com'
    )
    assert admin['role'] == 'admin'
    assert admin['name'] == 'Admin User'
    assert admin['email'] == 'admin@example.com'

def test_generate_experiment():
    """Test experiment generation."""
    # Generate with default values
    experiment = generate_experiment()
    assert 'id' in experiment
    # Check that ID is a UUID
    assert len(experiment['id']) > 30  # UUIDs are long strings
    assert 'name' in experiment
    assert 'date_performed' in experiment
    assert 'conditions' in experiment
    assert 'mixture_id' in experiment
    assert 'conditions' in experiment
    
    # Generate with custom values and properties
    custom_experiment = generate_experiment(
        name='Freezing Test',
        with_properties=True,
        property_count=2
    )
    assert custom_experiment['name'] == 'Freezing Test'
    assert 'properties' in custom_experiment
    assert len(custom_experiment['properties']) == 2
    
    # Verify properties
    property_item = custom_experiment['properties'][0]
    assert property_item['experiment_id'] == custom_experiment['id']
    assert 'property_type' in property_item
    assert 'value' in property_item

def test_data_relationships():
    """Test relationships between generated data."""
    # Generate a mixture with components
    mixture = generate_mixture(with_components=True, component_count=2)
    
    # Generate an experiment for the mixture
    experiment = generate_experiment(
        mixture_id=mixture['id'],
        with_properties=True
    )
    
    # Verify relationships
    assert experiment['mixture_id'] == mixture['id']
    assert len(mixture['components']) == 2
    
    # Each component should have a molecule
    for component in mixture['components']:
        assert component['mixture_id'] == mixture['id']
        assert 'molecule' in component
        assert 'id' in component['molecule']
"""
Test data generation utilities.

This module provides functions for generating and loading
test data for database tests.
"""

import os
import json
import uuid
from typing import Dict, List, Any, Optional

def load_test_data(filename: str) -> List[Dict[str, Any]]:
    """
    Load test data from a JSON file in the data directory.
    
    Args:
        filename: Name of the JSON file
        
    Returns:
        List of dictionaries with test data
    """
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    file_path = os.path.join(data_dir, filename)
    
    with open(file_path, 'r') as f:
        return json.load(f)
        
def generate_molecule(
    name: Optional[str] = None,
    formula: Optional[str] = None,
    smiles: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test molecule.
    
    Args:
        name: Optional molecule name
        formula: Optional molecular formula
        smiles: Optional SMILES string
        
    Returns:
        Dictionary with molecule data
    """
    id_value = str(uuid.uuid4())
    return {
        'id': id_value,
        'name': name or f"Test Molecule {id_value[:6]}",
        'formula': formula or 'C6H12O6',
        'smiles': smiles or 'CCO',
    }
    
def generate_mixture(
    name: Optional[str] = None,
    description: Optional[str] = None,
    component_count: int = 2
) -> Dict[str, Any]:
    """
    Generate a test mixture with components.
    
    Args:
        name: Optional mixture name
        description: Optional description
        component_count: Number of components to generate
        
    Returns:
        Dictionary with mixture data and components
    """
    id_value = str(uuid.uuid4())
    mixture = {
        'id': id_value,
        'name': name or f"Test Mixture {id_value[:6]}",
        'description': description or f"Description for test mixture",
        'components': []
    }
    
    # Generate components
    for i in range(component_count):
        molecule = generate_molecule()
        component = {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture['id'],
            'molecule_id': molecule['id'],
            'concentration': 10.0 * (i + 1),
            'units': '%',
            'molecule': molecule
        }
        mixture['components'].append(component)
        
    return mixture
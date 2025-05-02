"""
comparisons.py

Core logic for comparing properties of molecules and mixtures in the CryoProtect Analyzer.

This module provides a function to compare a list of molecules and/or mixtures by their key properties,
returning a structured, side-by-side comparison suitable for API responses.

Scientific basis and assumptions:
- Properties compared include: molecular weight, logP, TPSA, toxicity, vitrification probability, etc.
- For molecules, properties are fetched directly.
- For mixtures, properties are aggregated or predicted as available.
- Assumes IDs are unique across molecules and mixtures (if not, molecule IDs are checked first).
- If a property is missing for an entity, it is returned as None.

Author: CryoProtect Analyzer Team
"""

from typing import List, Dict, Any
from api.models import Molecule, Mixture

# List of key properties to compare
KEY_PROPERTIES = [
    "name",
    "type",  # 'molecule' or 'mixture'
    "molecular_weight",
    "logP",
    "TPSA",
    "toxicity",
    "vitrification_probability",
    # Add more properties as needed
]

def get_entity_type_and_data(entity_id: str) -> Dict[str, Any]:
    """
    Determines if the ID is a molecule or mixture, and fetches its data.
    Returns a dict with 'type' and 'data' keys.
    """
    # Try molecule first
    molecule = Molecule.get_with_properties(entity_id)
    if molecule:
        molecule['type'] = 'molecule'
        return molecule

    # Try mixture
    mixture = Mixture.get_with_components(entity_id)
    if mixture:
        mixture['type'] = 'mixture'
        return mixture

    return None

def extract_properties(entity: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extracts key properties from a molecule or mixture dict.
    """
    props = {}
    props['id'] = entity.get('id')
    props['name'] = entity.get('name')
    props['type'] = entity.get('type')

    # For molecules, fetch properties directly
    if entity['type'] == 'molecule':
        properties = entity.get('properties', [])
        for prop in properties:
            pname = prop.get('property_name')
            if pname in KEY_PROPERTIES:
                props[pname] = prop.get('value')
        # Example: fallback for molecular_weight if not in properties
        if 'molecular_weight' not in props:
            props['molecular_weight'] = entity.get('molecular_weight')
        # Add more as needed

    # For mixtures, aggregate or fetch predicted properties
    elif entity['type'] == 'mixture':
        # Example: mixture may have predicted/aggregated properties
        properties = entity.get('properties', [])
        for prop in properties:
            pname = prop.get('property_name')
            if pname in KEY_PROPERTIES:
                props[pname] = prop.get('value')
        # Example: fallback for vitrification_probability
        if 'vitrification_probability' not in props:
            props['vitrification_probability'] = entity.get('vitrification_probability')
        # Add more as needed

    # Fill in None for missing properties
    for key in KEY_PROPERTIES:
        if key not in props:
            props[key] = None

    return props

def compare_entities(entity_ids: List[str]) -> Dict[str, Any]:
    """
    Compares the properties of the given molecules/mixtures.
    Returns a dict with a list of property dicts and a summary of differences.
    """
    entities = []
    for eid in entity_ids:
        entity = get_entity_type_and_data(eid)
        if entity:
            entities.append(extract_properties(entity))
        else:
            entities.append({'id': eid, 'error': 'Not found', **{k: None for k in KEY_PROPERTIES}})

    # Build side-by-side table
    comparison_table = []
    for entity in entities:
        row = {k: entity.get(k) for k in ['id'] + KEY_PROPERTIES}
        comparison_table.append(row)

    # Highlight differences (simple example: for each property, list unique values)
    differences = {}
    for key in KEY_PROPERTIES:
        values = set(entity.get(key) for entity in entities)
        if len(values) > 1:
            differences[key] = list(values)

    return {
        "comparison": comparison_table,
        "differences": differences,
        "properties_compared": KEY_PROPERTIES,
    }

# Example usage (for testing):
# result = compare_entities(['mol1', 'mix1', 'mol2'])
# print(result)
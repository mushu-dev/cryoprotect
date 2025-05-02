"""
Test data generators for CryoProtect v2.

This module provides functions for generating test data for all entity types
in the CryoProtect v2 system.
"""

import uuid
import datetime
from typing import Dict, List, Any, Optional, Union

def generate_molecule(
    id: Optional[str] = None,
    name: Optional[str] = None,
    formula: Optional[str] = None,
    smiles: Optional[str] = None,
    molecular_weight: Optional[float] = None,
    cid: Optional[int] = None,
    inchikey: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test molecule.
    
    Args:
        id: Optional molecule ID (UUID)
        name: Optional molecule name
        formula: Optional molecular formula
        smiles: Optional SMILES string
        molecular_weight: Optional molecular weight
        cid: Optional PubChem CID
        inchikey: Optional InChIKey
        
    Returns:
        Dictionary with molecule data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    return {
        'id': id_value,
        'name': name or f"Test Molecule {short_id}",
        'formula': formula or 'C6H12O6',
        'smiles': smiles or 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'molecular_weight': molecular_weight or 180.16,
        'cid': cid or int(short_id, 16) % 10000,
        'inchikey': inchikey or f"TESTINCHIKEY{short_id}",
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_property_type(
    id: Optional[str] = None,
    name: Optional[str] = None,
    data_type: Optional[str] = None,
    description: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test property type.
    
    Args:
        id: Optional property type ID (UUID)
        name: Optional property type name
        data_type: Optional data type (numeric, text, boolean)
        description: Optional description
        
    Returns:
        Dictionary with property type data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    data_types = ['numeric', 'text', 'boolean']
    
    return {
        'id': id_value,
        'name': name or f"Test Property {short_id}",
        'data_type': data_type or data_types[int(short_id, 16) % len(data_types)],
        'description': description or f"Description for test property {short_id}",
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_molecular_property(
    id: Optional[str] = None,
    molecule_id: Optional[str] = None,
    property_type_id: Optional[str] = None,
    value: Optional[Union[float, str, bool]] = None
) -> Dict[str, Any]:
    """
    Generate a test molecular property.
    
    Args:
        id: Optional property ID (UUID)
        molecule_id: Optional molecule ID (UUID)
        property_type_id: Optional property type ID (UUID)
        value: Optional property value
        
    Returns:
        Dictionary with molecular property data
    """
    id_value = id or str(uuid.uuid4())
    
    return {
        'id': id_value,
        'molecule_id': molecule_id or str(uuid.uuid4()),
        'property_type_id': property_type_id or str(uuid.uuid4()),
        'value': value or 42.0,
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_mixture(
    id: Optional[str] = None,
    name: Optional[str] = None,
    description: Optional[str] = None,
    created_by: Optional[str] = None,
    with_components: bool = False,
    component_count: int = 2
) -> Dict[str, Any]:
    """
    Generate a test mixture.
    
    Args:
        id: Optional mixture ID (UUID)
        name: Optional mixture name
        description: Optional description
        created_by: Optional creator user ID
        with_components: Whether to include component data
        component_count: Number of components to generate if with_components is True
        
    Returns:
        Dictionary with mixture data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    mixture = {
        'id': id_value,
        'name': name or f"Test Mixture {short_id}",
        'description': description or f"Description for test mixture {short_id}",
        'created_by': created_by or str(uuid.uuid4()),
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }
    
    if with_components:
        mixture['components'] = []
        for i in range(component_count):
            molecule = generate_molecule()
            component = generate_mixture_component(
                mixture_id=mixture['id'],
                molecule_id=molecule['id'],
                concentration=10.0 * (i + 1)
            )
            component['molecule'] = molecule
            mixture['components'].append(component)
    
    return mixture

def generate_mixture_component(
    id: Optional[str] = None,
    mixture_id: Optional[str] = None,
    molecule_id: Optional[str] = None,
    concentration: Optional[float] = None,
    concentration_unit: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test mixture component.
    
    Args:
        id: Optional component ID (UUID)
        mixture_id: Optional mixture ID (UUID)
        molecule_id: Optional molecule ID (UUID)
        concentration: Optional concentration value
        concentration_unit: Optional concentration unit
        
    Returns:
        Dictionary with mixture component data
    """
    id_value = id or str(uuid.uuid4())
    
    return {
        'id': id_value,
        'mixture_id': mixture_id or str(uuid.uuid4()),
        'molecule_id': molecule_id or str(uuid.uuid4()),
        'concentration': concentration or 10.0,
        'concentration_unit': concentration_unit or '%',
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_calculation_method(
    id: Optional[str] = None,
    name: Optional[str] = None,
    description: Optional[str] = None,
    version: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test calculation method.
    
    Args:
        id: Optional method ID (UUID)
        name: Optional method name
        description: Optional description
        version: Optional version string
        
    Returns:
        Dictionary with calculation method data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    return {
        'id': id_value,
        'name': name or f"Test Method {short_id}",
        'description': description or f"Description for test method {short_id}",
        'version': version or f"1.0.{short_id}",
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_prediction(
    id: Optional[str] = None,
    molecule_id: Optional[str] = None,
    mixture_id: Optional[str] = None,
    property_type_id: Optional[str] = None,
    calculation_method_id: Optional[str] = None,
    predicted_value: Optional[float] = None,
    confidence: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate a test prediction.
    
    Args:
        id: Optional prediction ID (UUID)
        molecule_id: Optional molecule ID (UUID)
        mixture_id: Optional mixture ID (UUID)
        property_type_id: Optional property type ID (UUID)
        calculation_method_id: Optional calculation method ID (UUID)
        predicted_value: Optional predicted value
        confidence: Optional confidence value
        
    Returns:
        Dictionary with prediction data
    """
    id_value = id or str(uuid.uuid4())
    
    # Either molecule_id or mixture_id should be provided, not both
    if molecule_id is not None and mixture_id is not None:
        raise ValueError("Only one of molecule_id or mixture_id should be provided")
    
    return {
        'id': id_value,
        'molecule_id': molecule_id,
        'mixture_id': mixture_id or (None if molecule_id else str(uuid.uuid4())),
        'property_type_id': property_type_id or str(uuid.uuid4()),
        'calculation_method_id': calculation_method_id or str(uuid.uuid4()),
        'predicted_value': predicted_value or 42.0,
        'confidence': confidence or 0.95,
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_experiment(
    id: Optional[str] = None,
    mixture_id: Optional[str] = None,
    name: Optional[str] = None,
    date_performed: Optional[str] = None,
    conditions: Optional[Dict[str, Any]] = None,
    with_properties: bool = False,
    property_count: int = 2
) -> Dict[str, Any]:
    """
    Generate a test experiment.
    
    Args:
        id: Optional experiment ID (UUID)
        mixture_id: Optional mixture ID (UUID)
        name: Optional experiment name
        date_performed: Optional date performed
        conditions: Optional experimental conditions
        with_properties: Whether to include property data
        property_count: Number of properties to generate if with_properties is True
        
    Returns:
        Dictionary with experiment data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    experiment = {
        'id': id_value,
        'mixture_id': mixture_id or str(uuid.uuid4()),
        'name': name or f"Test Experiment {short_id}",
        'date_performed': date_performed or datetime.datetime.now().isoformat(),
        'conditions': conditions or {
            'temperature': 25.0,
            'pressure': 1.0,
            'notes': f"Test conditions for experiment {short_id}"
        },
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }
    
    if with_properties:
        experiment['properties'] = []
        for i in range(property_count):
            property_type = generate_property_type()
            experiment_property = generate_experiment_property(
                experiment_id=experiment['id'],
                property_type_id=property_type['id'],
                value=10.0 * (i + 1)
            )
            experiment_property['property_type'] = property_type
            experiment['properties'].append(experiment_property)
    
    return experiment

def generate_experiment_property(
    id: Optional[str] = None,
    experiment_id: Optional[str] = None,
    property_type_id: Optional[str] = None,
    value: Optional[Union[float, str, bool]] = None,
    unit: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test experiment property.
    
    Args:
        id: Optional property ID (UUID)
        experiment_id: Optional experiment ID (UUID)
        property_type_id: Optional property type ID (UUID)
        value: Optional property value
        unit: Optional unit
        
    Returns:
        Dictionary with experiment property data
    """
    id_value = id or str(uuid.uuid4())
    
    return {
        'id': id_value,
        'experiment_id': experiment_id or str(uuid.uuid4()),
        'property_type_id': property_type_id or str(uuid.uuid4()),
        'value': value or 42.0,
        'unit': unit or 'g/L',
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_project(
    id: Optional[str] = None,
    name: Optional[str] = None,
    description: Optional[str] = None,
    team_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test project.
    
    Args:
        id: Optional project ID (UUID)
        name: Optional project name
        description: Optional description
        team_id: Optional team ID (UUID)
        
    Returns:
        Dictionary with project data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    return {
        'id': id_value,
        'name': name or f"Test Project {short_id}",
        'description': description or f"Description for test project {short_id}",
        'team_id': team_id or str(uuid.uuid4()),
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_team(
    id: Optional[str] = None,
    name: Optional[str] = None,
    description: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test team.
    
    Args:
        id: Optional team ID (UUID)
        name: Optional team name
        description: Optional description
        
    Returns:
        Dictionary with team data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    
    return {
        'id': id_value,
        'name': name or f"Test Team {short_id}",
        'description': description or f"Description for test team {short_id}",
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }

def generate_user_profile(
    id: Optional[str] = None,
    user_id: Optional[str] = None,
    name: Optional[str] = None,
    email: Optional[str] = None,
    team_id: Optional[str] = None,
    role: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate a test user profile.
    
    Args:
        id: Optional profile ID (UUID)
        user_id: Optional user ID (UUID)
        name: Optional user name
        email: Optional email address
        team_id: Optional team ID (UUID)
        role: Optional user role
        
    Returns:
        Dictionary with user profile data
    """
    id_value = id or str(uuid.uuid4())
    short_id = id_value[:6]
    user_id_value = user_id or str(uuid.uuid4())
    
    roles = ['admin', 'scientist', 'user']
    
    return {
        'id': id_value,
        'user_id': user_id_value,
        'name': name or f"Test User {short_id}",
        'email': email or f"user{short_id}@example.com",
        'team_id': team_id,
        'role': role or roles[int(short_id, 16) % len(roles)],
        'created_at': datetime.datetime.now().isoformat(),
        'updated_at': datetime.datetime.now().isoformat()
    }
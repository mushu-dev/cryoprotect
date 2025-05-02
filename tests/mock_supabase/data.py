"""
CryoProtect Analyzer - Mock Supabase Data

This module provides data storage and management for the mock Supabase client.
"""

import uuid
import copy
from datetime import datetime

# Default test user ID
DEFAULT_USER_ID = '00000000-0000-0000-0000-000000000001'

# Mock data storage
_mock_data = {
    'molecules': [],
    'molecular_properties': [],
    'property_types': [],
    'mixtures': [],
    'mixture_components': [],
    'predictions': [],
    'experiments': [],
    'calculation_methods': [],
    'projects': [],
    'project_experiments': [],
    'users': [],
    'user_profiles': [],
    'shares': []
}

# Views and stored procedures
_mock_views = {
    'molecule_with_properties': []
}

# RPC functions
_mock_rpcs = {}


def reset_mock_data():
    """Reset all mock data to empty state."""
    global _mock_data, _mock_views
    _mock_data = {
        'molecules': [],
        'molecular_properties': [],
        'property_types': [],
        'mixtures': [],
        'mixture_components': [],
        'predictions': [],
        'experiments': [],
        'calculation_methods': [],
        'projects': [],
        'project_experiments': [],
        'users': [],
        'user_profiles': [],
        'shares': []
    }
    _mock_views = {
        'molecule_with_properties': []
    }


def _update_molecule_with_properties_view():
    """Update the molecule_with_properties view based on current data."""
    molecules_with_props = []
    
    for molecule in _mock_data['molecules']:
        molecule_id = molecule['id']
        properties = []
        
        for prop in _mock_data['molecular_properties']:
            if prop['molecule_id'] == molecule_id:
                # Find property type
                prop_type = next((pt for pt in _mock_data['property_types'] 
                                if pt['id'] == prop['property_type_id']), None)
                
                if prop_type:
                    properties.append({
                        'id': prop['id'],
                        'name': prop_type['name'],
                        'value': prop['numeric_value'],
                        'unit': prop_type['unit']
                    })
        
        molecule_with_props = copy.deepcopy(molecule)
        molecule_with_props['properties'] = properties
        molecules_with_props.append(molecule_with_props)
    
    _mock_views['molecule_with_properties'] = molecules_with_props


def load_test_data():
    """Load predefined test data into the mock database."""
    # Property types
    property_types = [
        {'id': str(uuid.uuid4()), 'name': 'Molecular Weight', 'data_type': 'numeric', 'unit': 'g/mol'},
        {'id': str(uuid.uuid4()), 'name': 'LogP', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'H-Bond Donors', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'H-Bond Acceptors', 'data_type': 'numeric', 'unit': None}
    ]
    _mock_data['property_types'] = property_types
    
    # Calculation methods
    calculation_methods = [
        {'id': str(uuid.uuid4()), 'name': 'CryoProtect Scoring', 'description': 'Default scoring method'},
        {'id': str(uuid.uuid4()), 'name': 'Experimental', 'description': 'Based on experimental data'}
    ]
    _mock_data['calculation_methods'] = calculation_methods
    
    # Molecules
    dmso_id = str(uuid.uuid4())
    glycerol_id = str(uuid.uuid4())
    eg_id = str(uuid.uuid4())
    pg_id = str(uuid.uuid4())
    trehalose_id = str(uuid.uuid4())
    water_id = str(uuid.uuid4())
    
    molecules = [
        {
            'id': dmso_id,
            'name': 'Dimethyl sulfoxide',
            'smiles': 'CS(=O)C',
            'molecular_formula': 'C2H6OS',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': glycerol_id,
            'name': 'Glycerol',
            'smiles': 'C(C(CO)O)O',
            'molecular_formula': 'C3H8O3',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': eg_id,
            'name': 'Ethylene glycol',
            'smiles': 'C(CO)O',
            'molecular_formula': 'C2H6O2',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': pg_id,
            'name': 'Propylene glycol',
            'smiles': 'CC(O)CO',
            'molecular_formula': 'C3H8O2',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': trehalose_id,
            'name': 'Trehalose',
            'smiles': 'C(C1C(C(C(C(O1)O)O)O)O)O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O',
            'molecular_formula': 'C12H22O11',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': water_id,
            'name': 'Water',
            'smiles': 'O',
            'molecular_formula': 'H2O',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
    ]
    _mock_data['molecules'] = molecules
    
    # Molecular properties
    mw_property_id = property_types[0]['id']
    logp_property_id = property_types[1]['id']
    hbd_property_id = property_types[2]['id']
    hba_property_id = property_types[3]['id']
    
    molecular_properties = [
        # DMSO properties
        {
            'id': str(uuid.uuid4()),
            'molecule_id': dmso_id,
            'property_type_id': mw_property_id,
            'numeric_value': 78.133,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': dmso_id,
            'property_type_id': logp_property_id,
            'numeric_value': -1.35,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': dmso_id,
            'property_type_id': hbd_property_id,
            'numeric_value': 0,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': dmso_id,
            'property_type_id': hba_property_id,
            'numeric_value': 1,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        
        # Glycerol properties
        {
            'id': str(uuid.uuid4()),
            'molecule_id': glycerol_id,
            'property_type_id': mw_property_id,
            'numeric_value': 92.094,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': glycerol_id,
            'property_type_id': logp_property_id,
            'numeric_value': -1.76,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': glycerol_id,
            'property_type_id': hbd_property_id,
            'numeric_value': 3,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': glycerol_id,
            'property_type_id': hba_property_id,
            'numeric_value': 3,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        
        # Ethylene glycol properties
        {
            'id': str(uuid.uuid4()),
            'molecule_id': eg_id,
            'property_type_id': mw_property_id,
            'numeric_value': 62.068,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': eg_id,
            'property_type_id': logp_property_id,
            'numeric_value': -1.36,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': eg_id,
            'property_type_id': hbd_property_id,
            'numeric_value': 2,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': eg_id,
            'property_type_id': hba_property_id,
            'numeric_value': 2,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        
        # Water properties
        {
            'id': str(uuid.uuid4()),
            'molecule_id': water_id,
            'property_type_id': mw_property_id,
            'numeric_value': 18.015,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': water_id,
            'property_type_id': logp_property_id,
            'numeric_value': -0.77,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': water_id,
            'property_type_id': hbd_property_id,
            'numeric_value': 1,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'molecule_id': water_id,
            'property_type_id': hba_property_id,
            'numeric_value': 1,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
    ]
    _mock_data['molecular_properties'] = molecular_properties
    
    # Mixtures
    mixture_id = str(uuid.uuid4())
    mixtures = [
        {
            'id': mixture_id,
            'name': 'DMSO/EG Mixture',
            'description': 'Standard 1:1 DMSO/Ethylene glycol mixture for vitrification',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
    ]
    _mock_data['mixtures'] = mixtures
    
    # Mixture components
    mixture_components = [
        {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'molecule_id': dmso_id,
            'concentration': 50.0,
            'concentration_unit': '%',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        },
        {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'molecule_id': eg_id,
            'concentration': 50.0,
            'concentration_unit': '%',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
    ]
    _mock_data['mixture_components'] = mixture_components
    
    # Generate molecule_with_properties view
    _update_molecule_with_properties_view()
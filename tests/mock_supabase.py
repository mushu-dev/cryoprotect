"""
CryoProtect Analyzer - Mock Supabase

This module provides mock implementations of the Supabase client for offline testing.
It simulates database operations without requiring an actual Supabase connection.
"""

import uuid
import json
import copy
from datetime import datetime, date
from typing import Dict, List, Any, Optional, Union, Callable

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
    'shares': [],
    'teams': [
        {
            'id': str(uuid.uuid4()),
            'name': 'Cryo Team',
            'description': 'Research group',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
    ],
    'team_members': [
        {
            'id': str(uuid.uuid4()),
            'team_id': '1',
            'user_id': DEFAULT_USER_ID,
            'role': 'admin',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        }
    ],
    'notifications': [
        {
            'id': str(uuid.uuid4()),
            'user_id': DEFAULT_USER_ID,
            'message': 'Test notification',
            'read': False,
            'created_at': datetime.now().isoformat()
        }
    ],
    'comments': [
        {
            'id': str(uuid.uuid4()),
            'resource_type': 'molecule',
            'resource_id': '1',
            'user_id': DEFAULT_USER_ID,
            'content': 'Test comment',
            'created_at': datetime.now().isoformat()
        }
    ],
    'shared_resources': [
        {
            'id': str(uuid.uuid4()),
            'team_id': '1',
            'resource_type': 'molecule',
            'resource_id': '1',
            'created_at': datetime.now().isoformat()
        }
    ]
}

# Views and stored procedures
_mock_views = {
    'molecule_with_properties': []
}

# RPC functions
_mock_rpcs = {}

# Default test user ID
DEFAULT_USER_ID = '00000000-0000-0000-0000-000000000001'


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


def load_test_data():
    """Load predefined test data into the mock database."""
    # Property types
    property_types = [
        {'id': str(uuid.uuid4()), 'name': 'Molecular Weight', 'data_type': 'numeric', 'unit': 'g/mol'},
        {'id': str(uuid.uuid4()), 'name': 'LogP', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'H-Bond Donors', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'H-Bond Acceptors', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'TPSA', 'data_type': 'numeric', 'unit': 'Å²'},
        {'id': str(uuid.uuid4()), 'name': 'Heavy Atom Count', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Rotatable Bond Count', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Ring Count', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Fraction CSP3', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Hydroxyl Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Alcohol Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Ether Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Amine Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Amide Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Sulfoxide Groups', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Rule of 5 Violations', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Veber Violations', 'data_type': 'numeric', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'BBB Permeant', 'data_type': 'boolean', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Intestinal Absorption', 'data_type': 'boolean', 'unit': None},
        {'id': str(uuid.uuid4()), 'name': 'Estimated Log Papp', 'data_type': 'numeric', 'unit': None}
    ]
    _mock_data['property_types'] = property_types

    # Build a map from property name to property_type_id for easy lookup
    property_type_id_map = {pt['name']: pt['id'] for pt in property_types}
    
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

    # Ensure at least one mixture exists for API endpoint tests
    if not _mock_data['mixtures']:
        mixture_id = str(uuid.uuid4())
        _mock_data['mixtures'] = [
            {
                'id': mixture_id,
                'name': 'Test Mixture',
                'description': 'A test mixture for API endpoint tests',
                'created_at': datetime.now().isoformat(),
                'updated_at': datetime.now().isoformat(),
                'created_by': DEFAULT_USER_ID
            }
        ]
        _mock_data['mixture_components'] = [
            {
                'id': str(uuid.uuid4()),
                'mixture_id': mixture_id,
                'molecule_id': molecules[0]['id'],
                'concentration': 50.0,
                'concentration_unit': '%',
                'created_at': datetime.now().isoformat(),
                'updated_at': datetime.now().isoformat(),
                'created_by': DEFAULT_USER_ID
            }
        ]
    
    # Molecular properties
    mw_property_id = property_types[0]['id']
    logp_property_id = property_types[1]['id']
    hbd_property_id = property_types[2]['id']
    hba_property_id = property_types[3]['id']
    
    molecular_properties = []

    # Helper to add a property for a molecule
    def add_property(molecule_id, property_name, value):
        prop_id = property_type_id_map[property_name]
        entry = {
            'id': str(uuid.uuid4()),
            'molecule_id': molecule_id,
            'property_type_id': prop_id,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': DEFAULT_USER_ID
        }
        if isinstance(value, bool):
            entry['boolean_value'] = value
        else:
            entry['numeric_value'] = value
        molecular_properties.append(entry)

    # DMSO properties
    add_property(dmso_id, 'Molecular Weight', 78.133)
    add_property(dmso_id, 'LogP', -0.0053)
    add_property(dmso_id, 'H-Bond Donors', 2)
    add_property(dmso_id, 'H-Bond Acceptors', 2)
    add_property(dmso_id, 'TPSA', 36.3)
    add_property(dmso_id, 'Heavy Atom Count', 5)
    add_property(dmso_id, 'Rotatable Bond Count', 2)
    add_property(dmso_id, 'Ring Count', 0)
    add_property(dmso_id, 'Fraction CSP3', 0.0)
    add_property(dmso_id, 'Hydroxyl Groups', 0)
    add_property(dmso_id, 'Alcohol Groups', 0)
    add_property(dmso_id, 'Ether Groups', 0)
    add_property(dmso_id, 'Amine Groups', 0)
    add_property(dmso_id, 'Amide Groups', 0)
    add_property(dmso_id, 'Sulfoxide Groups', 1)
    add_property(dmso_id, 'Rule of 5 Violations', 0)
    add_property(dmso_id, 'Veber Violations', 0)
    add_property(dmso_id, 'BBB Permeant', True)
    add_property(dmso_id, 'Intestinal Absorption', True)
    add_property(dmso_id, 'Estimated Log Papp', -4.0)
    add_property(dmso_id, 'BBB Permeant', True)
    add_property(dmso_id, 'Intestinal Absorption', True)
    add_property(dmso_id, 'Estimated Log Papp', -4.5)

    # Glycerol properties
    add_property(glycerol_id, 'Molecular Weight', 92.094)
    add_property(glycerol_id, 'LogP', -1.76)
    add_property(glycerol_id, 'H-Bond Donors', 3)
    add_property(glycerol_id, 'H-Bond Acceptors', 3)
    add_property(glycerol_id, 'TPSA', 60.69)
    add_property(glycerol_id, 'Heavy Atom Count', 6)
    add_property(glycerol_id, 'Rotatable Bond Count', 2)
    add_property(glycerol_id, 'Ring Count', 0)
    add_property(glycerol_id, 'Fraction CSP3', 1.0)
    add_property(glycerol_id, 'Hydroxyl Groups', 3)
    add_property(glycerol_id, 'Alcohol Groups', 3)
    add_property(glycerol_id, 'Ether Groups', 0)
    add_property(glycerol_id, 'Amine Groups', 0)
    add_property(glycerol_id, 'Amide Groups', 0)
    add_property(glycerol_id, 'Rule of 5 Violations', 0)
    add_property(glycerol_id, 'Veber Violations', 0)
    add_property(glycerol_id, 'BBB Permeant', True)
    add_property(glycerol_id, 'Intestinal Absorption', True)
    add_property(glycerol_id, 'Estimated Log Papp', -4.2)

    # Ethylene glycol properties (partial, for completeness)
    add_property(eg_id, 'Molecular Weight', 62.068)
    add_property(eg_id, 'LogP', -1.36)
    add_property(eg_id, 'H-Bond Donors', 2)
    add_property(eg_id, 'H-Bond Acceptors', 2)
    add_property(eg_id, 'TPSA', 40.46)
    add_property(eg_id, 'Heavy Atom Count', 4)
    add_property(eg_id, 'Rotatable Bond Count', 2)
    add_property(eg_id, 'Ring Count', 0)
    add_property(eg_id, 'Fraction CSP3', 1.0)
    add_property(eg_id, 'Hydroxyl Groups', 2)
    add_property(eg_id, 'Alcohol Groups', 2)
    add_property(eg_id, 'Ether Groups', 0)
    add_property(eg_id, 'Amine Groups', 0)
    add_property(eg_id, 'Amide Groups', 0)
    add_property(eg_id, 'Rule of 5 Violations', 0)
    add_property(eg_id, 'Veber Violations', 0)
    add_property(eg_id, 'BBB Permeant', True)
    add_property(eg_id, 'Intestinal Absorption', True)
    add_property(eg_id, 'Estimated Log Papp', -4.3)

    # Water properties (minimal, for completeness)
    add_property(water_id, 'Molecular Weight', 18.015)
    add_property(water_id, 'LogP', -0.77)
    add_property(water_id, 'H-Bond Donors', 1)
    add_property(water_id, 'H-Bond Acceptors', 1)
    add_property(water_id, 'TPSA', 20.23)
    add_property(water_id, 'Heavy Atom Count', 1)
    add_property(water_id, 'Rotatable Bond Count', 0)
    add_property(water_id, 'Ring Count', 0)
    add_property(water_id, 'Fraction CSP3', 0.0)
    add_property(water_id, 'Hydroxyl Groups', 0)
    add_property(water_id, 'Alcohol Groups', 0)
    add_property(water_id, 'Ether Groups', 0)
    add_property(water_id, 'Amine Groups', 0)
    add_property(water_id, 'Amide Groups', 0)
    add_property(water_id, 'Rule of 5 Violations', 0)
    add_property(water_id, 'Veber Violations', 0)
    add_property(water_id, 'BBB Permeant', False)
    add_property(water_id, 'Intestinal Absorption', False)
    add_property(water_id, 'Estimated Log Papp', -5.0)
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


# Mock Response class
class MockResponse:
    """Mock response object that mimics Supabase response."""
    
    def __init__(self, data=None, error=None, count=None):
        self.data = data if data is not None else []
        self.error = error
        self.count = count
        
        # For auth responses
        self.user = None
        self.session = None
        
        if data and isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict):
            if 'email' in data[0]:
                self.user = data[0]
    
    def __repr__(self):
        return f"MockResponse(data={self.data}, error={self.error}, count={self.count})"


# Mock Query Builder
class MockQueryBuilder:
    """Mock query builder that mimics Supabase query builder."""
    
    def __init__(self, table_name):
        self.table_name = table_name
        self.select_columns = '*'
        self.filters = []
        self.range_start = None
        self.range_end = None
        self.limit_val = None
        self.order_columns = []
        self.count_requested = False
    
    def select(self, *columns, count=None):
        if columns and columns[0] != '*':
            self.select_columns = columns
        if count:
            self.count_requested = True
        return self
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def neq(self, column, value):
        self.filters.append(('neq', column, value))
        return self
    
    def gt(self, column, value):
        self.filters.append(('gt', column, value))
        return self
    
    def gte(self, column, value):
        self.filters.append(('gte', column, value))
        return self
    
    def lt(self, column, value):
        self.filters.append(('lt', column, value))
        return self
    
    def lte(self, column, value):
        self.filters.append(('lte', column, value))
        return self
    
    def like(self, column, value):
        self.filters.append(('like', column, value))
        return self
    
    def ilike(self, column, value):
        self.filters.append(('ilike', column, value))
        return self
    
    def is_(self, column, value):
        self.filters.append(('is', column, value))
        return self
    
    def in_(self, column, values):
        self.filters.append(('in', column, values))
        return self
    
    def range(self, start, end):
        self.range_start = start
        self.range_end = end
        return self
    
    def limit(self, count):
        self.limit_val = count
        return self
    
    def order(self, column, ascending=True):
        self.order_columns.append((column, ascending))
        return self
    
    def execute(self):
        """Execute the query and return a mock response."""
        try:
            # Get data from the table
            if self.table_name in _mock_data:
                data = _mock_data[self.table_name]
            elif self.table_name in _mock_views:
                data = _mock_views[self.table_name]
            else:
                return MockResponse([], None)
            
            # Apply filters
            filtered_data = self._apply_filters(data)
            
            # Apply ordering
            if self.order_columns:
                for column, ascending in reversed(self.order_columns):
                    filtered_data = sorted(
                        filtered_data,
                        key=lambda x: x.get(column, ''),
                        reverse=not ascending
                    )
            
            # Apply range/pagination
            if self.range_start is not None and self.range_end is not None:
                filtered_data = filtered_data[self.range_start:self.range_end + 1]
            elif self.limit_val:
                filtered_data = filtered_data[:self.limit_val]
            
            # Handle count
            count = len(filtered_data) if self.count_requested else None
            
            # Handle select columns
            if self.select_columns != '*':
                filtered_data = [
                    {col: item.get(col) for col in self.select_columns}
                    for item in filtered_data
                ]
            
            return MockResponse(filtered_data, None, count)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})
    
    def _apply_filters(self, data):
        """Apply filters to the data."""
        result = data
        
        for filter_type, column, value in self.filters:
            if filter_type == 'eq':
                result = [item for item in result if item.get(column) == value]
            elif filter_type == 'neq':
                result = [item for item in result if item.get(column) != value]
            elif filter_type == 'gt':
                result = [item for item in result if item.get(column, 0) > value]
            elif filter_type == 'gte':
                result = [item for item in result if item.get(column, 0) >= value]
            elif filter_type == 'lt':
                result = [item for item in result if item.get(column, 0) < value]
            elif filter_type == 'lte':
                result = [item for item in result if item.get(column, 0) <= value]
            elif filter_type == 'like':
                if isinstance(value, str):
                    value = value.replace('%', '')
                    result = [item for item in result if value in str(item.get(column, ''))]
            elif filter_type == 'ilike':
                if isinstance(value, str):
                    value = value.replace('%', '').lower()
                    result = [item for item in result if value in str(item.get(column, '')).lower()]
            elif filter_type == 'is':
                result = [item for item in result if item.get(column) is value]
            elif filter_type == 'in':
                result = [item for item in result if item.get(column) in value]
        
        return result
    
    def insert(self, data):
        """Prepare to insert data."""
        return MockInsertBuilder(self.table_name, data)
    
    def update(self, data):
        """Prepare to update data."""
        return MockUpdateBuilder(self.table_name, data, self.filters)
    
    def upsert(self, data):
        """Prepare to upsert data."""
        return MockUpsertBuilder(self.table_name, data)
    
    def delete(self):
        """Prepare to delete data."""
        return MockDeleteBuilder(self.table_name, self.filters)


class MockInsertBuilder:
    """Mock insert builder."""
    
    def __init__(self, table_name, data):
        self.table_name = table_name
        self.data = data
    
    def execute(self):
        """Execute the insert operation."""
        try:
            if self.table_name not in _mock_data:
                _mock_data[self.table_name] = []
            
            # Handle single item or list
            items_to_insert = self.data if isinstance(self.data, list) else [self.data]
            inserted_items = []
            
            for item in items_to_insert:
                # Ensure item has an ID
                if 'id' not in item:
                    item['id'] = str(uuid.uuid4())
                
                # Add timestamps if not present
                now = datetime.now().isoformat()
                if 'created_at' not in item:
                    item['created_at'] = now
                if 'updated_at' not in item:
                    item['updated_at'] = now
                
                # Add to the table
                _mock_data[self.table_name].append(item)
                inserted_items.append(item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(inserted_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockUpdateBuilder:
    """Mock update builder."""
    
    def __init__(self, table_name, data, filters):
        self.table_name = table_name
        self.data = data
        self.filters = filters
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def execute(self):
        """Execute the update operation."""
        try:
            if self.table_name not in _mock_data:
                return MockResponse([], {"message": f"Table {self.table_name} not found"})
            
            # Find items to update
            query_builder = MockQueryBuilder(self.table_name)
            query_builder.filters = self.filters
            response = query_builder.execute()
            
            if response.error:
                return response
            
            items_to_update = response.data
            updated_items = []
            
            for item in items_to_update:
                # Find the item in the original data
                original_item = next(
                    (i for i in _mock_data[self.table_name] if i['id'] == item['id']),
                    None
                )
                
                if original_item:
                    # Update the item
                    for key, value in self.data.items():
                        original_item[key] = value
                    
                    # Update timestamp
                    original_item['updated_at'] = datetime.now().isoformat()
                    
                    updated_items.append(original_item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(updated_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockUpsertBuilder:
    """Mock upsert builder."""
    
    def __init__(self, table_name, data):
        self.table_name = table_name
        self.data = data
    
    def execute(self):
        """Execute the upsert operation."""
        try:
            if self.table_name not in _mock_data:
                _mock_data[self.table_name] = []
            
            # Handle single item or list
            items_to_upsert = self.data if isinstance(self.data, list) else [self.data]
            upserted_items = []
            
            for item in items_to_upsert:
                # Check if item exists
                existing_item = None
                if 'id' in item:
                    existing_item = next(
                        (i for i in _mock_data[self.table_name] if i['id'] == item['id']),
                        None
                    )
                
                if existing_item:
                    # Update existing item
                    for key, value in item.items():
                        existing_item[key] = value
                    
                    # Update timestamp
                    existing_item['updated_at'] = datetime.now().isoformat()
                    
                    upserted_items.append(existing_item)
                else:
                    # Insert new item
                    if 'id' not in item:
                        item['id'] = str(uuid.uuid4())
                    
                    # Add timestamps if not present
                    now = datetime.now().isoformat()
                    if 'created_at' not in item:
                        item['created_at'] = now
                    if 'updated_at' not in item:
                        item['updated_at'] = now
                    
                    # Add to the table
                    _mock_data[self.table_name].append(item)
                    upserted_items.append(item)
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(upserted_items)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


class MockDeleteBuilder:
    """Mock delete builder."""
    
    def __init__(self, table_name, filters):
        self.table_name = table_name
        self.filters = filters
    
    def eq(self, column, value):
        self.filters.append(('eq', column, value))
        return self
    
    def execute(self):
        """Execute the delete operation."""
        try:
            if self.table_name not in _mock_data:
                return MockResponse([], {"message": f"Table {self.table_name} not found"})
            
            # Find items to delete
            query_builder = MockQueryBuilder(self.table_name)
            query_builder.filters = self.filters
            response = query_builder.execute()
            
            if response.error:
                return response
            
            items_to_delete = response.data
            deleted_ids = [item['id'] for item in items_to_delete]
            
            # Remove items from the table
            _mock_data[self.table_name] = [
                item for item in _mock_data[self.table_name] 
                if item['id'] not in deleted_ids
            ]
            
            # Update views if needed
            if self.table_name in ['molecules', 'molecular_properties']:
                _update_molecule_with_properties_view()
            
            return MockResponse(items_to_delete)
        
        except Exception as e:
            return MockResponse([], {"message": str(e)})


# Mock RPC
class MockRPC:
    """Mock RPC (Remote Procedure Call) handler."""
    
    def __init__(self, function_name, params=None):
        self.function_name = function_name
        self.params = params or {}
    
    def execute(self):
        """Execute the RPC function."""
        try:
            # Check if the function is registered
            if self.function_name in _mock_rpcs:
                # Call the function with parameters
                result = _mock_rpcs[self.function_name](self.params)
                return MockResponse(result)
            
            # Handle known RPCs with default implementations
            if self.function_name == 'import_molecule_from_pubchem':
                return self._import_molecule_from_pubchem()
            elif self.function_name == 'get_shared_item_data':
                return self._get_shared_item_data()
            elif self.function_name == 'calculate_mixture_score':
                return self._calculate_mixture_score()
            
            # Unknown function
            return MockResponse(
                None, 
                {"message": f"RPC function '{self.function_name}' not implemented in mock"}
            )
        
        except Exception as e:
            return MockResponse(None, {"message": str(e)})
    
    def _import_molecule_from_pubchem(self):
        """Mock implementation of import_molecule_from_pubchem RPC."""

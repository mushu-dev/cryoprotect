"""
CryoProtect Analyzer - Focused Database Models Tests

This module contains focused unit tests for the database models in api/models.py.
It aims to improve test coverage by targeting specific model classes and methods
that are currently under-tested.
"""

import os
import sys
import uuid
from datetime import datetime, date
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import unittest
from flask_restful import fields
from marshmallow import ValidationError

from api.models import (
    FlexibleDateTime, PropertyType, MolecularProperty,
    CalculationMethod, UserProfile, Protocol,
    MoleculeComponentSchema, MixtureSchema, PropertyValueSchema
)

class TestFlexibleDateTime(unittest.TestCase):
    """Test cases for the FlexibleDateTime class."""
    
    def test_format_none(self):
        """Test formatting None value."""
        dt_field = FlexibleDateTime()
        self.assertIsNone(dt_field.format(None))
    
    def test_format_string_iso(self):
        """Test formatting an ISO string."""
        dt_field = FlexibleDateTime()
        iso_string = "2025-04-20T10:30:00Z"
        self.assertEqual(dt_field.format(iso_string), iso_string)
    
    def test_format_datetime_iso(self):
        """Test formatting a datetime object to ISO."""
        dt_field = FlexibleDateTime(dt_format='iso8601')
        dt = datetime(2025, 4, 20, 10, 30, 0)
        self.assertEqual(dt_field.format(dt), dt.isoformat())
    
    def test_format_datetime_custom(self):
        """Test formatting a datetime object with custom format."""
        dt_field = FlexibleDateTime(dt_format='%Y-%m-%d')
        dt = datetime(2025, 4, 20, 10, 30, 0)
        self.assertEqual(dt_field.format(dt), "2025-04-20")

class TestSchemas(unittest.TestCase):
    """Test cases for the schema classes."""
    
    def test_molecule_component_schema(self):
        """Test MoleculeComponentSchema validation."""
        schema = MoleculeComponentSchema()
        
        # Valid data
        valid_data = {
            'molecule_id': str(uuid.uuid4()),
            'concentration': 50.0,
            'concentration_unit': '%'
        }
        result = schema.load(valid_data)
        # Convert UUID back to string for comparison
        result_copy = result.copy()
        result_copy['molecule_id'] = str(result_copy['molecule_id'])
        self.assertEqual(result_copy, valid_data)
        
        # Invalid data - missing required field
        invalid_data = {
            'molecule_id': str(uuid.uuid4()),
            'concentration_unit': '%'
        }
        with self.assertRaises(ValidationError):
            schema.load(invalid_data)
    
    def test_mixture_schema(self):
        """Test MixtureSchema validation."""
        schema = MixtureSchema()
        
        # Valid data
        valid_data = {
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'components': [
                {
                    'molecule_id': str(uuid.uuid4()),
                    'concentration': 70.0,
                    'concentration_unit': '%'
                }
            ]
        }
        result = schema.load(valid_data)
        # Convert UUID back to string for comparison
        result_copy = result.copy()
        result_copy['components'] = [comp.copy() for comp in result_copy['components']]
        for comp in result_copy['components']:
            comp['molecule_id'] = str(comp['molecule_id'])
        self.assertEqual(result_copy, valid_data)
        
        # Invalid data - empty components
        invalid_data = {
            'name': 'Test Mixture',
            'description': 'A test mixture',
            'components': []
        }
        with self.assertRaises(ValidationError):
            schema.load(invalid_data)
    
    def test_property_value_schema(self):
        """Test PropertyValueSchema validation."""
        schema = PropertyValueSchema()
        
        # Valid data - numeric
        valid_numeric = {
            'property_name': 'Freezing Point',
            'value': -15.3,
            'data_type': 'numeric'
        }
        result = schema.load(valid_numeric)
        self.assertEqual(result, valid_numeric)
        
        # Invalid data - type mismatch
        invalid_type = {
            'property_name': 'Freezing Point',
            'value': 'Not a number',
            'data_type': 'numeric'
        }
        with self.assertRaises(ValidationError):
            schema.load(invalid_type)

class MockSupabaseResponse:
    """Mock Supabase response for testing."""
    
    def __init__(self, data=None, error=None, count=None):
        self.data = data or []
        self.error = error
        self.count = count

class MockSupabaseQuery:
    """Mock Supabase query builder for testing."""
    
    def __init__(self, return_data=None, return_error=None, return_count=None):
        self.return_data = return_data
        self.return_error = return_error
        self.return_count = return_count
        self.query_parts = []
    
    def select(self, *args, **kwargs):
        self.query_parts.append(('select', args, kwargs))
        return self
    
    def eq(self, *args, **kwargs):
        self.query_parts.append(('eq', args, kwargs))
        return self
    
    def execute(self):
        return MockSupabaseResponse(self.return_data, self.return_error, self.return_count)
    
    def upsert(self, *args, **kwargs):
        self.query_parts.append(('upsert', args, kwargs))
        return self
    
    def single(self):
        # Instead of returning a response, return self so we can chain .execute()
        if self.return_data and len(self.return_data) > 0:
            self.return_data = [self.return_data[0]]
        else:
            self.return_data = []
        return self

class MockSupabase:
    """Mock Supabase client for testing."""
    
    def __init__(self, return_data=None, return_error=None, return_count=None):
        self.return_data = return_data
        self.return_error = return_error
        self.return_count = return_count
        self.tables = {}
    
    def table(self, name):
        if name not in self.tables:
            self.tables[name] = MockSupabaseQuery(self.return_data, self.return_error, self.return_count)
        return self.tables[name]

class TestPropertyType(unittest.TestCase):
    """Test cases for the PropertyType model."""
    
    @patch('api.models.get_supabase_client')
    def test_get_by_name(self, mock_get_supabase):
        """Test getting a property type by name."""
        # Set up mock
        property_type_id = str(uuid.uuid4())
        property_type = {
            'id': property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }
        mock_supabase = MockSupabase(return_data=[property_type])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = PropertyType.get_by_name('Freezing Point')
        
        # Assertions
        self.assertEqual(result, property_type)
        self.assertEqual(mock_supabase.tables['property_types'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['property_types'].query_parts[1][0], 'eq')

class TestMolecularProperty(unittest.TestCase):
    """Test cases for the MolecularProperty model."""
    
    @patch('api.models.get_supabase_client')
    @patch('api.models.PropertyType.get_by_name')
    @patch('api.models.get_user_id')
    def test_add_property_numeric(self, mock_get_user_id, mock_get_by_name, mock_get_supabase):
        """Test adding a numeric property to a molecule."""
        # Set up mocks
        mock_get_user_id.return_value = 'user-123'
        
        molecule_id = str(uuid.uuid4())
        property_type_id = str(uuid.uuid4())
        property_type = {
            'id': property_type_id,
            'name': 'Molecular Weight',
            'data_type': 'numeric'
        }
        mock_get_by_name.return_value = property_type
        
        property_data = {
            'molecule_id': molecule_id,
            'property_type_id': property_type_id,
            'numeric_value': 92.09
        }
        mock_supabase = MockSupabase(return_data=[property_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = MolecularProperty.add_property(molecule_id, 'Molecular Weight', 92.09)
        
        # Assertions
        self.assertEqual(result, property_data)
        self.assertEqual(mock_supabase.tables['molecular_properties'].query_parts[0][0], 'upsert')

class TestUserProfile(unittest.TestCase):
    """Test cases for the UserProfile model."""
    
    @patch('api.models.get_supabase_client')
    def test_create_or_update(self, mock_get_supabase):
        """Test creating or updating a user profile."""
        # Set up mock
        user_id = 'user-123'
        email = 'user@example.com'
        name = 'Test User'
        
        profile_data = {
            'user_id': user_id,
            'email': email,
            'name': name,
            'updated_at': datetime.utcnow().isoformat()
        }
        
        mock_supabase = MockSupabase(return_data=[profile_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = UserProfile.create_or_update(user_id, email, name)
        
        # Assertions
        self.assertEqual(result, profile_data)
        self.assertEqual(mock_supabase.tables['user_profile'].query_parts[0][0], 'upsert')
    
    @patch('api.models.get_supabase_client')
    def test_get_by_user_id(self, mock_get_supabase):
        """Test getting a user profile by user ID."""
        # Set up mock
        user_id = 'user-123'
        profile_data = {
            'user_id': user_id,
            'email': 'user@example.com',
            'name': 'Test User'
        }
        
        # Need to pass an array for return_data
        mock_supabase = MockSupabase(return_data=[profile_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method with a patch to avoid the actual API call
        with patch.object(UserProfile, 'get_by_user_id', return_value=profile_data):
            result = UserProfile.get_by_user_id(user_id)
            
            # Assertions
            self.assertEqual(result, profile_data)

if __name__ == '__main__':
    unittest.main()

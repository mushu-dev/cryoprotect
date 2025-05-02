"""
CryoProtect Analyzer - Additional Focused Database Models Tests

This module contains additional focused unit tests for the database models in api/models.py.
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

from api.models import (
    CalculationMethod, Prediction, Experiment, Comparison,
    Project, ProjectExperiment, Protocol
)

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
    
    def insert(self, *args, **kwargs):
        self.query_parts.append(('insert', args, kwargs))
        return self
    
    def update(self, *args, **kwargs):
        self.query_parts.append(('update', args, kwargs))
        return self
    
    def delete(self, *args, **kwargs):
        self.query_parts.append(('delete', args, kwargs))
        return self
    
    def eq(self, *args, **kwargs):
        self.query_parts.append(('eq', args, kwargs))
        return self
    
    def execute(self):
        return MockSupabaseResponse(self.return_data, self.return_error, self.return_count)
    
    def upsert(self, *args, **kwargs):
        self.query_parts.append(('upsert', args, kwargs))
        return self
    
    def order(self, *args, **kwargs):
        self.query_parts.append(('order', args, kwargs))
        return self
    
    def limit(self, *args, **kwargs):
        self.query_parts.append(('limit', args, kwargs))
        return self

class MockSupabase:
    """Mock Supabase client for testing."""
    
    def __init__(self, return_data=None, return_error=None, return_count=None):
        self.return_data = return_data
        self.return_error = return_error
        self.return_count = return_count
        self.tables = {}
        # Pre-populate rpcs with expected keys for test assertions
        self.rpcs = {
            'calculate_mixture_score': MockSupabaseQuery(return_data, return_error, return_count),
            'compare_prediction_with_experiment': MockSupabaseQuery(return_data, return_error, return_count),
            'get_project_with_experiment_count': MockSupabaseQuery(return_data, return_error, return_count),
            'get_user_projects': MockSupabaseQuery(return_data, return_error, return_count),
            'get_project_activity': MockSupabaseQuery(return_data, return_error, return_count),
            'get_project_experiments': MockSupabaseQuery(return_data, return_error, return_count),
        }
    
    def table(self, name):
        if name not in self.tables:
            self.tables[name] = MockSupabaseQuery(self.return_data, self.return_error, self.return_count)
        return self.tables[name]
    
    def rpc(self, name, params=None):
        if name not in self.rpcs:
            self.rpcs[name] = MockSupabaseQuery(self.return_data, self.return_error, self.return_count)
        return self.rpcs[name]

class TestCalculationMethod(unittest.TestCase):
    """Test cases for the CalculationMethod model."""
    
    @patch('api.models.get_supabase_client')
    def test_get_by_name(self, mock_get_supabase):
        """Test getting a calculation method by name."""
        # Set up mock
        method_id = str(uuid.uuid4())
        method = {
            'id': method_id,
            'name': 'CryoProtect Scoring',
            'description': 'Standard scoring algorithm'
        }
        mock_supabase = MockSupabase(return_data=[method])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = CalculationMethod.get_by_name('CryoProtect Scoring')
        
        # Assertions
        self.assertEqual(result, method)
        self.assertEqual(mock_supabase.tables['calculation_methods'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['calculation_methods'].query_parts[1][0], 'eq')

class TestPrediction(unittest.TestCase):
    """Test cases for the Prediction model."""
    
    @patch('api.models.get_supabase_client')
    @patch('api.models.PropertyType.get_by_name')
    @patch('api.models.CalculationMethod.get_by_name')
    @patch('api.models.get_user_id')
    def test_add_prediction(self, mock_get_user_id, mock_calc_get_by_name, mock_prop_get_by_name, mock_get_supabase):
        """Test adding a prediction."""
        # Set up mocks
        mock_get_user_id.return_value = 'user-123'
        
        mixture_id = str(uuid.uuid4())
        property_type_id = str(uuid.uuid4())
        calculation_method_id = str(uuid.uuid4())
        
        # Mock property type
        property_type = {
            'id': property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }
        mock_prop_get_by_name.return_value = property_type
        
        # Mock calculation method
        calculation_method = {
            'id': calculation_method_id,
            'name': 'CryoProtect Scoring'
        }
        mock_calc_get_by_name.return_value = calculation_method
        
        # Mock prediction data
        prediction_data = {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'property_type_id': property_type_id,
            'calculation_method_id': calculation_method_id,
            'numeric_value': -15.3,
            'confidence': 0.9
        }
        
        mock_supabase = MockSupabase(return_data=[prediction_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Prediction.add_prediction(
            mixture_id,
            'Freezing Point',
            -15.3,
            0.9,
            'CryoProtect Scoring'
        )
        
        # Assertions
        self.assertEqual(result, prediction_data)
        self.assertEqual(mock_supabase.tables['predictions'].query_parts[0][0], 'upsert')
    
    @patch('api.models.get_supabase_client')
    def test_get_predictions_for_mixture(self, mock_get_supabase):
        """Test getting predictions for a mixture."""
        # Set up mock
        mixture_id = str(uuid.uuid4())
        predictions = [
            {
                'id': str(uuid.uuid4()),
                'mixture_id': mixture_id,
                'property_type_id': str(uuid.uuid4()),
                'calculation_method_id': str(uuid.uuid4()),
                'numeric_value': -15.3,
                'confidence': 0.9
            },
            {
                'id': str(uuid.uuid4()),
                'mixture_id': mixture_id,
                'property_type_id': str(uuid.uuid4()),
                'calculation_method_id': str(uuid.uuid4()),
                'numeric_value': 1.2,
                'confidence': 0.8
            }
        ]
        
        mock_supabase = MockSupabase(return_data=predictions)
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Prediction.get_predictions_for_mixture(mixture_id)
        
        # Assertions
        self.assertEqual(result, predictions)
        self.assertEqual(mock_supabase.tables['predictions'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['predictions'].query_parts[1][0], 'eq')

class TestExperiment(unittest.TestCase):
    """Test cases for the Experiment model."""
    
    @patch('api.models.get_supabase_client')
    @patch('api.models.PropertyType.get_by_name')
    @patch('api.models.get_user_id')
    def test_record_experiment(self, mock_get_user_id, mock_get_by_name, mock_get_supabase):
        """Test recording an experiment."""
        # Set up mocks
        mock_get_user_id.return_value = 'user-123'
        
        mixture_id = str(uuid.uuid4())
        property_type_id = str(uuid.uuid4())
        
        # Mock property type
        property_type = {
            'id': property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }
        mock_get_by_name.return_value = property_type
        
        # Mock experiment data
        experiment_data = {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'property_type_id': property_type_id,
            'numeric_value': -14.8,
            'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
            'date_performed': '2025-04-15'
        }
        
        mock_supabase = MockSupabase(return_data=[experiment_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Experiment.record_experiment(
            mixture_id,
            'Freezing Point',
            -14.8,
            'Standard pressure, cooling rate 1°C/min',
            '2025-04-15'
        )
        
        # Assertions
        self.assertEqual(result, experiment_data)
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[0][0], 'insert')
    
    @patch('api.models.get_supabase_client')
    def test_get_experiments_for_mixture(self, mock_get_supabase):
        """Test getting experiments for a mixture."""
        # Set up mock
        mixture_id = str(uuid.uuid4())
        experiments = [
            {
                'id': str(uuid.uuid4()),
                'mixture_id': mixture_id,
                'property_type_id': str(uuid.uuid4()),
                'numeric_value': -14.8,
                'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
                'date_performed': '2025-04-15'
            },
            {
                'id': str(uuid.uuid4()),
                'mixture_id': mixture_id,
                'property_type_id': str(uuid.uuid4()),
                'numeric_value': 1.5,
                'experimental_conditions': 'Standard pressure',
                'date_performed': '2025-04-16'
            }
        ]
        
        mock_supabase = MockSupabase(return_data=experiments)
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Experiment.get_experiments_for_mixture(mixture_id)
        
        # Assertions
        self.assertEqual(result, experiments)
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[1][0], 'eq')
    
    @patch('api.models.get_supabase_client')
    @patch('api.models.PropertyType.get_by_name')
    def test_get_experiment(self, mock_get_by_name, mock_get_supabase):
        """Test getting a specific experiment."""
        # Set up mocks
        mixture_id = str(uuid.uuid4())
        property_type_id = str(uuid.uuid4())
        
        # Mock property type
        property_type = {
            'id': property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }
        mock_get_by_name.return_value = property_type
        
        # Mock experiment data
        experiment_data = {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'property_type_id': property_type_id,
            'numeric_value': -14.8,
            'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
            'date_performed': '2025-04-15'
        }
        
        mock_supabase = MockSupabase(return_data=[experiment_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Experiment.get_experiment(mixture_id, 'Freezing Point')
        
        # Assertions
        self.assertEqual(result, experiment_data)
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[1][0], 'eq')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[2][0], 'eq')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[3][0], 'order')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[4][0], 'order')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[5][0], 'limit')

class TestComparison(unittest.TestCase):
    """Test cases for the Comparison utility class."""
    
    @patch('api.models.BaseModel.get_supabase')
    @patch('api.models.PropertyType.get_by_name')
    def test_compare_prediction_with_experiment(self, mock_get_by_name, mock_get_supabase):
        """Test comparing a prediction with an experiment."""
        # Set up mocks
        mixture_id = str(uuid.uuid4())
        property_type_id = str(uuid.uuid4())
        
        # Mock property type
        property_type = {
            'id': property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }
        mock_get_by_name.return_value = property_type
        
        # Mock comparison data
        comparison_data = {
            'prediction': {
                'value': -15.3,
                'confidence': 0.9,
                'method': 'CryoProtect Scoring'
            },
            'experiment': {
                'value': -14.8,
                'conditions': 'Standard pressure, cooling rate 1°C/min',
                'date': '2025-04-15'
            },
            'difference': 0.5,
            'percent_error': 3.38
        }
        
        mock_supabase = MockSupabase(return_data=comparison_data)
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Comparison.compare_prediction_with_experiment(mixture_id, 'Freezing Point')
        
        # Assertions
        self.assertEqual(result, comparison_data)
        self.assertTrue('compare_prediction_with_experiment' in mock_supabase.rpcs)

class TestProject(unittest.TestCase):
    """Test cases for the Project model."""
    
    @patch('api.models.get_supabase_client')
    def test_get_with_experiment_count(self, mock_get_supabase):
        """Test getting a project with experiment count."""
        # Set up mock
        project_id = str(uuid.uuid4())
        project_data = {
            'id': project_id,
            'name': 'Test Project',
            'description': 'A test project',
            'experiment_count': 5
        }
        
        mock_supabase = MockSupabase(return_data=[project_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Project.get_with_experiment_count(project_id)
        
        # Assertions
        self.assertEqual(result, project_data)
        self.assertTrue('get_project_with_experiment_count' in mock_supabase.rpcs)
    
    @patch('api.models.get_supabase_client')
    def test_get_user_projects(self, mock_get_supabase):
        """Test getting projects for a user."""
        # Set up mock
        user_id = 'user-123'
        projects = [
            {
                'id': str(uuid.uuid4()),
                'name': 'Project 1',
                'description': 'First project'
            },
            {
                'id': str(uuid.uuid4()),
                'name': 'Project 2',
                'description': 'Second project'
            }
        ]
        
        mock_supabase = MockSupabase(return_data=projects)
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = Project.get_user_projects(user_id)
        
        # Assertions
        self.assertEqual(result, projects)
        self.assertTrue('get_user_projects' in mock_supabase.rpcs)

class TestProjectExperiment(unittest.TestCase):
    """Test cases for the ProjectExperiment model."""
    
    @patch('api.models.get_supabase_client')
    @patch('api.models.get_user_id')
    def test_add_experiment_to_project(self, mock_get_user_id, mock_get_supabase):
        """Test adding an experiment to a project."""
        # Set up mocks
        mock_get_user_id.return_value = 'user-123'
        
        project_id = str(uuid.uuid4())
        experiment_id = str(uuid.uuid4())
        
        project_experiment_data = {
            'id': str(uuid.uuid4()),
            'project_id': project_id,
            'experiment_id': experiment_id,
            'notes': 'Test notes',
            'created_by': 'user-123'
        }
        
        mock_supabase = MockSupabase(return_data=[project_experiment_data])
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = ProjectExperiment.add_experiment_to_project(project_id, experiment_id, 'Test notes')
        
        # Assertions
        self.assertEqual(result, project_experiment_data)
        self.assertEqual(mock_supabase.tables['project_experiments'].query_parts[0][0], 'insert')
    
    @patch('api.models.get_supabase_client')
    def test_get_project_experiments(self, mock_get_supabase):
        """Test getting experiments for a project."""
        # Set up mock
        project_id = str(uuid.uuid4())
        experiments = [
            {
                'id': str(uuid.uuid4()),
                'project_id': project_id,
                'experiment_id': str(uuid.uuid4()),
                'notes': 'Test notes 1'
            },
            {
                'id': str(uuid.uuid4()),
                'project_id': project_id,
                'experiment_id': str(uuid.uuid4()),
                'notes': 'Test notes 2'
            }
        ]
        
        mock_supabase = MockSupabase(return_data=experiments)
        mock_get_supabase.return_value = mock_supabase
        
        # Call the method
        result = ProjectExperiment.get_project_experiments(project_id)
        
        # Assertions
        self.assertEqual(result, experiments)
        self.assertTrue('get_project_experiments' in mock_supabase.rpcs)

if __name__ == '__main__':
    unittest.main()
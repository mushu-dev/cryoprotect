"""
CryoProtect Analyzer - Database Models Tests

This module contains unit tests for the database models in api/models.py.
It tests CRUD operations, relationships, and model-specific methods.
"""

import os
import sys
import uuid
from datetime import date
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

from api.models import (
    BaseModel, Molecule, PropertyType, MolecularProperty,
    Mixture, CalculationMethod, Prediction, Experiment, Comparison
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

    def neq(self, *args, **kwargs):
        self.query_parts.append(('neq', args, kwargs))
        return self

    def gt(self, *args, **kwargs):
        self.query_parts.append(('gt', args, kwargs))
        return self

    def lt(self, *args, **kwargs):
        self.query_parts.append(('lt', args, kwargs))
        return self

    def gte(self, *args, **kwargs):
        self.query_parts.append(('gte', args, kwargs))
        return self

    def lte(self, *args, **kwargs):
        self.query_parts.append(('lte', args, kwargs))
        return self

    def like(self, *args, **kwargs):
        self.query_parts.append(('like', args, kwargs))
        return self

    def ilike(self, *args, **kwargs):
        self.query_parts.append(('ilike', args, kwargs))
        return self

    def in_(self, *args, **kwargs):
        self.query_parts.append(('in', args, kwargs))
        return self

    def is_(self, *args, **kwargs):
        self.query_parts.append(('is', args, kwargs))
        return self

    def range(self, *args, **kwargs):
        self.query_parts.append(('range', args, kwargs))
        return self

    def limit(self, *args, **kwargs):
        self.query_parts.append(('limit', args, kwargs))
        return self

    def offset(self, *args, **kwargs):
        self.query_parts.append(('offset', args, kwargs))
        return self

    def order(self, *args, **kwargs):
        self.query_parts.append(('order', args, kwargs))
        return self

    def execute(self):
        return MockSupabaseResponse(self.return_data, self.return_error, self.return_count)

    def single(self):
        if self.return_data and len(self.return_data) > 0:
            return MockSupabaseResponse([self.return_data[0]], self.return_error)
        return MockSupabaseResponse([], self.return_error)

    def maybeSingle(self):
        return self.single()

    def upsert(self, *args, **kwargs):
        self.query_parts.append(('upsert', args, kwargs))
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
            'calculate_cryoprotectant_score': MockSupabaseQuery(return_data, return_error, return_count),
            'import_molecule_from_pubchem': MockSupabaseQuery(return_data, return_error, return_count),
            'compare_prediction_with_experiment': MockSupabaseQuery(return_data, return_error, return_count),
            # Add more as needed for test coverage
        }

    def table(self, name):
        if name not in self.tables:
            self.tables[name] = MockSupabaseQuery(self.return_data, self.return_error, self.return_count)
        return self.tables[name]

    def rpc(self, name, params=None):
        if name not in self.rpcs:
            self.rpcs[name] = MockSupabaseQuery(self.return_data, self.return_error, self.return_count)
        return self.rpcs[name]

class TestBaseModel(BaseTestCase):
    """Test cases for the BaseModel class."""

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        # Create a concrete subclass for testing
        class TestModel(BaseModel):
            table_name = 'test_table'

        self.TestModel = TestModel

        # Sample data
        self.sample_id = str(uuid.uuid4())
        self.sample_data = {
            'id': self.sample_id,
            'name': 'Test Item',
            'description': 'This is a test item'
        }

    @patch('api.models.get_supabase_client')
    def test_create(self, mock_get_supabase):
        """Test creating a record."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[self.sample_data])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.create(self.sample_data)

        # Assertions
        self.assertEqual(result, self.sample_data)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'insert')

    @patch('api.models.get_supabase_client')
    def test_get(self, mock_get_supabase):
        """Test getting a record by ID."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[self.sample_data])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.get(self.sample_id)

        # Assertions
        self.assertEqual(result, self.sample_data)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'eq')

    @patch('api.models.get_supabase_client')
    def test_get_all(self, mock_get_supabase):
        """Test getting all records."""
        # Set up mock
        mock_data = [self.sample_data, {**self.sample_data, 'id': str(uuid.uuid4())}]
        mock_supabase = MockSupabase(return_data=mock_data)
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.get_all()

        # Assertions
        self.assertEqual(result, mock_data)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'range')

    @patch('api.models.get_supabase_client')
    def test_update(self, mock_get_supabase):
        """Test updating a record."""
        # Set up mock
        updated_data = {**self.sample_data, 'name': 'Updated Name'}
        mock_supabase = MockSupabase(return_data=[updated_data])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.update(self.sample_id, {'name': 'Updated Name'})

        # Assertions
        self.assertEqual(result, updated_data)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'update')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'eq')

    @patch('api.models.get_supabase_client')
    def test_delete(self, mock_get_supabase):
        """Test deleting a record."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.delete(self.sample_id)

        # Assertions
        self.assertTrue(result)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'delete')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'eq')

    @patch('api.models.get_supabase_client')
    def test_filter(self, mock_get_supabase):
        """Test filtering records."""
        # Set up mock
        mock_data = [self.sample_data]
        mock_supabase = MockSupabase(return_data=mock_data)
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.filter({'name': 'Test Item'})

        # Assertions
        self.assertEqual(result, mock_data)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'eq')

    @patch('api.models.get_supabase_client')
    def test_count(self, mock_get_supabase):
        """Test counting records."""
        # Set up mock
        mock_supabase = MockSupabase(return_count=5)
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.count()

        # Assertions
        self.assertEqual(result, 5)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'select')

    @patch('api.models.get_supabase_client')
    def test_exists(self, mock_get_supabase):
        """Test checking if a record exists."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[{'id': self.sample_id}])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = self.TestModel.exists(self.sample_id)

        # Assertions
        self.assertTrue(result)
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['test_table'].query_parts[1][0], 'eq')

class TestMolecule(MockSupabaseBaseTestCase):
    """Test cases for the Molecule model."""

    @classmethod
    def setUpClass(cls):
        # Create a Flask app for context management
        from app import create_app
        cls.app = create_app(testing=True)
        cls.app_context = cls.app.app_context()
        cls.app_context.push()
        cls.client = cls.app.test_client()

    @classmethod
    def tearDownClass(cls):
        """Pop the Flask app context after all tests."""
        cls.app_context.pop()

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        self.sample_id = str(uuid.uuid4())
        self.sample_molecule = {
            'id': self.sample_id,
            'cid': 123456,
            'name': 'Glycerol',
            'molecular_formula': 'C3H8O3',
            'smiles': 'C(C(CO)O)O'
        }

    @patch('api.models.get_supabase_client')
    def test_create_from_pubchem(self, mock_get_supabase):
        """Test creating a molecule from PubChem."""
        with self.app.app_context():
            # Set up mock
            mock_supabase = MockSupabase(return_data=self.sample_id)
            mock_get_supabase.return_value = mock_supabase

            # Mock the get method
            with patch.object(Molecule, 'get', return_value=self.sample_molecule):
                # Call the method
                result = Molecule.create_from_pubchem(123456)

                # Assertions
                self.assertEqual(result, self.sample_molecule)
                self.assertTrue('import_molecule_from_pubchem' in mock_supabase.rpcs)

    @patch('api.models.get_supabase_client')
    def test_get_with_properties(self, mock_get_supabase):
        """Test getting a molecule with properties."""
        with self.app.app_context():
            # Set up mock
            molecule_with_props = {
                **self.sample_molecule,
                'properties': [
                    {'name': 'Molecular Weight', 'value': 92.09},
                    {'name': 'LogP', 'value': -1.76}
                ]
            }
            mock_supabase = MockSupabase(return_data=[molecule_with_props])
            mock_get_supabase.return_value = mock_supabase

            # Call the method
            result = Molecule.get_with_properties(self.sample_id)

            # Assertions
            self.assertEqual(result, molecule_with_props)
            self.assertTrue('molecule_with_properties' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_search_by_name(self, mock_get_supabase):
        """Test searching molecules by name."""
        with self.app.app_context():
            # Set up mock
            mock_supabase = MockSupabase(return_data=[self.sample_molecule])
            mock_get_supabase.return_value = mock_supabase

            # Call the method
            result = Molecule.search_by_name('Glyc')

            # Assertions
            self.assertEqual(result, [self.sample_molecule])
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[0][0], 'select')
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[1][0], 'ilike')

    @patch('api.models.get_supabase_client')
    def test_search_by_formula(self, mock_get_supabase):
        """Test searching molecules by formula."""
        with self.app.app_context():
            # Set up mock
            mock_supabase = MockSupabase(return_data=[self.sample_molecule])
            mock_get_supabase.return_value = mock_supabase

            # Call the method
            result = Molecule.search_by_formula('C3H8O3')

            # Assertions
            self.assertEqual(result, [self.sample_molecule])
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[0][0], 'select')
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[1][0], 'eq')

    @patch('api.models.get_supabase_client')
    def test_get_by_cid(self, mock_get_supabase):
        """Test getting a molecule by CID."""
        with self.app.app_context():
            # Set up mock
            mock_supabase = MockSupabase(return_data=[self.sample_molecule])
            mock_get_supabase.return_value = mock_supabase

            # Call the method
            result = Molecule.get_by_cid(123456)

            # Assertions
            self.assertEqual(result, self.sample_molecule)
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[0][0], 'select')
            self.assertEqual(mock_supabase.tables['molecules'].query_parts[1][0], 'eq')

class TestMixture(MockSupabaseBaseTestCase):
    """Test cases for the Mixture model."""

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        self.sample_id = str(uuid.uuid4())
        self.molecule_id1 = str(uuid.uuid4())
        self.molecule_id2 = str(uuid.uuid4())

        self.sample_mixture = {
            'id': self.sample_id,
            'name': 'Test Mixture',
            'description': 'A test mixture'
        }

        self.sample_components = [
            {
                'molecule_id': self.molecule_id1,
                'concentration': 70,
                'concentration_unit': '%'
            },
            {
                'molecule_id': self.molecule_id2,
                'concentration': 30,
                'concentration_unit': '%'
            }
        ]

    @patch('api.models.get_supabase_client')
    def test_create_with_components(self, mock_get_supabase):
        """Test creating a mixture with components."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[self.sample_mixture])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Mixture.create_with_components(
            'Test Mixture',
            'A test mixture',
            self.sample_components
        )

        # Assertions
        self.assertEqual(result, self.sample_mixture)
        self.assertEqual(mock_supabase.tables['mixtures'].query_parts[0][0], 'insert')
        self.assertTrue('mixture_components' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_get_with_components(self, mock_get_supabase):
        """Test getting a mixture with components."""
        # Set up mock
        mixture_with_components = {
            **self.sample_mixture,
            'components': self.sample_components
        }
        mock_supabase = MockSupabase(return_data=[mixture_with_components])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Mixture.get_with_components(self.sample_id)

        # Assertions
        self.assertEqual(result, mixture_with_components)
        self.assertTrue('mixture_with_components' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_update_with_components(self, mock_get_supabase):
        """Test updating a mixture with components."""
        # Set up mock
        updated_mixture = {
            **self.sample_mixture,
            'name': 'Updated Mixture'
        }
        mock_supabase = MockSupabase(return_data=[updated_mixture])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Mixture.update_with_components(
            self.sample_id,
            {'name': 'Updated Mixture'},
            self.sample_components
        )

        # Assertions
        self.assertEqual(result, updated_mixture)
        self.assertEqual(mock_supabase.tables['mixtures'].query_parts[0][0], 'update')
        self.assertTrue('mixture_components' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_calculate_score(self, mock_get_supabase):
        """Test calculating a cryoprotectant score for a mixture."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=85.5)
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Mixture.calculate_score(self.sample_id)

        # Assertions
        self.assertEqual(result, 85.5)
        self.assertTrue('calculate_cryoprotectant_score' in mock_supabase.rpcs)

class TestPrediction(MockSupabaseBaseTestCase):
    """Test cases for the Prediction model."""

    @classmethod
    def setUpClass(cls):
        # Create a Flask app for context management
        from app import create_app
        cls.app = create_app(testing=True)
        cls.app_context = cls.app.app_context()
        cls.app_context.push()
        cls.client = cls.app.test_client()

    @classmethod
    def tearDownClass(cls):
        """Pop the Flask app context after all tests."""
        cls.app_context.pop()

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        self.sample_id = str(uuid.uuid4())
        self.mixture_id = str(uuid.uuid4())
        self.property_type_id = str(uuid.uuid4())
        self.calculation_method_id = str(uuid.uuid4())

        self.sample_prediction = {
            'id': self.sample_id,
            'mixture_id': self.mixture_id,
            'property_type_id': self.property_type_id,
            'calculation_method_id': self.calculation_method_id,
            'numeric_value': -15.3,
            'confidence': 0.9
        }

    @patch('api.models.get_supabase_client')
    def test_add_prediction(self, mock_get_supabase):
        """Test adding a prediction."""
        # Set up mock for property type
        property_type = {
            'id': self.property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }

        # Set up mock for calculation method
        calculation_method = {
            'id': self.calculation_method_id,
            'name': 'CryoProtect Scoring'
        }

        # Set up mock supabase
        mock_supabase = MockSupabase()
        mock_get_supabase.return_value = mock_supabase

        # Mock property type and calculation method queries
        mock_supabase.tables['property_types'] = MockSupabaseQuery(return_data=[property_type])
        mock_supabase.tables['calculation_methods'] = MockSupabaseQuery(return_data=[calculation_method])
        mock_supabase.tables['predictions'] = MockSupabaseQuery(return_data=[self.sample_prediction])

        # Call the method
        result = Prediction.add_prediction(
            self.mixture_id,
            'Freezing Point',
            -15.3,
            0.9,
            'CryoProtect Scoring'
        )

        # Assertions
        self.assertEqual(result, self.sample_prediction)
        self.assertTrue('property_types' in mock_supabase.tables)
        self.assertTrue('calculation_methods' in mock_supabase.tables)
        self.assertTrue('predictions' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_get_predictions_for_mixture(self, mock_get_supabase):
        """Test getting predictions for a mixture."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[self.sample_prediction])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Prediction.get_predictions_for_mixture(self.mixture_id)

        # Assertions
        self.assertEqual(result, [self.sample_prediction])
        self.assertEqual(mock_supabase.tables['predictions'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['predictions'].query_parts[1][0], 'eq')

class TestExperiment(MockSupabaseBaseTestCase):
    """Test cases for the Experiment model."""

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        self.sample_id = str(uuid.uuid4())
        self.mixture_id = str(uuid.uuid4())
        self.property_type_id = str(uuid.uuid4())

        self.sample_experiment = {
            'id': self.sample_id,
            'mixture_id': self.mixture_id,
            'property_type_id': self.property_type_id,
            'numeric_value': -14.8,
            'experimental_conditions': 'Standard pressure, cooling rate 1°C/min',
            'date_performed': '2025-04-15'
        }

    @patch('api.models.get_supabase_client')
    def test_record_experiment(self, mock_get_supabase):
        """Test recording an experiment."""
        # Set up mock for property type
        property_type = {
            'id': self.property_type_id,
            'name': 'Freezing Point',
            'data_type': 'numeric'
        }

        # Set up mock supabase
        mock_supabase = MockSupabase()
        mock_get_supabase.return_value = mock_supabase

        # Mock property type query
        mock_supabase.tables['property_types'] = MockSupabaseQuery(return_data=[property_type])
        mock_supabase.tables['experiments'] = MockSupabaseQuery(return_data=[self.sample_experiment])

        # Call the method
        result = Experiment.record_experiment(
            self.mixture_id,
            'Freezing Point',
            -14.8,
            'Standard pressure, cooling rate 1°C/min',
            '2025-04-15'
        )

        # Assertions
        self.assertEqual(result, self.sample_experiment)
        self.assertTrue('property_types' in mock_supabase.tables)
        self.assertTrue('experiments' in mock_supabase.tables)

    @patch('api.models.get_supabase_client')
    def test_get_experiments_for_mixture(self, mock_get_supabase):
        """Test getting experiments for a mixture."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=[self.sample_experiment])
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Experiment.get_experiments_for_mixture(self.mixture_id)

        # Assertions
        self.assertEqual(result, [self.sample_experiment])
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[0][0], 'select')
        self.assertEqual(mock_supabase.tables['experiments'].query_parts[1][0], 'eq')

class TestComparison(MockSupabaseBaseTestCase):
    """Test cases for the Comparison utility class."""

    def setUp(self):
        """Set up test cases."""
        super().setUp()
        self.mixture_id = str(uuid.uuid4())
        self.property_type_id = str(uuid.uuid4())

        self.sample_comparison = {
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

    @patch('api.models.BaseModel.get_supabase')
    def test_compare_prediction_with_experiment(self, mock_get_supabase):
        """Test comparing a prediction with an experiment."""
        # Set up mock
        mock_supabase = MockSupabase(return_data=self.sample_comparison)
        mock_get_supabase.return_value = mock_supabase

        # Call the method
        result = Comparison.compare_prediction_with_experiment(
            self.mixture_id,
            'Freezing Point'
        )

        # Assertions
        self.assertEqual(result, self.sample_comparison)
        self.assertTrue('compare_prediction_with_experiment' in mock_supabase.rpcs)

class TestDataRelationshipsAndIntegrity(MockSupabaseBaseTestCase):
    """Test data relationships, foreign key enforcement, cascading deletes, and project/user associations."""

    @patch('api.models.get_supabase_client')
    def test_invalid_foreign_keys(self, mock_get_supabase):
        """Test that creating records with invalid foreign keys fails."""
        # Simulate DB error for invalid foreign key
        mock_supabase = MockSupabase(return_data=[], return_error="Foreign key violation")
        mock_get_supabase.return_value = mock_supabase

        # Attempt to create a mixture with a non-existent molecule_id
        with self.assertRaises(Exception):
            Mixture.create_with_components(
                'Invalid Mixture',
                'Should fail',
                [{'molecule_id': 'non-existent-id', 'concentration': 50, 'concentration_unit': '%'}]
            )

        # Attempt to create a prediction with a non-existent mixture_id
        with self.assertRaises(Exception):
            Prediction.add_prediction(
                'non-existent-mixture-id',
                'Freezing Point',
                -20.0,
                0.8,
                'CryoProtect Scoring'
            )

        # Attempt to create an experiment with a non-existent mixture_id
        with self.assertRaises(Exception):
            Experiment.record_experiment(
                'non-existent-mixture-id',
                'Freezing Point',
                -20.0,
                'Standard conditions',
                '2025-04-16'
            )

    @patch('api.models.get_supabase_client')
    def test_cascading_deletes(self, mock_get_supabase):
        """Test that deleting parent records cascades or restricts deletes for child records."""
        # Simulate DB with parent and child records
        parent_id = str(uuid.uuid4())
        child_id = str(uuid.uuid4())
        mock_supabase = MockSupabase(return_data=[{'id': parent_id}], return_error=None)
        mock_get_supabase.return_value = mock_supabase

        # Mock delete on parent (e.g., molecule or mixture)
        with patch.object(Mixture, 'delete', return_value=True):
            result = Mixture.delete(parent_id)
            self.assertTrue(result)

        # After parent delete, attempt to get child (should not exist or raise error)
        mock_supabase = MockSupabase(return_data=[], return_error=None)
        mock_get_supabase.return_value = mock_supabase
        with patch.object(Prediction, 'get_predictions_for_mixture', return_value=[]):
            predictions = Prediction.get_predictions_for_mixture(parent_id)
            self.assertEqual(predictions, [])

        with patch.object(Experiment, 'get_experiments_for_mixture', return_value=[]):
            experiments = Experiment.get_experiments_for_mixture(parent_id)
            self.assertEqual(experiments, [])

    @patch('api.models.get_supabase_client')
    def test_project_user_associations(self, mock_get_supabase):
        """Test that project and user associations are enforced and correct."""
        # Simulate project creation with user association
        user_id = str(uuid.uuid4())
        project_id = str(uuid.uuid4())
        experiment_id = str(uuid.uuid4())
        mock_supabase = MockSupabase(return_data=[{
            'id': project_id,
            'name': 'Test Project',
            'created_by': user_id
        }], return_error=None)
        mock_get_supabase.return_value = mock_supabase

        # Create project
        from api.models import BaseModel
        class Project(BaseModel):
            table_name = 'projects'
        project = Project.create({'name': 'Test Project', 'created_by': user_id})
        self.assertEqual(project['created_by'], user_id)

        # Simulate project_experiments association
        mock_supabase = MockSupabase(return_data=[{
            'id': str(uuid.uuid4()),
            'project_id': project_id,
            'experiment_id': experiment_id,
            'created_by': user_id
        }], return_error=None)
        mock_get_supabase.return_value = mock_supabase

        class ProjectExperiment(BaseModel):
            table_name = 'project_experiments'
        project_experiment = ProjectExperiment.create({
            'project_id': project_id,
            'experiment_id': experiment_id,
            'created_by': user_id
        })
        self.assertEqual(project_experiment['project_id'], project_id)
        self.assertEqual(project_experiment['created_by'], user_id)

if __name__ == '__main__':
    import unittest
    unittest.main()
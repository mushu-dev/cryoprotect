"""
CryoProtect Analyzer - Database Health Check Tests

This module contains tests for the database health check utilities and API endpoints.
It tests the SchemaValidator, IntegrityChecker, PerformanceAnalyzer, and DatabaseHealthCheck
classes, as well as the health check API endpoints.
"""

import os
import sys
import json
import uuid
import pytest
from unittest.mock import patch, MagicMock, PropertyMock
from datetime import datetime, timedelta

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase, mock_rpc_function

# Import the health check modules
from database.utils.health_check import DatabaseHealthCheck
from database.utils.schema_validator import SchemaValidator
from database.utils.integrity_checker import IntegrityChecker
from database.utils.performance_analyzer import PerformanceAnalyzer

# Import API resources
from api.system_resources import (
    HealthResource, HealthDatabaseResource, HealthPerformanceResource
)

class TestSchemaValidator(MockSupabaseBaseTestCase):
    """Test cases for the SchemaValidator class."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Create a SchemaValidator instance
        self.validator = SchemaValidator()
        
        # Mock the connection pool
        self.mock_conn = MagicMock()
        
        # Set up mock responses for common queries
        self.setup_mock_responses()
    
    def setup_mock_responses(self):
        """Set up mock responses for database queries."""
        # Mock response for table existence check
        table_exists_response = MagicMock()
        table_exists_response.data = [{'exists': True}]
        table_exists_response.error = None
        self.mock_conn.rpc.return_value.execute.return_value = table_exists_response
        
        # Register mock RPC functions
        mock_rpc_function('exec_sql', [{'exists': True}])
        
        # Mock response for column definitions
        columns_response = MagicMock()
        columns_response.data = [
            {'column_name': 'id', 'data_type': 'uuid', 'is_nullable': 'NO'},
            {'column_name': 'name', 'data_type': 'character varying', 'is_nullable': 'NO'},
            {'column_name': 'smiles', 'data_type': 'character varying', 'is_nullable': 'YES'},
            {'column_name': 'created_at', 'data_type': 'timestamp with time zone', 'is_nullable': 'YES'},
            {'column_name': 'updated_at', 'data_type': 'timestamp with time zone', 'is_nullable': 'YES'}
        ]
        columns_response.error = None
        
        # Mock response for constraints
        constraints_response = MagicMock()
        constraints_response.data = [
            {'constraint_name': 'pk_molecules', 'constraint_type': 'PRIMARY KEY', 'table_name': 'molecules', 'column_name': 'id'},
            {'constraint_name': 'fk_mixture_components_molecule', 'constraint_type': 'FOREIGN KEY', 'table_name': 'mixture_components', 'column_name': 'molecule_id', 'referenced_table': 'molecules', 'referenced_column': 'id'}
        ]
        constraints_response.error = None
        
        # Mock response for views
        views_response = MagicMock()
        views_response.data = [
            {'view_name': 'molecule_details', 'view_definition': 'SELECT m.id, m.name, m.smiles, mp.property_name, mp.numeric_value FROM molecules m JOIN molecular_properties mp ON m.id = mp.molecule_id'}
        ]
        views_response.error = None
        
        # Mock response for RLS policies
        rls_response = MagicMock()
        rls_response.data = [
            {'table_name': 'molecules', 'rls_enabled': True},
            {'table_name': 'mixtures', 'rls_enabled': True},
            {'table_name': 'mixture_components', 'rls_enabled': True}
        ]
        rls_response.error = None
    
    @patch_supabase(load_data=True)
    def test_run_checks(self, mock_client):
        """Test the run_checks method of SchemaValidator."""
        # Run the checks
        results = self.validator.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIsInstance(results, dict)
        self.assertIn('status', results)
        self.assertIn('issues_found', results)
        self.assertIn('details', results)
        self.assertIn('recommendations', results)
    
    @patch_supabase(load_data=True)
    def test_detect_missing_table(self, mock_client):
        """Test detection of missing tables."""
        # Mock the table existence check to return False
        table_exists_response = MagicMock()
        table_exists_response.data = [{'exists': False}]
        table_exists_response.error = None
        self.mock_conn.rpc.return_value.execute.return_value = table_exists_response
        
        # Run the checks
        results = self.validator.run_checks(self.mock_conn)
        
        # Assertions
        self.assertEqual(results['status'], 'failed')
        self.assertGreater(results['issues_found'], 0)
    
    @patch_supabase(load_data=True)
    def test_detect_column_type_mismatch(self, mock_client):
        """Test detection of column type mismatches."""
        # Mock the column definitions to have a type mismatch
        def mock_execute(*args, **kwargs):
            if 'information_schema.tables' in args[0]['query']:
                response = MagicMock()
                response.data = [{'exists': True}]
                response.error = None
                return response
            elif 'information_schema.columns' in args[0]['query']:
                response = MagicMock()
                response.data = [
                    {'column_name': 'id', 'data_type': 'uuid', 'is_nullable': 'NO'},
                    {'column_name': 'name', 'data_type': 'integer', 'is_nullable': 'NO'},  # Type mismatch
                    {'column_name': 'smiles', 'data_type': 'character varying', 'is_nullable': 'YES'},
                    {'column_name': 'created_at', 'data_type': 'timestamp with time zone', 'is_nullable': 'YES'},
                    {'column_name': 'updated_at', 'data_type': 'timestamp with time zone', 'is_nullable': 'YES'}
                ]
                response.error = None
                return response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
        
        # Run the checks
        results = self.validator.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIn('status', results)
        self.assertGreater(results['issues_found'], 0)


class TestIntegrityChecker(MockSupabaseBaseTestCase):
    """Test cases for the IntegrityChecker class."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Create an IntegrityChecker instance
        self.checker = IntegrityChecker()
        
        # Mock the connection pool
        self.mock_conn = MagicMock()
        
        # Set up mock responses for common queries
        self.setup_mock_responses()
    
    def setup_mock_responses(self):
        """Set up mock responses for database queries."""
        # Mock response for orphaned records check
        orphaned_response = MagicMock()
        orphaned_response.data = []  # No orphaned records
        orphaned_response.error = None
        
        # Mock response for unique constraint check
        unique_response = MagicMock()
        unique_response.data = []  # No duplicate records
        unique_response.error = None
        
        # Mock response for required fields check
        required_response = MagicMock()
        required_response.data = []  # No null required fields
        required_response.error = None
        
        # Mock response for data consistency check
        consistency_response = MagicMock()
        consistency_response.data = []  # No consistency issues
        consistency_response.error = None
        
        # Set up the mock to return different responses based on the query
        def mock_execute(*args, **kwargs):
            if 'LEFT JOIN' in args[0]['query'] and 'IS NULL' in args[0]['query']:
                return orphaned_response
            elif 'GROUP BY' in args[0]['query'] and 'HAVING COUNT' in args[0]['query']:
                return unique_response
            elif 'IS NULL' in args[0]['query'] and 'WHERE' in args[0]['query']:
                return required_response
            elif 'concentration <= 0' in args[0]['query'] or 'numeric_value IS NOT NULL' in args[0]['query']:
                return consistency_response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
    
    @patch_supabase(load_data=True)
    def test_run_checks(self, mock_client):
        """Test the run_checks method of IntegrityChecker."""
        # Run the checks
        results = self.checker.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIsInstance(results, dict)
        self.assertIn('status', results)
        self.assertIn('issues_found', results)
        self.assertIn('details', results)
        self.assertIn('recommendations', results)
    
    @patch_supabase(load_data=True)
    def test_detect_orphaned_records(self, mock_client):
        """Test detection of orphaned records."""
        # Mock the orphaned records check to return orphaned records
        def mock_execute(*args, **kwargs):
            if 'LEFT JOIN' in args[0]['query'] and 'IS NULL' in args[0]['query']:
                response = MagicMock()
                response.data = [
                    {'id': str(uuid.uuid4()), 'mixture_id': str(uuid.uuid4())}
                ]  # Orphaned record
                response.error = None
                return response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
        
        # Run the checks
        results = self.checker.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIn('status', results)
        self.assertGreater(results['issues_found'], 0)
    
    @patch_supabase(load_data=True)
    def test_detect_duplicate_records(self, mock_client):
        """Test detection of duplicate records (unique constraint violations)."""
        # Mock the unique constraint check to return duplicate records
        def mock_execute(*args, **kwargs):
            if 'GROUP BY' in args[0]['query'] and 'HAVING COUNT' in args[0]['query']:
                response = MagicMock()
                response.data = [
                    {'smiles': 'C(C(CO)O)O', 'count': 2}
                ]  # Duplicate record
                response.error = None
                return response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
        
        # Run the checks
        results = self.checker.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIn('status', results)
        self.assertGreater(results['issues_found'], 0)


class TestPerformanceAnalyzer(MockSupabaseBaseTestCase):
    """Test cases for the PerformanceAnalyzer class."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Create a PerformanceAnalyzer instance
        self.analyzer = PerformanceAnalyzer()
        
        # Mock the connection pool
        self.mock_conn = MagicMock()
        
        # Set up mock responses for common queries
        self.setup_mock_responses()
    
    def setup_mock_responses(self):
        """Set up mock responses for database queries."""
        # Mock response for query execution time
        query_time_response = MagicMock()
        query_time_response.data = [{'result': 'success', 'execution_time_ms': 50}]
        query_time_response.error = None
        
        # Mock response for index usage
        index_usage_response = MagicMock()
        index_usage_response.data = [
            {'table_name': 'molecules', 'index_name': 'idx_molecules_id', 'column_name': 'id', 'usage_count': 100},
            {'table_name': 'molecules', 'index_name': 'idx_molecules_name', 'column_name': 'name', 'usage_count': 50}
        ]
        index_usage_response.error = None
        
        # Mock response for table statistics
        table_stats_response = MagicMock()
        table_stats_response.data = [
            {'table_name': 'molecules', 'row_count': 1000, 'total_size': 1024000, 'index_size': 512000}
        ]
        table_stats_response.error = None
        
        # Set up the mock to return different responses based on the query
        def mock_execute(*args, **kwargs):
            if any(q['name'] in args[0]['query'] for q in self.analyzer.test_queries):
                return query_time_response
            elif 'pg_stat_user_indexes' in args[0]['query']:
                return index_usage_response
            elif 'pg_stat_user_tables' in args[0]['query'] or 'pg_class' in args[0]['query']:
                return table_stats_response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
    
    @patch_supabase(load_data=True)
    def test_run_checks(self, mock_client):
        """Test the run_checks method of PerformanceAnalyzer."""
        # Run the checks
        results = self.analyzer.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIsInstance(results, dict)
        self.assertIn('status', results)
        self.assertIn('issues_found', results)
        self.assertIn('details', results)
        self.assertIn('recommendations', results)
    
    @patch_supabase(load_data=True)
    def test_detect_slow_queries(self, mock_client):
        """Test detection of slow queries."""
        # Mock the query execution time to be slow
        def mock_execute(*args, **kwargs):
            if any(q['name'] in args[0]['query'] for q in self.analyzer.test_queries):
                response = MagicMock()
                response.data = [{'result': 'success', 'execution_time_ms': 500}]  # Slow query
                response.error = None
                return response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
        
        # Run the checks
        results = self.analyzer.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIn('status', results)
        self.assertGreater(results['issues_found'], 0)


class TestDatabaseHealthCheck(MockSupabaseBaseTestCase):
    """Test cases for the DatabaseHealthCheck class."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Mock the validator classes
        self.mock_schema_validator = MagicMock()
        self.mock_integrity_checker = MagicMock()
        self.mock_performance_analyzer = MagicMock()
        
        # Set up mock responses for validators
        self.setup_mock_responses()
        
        # Create a DatabaseHealthCheck instance with mocked validators
        with patch('database.utils.health_check.SchemaValidator', return_value=self.mock_schema_validator), \
             patch('database.utils.health_check.IntegrityChecker', return_value=self.mock_integrity_checker), \
             patch('database.utils.health_check.PerformanceAnalyzer', return_value=self.mock_performance_analyzer):
            self.health_check = DatabaseHealthCheck(MagicMock())
    
    def setup_mock_responses(self):
        """Set up mock responses for validators."""
        # Mock schema validator response
        self.mock_schema_validator.run_checks.return_value = {
            'status': 'passed',
            'issues_found': 0,
            'details': {'tables': {'status': 'passed'}, 'columns': {'status': 'passed'}},
            'recommendations': []
        }
        
        # Mock integrity checker response
        self.mock_integrity_checker.run_checks.return_value = {
            'status': 'passed',
            'issues_found': 0,
            'details': {'orphaned_records': {'status': 'passed'}, 'unique_constraints': {'status': 'passed'}},
            'recommendations': []
        }
        
        # Mock performance analyzer response
        self.mock_performance_analyzer.run_checks.return_value = {
            'status': 'passed',
            'issues_found': 0,
            'details': {'query_times': {'status': 'passed'}, 'indexes': {'status': 'passed'}},
            'recommendations': []
        }
    
    def test_run_health_check_all_categories(self):
        """Test running health check for all categories."""
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Assertions
        self.assertIsInstance(results, dict)
        self.assertIn('overall_status', results)
        self.assertIn('categories_checked', results)
        self.assertIn('total_issues_found', results)
        self.assertIn('results_by_category', results)
        self.assertIn('recommendations', results)
        
        # Check that all validators were called
        self.mock_schema_validator.run_checks.assert_called_once()
        self.mock_integrity_checker.run_checks.assert_called_once()
        self.mock_performance_analyzer.run_checks.assert_called_once()
    
    def test_run_health_check_specific_category(self):
        """Test running health check for a specific category."""
        # Run the health check for schema only
        results = self.health_check.run_health_check(categories=['schema'])
        
        # Assertions
        self.assertIsInstance(results, dict)
        self.assertIn('overall_status', results)
        self.assertEqual(results['categories_checked'], ['schema'])
        
        # Check that only the schema validator was called
        self.mock_schema_validator.run_checks.assert_called_once()
        self.mock_integrity_checker.run_checks.assert_not_called()
        self.mock_performance_analyzer.run_checks.assert_not_called()
    
    def test_calculate_status_all_passed(self):
        """Test status calculation when all checks pass."""
        # Run the health check with all passed
        results = self.health_check.run_health_check()
        
        # Assertions
        self.assertEqual(results['overall_status'], 'passed')
    
    def test_calculate_status_with_warning(self):
        """Test status calculation when there's a warning."""
        # Set up mock responses with a warning
        self.mock_integrity_checker.run_checks.return_value = {
            'status': 'warning',
            'issues_found': 1,
            'details': {'orphaned_records': {'status': 'warning'}},
            'recommendations': ['Fix orphaned records']
        }
        
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Assertions
        self.assertEqual(results['overall_status'], 'warning')
        self.assertEqual(results['total_issues_found'], 1)
        self.assertIn('Fix orphaned records', results['recommendations'])
    
    def test_calculate_status_with_failure(self):
        """Test status calculation when there's a failure."""
        # Set up mock responses with a failure
        self.mock_schema_validator.run_checks.return_value = {
            'status': 'failed',
            'issues_found': 2,
            'details': {'tables': {'status': 'failed'}},
            'recommendations': ['Create missing tables']
        }
        
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Assertions
        self.assertEqual(results['overall_status'], 'failed')
        self.assertEqual(results['total_issues_found'], 2)
        self.assertIn('Create missing tables', results['recommendations'])
    
    def test_generate_report_json(self):
        """Test generating a JSON report."""
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Generate a JSON report
        report = self.health_check.generate_report(results, format='json')
        
        # Assertions
        self.assertIsInstance(report, dict)
        self.assertIn('overall_status', report)
    
    def test_generate_report_markdown(self):
        """Test generating a Markdown report."""
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Generate a Markdown report
        report = self.health_check.generate_report(results, format='markdown')
        
        # Assertions
        self.assertIsInstance(report, str)
        self.assertIn('# Database Health Check Report', report)
    
    def test_generate_report_html(self):
        """Test generating an HTML report."""
        # Run the health check
        results = self.health_check.run_health_check()
        
        # Generate an HTML report
        report = self.health_check.generate_report(results, format='html')
        
        # Assertions
        self.assertIsInstance(report, str)
        self.assertIn('<!DOCTYPE html>', report)
        self.assertIn('<title>Database Health Check Report</title>', report)


class TestHealthCheckAPI(MockSupabaseBaseTestCase):
    """Test cases for the health check API endpoints."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()
        
        # Mock the DatabaseHealthCheck class
        self.mock_health_check_patcher = patch('api.system_resources.DatabaseHealthCheck')
        self.mock_health_check_class = self.mock_health_check_patcher.start()
        self.mock_health_check = self.mock_health_check_class.return_value
        
        # Set up mock responses
        self.setup_mock_responses()
    
    def tearDown(self):
        """Tear down test case."""
        self.mock_health_check_patcher.stop()
        super().tearDown()
    
    def setup_mock_responses(self):
        """Set up mock responses for health check."""
        # Mock health check response
        self.mock_health_check.run_health_check.return_value = {
            'timestamp': datetime.now().isoformat(),
            'overall_status': 'passed',
            'categories_checked': ['schema', 'integrity', 'performance'],
            'total_issues_found': 0,
            'results_by_category': {
                'schema': {'status': 'passed', 'issues_found': 0},
                'integrity': {'status': 'passed', 'issues_found': 0},
                'performance': {'status': 'passed', 'issues_found': 0}
            },
            'recommendations': []
        }
        
        # Mock report generation
        self.mock_health_check.generate_report.return_value = "# Database Health Check Report\n\nAll checks passed."
    
    def test_health_endpoint(self):
        """Test the /health endpoint."""
        # Make request
        response = self.client.get('/health')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('status', data)
        self.assertIn('version', data)
        self.assertIn('uptime', data)
        self.assertIn('database_status', data)
    
    @patch('api.system_resources.token_required', lambda f: f)  # Mock token_required decorator
    @patch('api.system_resources.require_admin', lambda f: f)   # Mock require_admin decorator
    def test_health_database_endpoint(self):
        """Test the /health/database endpoint."""
        # Make request
        response = self.client.get('/health/database')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('status', data)
        self.assertIn('overall_health', data)
        self.assertIn('schema_status', data)
        self.assertIn('integrity_status', data)
        self.assertIn('performance_status', data)
        
        # Check that DatabaseHealthCheck.run_health_check was called
        self.mock_health_check.run_health_check.assert_called_once()
    
    @patch('api.system_resources.token_required', lambda f: f)  # Mock token_required decorator
    @patch('api.system_resources.require_admin', lambda f: f)   # Mock require_admin decorator
    def test_health_database_endpoint_with_verbosity(self):
        """Test the /health/database endpoint with verbosity parameter."""
        # Make request with minimal verbosity
        response = self.client.get('/health/database?verbosity=minimal')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        
        # Check that DatabaseHealthCheck.run_health_check was called with correct parameters
        self.mock_health_check.run_health_check.assert_called_with(categories=['schema'])
        
        # Reset mock
        self.mock_health_check.run_health_check.reset_mock()
        
        # Make request with detailed verbosity
        response = self.client.get('/health/database?verbosity=detailed')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        
        # Check that DatabaseHealthCheck.run_health_check was called with no parameters
        self.mock_health_check.run_health_check.assert_called_with()
    
    @patch('api.system_resources.token_required', lambda f: f)  # Mock token_required decorator
    @patch('api.system_resources.require_admin', lambda f: f)   # Mock require_admin decorator
    def test_health_database_endpoint_with_format(self):
        """Test the /health/database endpoint with format parameter."""
        # Make request with markdown format
        response = self.client.get('/health/database?format=markdown')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.mimetype, 'text/markdown')
        
        # Check that DatabaseHealthCheck.generate_report was called with correct parameters
        self.mock_health_check.generate_report.assert_called_once()
        
        # Reset mocks
        self.mock_health_check.generate_report.reset_mock()
        
        # Make request with HTML format
        response = self.client.get('/health/database?format=html')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.mimetype, 'text/html')
        
        # Check that DatabaseHealthCheck.generate_report was called with correct parameters
        self.mock_health_check.generate_report.assert_called_once()
    
    @patch('api.system_resources.token_required', lambda f: f)  # Mock token_required decorator
    @patch('api.system_resources.require_admin', lambda f: f)   # Mock require_admin decorator
    @patch('api.system_resources.psutil.cpu_percent', return_value=25.0)
    @patch('api.system_resources.psutil.virtual_memory')
    def test_health_performance_endpoint(self, mock_memory, mock_cpu):
        """Test the /health/performance endpoint."""
        # Set up mock memory
        mock_memory.return_value = MagicMock(
            total=16000000000,
            available=8000000000,
            used=8000000000,
            percent=50.0
        )
        
        # Make request
        response = self.client.get('/health/performance')
        
        # Assertions
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('status', data)
        self.assertIn('metrics', data)
        self.assertIn('database_metrics', data)
        self.assertIn('api_metrics', data)
        
        # Check that DatabaseHealthCheck.run_health_check was called with correct parameters
        self.mock_health_check.run_health_check.assert_called_with(categories=['performance'])
    
    @patch('api.system_resources.token_required', lambda f: f)  # Mock token_required decorator
    @patch('api.system_resources.require_admin', lambda f: f)   # Mock require_admin decorator
    def test_authentication_required(self):
        """Test that authentication is required for protected endpoints."""
        # Remove the mocks for token_required and require_admin
        with patch('api.system_resources.token_required', wraps=lambda f: f) as mock_token_required, \
             patch('api.system_resources.require_admin', wraps=lambda f: f) as mock_require_admin:
            
            # Make request to protected endpoint
            response = self.client.get('/health/database')
            
            # Check that the decorators were called
            self.assertTrue(mock_token_required.called)
            self.assertTrue(mock_require_admin.called)
    
    def test_error_handling(self):
        """Test error handling in health check endpoints."""
        # Mock the health check to raise an exception
        self.mock_health_check.run_health_check.side_effect = Exception("Test error")
        
        # Make request to health endpoint (which doesn't use DatabaseHealthCheck)
        response = self.client.get('/health')
        
        # Assertions - should still return 200 with degraded status
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.data)
        self.assertIn('status', data)
        
        # Make request to database health endpoint
        with patch('api.system_resources.token_required', lambda f: f), \
             patch('api.system_resources.require_admin', lambda f: f), \
             patch('api.system_resources.handle_error') as mock_handle_error:
            
            # Mock handle_error to return a response
            mock_handle_error.return_value = ({'error': 'Test error'}, 500)
            
            # Make request
            response = self.client.get('/health/database')
            
            # Check that handle_error was called
            self.assertTrue(mock_handle_error.called)


if __name__ == '__main__':
    unittest.main()
    
    @patch_supabase(load_data=True)
    def test_detect_missing_indexes(self, mock_client):
        """Test detection of missing indexes."""
        # Mock the index usage check to return missing indexes
        def mock_execute(*args, **kwargs):
            if 'pg_stat_user_indexes' in args[0]['query']:
                response = MagicMock()
                response.data = [
                    {'table_name': 'molecules', 'index_name': 'idx_molecules_id', 'column_name': 'id', 'usage_count': 100}
                    # Missing index for 'name'
                ]
                response.error = None
                return response
            else:
                response = MagicMock()
                response.data = []
                response.error = None
                return response
        
        self.mock_conn.rpc.return_value.execute.side_effect = mock_execute
        
        # Run the checks
        results = self.analyzer.run_checks(self.mock_conn)
        
        # Assertions
        self.assertIn('status', results)
        self.assertGreater(results['issues_found'], 0)

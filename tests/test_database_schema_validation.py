"""
CryoProtect Analyzer - Database Schema Validation Tests

This module contains tests to validate the database schema against the requirements.
It verifies table names, relationships, constraints, and RLS policies.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the base test case
from tests.base_test_case import MockSupabaseBaseTestCase
from tests.mock_supabase.helpers import patch_supabase

class TestDatabaseSchemaValidation(MockSupabaseBaseTestCase):
    """Test cases for database schema validation."""

    def setUp(self):
        """Set up test data for each test."""
        # Call the parent setUp to initialize the mock Supabase
        super().setUp()

    @patch_supabase(load_data=True)
    def test_table_names_are_plural(self, mock_client):
        """Test that all table names use plural form."""
        # Get all tables from the database
        response = mock_client.rpc('get_all_tables').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting tables: {response.error}")
        
        # Get the table names
        tables = response.data if response.data else []
        
        # Define expected plural table names
        expected_plural_tables = [
            'molecules',
            'property_types',
            'molecular_properties',
            'mixtures',
            'mixture_components',
            'calculation_methods',
            'predictions',
            'experiments',
            'projects',
            'teams',
            'user_profiles'
        ]
        
        # Check that all expected tables exist
        for table in expected_plural_tables:
            self.assertIn(table, tables, f"Table '{table}' not found in database")
        
        # Check that all table names are plural
        singular_tables = []
        for table in tables:
            # Skip system tables and views
            if table.startswith('pg_') or table.startswith('information_schema.'):
                continue
            
            # Check if the table name is plural
            if table.endswith('s') or table.endswith('es'):
                continue
            
            singular_tables.append(table)
        
        self.assertEqual(len(singular_tables), 0, 
                         f"Found tables with singular names: {singular_tables}")

    @patch_supabase(load_data=True)
    def test_foreign_key_constraints(self, mock_client):
        """Test that all relationships have proper FK constraints."""
        # Get all foreign key constraints from the database
        response = mock_client.rpc('get_all_foreign_keys').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting foreign keys: {response.error}")
        
        # Get the foreign keys
        foreign_keys = response.data if response.data else []
        
        # Define expected foreign key relationships
        expected_fk_relationships = [
            ('molecular_properties', 'molecule_id', 'molecules', 'id'),
            ('molecular_properties', 'property_type_id', 'property_types', 'id'),
            ('mixture_components', 'mixture_id', 'mixtures', 'id'),
            ('mixture_components', 'molecule_id', 'molecules', 'id'),
            ('predictions', 'molecule_id', 'molecules', 'id'),
            ('predictions', 'mixture_id', 'mixtures', 'id'),
            ('predictions', 'property_type_id', 'property_types', 'id'),
            ('predictions', 'calculation_method_id', 'calculation_methods', 'id'),
            ('experiments', 'mixture_id', 'mixtures', 'id'),
            ('experiment_properties', 'experiment_id', 'experiments', 'id'),
            ('experiment_properties', 'property_type_id', 'property_types', 'id'),
            ('projects', 'team_id', 'teams', 'id'),
            ('user_profiles', 'team_id', 'teams', 'id')
        ]
        
        # Check that all expected foreign key relationships exist
        missing_fks = []
        for fk in expected_fk_relationships:
            child_table, child_column, parent_table, parent_column = fk
            
            # Find the foreign key in the list
            found = False
            for foreign_key in foreign_keys:
                if (foreign_key['child_table'] == child_table and
                    foreign_key['child_column'] == child_column and
                    foreign_key['parent_table'] == parent_table and
                    foreign_key['parent_column'] == parent_column):
                    found = True
                    break
            
            if not found:
                missing_fks.append(fk)
        
        self.assertEqual(len(missing_fks), 0, 
                         f"Missing foreign key constraints: {missing_fks}")

    @patch_supabase(load_data=True)
    def test_rls_enabled_on_all_tables(self, mock_client):
        """Test that RLS is enabled on all tables."""
        # Get all tables with RLS status from the database
        response = mock_client.rpc('get_tables_with_rls_status').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting RLS status: {response.error}")
        
        # Get the tables with RLS status
        tables_with_rls = response.data if response.data else []
        
        # Define tables that should have RLS enabled
        tables_requiring_rls = [
            'molecules',
            'property_types',
            'molecular_properties',
            'mixtures',
            'mixture_components',
            'calculation_methods',
            'predictions',
            'experiments',
            'projects',
            'teams',
            'user_profiles'
        ]
        
        # Check that RLS is enabled on all required tables
        tables_without_rls = []
        for table in tables_requiring_rls:
            found = False
            for table_with_rls in tables_with_rls:
                if table_with_rls['table_name'] == table:
                    if not table_with_rls['rls_enabled']:
                        tables_without_rls.append(table)
                    found = True
                    break
            
            if not found:
                tables_without_rls.append(table)
        
        self.assertEqual(len(tables_without_rls), 0, 
                         f"Tables without RLS enabled: {tables_without_rls}")

    @patch_supabase(load_data=True)
    def test_rls_policies_for_different_roles(self, mock_client):
        """Test RLS policies with different user roles."""
        # Get all RLS policies from the database
        response = mock_client.rpc('get_all_rls_policies').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting RLS policies: {response.error}")
        
        # Get the RLS policies
        rls_policies = response.data if response.data else []
        
        # Define expected RLS policies for different roles
        expected_policies = [
            # Anonymous users can view data
            {'table': 'molecules', 'operation': 'SELECT', 'role': 'anon', 'using': True},
            {'table': 'mixtures', 'operation': 'SELECT', 'role': 'anon', 'using': True},
            
            # Authenticated users can insert data
            {'table': 'molecules', 'operation': 'INSERT', 'role': 'authenticated', 'using': True},
            {'table': 'mixtures', 'operation': 'INSERT', 'role': 'authenticated', 'using': True},
            
            # Users can only update/delete their own data
            {'table': 'molecules', 'operation': 'UPDATE', 'role': 'authenticated', 'using': 'created_by = auth.uid()'},
            {'table': 'mixtures', 'operation': 'UPDATE', 'role': 'authenticated', 'using': 'created_by = auth.uid()'},
            {'table': 'molecules', 'operation': 'DELETE', 'role': 'authenticated', 'using': 'created_by = auth.uid()'},
            {'table': 'mixtures', 'operation': 'DELETE', 'role': 'authenticated', 'using': 'created_by = auth.uid()'}
        ]
        
        # Check that all expected RLS policies exist
        missing_policies = []
        for policy in expected_policies:
            found = False
            for rls_policy in rls_policies:
                if (rls_policy['table_name'] == policy['table'] and
                    rls_policy['operation'] == policy['operation'] and
                    rls_policy['role'] == policy['role']):
                    # For simple policies, just check that they exist
                    if policy['using'] is True:
                        found = True
                        break
                    # For complex policies, check the using expression
                    elif policy['using'] in rls_policy['using_expression']:
                        found = True
                        break
            
            if not found:
                missing_policies.append(policy)
        
        self.assertEqual(len(missing_policies), 0, 
                         f"Missing RLS policies: {missing_policies}")

    @patch_supabase(load_data=True)
    def test_views_exist(self, mock_client):
        """Test that required views exist."""
        # Get all views from the database
        response = mock_client.rpc('get_all_views').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting views: {response.error}")
        
        # Get the views
        views = response.data if response.data else []
        
        # Define expected views
        expected_views = [
            'molecule_with_properties',
            'mixture_with_components'
        ]
        
        # Check that all expected views exist
        for view in expected_views:
            found = False
            for db_view in views:
                if db_view['view_name'] == view:
                    found = True
                    break
            
            self.assertTrue(found, f"View '{view}' not found in database")

    @patch_supabase(load_data=True)
    def test_functions_exist(self, mock_client):
        """Test that required functions exist."""
        # Get all functions from the database
        response = mock_client.rpc('get_all_functions').execute()
        
        # Check for errors
        self.assertFalse(response.error, f"Error getting functions: {response.error}")
        
        # Get the functions
        functions = response.data if response.data else []
        
        # Define expected functions
        expected_functions = [
            'import_molecule_from_pubchem',
            'calculate_mixture_score',
            'compare_prediction_with_experiment'
        ]
        
        # Check that all expected functions exist
        for function in expected_functions:
            found = False
            for db_function in functions:
                if db_function['function_name'] == function:
                    found = True
                    break
            
            self.assertTrue(found, f"Function '{function}' not found in database")

if __name__ == '__main__':
    unittest.main()
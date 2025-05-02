#!/usr/bin/env python3
"""
CryoProtect v2 - Test Database Audit Functionality

This script tests the database resource verification functionality in supabase_database_audit.py.
It verifies that the script correctly checks for required tables, foreign key relationships,
and database quotas, and that it aborts execution when critical resources are missing.

Usage:
    python -m tests.test_database_audit
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import json
from pathlib import Path

# Add parent directory to path to import supabase_database_audit
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the module to test
import supabase_database_audit

class TestDatabaseAudit(unittest.TestCase):
    """Test cases for the database audit functionality."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock config
        self.mock_config = MagicMock()
        self.mock_config.SUPABASE_URL = "https://example.supabase.co"
        self.mock_config.SUPABASE_KEY = "mock-api-key"
        self.mock_config.SUPABASE_USER = "test@example.com"
        self.mock_config.SUPABASE_PASSWORD = "password"
        self.mock_config.__class__.__name__ = "TestingConfig"
        self.mock_config.as_dict.return_value = {
            "SUPABASE_URL": "https://example.supabase.co",
            "SUPABASE_KEY": "[REDACTED]",
            "LOG_LEVEL": "INFO"
        }

        # Create mock supabase client
        self.mock_supabase = MagicMock()
        
        # Ensure directories exist
        Path("logs").mkdir(exist_ok=True)
        Path("reports").mkdir(exist_ok=True)

    @patch('supabase_database_audit.BaseConfig.from_env')
    @patch('supabase_database_audit.create_client')
    def test_connect_to_supabase_success(self, mock_create_client, mock_from_env):
        """Test successful connection to Supabase."""
        # Configure mocks
        mock_from_env.return_value = self.mock_config
        mock_create_client.return_value = self.mock_supabase
        
        # Mock successful authentication
        auth_response = MagicMock()
        auth_response.error = None
        self.mock_supabase.auth.sign_in_with_password.return_value = auth_response
        
        # Mock successful connection test
        self.mock_supabase.rpc.return_value.execute.return_value = MagicMock()
        
        # Call the function
        result = supabase_database_audit.connect_to_supabase()
        
        # Assertions
        self.assertEqual(result, self.mock_supabase)
        mock_create_client.assert_called_once_with(
            self.mock_config.SUPABASE_URL, 
            self.mock_config.SUPABASE_KEY
        )
        self.mock_supabase.auth.sign_in_with_password.assert_called_once()
        self.mock_supabase.rpc.assert_called_once_with('get_service_role')

    @patch('supabase_database_audit.BaseConfig.from_env')
    @patch('supabase_database_audit.create_client')
    @patch('supabase_database_audit.sys.exit')
    def test_connect_to_supabase_failure(self, mock_exit, mock_create_client, mock_from_env):
        """Test connection failure to Supabase."""
        # Configure mocks
        mock_from_env.return_value = self.mock_config
        mock_create_client.return_value = self.mock_supabase
        
        # Mock connection test failure
        self.mock_supabase.rpc.return_value.execute.side_effect = Exception("Connection failed")
        
        # Call the function
        supabase_database_audit.connect_to_supabase()
        
        # Assertions
        mock_exit.assert_called_once_with(1)

    @patch('supabase_database_audit.logger')
    def test_verify_required_tables_all_exist(self, mock_logger):
        """Test verification when all required tables exist."""
        # Configure mock responses for each table
        self.mock_supabase.table.return_value.select.return_value.limit.return_value.execute.return_value = MagicMock(data=[{"id": 1}])
        
        # Call the function
        result, details = supabase_database_audit.verify_required_tables(self.mock_supabase)
        
        # Assertions
        self.assertTrue(result)
        self.assertEqual(len(details["existing"]), len(supabase_database_audit.REQUIRED_TABLES))
        self.assertEqual(len(details["missing"]), 0)
        mock_logger.info.assert_any_call("All required tables exist in the database")

    @patch('supabase_database_audit.logger')
    def test_verify_required_tables_missing(self, mock_logger):
        """Test verification when some required tables are missing."""
        # Configure mock to simulate missing tables
        def mock_table_response(table_name):
            mock_select = MagicMock()
            mock_limit = MagicMock()
            mock_execute = MagicMock()
            
            if table_name in ["molecules", "molecular_properties"]:
                # These tables exist
                mock_execute.return_value = MagicMock(data=[{"id": 1}])
            else:
                # These tables don't exist
                mock_execute.side_effect = Exception(f"Table {table_name} does not exist")
            
            mock_limit.return_value.execute = mock_execute
            mock_select.return_value.limit = mock_limit
            return mock_select
        
        self.mock_supabase.table.side_effect = lambda name: mock_table_response(name)
        
        # Call the function
        result, details = supabase_database_audit.verify_required_tables(self.mock_supabase)
        
        # Assertions
        self.assertFalse(result)
        self.assertEqual(len(details["existing"]), 2)  # Only molecules and molecular_properties exist
        self.assertEqual(len(details["missing"]), len(supabase_database_audit.REQUIRED_TABLES) - 2)
        mock_logger.critical.assert_any_call("Database schema verification failed. Required tables are missing.")

    @patch('supabase_database_audit.logger')
    def test_check_foreign_key_integrity_success(self, mock_logger):
        """Test foreign key integrity check when all relationships are valid."""
        # Configure mock responses for foreign key checks
        def mock_table_response(table_name):
            mock_select = MagicMock()
            mock_eq = MagicMock()
            mock_execute = MagicMock()
            
            # For foreign key value queries
            if table_name in supabase_database_audit.TABLE_RELATIONSHIPS:
                mock_execute.return_value = MagicMock(data=[{"molecule_id": "123", "mixture_id": "456"}])
            
            # For parent table existence checks
            mock_eq.return_value.execute.return_value = MagicMock(data=[{"id": "123"}])
            
            mock_select.return_value.execute = mock_execute
            mock_select.return_value.eq = mock_eq
            return mock_select
        
        self.mock_supabase.table.side_effect = lambda name: mock_table_response(name)
        
        # Call the function
        result, details = supabase_database_audit.check_foreign_key_integrity(self.mock_supabase)
        
        # Assertions
        self.assertTrue(result)
        mock_logger.info.assert_any_call("Foreign key integrity checks passed")

    @patch('supabase_database_audit.logger')
    def test_check_foreign_key_integrity_failure(self, mock_logger):
        """Test foreign key integrity check when some relationships are invalid."""
        # Configure mock responses for foreign key checks
        def mock_table_response(table_name):
            mock_select = MagicMock()
            mock_eq = MagicMock()
            mock_execute = MagicMock()
            
            # For foreign key value queries
            if table_name in supabase_database_audit.TABLE_RELATIONSHIPS:
                mock_execute.return_value = MagicMock(data=[{"molecule_id": "123", "mixture_id": "456"}])
            
            # For parent table existence checks - simulate failure for mixture table
            if table_name == "mixture":
                mock_eq.return_value.execute.return_value = MagicMock(data=[])  # Empty result = FK violation
            else:
                mock_eq.return_value.execute.return_value = MagicMock(data=[{"id": "123"}])
            
            mock_select.return_value.execute = mock_execute
            mock_select.return_value.eq = mock_eq
            return mock_select
        
        self.mock_supabase.table.side_effect = lambda name: mock_table_response(name)
        
        # Call the function
        result, details = supabase_database_audit.check_foreign_key_integrity(self.mock_supabase)
        
        # Assertions
        self.assertFalse(result)
        mock_logger.critical.assert_any_call("Database relationship verification failed. Foreign key constraints are violated.")

    @patch('supabase_database_audit.BaseConfig.from_env')
    @patch('supabase_database_audit.connect_to_supabase')
    @patch('supabase_database_audit.verify_required_tables')
    @patch('supabase_database_audit.check_foreign_key_integrity')
    @patch('supabase_database_audit.check_database_quotas')
    @patch('supabase_database_audit.get_table_counts')
    @patch('supabase_database_audit.analyze_data_completeness')
    @patch('supabase_database_audit.generate_dependency_graph')
    @patch('supabase_database_audit.generate_audit_report')
    @patch('supabase_database_audit.sys.exit')
    def test_main_verification_failure(self, mock_exit, mock_report, mock_graph, mock_completeness, 
                                      mock_counts, mock_quotas, mock_integrity, mock_tables, 
                                      mock_connect, mock_config):
        """Test main function when verification fails."""
        # Configure mocks
        mock_config.return_value = self.mock_config
        mock_connect.return_value = self.mock_supabase
        mock_tables.return_value = (False, {"missing": ["mixture"], "existing": ["molecules"]})
        mock_integrity.return_value = (False, {})
        mock_quotas.return_value = {"status": "PASS", "details": {}, "warnings": []}
        mock_counts.return_value = {"molecules": 10}
        mock_completeness.return_value = {}
        mock_graph.return_value = {"nodes": [], "edges": []}
        mock_report.return_value = {
            "summary": {
                "total_tables": 1,
                "empty_tables": 0,
                "populated_tables": 1,
                "integrity_issues": 1,
                "completeness_issues": 0
            },
            "population_recommendations": []
        }
        
        # Call the function
        supabase_database_audit.main()
        
        # Assertions
        mock_exit.assert_called_once_with(1)

    @patch('supabase_database_audit.BaseConfig.from_env')
    @patch('supabase_database_audit.connect_to_supabase')
    @patch('supabase_database_audit.verify_required_tables')
    @patch('supabase_database_audit.check_foreign_key_integrity')
    @patch('supabase_database_audit.check_database_quotas')
    @patch('supabase_database_audit.get_table_counts')
    @patch('supabase_database_audit.analyze_data_completeness')
    @patch('supabase_database_audit.generate_dependency_graph')
    @patch('supabase_database_audit.generate_audit_report')
    @patch('supabase_database_audit.sys.exit')
    def test_main_verification_success(self, mock_exit, mock_report, mock_graph, mock_completeness, 
                                      mock_counts, mock_quotas, mock_integrity, mock_tables, 
                                      mock_connect, mock_config):
        """Test main function when verification succeeds."""
        # Configure mocks
        mock_config.return_value = self.mock_config
        mock_connect.return_value = self.mock_supabase
        mock_tables.return_value = (True, {"missing": [], "existing": ["molecules", "mixture"]})
        mock_integrity.return_value = (True, {})
        mock_quotas.return_value = {"status": "PASS", "details": {}, "warnings": []}
        mock_counts.return_value = {"molecules": 10, "mixture": 5}
        mock_completeness.return_value = {}
        mock_graph.return_value = {"nodes": [], "edges": []}
        mock_report.return_value = {
            "summary": {
                "total_tables": 2,
                "empty_tables": 0,
                "populated_tables": 2,
                "integrity_issues": 0,
                "completeness_issues": 0
            },
            "population_recommendations": []
        }
        
        # Call the function
        supabase_database_audit.main()
        
        # Assertions
        mock_exit.assert_not_called()  # Should not exit when verification succeeds

if __name__ == '__main__':
    unittest.main()
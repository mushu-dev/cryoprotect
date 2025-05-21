#!/usr/bin/env python3
"""
Test script for RLS verification and remediation utilities.

This script tests the functionality of the rls_utils.py module to ensure that
it correctly verifies and remediates RLS settings for the molecule and
molecular_property tables.
"""

import os
import sys
import json
import unittest
from pathlib import Path
from typing import Dict, Any, List

# Import the RLS utilities
from rls_utils import (
    verify_rls, remediate_rls, ensure_rls_restored,
    is_rls_enabled, get_table_policies, TABLES_TO_VERIFY,
    get_supabase_client, execute_ddl
)

class TestRLSUtils(unittest.TestCase):
    """Test cases for RLS verification and remediation utilities."""
    
    def setUp(self):
        """Set up the test environment."""
        # Get a Supabase client
        self.client = get_supabase_client()
        
        # Ensure the logs directory exists
        Path("logs").mkdir(exist_ok=True)
        
        # Clear the audit log for clean testing
        self.audit_log_path = Path("logs/rls_audit.jsonl")
        if self.audit_log_path.exists():
            # Backup the existing log
            backup_path = Path("logs/rls_audit.jsonl.bak")
            if backup_path.exists():
                backup_path.unlink()
            self.audit_log_path.rename(backup_path)
        
        # Create a fresh audit log
        with open(self.audit_log_path, "w") as f:
            f.write("")
    
    def tearDown(self):
        """Clean up after tests."""
        # Restore RLS settings to ensure they're correct after tests
        issues = verify_rls()
        if issues:
            remediate_rls(issues)
    
    def test_is_rls_enabled(self):
        """Test the is_rls_enabled function."""
        for table in TABLES_TO_VERIFY:
            # Ensure RLS is enabled
            execute_ddl(f"ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;")
            
            # Check if RLS is enabled
            self.assertTrue(is_rls_enabled(table), f"RLS should be enabled for {table}")
            
            # Disable RLS
            execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
            
            # Check if RLS is disabled
            self.assertFalse(is_rls_enabled(table), f"RLS should be disabled for {table}")
            
            # Re-enable RLS
            execute_ddl(f"ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;")
    
    def test_get_table_policies(self):
        """Test the get_table_policies function."""
        for table in TABLES_TO_VERIFY:
            # Get policies
            policies = get_table_policies(table)
            
            # Check that we got some policies
            self.assertIsInstance(policies, dict, f"Policies for {table} should be a dictionary")
            
            # Check that we have at least the project member policies
            expected_policies = [
                f"Select {table}s for project members",
                f"Insert {table}s for project members",
                f"Update {table}s for project members",
                f"Delete {table}s for project members"
            ]
            
            for policy_name in expected_policies:
                self.assertIn(policy_name, policies, f"Policy {policy_name} should exist for {table}")
    
    def test_verify_rls_with_all_correct(self):
        """Test verify_rls when all settings are correct."""
        # Ensure all settings are correct
        issues = verify_rls()
        if issues:
            remediate_rls(issues)
        
        # Verify again - should be no issues
        issues = verify_rls()
        self.assertEqual(issues, {}, "There should be no RLS issues when all settings are correct")
    
    def test_verify_rls_with_disabled_rls(self):
        """Test verify_rls when RLS is disabled."""
        for table in TABLES_TO_VERIFY:
            # Disable RLS
            execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
            
            # Verify RLS
            issues = verify_rls()
            
            # Check that the issue was detected
            self.assertIn(table, issues, f"verify_rls should detect disabled RLS for {table}")
            self.assertTrue(issues[table].get("rls_disabled", False), 
                           f"verify_rls should report rls_disabled for {table}")
            
            # Re-enable RLS
            execute_ddl(f"ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;")
    
    def test_verify_rls_with_missing_policy(self):
        """Test verify_rls when a policy is missing."""
        for table in TABLES_TO_VERIFY:
            # Drop a policy
            policy_name = f"Select {table}s for project members"
            execute_ddl(f"DROP POLICY IF EXISTS \"{policy_name}\" ON {table};")
            
            # Verify RLS
            issues = verify_rls()
            
            # Check that the issue was detected
            self.assertIn(table, issues, f"verify_rls should detect missing policy for {table}")
            self.assertIn("missing_policies", issues[table], 
                         f"verify_rls should report missing_policies for {table}")
            self.assertIn(policy_name, issues[table]["missing_policies"], 
                         f"verify_rls should report {policy_name} as missing for {table}")
            
            # Remediate
            remediate_rls(issues)
    
    def test_verify_rls_with_incorrect_policy(self):
        """Test verify_rls when a policy is incorrect."""
        for table in TABLES_TO_VERIFY:
            # Drop and recreate a policy with incorrect definition
            policy_name = f"Update {table}s for project members"
            execute_ddl(f"DROP POLICY IF EXISTS \"{policy_name}\" ON {table};")
            
            incorrect_policy = f"""
            CREATE POLICY "{policy_name}" ON {table}
            FOR UPDATE
            USING (true);
            """
            execute_ddl(incorrect_policy)
            
            # Verify RLS
            issues = verify_rls()
            
            # Check that the issue was detected
            self.assertIn(table, issues, f"verify_rls should detect incorrect policy for {table}")
            self.assertIn("incorrect_policies", issues[table], 
                         f"verify_rls should report incorrect_policies for {table}")
            
            # The policy name in incorrect_policies includes the issue type
            incorrect_policy_found = False
            for policy in issues[table]["incorrect_policies"]:
                if policy.startswith(policy_name):
                    incorrect_policy_found = True
                    break
            
            self.assertTrue(incorrect_policy_found, 
                           f"verify_rls should report {policy_name} as incorrect for {table}")
            
            # Remediate
            remediate_rls(issues)
    
    def test_verify_rls_with_unauthorized_policy(self):
        """Test verify_rls when an unauthorized policy exists."""
        for table in TABLES_TO_VERIFY:
            # Create an unauthorized policy
            policy_name = f"Unauthorized policy for {table}"
            unauthorized_policy = f"""
            CREATE POLICY "{policy_name}" ON {table}
            FOR SELECT
            USING (true);
            """
            execute_ddl(unauthorized_policy)
            
            # Verify RLS
            issues = verify_rls()
            
            # Check that the issue was detected
            self.assertIn(table, issues, f"verify_rls should detect unauthorized policy for {table}")
            self.assertIn("unauthorized_policies", issues[table], 
                         f"verify_rls should report unauthorized_policies for {table}")
            self.assertIn(policy_name, issues[table]["unauthorized_policies"], 
                         f"verify_rls should report {policy_name} as unauthorized for {table}")
            
            # Remediate
            remediate_rls(issues)
    
    def test_remediate_rls(self):
        """Test remediate_rls."""
        for table in TABLES_TO_VERIFY:
            # Create multiple issues
            # 1. Disable RLS
            execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
            
            # 2. Drop a policy
            policy_name = f"Select {table}s for project members"
            execute_ddl(f"DROP POLICY IF EXISTS \"{policy_name}\" ON {table};")
            
            # 3. Create an incorrect policy
            incorrect_policy_name = f"Update {table}s for project members"
            execute_ddl(f"DROP POLICY IF EXISTS \"{incorrect_policy_name}\" ON {table};")
            
            incorrect_policy = f"""
            CREATE POLICY "{incorrect_policy_name}" ON {table}
            FOR UPDATE
            USING (true);
            """
            execute_ddl(incorrect_policy)
            
            # 4. Create an unauthorized policy
            unauthorized_policy_name = f"Unauthorized policy for {table}"
            unauthorized_policy = f"""
            CREATE POLICY "{unauthorized_policy_name}" ON {table}
            FOR SELECT
            USING (true);
            """
            execute_ddl(unauthorized_policy)
            
            # Verify RLS to get issues
            issues = verify_rls()
            
            # Remediate issues
            success = remediate_rls(issues)
            self.assertTrue(success, f"remediate_rls should succeed for {table}")
            
            # Verify again - should be no issues
            issues_after = verify_rls()
            self.assertNotIn(table, issues_after, 
                            f"There should be no issues for {table} after remediation")
    
    def test_ensure_rls_restored_decorator(self):
        """Test the ensure_rls_restored decorator."""
        # Define a test function that disables RLS
        @ensure_rls_restored
        def test_function():
            for table in TABLES_TO_VERIFY:
                execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
            return "Test function executed"
        
        # Ensure RLS is enabled before the test
        issues = verify_rls()
        if issues:
            remediate_rls(issues)
        
        # Call the test function
        result = test_function()
        self.assertEqual(result, "Test function executed", "The decorated function should execute normally")
        
        # Verify that RLS is still enabled after the function call
        for table in TABLES_TO_VERIFY:
            self.assertTrue(is_rls_enabled(table), 
                           f"RLS should be enabled for {table} after decorated function call")
    
    def test_ensure_rls_restored_decorator_with_exception(self):
        """Test the ensure_rls_restored decorator when the function raises an exception."""
        # Define a test function that disables RLS and raises an exception
        @ensure_rls_restored
        def test_function_with_exception():
            for table in TABLES_TO_VERIFY:
                execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
            raise ValueError("Test exception")
        
        # Ensure RLS is enabled before the test
        issues = verify_rls()
        if issues:
            remediate_rls(issues)
        
        # Call the test function and expect an exception
        with self.assertRaises(ValueError):
            test_function_with_exception()
        
        # Verify that RLS is still enabled after the function call
        for table in TABLES_TO_VERIFY:
            self.assertTrue(is_rls_enabled(table), 
                           f"RLS should be enabled for {table} after decorated function exception")
    
    def test_audit_logging(self):
        """Test that audit logging works correctly."""
        # Create an issue
        table = TABLES_TO_VERIFY[0]
        execute_ddl(f"ALTER TABLE {table} DISABLE ROW LEVEL SECURITY;")
        
        # Verify RLS
        issues = verify_rls()
        
        # Remediate issues
        remediate_rls(issues)
        
        # Check that audit log entries were created
        self.assertTrue(self.audit_log_path.exists(), "Audit log file should exist")
        
        # Read the audit log
        audit_entries = []
        with open(self.audit_log_path, "r") as f:
            for line in f:
                if line.strip():
                    audit_entries.append(json.loads(line))
        
        # Check that we have at least two entries (verify and remediate)
        self.assertGreaterEqual(len(audit_entries), 2, 
                               "Audit log should have at least two entries")
        
        # Check that we have a verify entry
        verify_entries = [e for e in audit_entries if e["action"] == "verify"]
        self.assertGreaterEqual(len(verify_entries), 1, 
                               "Audit log should have at least one verify entry")
        
        # Check that we have a remediate entry
        remediate_entries = [e for e in audit_entries if e["action"] == "remediate"]
        self.assertGreaterEqual(len(remediate_entries), 1, 
                               "Audit log should have at least one remediate entry")

if __name__ == "__main__":
    unittest.main()
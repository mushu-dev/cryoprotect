#!/usr/bin/env python
"""
test_rls_policy_improvements.py: Validates that the RLS policy improvements have been
properly applied and are functioning correctly.

This script extends the original RLS verification tests to validate the improved
team-based access model and security definer functions.
"""

import os
import sys
import unittest
import psycopg2
import json
import uuid
from datetime import datetime
import logging

# Add the parent directory to the path to import from the project root
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the RLS test helper
from tests.rls_test_helper import RLSTestHelper

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create test UUIDs for consistency
admin_user_id = "cbaff98e-1e2a-4e9e-a654-6a89365db0d3"
member_user_id = "7f9c24e0-9c3b-4a28-8f89-3e5e4e9b7e5d"
non_member_user_id = "5a8e7d2c-1f3b-4d7a-9b8c-2e1d3f4c5a6b"

team1_id = "1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d"
team2_id = "2b3c4d5e-6f7a-8b9c-0d1e-2f3a4b5c6d7e"

project1_id = "3c4d5e6f-7a8b-9c0d-1e2f-3a4b5c6d7e8f"
project2_id = "4d5e6f7a-8b9c-0d1e-2f3a-4b5c6d7e8f9a"

molecule1_id = "5e6f7a8b-9c0d-1e2f-3a4b-5c6d7e8f9a0b"
molecule2_id = "6f7a8b9c-0d1e-2f3a-4b5c-6d7e8f9a0b1c"

experiment1_id = "7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2d"
experiment2_id = "8b9c0d1e-2f3a-4b5c-6d7e-8f9a0b1c2d3e"

mixture1_id = "9c0d1e2f-3a4b-5c6d-7e8f-9a0b1c2d3e4f"
mixture2_id = "0d1e2f3a-4b5c-6d7e-8f9a-0b1c2d3e4f5a"

class RLSSecurityDefinerFunctionTests(unittest.TestCase):
    """Tests for the security definer functions that support RLS policies."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment once for the whole class."""
        cls.helper = RLSTestHelper()
        cls.helper.create_test_data()
        
        # Store UUIDs for use in tests
        cls.admin_user_id = admin_user_id
        cls.member_user_id = member_user_id
        cls.non_member_user_id = non_member_user_id
        cls.team1_id = team1_id
        cls.team2_id = team2_id
        cls.project1_id = project1_id
        cls.project2_id = project2_id
    
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests have run."""
        cls.helper.cleanup_test_data()
    
    def test_is_team_member_function(self):
        """TC-RLSI-001: Test auth.is_team_member function."""
        # Test with admin user (should be true)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT auth.is_team_member('{self.team1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with member user (should be true)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT auth.is_team_member('{self.team1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with non-member user (should be false)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT auth.is_team_member('{self.team1_id}')"
        )
        self.assertFalse(result[0][0])
    
    def test_is_team_admin_function(self):
        """TC-RLSI-002: Test auth.is_team_admin function."""
        # Test with admin user (should be true)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT auth.is_team_admin('{self.team1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with member user (should be false - regular member, not admin)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT auth.is_team_admin('{self.team1_id}')"
        )
        self.assertFalse(result[0][0])
        
        # Test with non-member user (should be false)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT auth.is_team_admin('{self.team1_id}')"
        )
        self.assertFalse(result[0][0])
    
    def test_get_user_teams_function(self):
        """TC-RLSI-003: Test auth.get_user_teams function."""
        # Test with admin user (should return team IDs)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            "SELECT * FROM auth.get_user_teams()"
        )
        self.assertTrue(len(result) > 0)
        
        # Test with member user (should return team IDs)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            "SELECT * FROM auth.get_user_teams()"
        )
        self.assertTrue(len(result) > 0)
        
        # Test with non-member user (should return empty set)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            "SELECT * FROM auth.get_user_teams()"
        )
        self.assertEqual(len(result), 0)
    
    def test_can_access_project_function(self):
        """TC-RLSI-004: Test auth.can_access_project function."""
        # Test with admin user for project1 (should be true)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT auth.can_access_project('{self.project1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with member user for project1 (should be true)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT auth.can_access_project('{self.project1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with non-member user for project1 (should be false)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT auth.can_access_project('{self.project1_id}')"
        )
        self.assertFalse(result[0][0])
    
    def test_can_manage_project_function(self):
        """TC-RLSI-005: Test auth.can_manage_project function."""
        # Test with admin user for project1 (should be true)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT auth.can_manage_project('{self.project1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with member user for project1 (should be false - regular member, not admin)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT auth.can_manage_project('{self.project1_id}')"
        )
        self.assertFalse(result[0][0])
        
        # Test with non-member user for project1 (should be false)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT auth.can_manage_project('{self.project1_id}')"
        )
        self.assertFalse(result[0][0])
    
    def test_get_accessible_projects_function(self):
        """TC-RLSI-006: Test auth.get_accessible_projects function."""
        # Test with admin user (should return project IDs)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            "SELECT * FROM auth.get_accessible_projects()"
        )
        self.assertTrue(len(result) > 0)
        
        # Test with member user (should return project IDs)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            "SELECT * FROM auth.get_accessible_projects()"
        )
        self.assertTrue(len(result) > 0)
        
        # Test with non-member user (should return empty set)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            "SELECT * FROM auth.get_accessible_projects()"
        )
        self.assertEqual(len(result), 0)
    
    def test_can_access_molecule_function(self):
        """TC-RLSI-007: Test auth.can_access_molecule function."""
        # Test with admin user for molecule1 (should be true)
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT auth.can_access_molecule('{molecule1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with member user for molecule1 (should be true)
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT auth.can_access_molecule('{molecule1_id}')"
        )
        self.assertTrue(result[0][0])
        
        # Test with non-member user for molecule1 (should be false)
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT auth.can_access_molecule('{molecule1_id}')"
        )
        self.assertFalse(result[0][0])

class RLSTablePolicyTests(unittest.TestCase):
    """Tests for the RLS policies on tables."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment once for the whole class."""
        cls.helper = RLSTestHelper()
        cls.helper.create_test_data()
        
        # Store UUIDs for use in tests
        cls.admin_user_id = admin_user_id
        cls.member_user_id = member_user_id
        cls.non_member_user_id = non_member_user_id
        cls.team1_id = team1_id
        cls.team2_id = team2_id
        cls.project1_id = project1_id
        cls.project2_id = project2_id
    
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests have run."""
        cls.helper.cleanup_test_data()
    
    def test_projects_select_policy(self):
        """TC-RLSP-001: Test SELECT access to projects table."""
        # Test admin user can select all projects in their team
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"SELECT COUNT(*) FROM projects WHERE team_id = '{self.team1_id}'"
        )
        self.assertGreater(result[0][0], 0)
        
        # Test member user can select all projects in their team
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"SELECT COUNT(*) FROM projects WHERE team_id = '{self.team1_id}'"
        )
        self.assertGreater(result[0][0], 0)
        
        # Test non-member user cannot select projects from team1
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"SELECT COUNT(*) FROM projects WHERE team_id = '{self.team1_id}'"
        )
        self.assertEqual(result[0][0], 0)
    
    def test_projects_insert_policy(self):
        """TC-RLSP-002: Test INSERT access to projects table."""
        # Test admin user can insert projects for their team
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            INSERT INTO projects (id, name, team_id, created_at) 
            VALUES ('{uuid.uuid4()}', 'Test Admin Project', '{self.team1_id}', NOW())
            RETURNING id
            """
        )
        self.assertIsNotNone(result)
        
        # Test member user can insert projects for their team
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"""
            INSERT INTO projects (id, name, team_id, created_at) 
            VALUES ('{uuid.uuid4()}', 'Test Member Project', '{self.team1_id}', NOW())
            RETURNING id
            """
        )
        self.assertIsNotNone(result)
        
        # Test non-member user cannot insert projects for team1
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.non_member_user_id,
                f"""
                INSERT INTO projects (id, name, team_id, created_at) 
                VALUES ('{uuid.uuid4()}', 'Test Non-Member Project', '{self.team1_id}', NOW())
                RETURNING id
                """
            )
    
    def test_projects_update_policy(self):
        """TC-RLSP-003: Test UPDATE access to projects table."""
        # Test admin user can update projects for their team
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            UPDATE projects 
            SET name = 'Updated by Admin' 
            WHERE id = '{self.project1_id}'
            RETURNING id
            """
        )
        self.assertIsNotNone(result)
        
        # Test member user cannot update projects (only admins can)
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.member_user_id,
                f"""
                UPDATE projects 
                SET name = 'Updated by Member' 
                WHERE id = '{self.project1_id}'
                RETURNING id
                """
            )
        
        # Test non-member user cannot update projects
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.non_member_user_id,
                f"""
                UPDATE projects 
                SET name = 'Updated by Non-Member' 
                WHERE id = '{self.project1_id}'
                RETURNING id
                """
            )
    
    def test_projects_delete_policy(self):
        """TC-RLSP-004: Test DELETE access to projects table."""
        # Create a temporary project to delete
        temp_project_id = str(uuid.uuid4())
        self.helper.execute_sql(
            f"""
            INSERT INTO projects (id, name, team_id, created_at) 
            VALUES ('{temp_project_id}', 'Temp Project', '{self.team1_id}', NOW())
            """
        )
        
        # Test admin user can delete projects for their team
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            DELETE FROM projects 
            WHERE id = '{temp_project_id}'
            RETURNING id
            """
        )
        self.assertIsNotNone(result)
        
        # Create another temporary project
        temp_project_id = str(uuid.uuid4())
        self.helper.execute_sql(
            f"""
            INSERT INTO projects (id, name, team_id, created_at) 
            VALUES ('{temp_project_id}', 'Temp Project', '{self.team1_id}', NOW())
            """
        )
        
        # Test member user cannot delete projects (only admins can)
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.member_user_id,
                f"""
                DELETE FROM projects 
                WHERE id = '{temp_project_id}'
                RETURNING id
                """
            )
        
        # Test non-member user cannot delete projects
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.non_member_user_id,
                f"""
                DELETE FROM projects 
                WHERE id = '{temp_project_id}'
                RETURNING id
                """
            )
        
        # Clean up
        self.helper.execute_sql(f"DELETE FROM projects WHERE id = '{temp_project_id}'")
    
    def test_molecules_access_through_experiments(self):
        """TC-RLSP-005: Test access to molecules through experiments."""
        # Test admin user can select molecules through experiments
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            SELECT m.* FROM molecules m
            JOIN experiments e ON e.molecule_id = m.id
            WHERE e.project_id = '{self.project1_id}'
            """
        )
        self.assertGreater(len(result), 0)
        
        # Test member user can select molecules through experiments
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"""
            SELECT m.* FROM molecules m
            JOIN experiments e ON e.molecule_id = m.id
            WHERE e.project_id = '{self.project1_id}'
            """
        )
        self.assertGreater(len(result), 0)
        
        # Test non-member user cannot see molecules through experiments
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"""
            SELECT m.* FROM molecules m
            JOIN experiments e ON e.molecule_id = m.id
            WHERE e.project_id = '{self.project1_id}'
            """
        )
        self.assertEqual(len(result), 0)
    
    def test_molecular_properties_policy(self):
        """TC-RLSP-006: Test access to molecular_properties table."""
        # Insert test molecular properties
        self.helper.execute_sql(
            f"""
            INSERT INTO molecular_properties (id, molecule_id, property_name, property_value, created_at)
            VALUES ('{uuid.uuid4()}', '{molecule1_id}', 'test_property', 'test_value', NOW())
            """
        )
        
        # Test admin user can select molecular properties
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            SELECT * FROM molecular_properties
            WHERE molecule_id = '{molecule1_id}'
            """
        )
        self.assertGreater(len(result), 0)
        
        # Test member user can select molecular properties
        result = self.helper.execute_sql_as_user(
            self.member_user_id,
            f"""
            SELECT * FROM molecular_properties
            WHERE molecule_id = '{molecule1_id}'
            """
        )
        self.assertGreater(len(result), 0)
        
        # Test non-member user cannot select molecular properties
        result = self.helper.execute_sql_as_user(
            self.non_member_user_id,
            f"""
            SELECT * FROM molecular_properties
            WHERE molecule_id = '{molecule1_id}'
            """
        )
        self.assertEqual(len(result), 0)
    
    def test_user_teams_policy(self):
        """TC-RLSP-007: Test access to user_teams table."""
        # Test admin user can select team members
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            SELECT * FROM user_teams
            WHERE team_id = '{self.team1_id}'
            """
        )
        self.assertGreater(len(result), 0)
        
        # Test admin user can add new team members
        new_user_id = str(uuid.uuid4())
        result = self.helper.execute_sql_as_user(
            self.admin_user_id,
            f"""
            INSERT INTO user_teams (id, user_id, team_id, role, created_at)
            VALUES ('{uuid.uuid4()}', '{new_user_id}', '{self.team1_id}', 'member', NOW())
            RETURNING id
            """
        )
        self.assertIsNotNone(result)
        
        # Test regular member cannot add team members
        with self.assertRaises(Exception):
            self.helper.execute_sql_as_user(
                self.member_user_id,
                f"""
                INSERT INTO user_teams (id, user_id, team_id, role, created_at)
                VALUES ('{uuid.uuid4()}', '{str(uuid.uuid4())}', '{self.team1_id}', 'member', NOW())
                RETURNING id
                """
            )

class RLSPerformanceTests(unittest.TestCase):
    """Tests for RLS policy performance with the new security definer functions."""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment once for the whole class."""
        cls.helper = RLSTestHelper()
        cls.helper.create_test_data()
        
        # Store UUIDs for use in tests
        cls.admin_user_id = admin_user_id
        cls.member_user_id = member_user_id
        cls.team1_id = team1_id
        cls.project1_id = project1_id
    
    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests have run."""
        cls.helper.cleanup_test_data()
    
    def test_query_performance(self):
        """TC-RLSPERF-001: Test query performance with the new RLS policies."""
        # Measure performance of a query with the new RLS policies
        start_time = datetime.now()
        
        self.helper.execute_sql_as_user(
            self.member_user_id,
            f"""
            EXPLAIN ANALYZE
            SELECT m.* FROM molecules m
            JOIN experiments e ON e.molecule_id = m.id
            WHERE e.project_id = '{self.project1_id}'
            """
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # For now, just log the duration - in a real test we would compare with a baseline
        logger.info(f"Query performance test duration: {duration} seconds")
        
        # Ensure the query completes in a reasonable time (adjust as needed)
        self.assertLess(duration, 2.0)  # Should complete in under 2 seconds

def run_tests():
    """Run the tests and return the results."""
    # Create a test suite
    suite = unittest.TestSuite()
    
    # Add test cases
    suite.addTest(unittest.makeSuite(RLSSecurityDefinerFunctionTests))
    suite.addTest(unittest.makeSuite(RLSTablePolicyTests))
    suite.addTest(unittest.makeSuite(RLSPerformanceTests))
    
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return the test result
    return result

if __name__ == "__main__":
    # Run the tests
    result = run_tests()
    
    # Exit with appropriate status code
    sys.exit(0 if result.wasSuccessful() else 1)
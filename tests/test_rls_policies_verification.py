"""
Comprehensive RLS Policy Verification Tests
Based on the test cases defined in DATABASE_VERIFICATION_PLAN.md
"""
import os
import sys
import unittest
import json
import logging
from dotenv import load_dotenv

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import RLS test helper
from tests.rls_test_helper import RLSTestHelper

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='rls_verification_tests.log'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Get database connection parameters
DB_CONNECTION_PARAMS = {
    'host': os.environ.get('SUPABASE_DB_HOST'),
    'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
    'dbname': os.environ.get('SUPABASE_DB_NAME'),
    'user': os.environ.get('SUPABASE_DB_USER'),
    'password': os.environ.get('SUPABASE_DB_PASSWORD'),
    'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
}

class RLSSecurityDefinerFunctionTests(unittest.TestCase):
    """Tests for Security Definer Functions (TC-RLS-001 to TC-RLS-007)"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for all tests"""
        cls.helper = RLSTestHelper(DB_CONNECTION_PARAMS)
        cls.test_data = cls.helper.create_test_data()
        logger.info(f"Test data created: {json.dumps(cls.test_data)}")

        # Store user and project IDs for easy access
        cls.user1_id = cls.test_data['users'][0]
        cls.user2_id = cls.test_data['users'][1]
        cls.project1_id = cls.test_data['projects'][0]
        cls.project2_id = cls.test_data['projects'][1]
        cls.team_id = cls.test_data['team']
        cls.molecule1_id = cls.test_data['molecules'][0] if cls.test_data.get('molecules') else None
        cls.molecule2_id = cls.test_data['molecules'][1] if len(cls.test_data.get('molecules', [])) > 1 else None
        cls.experiment1_id = cls.test_data['experiments'][0] if cls.test_data.get('experiments') else None
        cls.experiment2_id = cls.test_data['experiments'][1] if len(cls.test_data.get('experiments', [])) > 1 else None
        cls.mixture1_id = cls.test_data['mixtures'][0] if cls.test_data.get('mixtures') else None
        cls.mixture2_id = cls.test_data['mixtures'][1] if len(cls.test_data.get('mixtures', [])) > 1 else None

    @classmethod
    def tearDownClass(cls):
        """Clean up test data after all tests"""
        # Uncomment the following lines to clean up test data
        # Don't clean up by default to allow manual inspection
        # cls.helper.cleanup_test_data(cls.test_data)
        # logger.info("Test data cleaned up")
        pass
    
    def test_is_project_member_valid(self):
        """TC-RLS-001: Verify is_project_member function with valid project"""
        # Test user1 is a member of project1
        result = self.helper.test_security_definer_function(
            self.user1_id, 'is_project_member', [self.project1_id]
        )
        self.assertTrue(result)
        
        # Test user2 is a member of project2
        result = self.helper.test_security_definer_function(
            self.user2_id, 'is_project_member', [self.project2_id]
        )
        self.assertTrue(result)
    
    def test_is_project_member_invalid(self):
        """TC-RLS-002: Verify is_project_member function with invalid project"""
        # Test user1 is NOT a member of project2
        result = self.helper.test_security_definer_function(
            self.user1_id, 'is_project_member', [self.project2_id]
        )
        self.assertFalse(result)
        
        # Test user2 is NOT a member of project1
        result = self.helper.test_security_definer_function(
            self.user2_id, 'is_project_member', [self.project1_id]
        )
        self.assertFalse(result)
    
    def test_user_projects(self):
        """TC-RLS-003: Verify user_projects function"""
        # This test requires the user_projects function to exist
        # Check if it exists first
        sql = "SELECT COUNT(*) FROM pg_proc WHERE proname = 'user_projects'"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or result.get('data', [{}])[0].get('count', 0) == 0:
            self.skipTest("user_projects function does not exist")
            
        # Test user1's projects
        sql = f"SELECT user_projects()"
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        
        # If the function returns a set, we need to check it contains project1_id
        # The exact format depends on the function implementation
        if 'error' not in result:
            data_str = str(result.get('data', []))
            self.assertIn(self.project1_id, data_str)
            self.assertNotIn(self.project2_id, data_str)
    
    def test_is_team_member(self):
        """TC-RLS-004: Verify is_team_member function"""
        # This test requires the is_team_member function to exist
        # Check if it exists first
        sql = "SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_team_member'"
        result = self.helper.execute_sql(sql, role="service_role")

        if 'error' in result or result.get('data', [{}])[0].get('count', 0) == 0:
            self.skipTest("is_team_member function does not exist")

        # We already have a team created in setup with ID "33333333-3333-3333-3333-333333333333"
        # User1 and User2 are already added to the team in different roles
        team_id = self.team_id
        
        # Test user1 is a member of the team
        result = self.helper.test_security_definer_function(
            self.user1_id, 'is_team_member', [team_id]
        )
        self.assertTrue(result)
        
        # Test user2 is NOT a member of the team
        result = self.helper.test_security_definer_function(
            self.user2_id, 'is_team_member', [team_id]
        )
        self.assertFalse(result)
    
    def test_is_project_owner(self):
        """TC-RLS-005: Verify is_project_owner function"""
        # This test requires the is_project_owner function to exist
        # Check if it exists first
        sql = "SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_project_owner'"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or result.get('data', [{}])[0].get('count', 0) == 0:
            self.skipTest("is_project_owner function does not exist")
        
        # Test user1 is an owner of project1
        result = self.helper.test_security_definer_function(
            self.user1_id, 'is_project_owner', [self.project1_id]
        )
        self.assertTrue(result)
        
        # Test user2 is NOT an owner of project1
        result = self.helper.test_security_definer_function(
            self.user2_id, 'is_project_owner', [self.project1_id]
        )
        self.assertFalse(result)
    
    def test_molecule_in_user_project(self):
        """TC-RLS-006: Verify molecule_in_user_project function"""
        # This test requires the molecule_in_user_project function to exist
        # Check if it exists first
        sql = "SELECT COUNT(*) FROM pg_proc WHERE proname = 'molecule_in_user_project'"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or result.get('data', [{}])[0].get('count', 0) == 0:
            self.skipTest("molecule_in_user_project function does not exist")
        
        # Get a molecule from project1
        sql = f"SELECT id FROM molecule WHERE project_id = '{self.project1_id}' LIMIT 1"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or len(result.get('data', [])) == 0:
            self.skipTest("No molecules found in test project")
        
        molecule1_id = result['data'][0]['id']
        
        # Test user1 can access molecule1
        result = self.helper.test_security_definer_function(
            self.user1_id, 'molecule_in_user_project', [molecule1_id]
        )
        self.assertTrue(result)
        
        # Test user2 cannot access molecule1
        result = self.helper.test_security_definer_function(
            self.user2_id, 'molecule_in_user_project', [molecule1_id]
        )
        self.assertFalse(result)
    
    def test_mixture_in_user_project(self):
        """TC-RLS-007: Verify mixture_in_user_project function"""
        # This test requires the mixture_in_user_project function to exist
        # Check if it exists first
        sql = "SELECT COUNT(*) FROM pg_proc WHERE proname = 'mixture_in_user_project'"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or result.get('data', [{}])[0].get('count', 0) == 0:
            self.skipTest("mixture_in_user_project function does not exist")
        
        # Get a mixture from project1
        sql = f"SELECT id FROM mixture WHERE project_id = '{self.project1_id}' LIMIT 1"
        result = self.helper.execute_sql(sql, role="service_role")
        
        if 'error' in result or len(result.get('data', [])) == 0:
            self.skipTest("No mixtures found in test project")
        
        mixture1_id = result['data'][0]['id']
        
        # Test user1 can access mixture1
        result = self.helper.test_security_definer_function(
            self.user1_id, 'mixture_in_user_project', [mixture1_id]
        )
        self.assertTrue(result)
        
        # Test user2 cannot access mixture1
        result = self.helper.test_security_definer_function(
            self.user2_id, 'mixture_in_user_project', [mixture1_id]
        )
        self.assertFalse(result)

class RLSTableAccessPolicyTests(unittest.TestCase):
    """Tests for Table Access Policies (TC-RLS-010 to TC-RLS-017)"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for all tests"""
        cls.helper = RLSTestHelper(DB_CONNECTION_PARAMS)
        cls.test_data = cls.helper.create_test_data()
        logger.info(f"Test data created: {json.dumps(cls.test_data)}")

        # Store user and project IDs for easy access
        cls.user1_id = cls.test_data['users'][0]
        cls.user2_id = cls.test_data['users'][1]
        cls.project1_id = cls.test_data['projects'][0]
        cls.project2_id = cls.test_data['projects'][1]
        cls.team_id = cls.test_data['team']
        cls.molecule1_id = cls.test_data['molecules'][0] if cls.test_data.get('molecules') else None
        cls.molecule2_id = cls.test_data['molecules'][1] if len(cls.test_data.get('molecules', [])) > 1 else None
        cls.experiment1_id = cls.test_data['experiments'][0] if cls.test_data.get('experiments') else None
        cls.experiment2_id = cls.test_data['experiments'][1] if len(cls.test_data.get('experiments', [])) > 1 else None
        cls.mixture1_id = cls.test_data['mixtures'][0] if cls.test_data.get('mixtures') else None
        cls.mixture2_id = cls.test_data['mixtures'][1] if len(cls.test_data.get('mixtures', [])) > 1 else None

    @classmethod
    def tearDownClass(cls):
        """Clean up test data after all tests"""
        # Uncomment the following lines to clean up test data
        # Don't clean up by default to allow manual inspection
        # cls.helper.cleanup_test_data(cls.test_data)
        # logger.info("Test data cleaned up")
        pass
    
    def test_select_molecule_as_member(self):
        """TC-RLS-010: Test SELECT access to molecule table as project member"""
        # Test user1 can select molecules from project1 (via experiments)
        result = self.helper.test_user_access(
            self.user1_id, 'molecules', f"id IN (SELECT molecule_id FROM experiments WHERE project_id = '{self.project1_id}')"
        )
        self.assertGreater(result.get('count', 0), 0)
        self.assertNotIn('error', result)

    def test_select_molecule_as_non_member(self):
        """TC-RLS-011: Test SELECT access to molecule table as non-member"""
        # Test user1 cannot select molecules from project2 (via experiments)
        result = self.helper.test_user_access(
            self.user1_id, 'molecules', f"id IN (SELECT molecule_id FROM experiments WHERE project_id = '{self.project2_id}')"
        )
        self.assertEqual(result.get('count', -1), 0)

    def test_insert_molecule_as_member(self):
        """TC-RLS-012: Test INSERT into molecule table as project member"""
        # Test user1 can insert molecule
        molecule_id = str(uuid.uuid4())

        # First insert the molecule
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Test Insert Molecule', 'CCO', '{self.user1_id}')
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        self.assertNotIn('error', result)

        # Then link it to project1 via experiment
        if 'error' not in result:
            experiment_id = str(uuid.uuid4())
            sql = f"""
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                VALUES ('{experiment_id}', '{molecule_id}', '{self.project1_id}', 'Test Insert Experiment', '{self.user1_id}')
                RETURNING id
            """
            result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
            self.assertNotIn('error', result)

    def test_insert_molecule_as_non_member(self):
        """TC-RLS-013: Test INSERT into molecule table as non-member"""
        # Test user1 cannot insert experiment into project2
        molecule_id = str(uuid.uuid4())

        # First insert the molecule (this should work)
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Test Insert Molecule Non-Member', 'CCCO', '{self.user1_id}')
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)

        # Try to link to project2 (this should fail)
        if 'error' not in result:
            experiment_id = str(uuid.uuid4())
            sql = f"""
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                VALUES ('{experiment_id}', '{molecule_id}', '{self.project2_id}', 'Test Insert Experiment Fail', '{self.user1_id}')
                RETURNING id
            """
            result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
            self.assertIn('error', result)

    def test_update_molecule_as_member(self):
        """TC-RLS-014: Test UPDATE on molecule table as project member"""
        # Get a molecule linked to project1 via experiment
        sql = f"""
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.molecule_id
            WHERE e.project_id = '{self.project1_id}'
            LIMIT 1
        """
        result = self.helper.execute_sql(sql, role="service_role")

        if 'error' in result or len(result.get('data', [])) == 0:
            self.skipTest("No molecules found linked to test project")

        molecule1_id = result['data'][0]['id']

        # Test user1 can update molecule1
        sql = f"""
            UPDATE molecules
            SET name = 'Updated Test Molecule'
            WHERE id = '{molecule1_id}'
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        self.assertNotIn('error', result)

    def test_update_molecule_as_non_member(self):
        """TC-RLS-015: Test UPDATE on molecule table as non-member"""
        # Get a molecule linked to project2 via experiment
        sql = f"""
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.molecule_id
            WHERE e.project_id = '{self.project2_id}'
            LIMIT 1
        """
        result = self.helper.execute_sql(sql, role="service_role")

        if 'error' in result or len(result.get('data', [])) == 0:
            self.skipTest("No molecules found linked to test project")

        molecule2_id = result['data'][0]['id']

        # Test user1 cannot update molecule2
        sql = f"""
            UPDATE molecules
            SET name = 'Should Fail Update'
            WHERE id = '{molecule2_id}'
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        self.assertIn('error', result)

    def test_delete_molecule_as_member(self):
        """TC-RLS-016: Test DELETE from molecule table as project member"""
        # Insert a molecule to delete
        molecule_id = str(uuid.uuid4())

        # First create the molecule
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Test Delete Molecule', 'CO', '{self.user1_id}')
            RETURNING id
        """
        molecule_result = self.helper.execute_sql(sql, role="service_role")

        # Link it to project1 via experiment
        if 'error' not in molecule_result:
            experiment_id = str(uuid.uuid4())
            sql = f"""
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                VALUES ('{experiment_id}', '{molecule_id}', '{self.project1_id}', 'Test Delete Experiment', '{self.user1_id}')
                RETURNING id
            """
            self.helper.execute_sql(sql, role="service_role")

        # Clean up the experiment first (it has a foreign key to molecule)
        sql = f"""
            DELETE FROM experiments
            WHERE molecule_id = '{molecule_id}'
        """
        self.helper.execute_sql(sql, role="service_role")

        # Test user1 can delete molecule
        sql = f"""
            DELETE FROM molecules
            WHERE id = '{molecule_id}'
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        self.assertNotIn('error', result)

    def test_delete_molecule_as_non_member(self):
        """TC-RLS-017: Test DELETE from molecule table as non-member"""
        # Insert a molecule to delete
        molecule_id = str(uuid.uuid4())

        # First create the molecule
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Test Delete Molecule Fail', 'COO', '{self.user2_id}')
            RETURNING id
        """
        molecule_result = self.helper.execute_sql(sql, role="service_role")

        # Link it to project2 via experiment
        if 'error' not in molecule_result:
            experiment_id = str(uuid.uuid4())
            sql = f"""
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                VALUES ('{experiment_id}', '{molecule_id}', '{self.project2_id}', 'Test Delete Experiment Fail', '{self.user2_id}')
                RETURNING id
            """
            self.helper.execute_sql(sql, role="service_role")

        # Test user1 cannot delete molecule from project2
        sql = f"""
            DELETE FROM molecules
            WHERE id = '{molecule_id}'
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="authenticated", user_id=self.user1_id)
        self.assertIn('error', result)

        # Clean up experiment for the next test
        sql = f"DELETE FROM experiments WHERE molecule_id = '{molecule_id}'"
        self.helper.execute_sql(sql, role="service_role")

        # Clean up molecule for the next test
        sql = f"DELETE FROM molecules WHERE id = '{molecule_id}'"
        self.helper.execute_sql(sql, role="service_role")

class RLSRelationshipPolicyTests(unittest.TestCase):
    """Tests for Relationship Policies (TC-RLS-020 to TC-RLS-025)"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for all tests"""
        cls.helper = RLSTestHelper(DB_CONNECTION_PARAMS)
        cls.test_data = cls.helper.create_test_data()
        logger.info(f"Test data created: {json.dumps(cls.test_data)}")

        # Store user and project IDs for easy access
        cls.user1_id = cls.test_data['users'][0]
        cls.user2_id = cls.test_data['users'][1]
        cls.project1_id = cls.test_data['projects'][0]
        cls.project2_id = cls.test_data['projects'][1]
        cls.team_id = cls.test_data['team']
        cls.molecule1_id = cls.test_data['molecules'][0] if cls.test_data.get('molecules') else None
        cls.molecule2_id = cls.test_data['molecules'][1] if len(cls.test_data.get('molecules', [])) > 1 else None
        cls.experiment1_id = cls.test_data['experiments'][0] if cls.test_data.get('experiments') else None
        cls.experiment2_id = cls.test_data['experiments'][1] if len(cls.test_data.get('experiments', [])) > 1 else None
        cls.mixture1_id = cls.test_data['mixtures'][0] if cls.test_data.get('mixtures') else None
        cls.mixture2_id = cls.test_data['mixtures'][1] if len(cls.test_data.get('mixtures', [])) > 1 else None

        # Create related data for testing
        cls.setup_related_data()
    
    @classmethod
    def setup_related_data(cls):
        """Set up related data for relationship policy tests"""
        # Create molecular properties
        if cls.molecule1_id and cls.molecule2_id:
            # Check if molecular properties exist
            sql = f"SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = '{cls.molecule1_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO molecular_properties (id, molecule_id, property_name, property_value, unit, created_by)
                    VALUES
                    ('{str(uuid.uuid4())}', '{cls.molecule1_id}', 'boiling_point', '100.0', 'C', '{cls.user1_id}'),
                    ('{str(uuid.uuid4())}', '{cls.molecule1_id}', 'melting_point', '0.0', 'C', '{cls.user1_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")

            sql = f"SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = '{cls.molecule2_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO molecular_properties (id, molecule_id, property_name, property_value, unit, created_by)
                    VALUES
                    ('{str(uuid.uuid4())}', '{cls.molecule2_id}', 'boiling_point', '120.0', 'C', '{cls.user2_id}'),
                    ('{str(uuid.uuid4())}', '{cls.molecule2_id}', 'melting_point', '5.0', 'C', '{cls.user2_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")

        # Create mixture components
        if cls.mixture1_id and cls.mixture2_id and cls.molecule1_id and cls.molecule2_id:
            # Check if mixture components exist
            sql = f"SELECT COUNT(*) FROM mixture_components WHERE mixture_id = '{cls.mixture1_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO mixture_components (id, mixture_id, molecule_id, concentration, concentration_unit, created_by)
                    VALUES ('{str(uuid.uuid4())}', '{cls.mixture1_id}', '{cls.molecule1_id}', 50.0, 'percent', '{cls.user1_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")

            sql = f"SELECT COUNT(*) FROM mixture_components WHERE mixture_id = '{cls.mixture2_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO mixture_components (id, mixture_id, molecule_id, concentration, concentration_unit, created_by)
                    VALUES ('{str(uuid.uuid4())}', '{cls.mixture2_id}', '{cls.molecule2_id}', 75.0, 'percent', '{cls.user2_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")

        # Create experiment properties
        if cls.experiment1_id and cls.experiment2_id:
            # Check if experiment properties exist
            sql = f"SELECT COUNT(*) FROM experiment_properties WHERE experiment_id = '{cls.experiment1_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO experiment_properties (id, experiment_id, property_type_id, numeric_value, unit, created_by)
                    VALUES ('{str(uuid.uuid4())}', '{cls.experiment1_id}', '{str(uuid.uuid4())}', 25.0, 'C', '{cls.user1_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")

            sql = f"SELECT COUNT(*) FROM experiment_properties WHERE experiment_id = '{cls.experiment2_id}'"
            result = cls.helper.execute_sql(sql, role="service_role")
            if 'error' not in result and result.get('data', [{}])[0].get('count', 0) == 0:
                sql = f"""
                    INSERT INTO experiment_properties (id, experiment_id, property_type_id, numeric_value, unit, created_by)
                    VALUES ('{str(uuid.uuid4())}', '{cls.experiment2_id}', '{str(uuid.uuid4())}', 30.0, 'C', '{cls.user2_id}')
                """
                cls.helper.execute_sql(sql, role="service_role")
        # Create molecules with properties
        sql = f"""
            -- Get or create a molecule for project1
            WITH mol1 AS (
                SELECT id FROM molecules WHERE project_id = '{cls.project1_id}' LIMIT 1
            ), new_mol1 AS (
                INSERT INTO molecules (id, project_id, name, smiles)
                SELECT gen_random_uuid(), '{cls.project1_id}', 'Related Test Molecule 1', 'CC'
                WHERE NOT EXISTS (SELECT 1 FROM mol1)
                RETURNING id
            )
            SELECT id FROM mol1
            UNION ALL
            SELECT id FROM new_mol1;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.molecule1_id = result['data'][0]['id']

            # Create molecular property for molecule1
            sql = f"""
                INSERT INTO molecular_properties (molecule_id, property_name, property_value, property_unit)
                VALUES ('{cls.molecule1_id}', 'Test Property', '123.45', 'units')
                ON CONFLICT (molecule_id, property_name) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")

        # Get or create a molecule for project2
        sql = f"""
            WITH mol2 AS (
                SELECT id FROM molecules WHERE project_id = '{cls.project2_id}' LIMIT 1
            ), new_mol2 AS (
                INSERT INTO molecules (id, project_id, name, smiles)
                SELECT gen_random_uuid(), '{cls.project2_id}', 'Related Test Molecule 2', 'CCC'
                WHERE NOT EXISTS (SELECT 1 FROM mol2)
                RETURNING id
            )
            SELECT id FROM mol2
            UNION ALL
            SELECT id FROM new_mol2;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.molecule2_id = result['data'][0]['id']

            # Create molecular property for molecule2
            sql = f"""
                INSERT INTO molecular_properties (molecule_id, property_name, property_value, property_unit)
                VALUES ('{cls.molecule2_id}', 'Test Property 2', '67.89', 'units')
                ON CONFLICT (molecule_id, property_name) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")

        # Create mixtures with components
        sql = f"""
            -- Get or create a mixture for project1
            WITH mix1 AS (
                SELECT id FROM mixtures WHERE project_id = '{cls.project1_id}' LIMIT 1
            ), new_mix1 AS (
                INSERT INTO mixtures (id, project_id, name)
                SELECT gen_random_uuid(), '{cls.project1_id}', 'Related Test Mixture 1'
                WHERE NOT EXISTS (SELECT 1 FROM mix1)
                RETURNING id
            )
            SELECT id FROM mix1
            UNION ALL
            SELECT id FROM new_mix1;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.mixture1_id = result['data'][0]['id']

            # Create mixture component for mixture1
            sql = f"""
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
                VALUES ('{cls.mixture1_id}', '{cls.molecule1_id}', 50.0)
                ON CONFLICT (mixture_id, molecule_id) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")

        # Create mixtures with components for project2
        sql = f"""
            -- Get or create a mixture for project2
            WITH mix2 AS (
                SELECT id FROM mixtures WHERE project_id = '{cls.project2_id}' LIMIT 1
            ), new_mix2 AS (
                INSERT INTO mixtures (id, project_id, name)
                SELECT gen_random_uuid(), '{cls.project2_id}', 'Related Test Mixture 2'
                WHERE NOT EXISTS (SELECT 1 FROM mix2)
                RETURNING id
            )
            SELECT id FROM mix2
            UNION ALL
            SELECT id FROM new_mix2;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.mixture2_id = result['data'][0]['id']

            # Create mixture component for mixture2
            sql = f"""
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
                VALUES ('{cls.mixture2_id}', '{cls.molecule2_id}', 75.0)
                ON CONFLICT (mixture_id, molecule_id) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")

        # Create experiments with properties
        sql = f"""
            -- Get or create an experiment for project1
            WITH exp1 AS (
                SELECT id FROM experiments WHERE project_id = '{cls.project1_id}' LIMIT 1
            ), new_exp1 AS (
                INSERT INTO experiments (id, project_id, name, description)
                SELECT gen_random_uuid(), '{cls.project1_id}', 'Related Test Experiment 1', 'Test experiment'
                WHERE NOT EXISTS (SELECT 1 FROM exp1)
                RETURNING id
            )
            SELECT id FROM exp1
            UNION ALL
            SELECT id FROM new_exp1;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.experiment1_id = result['data'][0]['id']

            # Create experiment property for experiment1
            sql = f"""
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, property_unit)
                VALUES ('{cls.experiment1_id}', 'Test Exp Property', '42.0', 'units')
                ON CONFLICT (experiment_id, property_name) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")

        # Create experiments with properties for project2
        sql = f"""
            -- Get or create an experiment for project2
            WITH exp2 AS (
                SELECT id FROM experiments WHERE project_id = '{cls.project2_id}' LIMIT 1
            ), new_exp2 AS (
                INSERT INTO experiments (id, project_id, name, description)
                SELECT gen_random_uuid(), '{cls.project2_id}', 'Related Test Experiment 2', 'Test experiment 2'
                WHERE NOT EXISTS (SELECT 1 FROM exp2)
                RETURNING id
            )
            SELECT id FROM exp2
            UNION ALL
            SELECT id FROM new_exp2;
        """
        result = cls.helper.execute_sql(sql, role="service_role")
        if 'error' not in result and len(result.get('data', [])) > 0:
            cls.experiment2_id = result['data'][0]['id']

            # Create experiment property for experiment2
            sql = f"""
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, property_unit)
                VALUES ('{cls.experiment2_id}', 'Test Exp Property 2', '84.0', 'units')
                ON CONFLICT (experiment_id, property_name) DO NOTHING
            """
            cls.helper.execute_sql(sql, role="service_role")
    
    @classmethod
    def tearDownClass(cls):
        """Clean up test data after all tests"""
        # Uncomment the following lines to clean up test data
        # Don't clean up by default to allow manual inspection
        # cls.helper.cleanup_test_data(cls.test_data)
        # logger.info("Test data cleaned up")
        pass
    
    def test_molecular_property_user_project_access(self):
        """TC-RLS-020: Test access to molecular_property for molecule in user's project"""
        # Test user1 can select molecular properties for molecule1
        result = self.helper.test_user_access(
            self.user1_id, 'molecular_properties', f"molecule_id = '{self.molecule1_id}'"
        )
        self.assertGreater(result.get('count', 0), 0)
        self.assertNotIn('error', result)

    def test_molecular_property_non_user_project_access(self):
        """TC-RLS-021: Test access to molecular_property for molecule not in user's project"""
        # Test user1 cannot select molecular properties for molecule2
        result = self.helper.test_user_access(
            self.user1_id, 'molecular_properties', f"molecule_id = '{self.molecule2_id}'"
        )
        self.assertEqual(result.get('count', -1), 0)

    def test_mixture_component_user_project_access(self):
        """TC-RLS-022: Test access to mixture_component for mixture in user's project"""
        # Test user1 can select mixture components for mixture1
        result = self.helper.test_user_access(
            self.user1_id, 'mixture_components', f"mixture_id = '{self.mixture1_id}'"
        )
        self.assertGreater(result.get('count', 0), 0)
        self.assertNotIn('error', result)

    def test_mixture_component_non_user_project_access(self):
        """TC-RLS-023: Test access to mixture_component for mixture not in user's project"""
        # Test user1 cannot select mixture components for mixture2
        result = self.helper.test_user_access(
            self.user1_id, 'mixture_components', f"mixture_id = '{self.mixture2_id}'"
        )
        self.assertEqual(result.get('count', -1), 0)

    def test_experiment_property_user_project_access(self):
        """TC-RLS-024: Test access to experiment_property for experiment in user's project"""
        # Test user1 can select experiment properties for experiment1
        result = self.helper.test_user_access(
            self.user1_id, 'experiment_properties', f"experiment_id = '{self.experiment1_id}'"
        )
        self.assertGreater(result.get('count', 0), 0)
        self.assertNotIn('error', result)

    def test_experiment_property_non_user_project_access(self):
        """TC-RLS-025: Test access to experiment_property for experiment not in user's project"""
        # Test user1 cannot select experiment properties for experiment2
        result = self.helper.test_user_access(
            self.user1_id, 'experiment_properties', f"experiment_id = '{self.experiment2_id}'"
        )
        self.assertEqual(result.get('count', -1), 0)

class RLSServiceRoleAccessTests(unittest.TestCase):
    """Tests for Service Role Access (TC-RLS-030 to TC-RLS-033)"""

    @classmethod
    def setUpClass(cls):
        """Set up test data for all tests"""
        cls.helper = RLSTestHelper(DB_CONNECTION_PARAMS)
        cls.test_data = cls.helper.create_test_data()
        logger.info(f"Test data created: {json.dumps(cls.test_data)}")

        # Store user and project IDs for easy access
        cls.user1_id = cls.test_data['users'][0]
        cls.user2_id = cls.test_data['users'][1]
        cls.project1_id = cls.test_data['projects'][0]
        cls.project2_id = cls.test_data['projects'][1]
        cls.team_id = cls.test_data['team']
        cls.molecule1_id = cls.test_data['molecules'][0] if cls.test_data.get('molecules') else None
        cls.molecule2_id = cls.test_data['molecules'][1] if len(cls.test_data.get('molecules', [])) > 1 else None
        cls.experiment1_id = cls.test_data['experiments'][0] if cls.test_data.get('experiments') else None
        cls.experiment2_id = cls.test_data['experiments'][1] if len(cls.test_data.get('experiments', [])) > 1 else None
        cls.mixture1_id = cls.test_data['mixtures'][0] if cls.test_data.get('mixtures') else None
        cls.mixture2_id = cls.test_data['mixtures'][1] if len(cls.test_data.get('mixtures', [])) > 1 else None

    @classmethod
    def tearDownClass(cls):
        """Clean up test data after all tests"""
        # Uncomment the following lines to clean up test data
        # Don't clean up by default to allow manual inspection
        # cls.helper.cleanup_test_data(cls.test_data)
        # logger.info("Test data cleaned up")
        pass
    
    def test_service_role_select_access(self):
        """TC-RLS-030: Test SELECT access as service role to all tables"""
        # Test tables that should be accessible
        tables = [
            'molecules', 'mixtures', 'experiments', 'projects', 'user_profile',
            'molecular_properties', 'mixture_components', 'experiment_properties'
        ]

        for table in tables:
            result = self.helper.test_service_role_access(table)
            self.assertNotIn('error', result, f"Error accessing {table}: {result.get('error', '')}")

    def test_service_role_insert_access(self):
        """TC-RLS-031: Test INSERT access as service role to all tables"""
        # Test insert into molecule
        molecule_id = str(uuid.uuid4())
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Service Role Test Molecule', 'CCCC', '{self.user1_id}')
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="service_role")
        self.assertNotIn('error', result)

        # Test insert into mixture
        mixture_id = str(uuid.uuid4())
        sql = f"""
            INSERT INTO mixtures (id, project_id, name, created_by)
            VALUES ('{mixture_id}', '{self.project1_id}', 'Service Role Test Mixture', '{self.user1_id}')
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="service_role")
        self.assertNotIn('error', result)

    def test_service_role_update_access(self):
        """TC-RLS-032: Test UPDATE access as service role to all tables"""
        # Find a molecule to update
        sql = f"SELECT id FROM molecules LIMIT 1"
        result = self.helper.execute_sql(sql, role="service_role")

        if 'error' not in result and len(result.get('data', [])) > 0:
            molecule_id = result['data'][0]['id']

            # Test update molecule
            sql = f"""
                UPDATE molecules
                SET name = 'Service Role Updated Test Molecule'
                WHERE id = '{molecule_id}'
                RETURNING id
            """
            result = self.helper.execute_sql(sql, role="service_role")
            self.assertNotIn('error', result)

    def test_service_role_delete_access(self):
        """TC-RLS-033: Test DELETE access as service role to all tables"""
        # Create a molecule to delete
        molecule_id = str(uuid.uuid4())
        sql = f"""
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('{molecule_id}', 'Service Role Delete Test', 'CCCCC', '{self.user1_id}')
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="service_role")
        self.assertNotIn('error', result)

        # Test delete molecule
        sql = f"""
            DELETE FROM molecules
            WHERE id = '{molecule_id}'
            RETURNING id
        """
        result = self.helper.execute_sql(sql, role="service_role")
        self.assertNotIn('error', result)

# Import uuid for generating test IDs
import uuid

if __name__ == '__main__':
    unittest.main()
"""
RLS Test Helper - Utility to test RLS policies
"""
import psycopg2
from psycopg2.extras import RealDictCursor
import uuid
import json
import os
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='rls_tests.log'
)
logger = logging.getLogger(__name__)

class RLSTestHelper:
    def __init__(self, conn_params):
        """Initialize RLS test helper with database connection parameters"""
        self.conn_params = conn_params
        logger.info("RLSTestHelper initialized with connection parameters")
        
    def create_test_data(self):
        """Create test data for RLS testing"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            conn.autocommit = False
            cursor = conn.cursor()
            
            logger.info("Creating test data for RLS testing")
            
            # Create test users (using fixed UUIDs for consistent testing)
            user1_id = "00000000-0000-0000-0000-000000000001"
            user2_id = "00000000-0000-0000-0000-000000000002"

            # Check if users exist in the auth.users table
            cursor.execute(f"SELECT id FROM auth.users WHERE id = '{user1_id}'")
            if not cursor.fetchone():
                # Create entry in auth.users table for user1
                cursor.execute(f"""
                    INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, email_confirmed_at)
                    VALUES (
                        '{user1_id}',
                        'testuser1@example.com',
                        CURRENT_TIMESTAMP,
                        CURRENT_TIMESTAMP,
                        'authenticated',
                        'authenticated',
                        CURRENT_TIMESTAMP
                    )
                """)
                logger.info(f"Created test user1 in auth.users table with ID {user1_id}")

            cursor.execute(f"SELECT id FROM auth.users WHERE id = '{user2_id}'")
            if not cursor.fetchone():
                # Create entry in auth.users table for user2
                cursor.execute(f"""
                    INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, email_confirmed_at)
                    VALUES (
                        '{user2_id}',
                        'testuser2@example.com',
                        CURRENT_TIMESTAMP,
                        CURRENT_TIMESTAMP,
                        'authenticated',
                        'authenticated',
                        CURRENT_TIMESTAMP
                    )
                """)
                logger.info(f"Created test user2 in auth.users table with ID {user2_id}")

            # Create test projects
            project1_id = "11111111-1111-1111-1111-111111111111"
            project2_id = "22222222-2222-2222-2222-222222222222"
            
            # First check if projects exist
            cursor.execute(f"""
                SELECT COUNT(*) FROM projects
                WHERE id IN ('{project1_id}', '{project2_id}')
            """)
            project_count = cursor.fetchone()[0]

            if project_count < 2:
                # Create test projects if they don't exist
                cursor.execute(f"""
                    INSERT INTO projects (id, name, description)
                    VALUES
                      ('{project1_id}', 'Test Project 1', 'Project for RLS testing'),
                      ('{project2_id}', 'Test Project 2', 'Another project for RLS testing')
                    ON CONFLICT (id) DO NOTHING
                """)
            
            # Check if user profiles exist
            cursor.execute(f"""
                SELECT COUNT(*) FROM user_profile
                WHERE auth_user_id IN ('{user1_id}', '{user2_id}')
            """)
            profile_count = cursor.fetchone()[0]

            if profile_count < 2:
                # Check for existing profiles one by one to avoid ON CONFLICT
                cursor.execute(f"SELECT id FROM user_profile WHERE auth_user_id = '{user1_id}'")
                if not cursor.fetchone():
                    cursor.execute(f"""
                        INSERT INTO user_profile (auth_user_id, display_name, email)
                        VALUES ('{user1_id}', 'Test User 1', 'testuser1@example.com')
                    """)

                cursor.execute(f"SELECT id FROM user_profile WHERE auth_user_id = '{user2_id}'")
                if not cursor.fetchone():
                    cursor.execute(f"""
                        INSERT INTO user_profile (auth_user_id, display_name, email)
                        VALUES ('{user2_id}', 'Test User 2', 'testuser2@example.com')
                    """)

                # Re-check profile count
                cursor.execute(f"""
                    SELECT COUNT(*) FROM user_profile
                    WHERE auth_user_id IN ('{user1_id}', '{user2_id}')
                """)
                profile_count = cursor.fetchone()[0]
                logger.info(f"Created {profile_count} user profiles")

            # Create team for testing
            team_id = "33333333-3333-3333-3333-333333333333"

            # Check if team exists first
            cursor.execute(f"SELECT id FROM teams WHERE id = '{team_id}'")
            if not cursor.fetchone():
                cursor.execute(f"""
                    INSERT INTO teams (id, name, description)
                    VALUES ('{team_id}', 'Test Team', 'Team for RLS testing')
                """)
                logger.info(f"Created test team with ID {team_id}")

            # Link users to team
            # Check if user1 is already linked to team
            cursor.execute(f"SELECT id FROM team_members WHERE team_id = '{team_id}' AND user_id = '{user1_id}'")
            if not cursor.fetchone():
                cursor.execute(f"""
                    INSERT INTO team_members (team_id, user_id, role)
                    VALUES ('{team_id}', '{user1_id}', 'admin')
                """)
                logger.info(f"Linked user1 {user1_id} to team {team_id}")

            # Check if user2 is already linked to team
            cursor.execute(f"SELECT id FROM team_members WHERE team_id = '{team_id}' AND user_id = '{user2_id}'")
            if not cursor.fetchone():
                cursor.execute(f"""
                    INSERT INTO team_members (team_id, user_id, role)
                    VALUES ('{team_id}', '{user2_id}', 'viewer')
                """)
                logger.info(f"Linked user2 {user2_id} to team {team_id}")

            # Link projects to team
            cursor.execute(f"""
                UPDATE projects SET team_id = '{team_id}'
                WHERE id IN ('{project1_id}', '{project2_id}')
            """)
            
            # Create test molecules
            # Molecules don't have project_id, but we'll still use two distinct molecules
            cursor.execute(f"SELECT COUNT(*) FROM molecules WHERE name = 'Test Molecule 1' OR name = 'Test Molecule 2'")
            molecule_count = cursor.fetchone()[0]

            if molecule_count < 2:
                # Create separate molecules for each project
                molecule1_id = str(uuid.uuid4())
                molecule2_id = str(uuid.uuid4())

                # Create molecules (no project_id)
                cursor.execute(f"""
                    INSERT INTO molecules (id, name, smiles, created_by)
                    VALUES
                      ('{molecule1_id}', 'Test Molecule 1', 'C', '{user1_id}')
                """)
                logger.info(f"Created test molecule 1 with ID {molecule1_id}")

                cursor.execute(f"""
                    INSERT INTO molecules (id, name, smiles, created_by)
                    VALUES
                      ('{molecule2_id}', 'Test Molecule 2', 'CC', '{user2_id}')
                """)
                logger.info(f"Created test molecule 2 with ID {molecule2_id}")

                # Now create experiments that link molecules to projects
                experiment1_id = str(uuid.uuid4())
                cursor.execute(f"""
                    INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                    VALUES ('{experiment1_id}', '{molecule1_id}', '{project1_id}', 'Experiment 1', '{user1_id}')
                """)
                logger.info(f"Created experiment 1 linking molecule 1 to project 1")

                experiment2_id = str(uuid.uuid4())
                cursor.execute(f"""
                    INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
                    VALUES ('{experiment2_id}', '{molecule2_id}', '{project2_id}', 'Experiment 2', '{user2_id}')
                """)
                logger.info(f"Created experiment 2 linking molecule 2 to project 2")

            # Check if mixtures exist
            cursor.execute(f"""
                SELECT COUNT(*) FROM mixtures
                WHERE project_id IN ('{project1_id}', '{project2_id}')
            """)
            mixture_count = cursor.fetchone()[0]

            if mixture_count < 2:
                # Create test mixtures
                mixture1_id = str(uuid.uuid4())
                mixture2_id = str(uuid.uuid4())

                # Mixtures have project_id
                cursor.execute(f"""
                    INSERT INTO mixtures (id, project_id, name, created_by)
                    VALUES ('{mixture1_id}', '{project1_id}', 'Test Mixture 1', '{user1_id}')
                """)
                logger.info(f"Created test mixture 1 with ID {mixture1_id}")

                cursor.execute(f"""
                    INSERT INTO mixtures (id, project_id, name, created_by)
                    VALUES ('{mixture2_id}', '{project2_id}', 'Test Mixture 2', '{user2_id}')
                """)
                logger.info(f"Created test mixture 2 with ID {mixture2_id}")
            
            # Commit changes
            conn.commit()
            logger.info("Test data created successfully")
            
            # Return test data IDs
            # Store molecule and experiment IDs if they were created
            if molecule_count < 2 and 'molecule1_id' in locals() and 'molecule2_id' in locals():
                molecules = [molecule1_id, molecule2_id]
                experiments = [experiment1_id, experiment2_id]
            else:
                # Retrieve existing molecule IDs
                cursor.execute("SELECT id FROM molecules WHERE name = 'Test Molecule 1' LIMIT 1")
                result = cursor.fetchone()
                molecule1_id = result[0] if result else None

                cursor.execute("SELECT id FROM molecules WHERE name = 'Test Molecule 2' LIMIT 1")
                result = cursor.fetchone()
                molecule2_id = result[0] if result else None

                # Retrieve existing experiment IDs
                cursor.execute(f"SELECT id FROM experiments WHERE project_id = '{project1_id}' LIMIT 1")
                result = cursor.fetchone()
                experiment1_id = result[0] if result else None

                cursor.execute(f"SELECT id FROM experiments WHERE project_id = '{project2_id}' LIMIT 1")
                result = cursor.fetchone()
                experiment2_id = result[0] if result else None

                molecules = [molecule1_id, molecule2_id] if molecule1_id and molecule2_id else []
                experiments = [experiment1_id, experiment2_id] if experiment1_id and experiment2_id else []

            # Retrieve mixture IDs
            cursor.execute(f"SELECT id FROM mixtures WHERE project_id = '{project1_id}' LIMIT 1")
            result = cursor.fetchone()
            mixture1_id = result[0] if result else None

            cursor.execute(f"SELECT id FROM mixtures WHERE project_id = '{project2_id}' LIMIT 1")
            result = cursor.fetchone()
            mixture2_id = result[0] if result else None

            mixtures = [mixture1_id, mixture2_id] if mixture1_id and mixture2_id else []

            return {
                'users': [user1_id, user2_id],
                'projects': [project1_id, project2_id],
                'team': team_id,
                'molecules': molecules,
                'experiments': experiments,
                'mixtures': mixtures
            }
        except Exception as e:
            # Rollback on error
            if conn:
                conn.rollback()
            logger.error(f"Error creating test data: {str(e)}")
            raise e
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def test_user_access(self, user_id, table_name, condition=None):
        """Test user access to a table"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Set auth context
            cursor.execute(f"SET LOCAL ROLE authenticated")
            cursor.execute(f"SET LOCAL auth.uid = '{user_id}'")
            logger.info(f"Testing user {user_id} access to {table_name}")
            
            # Execute query
            query = f"SELECT * FROM {table_name}"
            if condition:
                query += f" WHERE {condition}"
            logger.info(f"Executing query: {query}")
            cursor.execute(query)
            
            # Get results
            results = cursor.fetchall()
            logger.info(f"Query returned {len(results)} results")
            
            # Return result count and data
            return {
                'count': len(results),
                'data': [dict(row) for row in results]
            }
        except Exception as e:
            # Handle error
            logger.error(f"Error testing user access: {str(e)}")
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def test_security_definer_function(self, user_id, function_name, args=None):
        """Test a security definer function"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor()
            
            # Set auth context
            cursor.execute(f"SET LOCAL ROLE authenticated")
            cursor.execute(f"SET LOCAL auth.uid = '{user_id}'")
            logger.info(f"Testing function {function_name} for user {user_id}")
            
            # Prepare function call
            if args:
                arg_list = ", ".join([f"'{arg}'" if isinstance(arg, str) else str(arg) for arg in args])
                query = f"SELECT {function_name}({arg_list})"
            else:
                query = f"SELECT {function_name}()"
                
            logger.info(f"Executing function: {query}")
            cursor.execute(query)
            
            # Get result
            result = cursor.fetchone()[0]
            logger.info(f"Function returned: {result}")
            
            return result
        except Exception as e:
            # Handle error
            logger.error(f"Error testing security definer function: {str(e)}")
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def test_service_role_access(self, table_name):
        """Test service role access to a table"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Set auth context
            cursor.execute(f"SET ROLE service_role")
            logger.info(f"Testing service role access to {table_name}")
            
            # Execute query
            cursor.execute(f"SELECT * FROM {table_name} LIMIT 10")
            
            # Get results
            results = cursor.fetchall()
            logger.info(f"Query returned {len(results)} results")
            
            # Return result count and data
            return {
                'count': len(results),
                'data': [dict(row) for row in results]
            }
        except Exception as e:
            # Handle error
            logger.error(f"Error testing service role access: {str(e)}")
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def execute_sql(self, sql, role="authenticated", user_id=None):
        """Execute SQL with specific role and user_id"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Set auth context
            cursor.execute(f"SET ROLE {role}")
            if user_id:
                cursor.execute(f"SET LOCAL auth.uid = '{user_id}'")
            logger.info(f"Executing SQL as {role} (user_id: {user_id}): {sql[:100]}...")
            
            # Execute query
            cursor.execute(sql)
            
            # Try to fetch results
            try:
                results = cursor.fetchall()
                logger.info(f"Query returned {len(results)} results")
                return {
                    'count': len(results),
                    'data': [dict(row) for row in results]
                }
            except psycopg2.ProgrammingError:
                # No results to fetch (for INSERT, UPDATE, etc.)
                rowcount = cursor.rowcount
                logger.info(f"Query affected {rowcount} rows")
                return {
                    'rowcount': rowcount
                }
                
        except Exception as e:
            # Handle error
            logger.error(f"Error executing SQL: {str(e)}")
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def cleanup_test_data(self, test_data):
        """Clean up test data"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            conn.autocommit = False
            cursor = conn.cursor()

            logger.info("Cleaning up test data")

            # Get the team ID (used for cleanup)
            team_id = "33333333-3333-3333-3333-333333333333"

            # Clean up data in reverse order of dependencies
            # Delete experiments
            for experiment_id in test_data.get('experiments', []):
                cursor.execute(f"DELETE FROM experiments WHERE id = '{experiment_id}'")
                logger.info(f"Deleted experiment {experiment_id}")

            # Delete mixtures
            for mixture_id in test_data.get('mixtures', []):
                cursor.execute(f"DELETE FROM mixtures WHERE id = '{mixture_id}'")
                logger.info(f"Deleted mixture {mixture_id}")

            # Delete molecules
            for molecule_id in test_data.get('molecules', []):
                cursor.execute(f"DELETE FROM molecules WHERE id = '{molecule_id}'")
                logger.info(f"Deleted molecule {molecule_id}")

            # Delete projects
            for project_id in test_data['projects']:
                cursor.execute(f"DELETE FROM projects WHERE id = '{project_id}'")
                logger.info(f"Deleted project {project_id}")

            # Delete team members
            for user_id in test_data['users']:
                # First clean up team_members which has FK to users
                cursor.execute(f"DELETE FROM team_members WHERE user_id = '{user_id}'")
                logger.info(f"Deleted team members for user {user_id}")

                # Delete user profiles
                cursor.execute(f"DELETE FROM user_profile WHERE auth_user_id = '{user_id}'")
                logger.info(f"Deleted user profile for user {user_id}")

                # Delete users
                cursor.execute(f"DELETE FROM auth.users WHERE id = '{user_id}'")
                logger.info(f"Deleted user {user_id} from auth.users")

            # Delete team
            cursor.execute(f"DELETE FROM teams WHERE id = '{team_id}'")
            logger.info(f"Deleted team {team_id}")

            # Commit changes
            conn.commit()
            logger.info("Test data cleanup completed successfully")

            return True
        except Exception as e:
            # Rollback on error
            if conn:
                conn.rollback()
            logger.error(f"Error cleaning up test data: {str(e)}")
            raise e
        finally:
            # Close connection
            if conn:
                conn.close()
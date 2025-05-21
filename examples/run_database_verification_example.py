"""
Example script to show how to use the database verification framework.

This script demonstrates how to run specific database verification tests 
and analyze the results. It can be used as a starting point for custom 
verification workflows.
"""
import os
import sys
import argparse
import logging
from dotenv import load_dotenv

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import verification modules
from tests.rls_test_helper import RLSTestHelper

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='database_verification_example.log'
)
logger = logging.getLogger(__name__)

def test_connection_and_setup():
    """Test database connection and create test data"""
    # Load environment variables
    load_dotenv()
    
    # Get database connection parameters
    db_params = {
        'host': os.environ.get('SUPABASE_DB_HOST'),
        'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
        'dbname': os.environ.get('SUPABASE_DB_NAME'),
        'user': os.environ.get('SUPABASE_DB_USER'),
        'password': os.environ.get('SUPABASE_DB_PASSWORD'),
        'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Print connection parameters (without password)
    print("\nDatabase Connection Parameters:")
    for key, value in db_params.items():
        if key != 'password':
            print(f"  {key}: {value}")
    
    # Create RLS test helper
    helper = RLSTestHelper(db_params)
    
    # Create test data
    print("\nCreating test data...")
    test_data = helper.create_test_data()
    print(f"Test data created: {test_data}")
    
    # Test security definer function
    print("\nTesting is_project_member function...")
    result = helper.test_security_definer_function(
        test_data['users'][0], 'is_project_member', [test_data['projects'][0]]
    )
    print(f"Result for user's own project: {result}")
    
    result = helper.test_security_definer_function(
        test_data['users'][0], 'is_project_member', [test_data['projects'][1]]
    )
    print(f"Result for another user's project: {result}")
    
    # Test user access to molecule table
    print("\nTesting user access to molecule table...")
    result = helper.test_user_access(
        test_data['users'][0], 'molecule', f"project_id = '{test_data['projects'][0]}'"
    )
    print(f"User can access own project's molecules: {result.get('count', 0)} found")
    
    result = helper.test_user_access(
        test_data['users'][0], 'molecule', f"project_id = '{test_data['projects'][1]}'"
    )
    print(f"User can access another project's molecules: {result.get('count', 0)} found")
    
    # Test service role access
    print("\nTesting service role access...")
    result = helper.test_service_role_access('molecule')
    print(f"Service role can access molecules: {result.get('count', 0)} found")
    
    return test_data, helper

def test_sql_queries(helper):
    """Test SQL queries with different roles"""
    print("\nTesting SQL queries with different roles...")
    
    # Test as authenticated user
    print("\nAs authenticated user:")
    result = helper.execute_sql(
        "SELECT COUNT(*) FROM molecule",
        role="authenticated",
        user_id="00000000-0000-0000-0000-000000000001"
    )
    print(f"Molecule count: {result.get('data', [{}])[0].get('count', 0)}")
    
    # Test as service role
    print("\nAs service role:")
    result = helper.execute_sql(
        "SELECT COUNT(*) FROM molecule",
        role="service_role"
    )
    print(f"Molecule count: {result.get('data', [{}])[0].get('count', 0)}")
    
    # Try to insert as authenticated user
    print("\nInserting as authenticated user:")
    result = helper.execute_sql(
        """
        INSERT INTO molecule (id, project_id, name, smiles)
        VALUES (gen_random_uuid(), '11111111-1111-1111-1111-111111111111', 'Test Insert', 'C')
        RETURNING id
        """,
        role="authenticated",
        user_id="00000000-0000-0000-0000-000000000001"
    )
    if 'error' in result:
        print(f"Error: {result['error']}")
    else:
        print(f"Inserted molecule with ID: {result.get('data', [{}])[0].get('id', 'unknown')}")

def analyze_rls_policies(helper):
    """Analyze RLS policies in the database"""
    print("\nAnalyzing RLS policies...")
    
    # Get RLS policies
    result = helper.execute_sql(
        """
        SELECT tablename, policyname, cmd, roles, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        ORDER BY tablename, policyname
        """,
        role="service_role"
    )
    
    if 'error' in result:
        print(f"Error: {result['error']}")
        return
    
    policies = result.get('data', [])
    print(f"Found {len(policies)} RLS policies")
    
    # Group policies by table
    policies_by_table = {}
    for policy in policies:
        table = policy['tablename']
        if table not in policies_by_table:
            policies_by_table[table] = []
        policies_by_table[table].append(policy)
    
    # Print policies by table
    for table, table_policies in policies_by_table.items():
        print(f"\nTable: {table}")
        for policy in table_policies:
            print(f"  Policy: {policy['policyname']}")
            print(f"  Command: {policy['cmd']}")
            print(f"  Roles: {policy['roles']}")
            
            # Check if policy uses security definer functions
            if policy['qual']:
                using_secdef = any(func in policy['qual'] for func in [
                    'is_project_member', 'is_project_owner', 'user_projects',
                    'is_team_member', 'molecule_in_user_project', 'mixture_in_user_project'
                ])
                print(f"  Using Security Definer: {'Yes' if using_secdef else 'No'}")
            
            print()

def main():
    """Main function to run the database verification example"""
    parser = argparse.ArgumentParser(description='Database Verification Example')
    parser.add_argument('--test', choices=['connection', 'rls', 'sql', 'all'], default='all',
                      help='Which test to run (default: all)')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Database Verification Example")
    print("=" * 60)
    
    try:
        if args.test in ['connection', 'all']:
            # Test connection and setup
            test_data, helper = test_connection_and_setup()
            
            if args.test == 'all':
                # Run other tests
                test_sql_queries(helper)
                analyze_rls_policies(helper)
        
        elif args.test == 'rls':
            # Load environment variables
            load_dotenv()
            
            # Get database connection parameters
            db_params = {
                'host': os.environ.get('SUPABASE_DB_HOST'),
                'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
                'dbname': os.environ.get('SUPABASE_DB_NAME'),
                'user': os.environ.get('SUPABASE_DB_USER'),
                'password': os.environ.get('SUPABASE_DB_PASSWORD'),
                'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
            }
            
            # Create RLS test helper
            helper = RLSTestHelper(db_params)
            
            # Create test data
            test_data = helper.create_test_data()
            
            # Analyze RLS policies
            analyze_rls_policies(helper)
        
        elif args.test == 'sql':
            # Load environment variables
            load_dotenv()
            
            # Get database connection parameters
            db_params = {
                'host': os.environ.get('SUPABASE_DB_HOST'),
                'port': os.environ.get('SUPABASE_DB_PORT', '5432'),
                'dbname': os.environ.get('SUPABASE_DB_NAME'),
                'user': os.environ.get('SUPABASE_DB_USER'),
                'password': os.environ.get('SUPABASE_DB_PASSWORD'),
                'sslmode': os.environ.get('SUPABASE_DB_SSLMODE', 'require')
            }
            
            # Create RLS test helper
            helper = RLSTestHelper(db_params)
            
            # Test SQL queries
            test_sql_queries(helper)
        
        print("\nExample completed successfully!")
    
    except Exception as e:
        print(f"Error running example: {str(e)}")
        logger.error(f"Error running example: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == '__main__':
    main()
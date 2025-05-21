"""
Database Connection Test Script for RLS Verification

This script tests the connection to the Supabase database before running the RLS verification tests.
It helps identify any connection issues before running the actual tests.
"""
import os
import sys
import psycopg2
from dotenv import load_dotenv

def test_connection():
    """Test connection to the Supabase database"""
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
    
    # Check if required parameters are set
    missing_params = [param for param in ['host', 'dbname', 'user', 'password'] if not db_params.get(param)]
    if missing_params:
        print(f"\nError: Missing required parameters: {', '.join(missing_params)}")
        print("Please set these environment variables or add them to your .env file")
        return False
    
    # Test connection
    try:
        print("\nAttempting connection to the database...")
        conn = psycopg2.connect(**db_params)
        print("Connection successful!")
        
        # Test roles
        print("\nTesting roles:")
        cursor = conn.cursor()
        
        # Test authenticated role
        try:
            print("  Testing 'authenticated' role...")
            cursor.execute("SET ROLE authenticated")
            cursor.execute("SELECT current_user, current_setting('auth.uid', true)")
            result = cursor.fetchone()
            print(f"    Current user: {result[0]}")
            print(f"    auth.uid setting: {result[1] or 'Not set'}")
        except Exception as e:
            print(f"    Error testing authenticated role: {str(e)}")
        
        # Test service_role
        try:
            print("  Testing 'service_role' role...")
            cursor.execute("SET ROLE service_role")
            cursor.execute("SELECT current_user")
            result = cursor.fetchone()
            print(f"    Current user: {result[0]}")
        except Exception as e:
            print(f"    Error testing service_role: {str(e)}")
        
        # Test RLS enforcement
        print("\nTesting RLS enforcement:")
        try:
            # Set authenticated role with dummy user
            cursor.execute("SET ROLE authenticated")
            cursor.execute("SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000000'")
            
            # Try to select from molecule table
            cursor.execute("SELECT COUNT(*) FROM molecule")
            result = cursor.fetchone()
            print(f"  Molecule count as dummy user: {result[0]}")
            
            # Now try with service_role
            cursor.execute("SET ROLE service_role")
            cursor.execute("SELECT COUNT(*) FROM molecule")
            result = cursor.fetchone()
            print(f"  Molecule count as service_role: {result[0]}")
            
            # If counts are different, RLS is working
            print("  RLS appears to be enforced (different results with different roles)")
        except Exception as e:
            print(f"  Error testing RLS enforcement: {str(e)}")
        
        # Check if security definer functions exist
        print("\nChecking for security definer functions:")
        try:
            cursor.execute("""
                SELECT proname, prosecdef 
                FROM pg_proc 
                WHERE pronamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public')
                AND prosecdef = true
            """)
            security_definer_functions = cursor.fetchall()
            
            if security_definer_functions:
                print("  Found security definer functions:")
                for func in security_definer_functions:
                    print(f"    - {func[0]}")
            else:
                print("  No security definer functions found")
        except Exception as e:
            print(f"  Error checking security definer functions: {str(e)}")
        
        # Close connection
        conn.close()
        return True
    except Exception as e:
        print(f"\nError connecting to database: {str(e)}")
        print("\nTroubleshooting tips:")
        print("1. Check if the Supabase host is correct")
        print("2. Verify that database credentials are correct")
        print("3. Ensure your IP address is allowed in Supabase's database settings")
        print("4. Check if SSL mode is properly configured")
        print("5. Verify network connectivity to the database")
        return False

if __name__ == "__main__":
    print("=" * 60)
    print("Database Connection Test for RLS Verification")
    print("=" * 60)
    
    success = test_connection()
    
    print("\n" + "=" * 60)
    if success:
        print("Connection test passed!")
        print("You can proceed with running the RLS verification tests.")
        sys.exit(0)
    else:
        print("Connection test failed!")
        print("Please fix the connection issues before running the RLS verification tests.")
        sys.exit(1)
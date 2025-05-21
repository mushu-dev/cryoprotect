"""
List Tables in the Supabase Database
"""
import os
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

def list_tables():
    """List all tables in the database"""
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
    
    conn = None
    try:
        # Connect to database
        print("\nConnecting to database...")
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # List schemas
        print("\nAvailable Schemas:")
        cursor.execute("""
            SELECT schema_name
            FROM information_schema.schemata
            WHERE schema_name NOT IN ('information_schema', 'pg_catalog', 'pg_toast')
            ORDER BY schema_name
        """)
        schemas = cursor.fetchall()
        for schema in schemas:
            print(f"  - {schema['schema_name']}")
        
        # List tables in public schema
        print("\nTables in public schema:")
        cursor.execute("""
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'public'
            AND table_type = 'BASE TABLE'
            ORDER BY table_name
        """)
        tables = cursor.fetchall()
        if tables:
            for table in tables:
                print(f"  - {table['table_name']}")
        else:
            print("  No tables found in public schema")
        
        # List RLS policies if any
        print("\nRLS Policies:")
        cursor.execute("""
            SELECT tablename, policyname
            FROM pg_policies
            WHERE schemaname = 'public'
            ORDER BY tablename, policyname
        """)
        policies = cursor.fetchall()
        if policies:
            for policy in policies:
                print(f"  - {policy['tablename']}: {policy['policyname']}")
        else:
            print("  No RLS policies found")
        
        # List security definer functions if any
        print("\nSecurity Definer Functions:")
        cursor.execute("""
            SELECT proname
            FROM pg_proc
            JOIN pg_namespace ON pg_proc.pronamespace = pg_namespace.oid
            WHERE nspname = 'public'
            AND prosecdef = true
            ORDER BY proname
        """)
        functions = cursor.fetchall()
        if functions:
            for func in functions:
                print(f"  - {func['proname']}")
        else:
            print("  No security definer functions found")
        
        conn.close()
        print("\nConnection closed.")
    except Exception as e:
        print(f"\nError: {str(e)}")
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    list_tables()
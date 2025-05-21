"""
Verify Optimized RLS Implementation

This script generates a verification plan specific to the optimized RLS implementation.
It analyzes the database schema, existing RLS policies, and security definer functions,
and generates a comprehensive report on the implementation status.
"""
import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import datetime
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='optimized_rls_verification.log'
)
logger = logging.getLogger(__name__)

def get_db_connection():
    """Get database connection based on environment variables"""
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
    
    # Connect to database
    conn = psycopg2.connect(**db_params)
    return conn

def get_security_definer_functions(conn):
    """Get all security definer functions in the public schema"""
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    cursor.execute("""
        SELECT 
            p.proname AS function_name,
            pg_get_functiondef(p.oid) AS function_definition,
            p.prosecdef AS security_definer,
            obj_description(p.oid, 'pg_proc') AS description,
            l.lanname AS language
        FROM pg_proc p
        JOIN pg_namespace n ON p.pronamespace = n.oid
        JOIN pg_language l ON p.prolang = l.oid
        WHERE n.nspname = 'public'
        AND p.prosecdef = true
        ORDER BY p.proname
    """)
    functions = cursor.fetchall()
    cursor.close()
    return functions

def get_rls_policies(conn):
    """Get all RLS policies in the database"""
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    cursor.execute("""
        SELECT 
            schemaname,
            tablename,
            policyname,
            roles,
            cmd,
            qual,
            with_check
        FROM pg_policies
        ORDER BY tablename, policyname
    """)
    policies = cursor.fetchall()
    cursor.close()
    return policies

def get_table_indexes(conn):
    """Get indexes for each table in the public schema"""
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    cursor.execute("""
        SELECT 
            t.relname AS table_name,
            i.relname AS index_name,
            a.attname AS column_name,
            ix.indisunique AS is_unique,
            ix.indisprimary AS is_primary
        FROM pg_class t
        JOIN pg_index ix ON t.oid = ix.indrelid
        JOIN pg_class i ON i.oid = ix.indexrelid
        JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(ix.indkey)
        JOIN pg_namespace n ON n.oid = t.relnamespace
        WHERE t.relkind = 'r'
        AND n.nspname = 'public'
        ORDER BY t.relname, i.relname, a.attnum
    """)
    indexes = cursor.fetchall()
    cursor.close()
    return indexes

def analyze_rls_implementation(conn):
    """Analyze the RLS implementation and generate a verification plan"""
    # Get security definer functions
    functions = get_security_definer_functions(conn)
    logger.info(f"Found {len(functions)} security definer functions")
    
    # Get RLS policies
    policies = get_rls_policies(conn)
    logger.info(f"Found {len(policies)} RLS policies")
    
    # Get table indexes
    indexes = get_table_indexes(conn)
    logger.info(f"Found {len(indexes)} indexes")
    
    # Analyze functions
    function_analysis = {}
    expected_functions = [
        'is_project_member', 'is_project_owner', 'user_projects',
        'is_team_member', 'molecule_in_user_project', 'mixture_in_user_project',
        'experiment_in_user_project'
    ]
    
    found_functions = [f['function_name'] for f in functions]
    
    for func_name in expected_functions:
        function_analysis[func_name] = {
            'exists': func_name in found_functions,
            'definition': next((f['function_definition'] for f in functions if f['function_name'] == func_name), None)
        }
    
    # Analyze policies
    policy_analysis = {}
    for policy in policies:
        table = policy['tablename']
        if table not in policy_analysis:
            policy_analysis[table] = []
        
        # Check if policy uses security definer functions
        using_sec_def = False
        for func_name in found_functions:
            if func_name in (policy['qual'] or '') or func_name in (policy['with_check'] or ''):
                using_sec_def = True
                break
        
        policy_analysis[table].append({
            'name': policy['policyname'],
            'cmd': policy['cmd'],
            'roles': policy['roles'],
            'using_security_definer': using_sec_def,
            'qual': policy['qual'],
            'with_check': policy['with_check']
        })
    
    # Analyze indexes for RLS performance
    index_analysis = {}
    rls_related_columns = ['project_id', 'user_id', 'team_id']
    
    for index in indexes:
        table = index['table_name']
        column = index['column_name']
        
        if column in rls_related_columns:
            if table not in index_analysis:
                index_analysis[table] = []
            
            index_analysis[table].append({
                'column': column,
                'index_name': index['index_name'],
                'is_unique': index['is_unique'],
                'is_primary': index['is_primary']
            })
    
    # Return the analysis
    return {
        'functions': function_analysis,
        'policies': policy_analysis,
        'indexes': index_analysis
    }

def generate_verification_report(analysis):
    """Generate a Markdown verification report"""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = f"reports/optimized_rls_verification_{timestamp}.md"
    
    # Ensure reports directory exists
    os.makedirs("reports", exist_ok=True)
    
    # Calculate statistics
    function_count = len(analysis['functions'])
    function_implemented = sum(1 for f in analysis['functions'].values() if f['exists'])
    function_pct = function_implemented / function_count * 100 if function_count > 0 else 0
    
    policy_count = sum(len(policies) for policies in analysis['policies'].values())
    policies_using_secdef = sum(1 for policies in analysis['policies'].values() 
                              for policy in policies if policy['using_security_definer'])
    policy_pct = policies_using_secdef / policy_count * 100 if policy_count > 0 else 0
    
    indexed_tables = len(analysis['indexes'])
    total_tables = len(analysis['policies'])
    index_pct = indexed_tables / total_tables * 100 if total_tables > 0 else 0
    
    # Generate the report
    with open(report_file, 'w') as f:
        f.write(f"# Optimized RLS Verification Report\n\n")
        f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Summary
        f.write("## Summary\n\n")
        f.write(f"- Security Definer Functions: {function_implemented}/{function_count} implemented ({function_pct:.1f}%)\n")
        f.write(f"- RLS Policies Using Security Definer Functions: {policies_using_secdef}/{policy_count} ({policy_pct:.1f}%)\n")
        f.write(f"- Tables with RLS-Related Indexes: {indexed_tables}/{total_tables} ({index_pct:.1f}%)\n\n")
        
        # Security Definer Functions
        f.write("## Security Definer Functions\n\n")
        f.write("| Function Name | Status | Notes |\n")
        f.write("|--------------|--------|-------|\n")
        
        for func_name, func_info in analysis['functions'].items():
            status = "✅ Implemented" if func_info['exists'] else "❌ Not Implemented"
            f.write(f"| {func_name} | {status} | |\n")
        
        # RLS Policies
        f.write("\n## RLS Policies\n\n")
        
        for table, policies in analysis['policies'].items():
            f.write(f"### {table}\n\n")
            f.write("| Policy Name | Command | Using Security Definer | Notes |\n")
            f.write("|-------------|---------|------------------------|-------|\n")
            
            for policy in policies:
                using_secdef = "✅ Yes" if policy['using_security_definer'] else "❌ No"
                f.write(f"| {policy['name']} | {policy['cmd']} | {using_secdef} | |\n")
            
            f.write("\n")
        
        # Indexes for RLS Performance
        f.write("## RLS Performance Indexes\n\n")
        f.write("| Table | Column | Index | Unique/Primary | Notes |\n")
        f.write("|-------|--------|-------|----------------|-------|\n")
        
        for table, indexes in analysis['indexes'].items():
            for idx in indexes:
                index_type = "Primary Key" if idx['is_primary'] else "Unique" if idx['is_unique'] else "Non-Unique"
                f.write(f"| {table} | {idx['column']} | {idx['index_name']} | {index_type} | |\n")
        
        # Verification Status
        f.write("\n## Verification Status\n\n")
        
        overall_status = "✅ Complete" if function_pct > 90 and policy_pct > 90 and index_pct > 90 else "🟡 Partial" if function_pct > 50 and policy_pct > 50 and index_pct > 50 else "❌ Incomplete"
        
        f.write(f"Overall Implementation Status: **{overall_status}**\n\n")
        
        if overall_status != "✅ Complete":
            f.write("### Action Items\n\n")
            
            if function_pct < 100:
                missing_functions = [f for f, info in analysis['functions'].items() if not info['exists']]
                f.write(f"1. Implement missing security definer functions: {', '.join(missing_functions)}\n")
            
            if policy_pct < 100:
                f.write("2. Update remaining RLS policies to use security definer functions\n")
            
            if index_pct < 100:
                missing_indexes = [t for t in analysis['policies'] if t not in analysis['indexes']]
                f.write(f"3. Add indexes for RLS-related columns on tables: {', '.join(missing_indexes)}\n")
        
        # Next Steps
        f.write("\n## Next Steps\n\n")
        f.write("1. Run the RLS policy verification tests\n")
        f.write("2. Fix any issues identified in the tests\n")
        f.write("3. Verify performance improvements with benchmarks\n")
        f.write("4. Document the optimized RLS implementation\n")
    
    return report_file

def main():
    """Main function to run the verification"""
    try:
        # Connect to database
        conn = get_db_connection()
        logger.info("Connected to database successfully")
        
        # Analyze RLS implementation
        analysis = analyze_rls_implementation(conn)
        logger.info("Analysis completed")
        
        # Generate report
        report_file = generate_verification_report(analysis)
        logger.info(f"Report generated: {report_file}")
        
        print(f"Optimized RLS verification completed. Report generated: {report_file}")
        
        # Close connection
        conn.close()
    except Exception as e:
        logger.error(f"Error verifying optimized RLS: {str(e)}")
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
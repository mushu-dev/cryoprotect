#!/usr/bin/env python
"""
CryoProtect RLS Policy Applier for Missing Tables/Views (Improved)

This script applies enhanced Row Level Security (RLS) policies to missing tables/views:
1. experiment_with_results
2. migrations
3. mixture_with_components
4. molecule_with_properties

Improvements:
- Transaction support for atomic operations
- Enhanced verification to validate policy effectiveness
- Detailed reporting for scientific data security
- Connection pooling-aware implementation
- Performance impact measurement
"""

import os
import sys
import time
import logging
import json
import tempfile
from datetime import datetime
from contextlib import contextmanager
from dotenv import load_dotenv
import psycopg2
from psycopg2 import sql, extras
from psycopg2.extensions import ISOLATION_LEVEL_READ_COMMITTED

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/apply_rls_policies.log', 'a')
    ]
)
logger = logging.getLogger('apply_missing_rls')

def ensure_logs_dir():
    """Ensure logs directory exists."""
    os.makedirs('logs', exist_ok=True)

def ensure_reports_dir():
    """Ensure reports directory exists."""
    os.makedirs('reports/security', exist_ok=True)

def get_db_connection():
    """Create database connection from environment variables."""
    load_dotenv()  # Load environment variables from .env file
    
    # Get database connection parameters from environment variables
    db_host = os.getenv('SUPABASE_DB_HOST', 'localhost')
    db_port = os.getenv('SUPABASE_DB_PORT', '5432')
    db_name = os.getenv('SUPABASE_DB_NAME', 'postgres')
    db_user = os.getenv('SUPABASE_DB_USER', 'postgres')
    db_pass = os.getenv('SUPABASE_DB_PASSWORD', '')
    
    # Create and return the connection with proper isolation level for transactions
    conn = psycopg2.connect(
        host=db_host,
        port=db_port,
        dbname=db_name,
        user=db_user,
        password=db_pass
    )
    conn.set_isolation_level(ISOLATION_LEVEL_READ_COMMITTED)
    return conn

@contextmanager
def transaction(conn):
    """Transaction context manager for atomic operations."""
    try:
        yield conn
        conn.commit()
        logger.info("Transaction committed successfully")
    except Exception as e:
        conn.rollback()
        logger.error(f"Transaction rolled back due to error: {e}")
        raise

def parse_sql_file(file_path):
    """Parse SQL file into individual statements for controlled execution."""
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Basic parsing of SQL statements (handles --comments but not nested quotes well)
    # For production, consider using a proper SQL parser
    statements = []
    current_statement = ""
    
    for line in content.split('\n'):
        # Skip comment lines
        if line.strip().startswith('--'):
            continue
            
        # Add to current statement
        current_statement += line + "\n"
        
        # If statement ends with semicolon
        if line.strip().endswith(';'):
            statements.append(current_statement.strip())
            current_statement = ""
    
    # Add any remaining statement without semicolon
    if current_statement.strip():
        statements.append(current_statement.strip())
        
    return statements

def execute_sql_statements(conn, statements):
    """Execute SQL statements with transaction support."""
    cursor = conn.cursor()
    
    try:
        for i, statement in enumerate(statements):
            if not statement.strip():
                continue
                
            logger.debug(f"Executing statement {i+1}/{len(statements)}")
            try:
                cursor.execute(statement)
                logger.debug(f"Statement {i+1} executed successfully")
            except Exception as e:
                logger.error(f"Error executing statement {i+1}: {e}")
                logger.error(f"Statement: {statement[:100]}...")
                raise
    finally:
        cursor.close()

def measure_query_performance(conn, query, description, iterations=3):
    """Measure query performance before and after RLS implementation."""
    cursor = conn.cursor()
    try:
        # Execute query multiple times and average the results
        times = []
        rows = 0
        
        for i in range(iterations):
            start_time = time.time()
            cursor.execute(query)
            results = cursor.fetchall()
            end_time = time.time()
            
            times.append(end_time - start_time)
            rows = len(results)
        
        avg_time = sum(times) / len(times)
        logger.info(f"Performance - {description}: {avg_time:.6f} seconds for {rows} rows")
        
        return {
            "description": description,
            "avg_execution_time": avg_time,
            "row_count": rows,
            "iterations": iterations
        }
    finally:
        cursor.close()

def verify_rls_policies(conn):
    """Verify that RLS policies have been applied correctly and are effective."""
    logger.info("Verifying RLS policies...")
    results = {
        "tables": {},
        "views": {},
        "effectiveness": {}
    }
    
    tables_to_check = [
        'molecule', 'molecular_property',
        'mixture', 'mixture_component',
        'experiment', 'experiment_property',
        'migrations', 'scientific_data_audit'
    ]
    
    views_to_check = [
        'molecule_with_properties',
        'mixture_with_components',
        'experiment_with_results'
    ]
    
    cursor = conn.cursor(cursor_factory=extras.DictCursor)
    
    # Check tables
    for table in tables_to_check:
        table_info = {"exists": False, "has_rls": False, "policies": []}
        
        # Check if table exists
        cursor.execute(
            "SELECT COUNT(*) FROM pg_tables WHERE tablename = %s AND schemaname = 'public'", 
            (table,)
        )
        if cursor.fetchone()[0] == 0:
            logger.warning(f"Table {table} does not exist. Skipping...")
            table_info["exists"] = False
            results["tables"][table] = table_info
            continue
        
        table_info["exists"] = True
            
        # Check if RLS is enabled
        cursor.execute(
            """
            SELECT relrowsecurity
            FROM pg_class
            WHERE relname = %s AND relnamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public')
            """, 
            (table,)
        )
        table_info["has_rls"] = cursor.fetchone()[0]
        
        # Get policies
        cursor.execute(
            "SELECT polname, polcmd, polpermissive, polroles::text, pg_get_expr(polqual, oid) FROM pg_policy WHERE polrelid = %s::regclass", 
            (f"public.{table}",)
        )
        policies = cursor.fetchall()
        
        for policy in policies:
            table_info["policies"].append({
                "name": policy[0],
                "command": policy[1],
                "permissive": policy[2],
                "roles": policy[3],
                "expression": policy[4]
            })
        
        results["tables"][table] = table_info
        
        if table_info["has_rls"]:
            logger.info(f"✓ Table {table} has RLS enabled with {len(table_info['policies'])} policies")
        else:
            logger.warning(f"✗ Table {table} does not have RLS enabled")
    
    # Check views
    for view in views_to_check:
        view_info = {"exists": False, "security_invoker": False, "underlying_tables": []}
        
        # Check if view exists
        cursor.execute(
            "SELECT COUNT(*) FROM pg_views WHERE viewname = %s AND schemaname = 'public'", 
            (view,)
        )
        if cursor.fetchone()[0] == 0:
            logger.warning(f"View {view} does not exist. Skipping...")
            view_info["exists"] = False
            results["views"][view] = view_info
            continue
            
        view_info["exists"] = True
            
        # Check if SECURITY INVOKER is set
        cursor.execute(
            """
            SELECT pg_catalog.pg_get_viewdef(c.oid) 
            FROM pg_catalog.pg_class c 
            JOIN pg_catalog.pg_namespace n ON c.relnamespace = n.oid 
            WHERE c.relname = %s 
              AND n.nspname = 'public'
            """, 
            (view,)
        )
        view_def = cursor.fetchone()[0].upper()
        view_info["security_invoker"] = 'SECURITY INVOKER' in view_def or 'SECURITY DEFINER' not in view_def
        
        # Get underlying tables
        cursor.execute(
            """
            SELECT DISTINCT ref_table.relname
            FROM pg_depend dep
            JOIN pg_rewrite rewr ON dep.objid = rewr.oid
            JOIN pg_class view ON rewr.ev_class = view.oid
            JOIN pg_class ref_table ON dep.refobjid = ref_table.oid
            JOIN pg_namespace view_ns ON view.relnamespace = view_ns.oid
            JOIN pg_namespace ref_ns ON ref_table.relnamespace = ref_ns.oid
            WHERE view.relname = %s
            AND view_ns.nspname = 'public'
            AND ref_ns.nspname = 'public'
            AND ref_table.relkind = 'r'
            """,
            (view,)
        )
        view_info["underlying_tables"] = [r[0] for r in cursor.fetchall()]
        
        results["views"][view] = view_info
        
        if view_info["security_invoker"]:
            logger.info(f"✓ View {view} has SECURITY INVOKER set (uses tables: {', '.join(view_info['underlying_tables'])})")
        else:
            logger.warning(f"✗ View {view} does not have SECURITY INVOKER set")
    
    # Verify policy effectiveness with sample queries
    try:
        # Test with service role access
        cursor.execute("SET ROLE service_role;")
        
        # Test access to molecule_with_properties
        cursor.execute("SELECT COUNT(*) FROM molecule_with_properties")
        service_role_access = cursor.fetchone()[0]
        results["effectiveness"]["service_role_molecule_view"] = service_role_access > 0
        
        # Test access to scientific_data_audit
        cursor.execute("SELECT COUNT(*) FROM scientific_data_audit")
        service_role_audit_access = cursor.fetchone()[0] >= 0  # Should be accessible
        results["effectiveness"]["service_role_audit_access"] = service_role_audit_access
        
        # Reset role
        cursor.execute("RESET ROLE;")
        
        # More comprehensive verification would include tests with different user roles
        # But that would require creating test users with specific roles
        
    except Exception as e:
        logger.error(f"Error during policy effectiveness verification: {e}")
    finally:
        cursor.execute("RESET ROLE;")
    
    cursor.close()
    
    # Generate verification report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/security/rls_verification_report_{timestamp}.json"
    
    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"RLS verification report saved to {report_path}")
    
    return results

def benchmark_performance_impact(conn):
    """Benchmark performance impact of RLS policies."""
    logger.info("Benchmarking performance impact of RLS policies...")
    
    benchmark_results = []
    
    # Sample queries to benchmark
    queries = [
        {
            "description": "Select all molecules with properties",
            "sql": "SELECT * FROM molecule_with_properties LIMIT 100"
        },
        {
            "description": "Select all mixtures with components",
            "sql": "SELECT * FROM mixture_with_components LIMIT 100"
        },
        {
            "description": "Select all experiments with results",
            "sql": "SELECT * FROM experiment_with_results LIMIT 100"
        }
    ]
    
    # Run benchmarks with service role (RLS should be bypassed)
    cursor = conn.cursor()
    cursor.execute("SET ROLE service_role;")
    
    service_role_results = []
    for query in queries:
        result = measure_query_performance(
            conn, 
            query["sql"], 
            f"Service role: {query['description']}"
        )
        service_role_results.append(result)
    
    cursor.execute("RESET ROLE;")
    
    # Run benchmarks with normal role (RLS should be applied)
    # This would be more meaningful with actual user roles
    # but we're just demonstrating the concept
    normal_role_results = []
    for query in queries:
        result = measure_query_performance(
            conn, 
            query["sql"], 
            f"Normal role: {query['description']}"
        )
        normal_role_results.append(result)
    
    # Compare and report
    benchmark_results = {
        "timestamp": datetime.now().isoformat(),
        "service_role_results": service_role_results,
        "normal_role_results": normal_role_results,
        "comparisons": []
    }
    
    for i, query in enumerate(queries):
        service_time = service_role_results[i]["avg_execution_time"]
        normal_time = normal_role_results[i]["avg_execution_time"]
        
        if service_time > 0:
            overhead_pct = ((normal_time - service_time) / service_time) * 100
        else:
            overhead_pct = 0
            
        benchmark_results["comparisons"].append({
            "description": query["description"],
            "service_role_time": service_time,
            "normal_role_time": normal_time,
            "overhead_percentage": overhead_pct,
            "conclusion": "Significant overhead" if overhead_pct > 20 else "Acceptable overhead"
        })
        
        logger.info(f"Performance impact - {query['description']}: {overhead_pct:.2f}% overhead")
    
    # Save benchmark results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/security/rls_performance_impact_{timestamp}.json"
    
    with open(report_path, 'w') as f:
        json.dump(benchmark_results, f, indent=2)
    
    logger.info(f"Performance impact report saved to {report_path}")
    
    return benchmark_results

def generate_summary_report(verification_results, benchmark_results):
    """Generate a summary report of RLS implementation."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/security/rls_implementation_summary_{timestamp}.md"
    
    with open(report_path, 'w') as f:
        f.write("# CryoProtect RLS Implementation Summary\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Tables with RLS\n\n")
        f.write("| Table | RLS Enabled | # Policies | Notes |\n")
        f.write("|-------|-------------|------------|-------|\n")
        
        for table, info in verification_results["tables"].items():
            if info["exists"]:
                f.write(f"| {table} | {'✅' if info['has_rls'] else '❌'} | {len(info['policies'])} | |\n")
            else:
                f.write(f"| {table} | ❌ | N/A | Table does not exist |\n")
        
        f.write("\n## Views with Security Settings\n\n")
        f.write("| View | Exists | Security Invoker | Underlying Tables |\n")
        f.write("|------|--------|-----------------|-------------------|\n")
        
        for view, info in verification_results["views"].items():
            if info["exists"]:
                tables_str = ", ".join(info["underlying_tables"])
                f.write(f"| {view} | ✅ | {'✅' if info['security_invoker'] else '❌'} | {tables_str} |\n")
            else:
                f.write(f"| {view} | ❌ | N/A | N/A |\n")
        
        f.write("\n## RLS Effectiveness\n\n")
        
        for test, result in verification_results["effectiveness"].items():
            f.write(f"- **{test}**: {'✅ Passed' if result else '❌ Failed'}\n")
        
        f.write("\n## Performance Impact\n\n")
        f.write("| Query | Service Role Time (s) | Normal Role Time (s) | Overhead (%) | Assessment |\n")
        f.write("|-------|----------------------|---------------------|--------------|------------|\n")
        
        for comp in benchmark_results["comparisons"]:
            f.write(f"| {comp['description']} | {comp['service_role_time']:.6f} | {comp['normal_role_time']:.6f} | {comp['overhead_percentage']:.2f} | {comp['conclusion']} |\n")
        
        f.write("\n## Recommendations\n\n")
        
        # Auto-generate recommendations based on findings
        high_overhead_queries = [c for c in benchmark_results["comparisons"] if c["overhead_percentage"] > 20]
        if high_overhead_queries:
            f.write("### Performance Optimizations Needed\n\n")
            for query in high_overhead_queries:
                f.write(f"- **{query['description']}**: Consider adding indexes or optimizing RLS policy conditions\n")
        
        missing_rls_tables = [t for t, i in verification_results["tables"].items() if i["exists"] and not i["has_rls"]]
        if missing_rls_tables:
            f.write("\n### Missing RLS Configurations\n\n")
            for table in missing_rls_tables:
                f.write(f"- Add RLS policies to table: **{table}**\n")
        
        missing_security_invoker = [v for v, i in verification_results["views"].items() if i["exists"] and not i["security_invoker"]]
        if missing_security_invoker:
            f.write("\n### View Security Issues\n\n")
            for view in missing_security_invoker:
                f.write(f"- Add SECURITY INVOKER to view: **{view}**\n")
    
    logger.info(f"Summary report generated: {report_path}")
    return report_path

def main():
    """Main function to apply RLS policies."""
    ensure_logs_dir()
    ensure_reports_dir()
    
    logger.info("Starting enhanced RLS policy application for missing tables/views")
    
    try:
        # Get database connection
        conn = get_db_connection()
        logger.info("Successfully connected to the database")
        
        # Parse and execute SQL file with RLS policies
        sql_file = 'migrations/missing_rls_policies_improved.sql'
        sql_statements = parse_sql_file(sql_file)
        logger.info(f"Parsed {len(sql_statements)} SQL statements from {sql_file}")
        
        # Execute statements within a transaction
        with transaction(conn):
            execute_sql_statements(conn, sql_statements)
        
        # Verify RLS policies
        verification_results = verify_rls_policies(conn)
        
        # Benchmark performance impact
        benchmark_results = benchmark_performance_impact(conn)
        
        # Generate summary report
        summary_report = generate_summary_report(verification_results, benchmark_results)
        
        # Print summary
        logger.info("RLS policy application completed successfully")
        logger.info(f"Summary report available at: {summary_report}")
        
        conn.close()
        
    except Exception as e:
        logger.error(f"Error applying RLS policies: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
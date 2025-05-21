#!/usr/bin/env python
"""
CryoProtect RLS Policy Verification Script

This script verifies the effectiveness of Row Level Security (RLS) policies by:
1. Simulating access as different user roles (service, admin, regular, anonymous)
2. Testing CRUD operations on key tables/views for each role
3. Measuring query performance with RLS enabled
4. Documenting results in a comprehensive verification report
"""

import os
import sys
import time
import json
import logging
import argparse
import statistics
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
        logging.FileHandler('logs/verify_rls_effectiveness.log', 'a')
    ]
)
logger = logging.getLogger('verify_rls')

def ensure_logs_dir():
    """Ensure logs directory exists."""
    os.makedirs('logs', exist_ok=True)

def ensure_reports_dir():
    """Ensure reports directory exists."""
    os.makedirs('reports/security', exist_ok=True)

def get_db_connection(project_id=None):
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

def measure_query_performance(conn, query, description, iterations=5):
    """Measure query performance with multiple iterations."""
    cursor = conn.cursor()
    try:
        # Execute query multiple times and measure performance
        times = []
        rows = 0
        
        for i in range(iterations):
            start_time = time.time()
            cursor.execute(query)
            results = cursor.fetchall()
            end_time = time.time()
            
            times.append(end_time - start_time)
            rows = len(results)
        
        # Calculate statistics
        avg_time = sum(times) / len(times)
        median_time = statistics.median(times)
        min_time = min(times)
        max_time = max(times)
        std_dev = statistics.stdev(times) if len(times) > 1 else 0
        
        logger.info(f"Performance - {description}: avg={avg_time:.6f}s, median={median_time:.6f}s, min={min_time:.6f}s, max={max_time:.6f}s for {rows} rows")
        
        return {
            "description": description,
            "avg_execution_time": avg_time,
            "median_execution_time": median_time,
            "min_execution_time": min_time,
            "max_execution_time": max_time,
            "std_deviation": std_dev,
            "row_count": rows,
            "iterations": iterations
        }
    finally:
        cursor.close()

def simulate_role_access(conn, role, operation, table, data=None):
    """Simulate access to a table as a specific role and perform the specified operation."""
    cursor = conn.cursor()
    result = {
        "role": role,
        "operation": operation,
        "table": table,
        "success": False,
        "rows_affected": 0,
        "error": None,
        "data": None
    }
    
    try:
        # Set role for this operation
        if role == "anonymous":
            # For anonymous, we need to reset the role and set auth.uid() to NULL
            cursor.execute("RESET ROLE;")
            cursor.execute("SET LOCAL auth.uid = NULL;")
        elif role == "service_role":
            cursor.execute("SET ROLE service_role;")
        elif role == "admin":
            # For admin, we need to set a specific admin user ID
            cursor.execute("RESET ROLE;")
            cursor.execute("SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';")
            # Also set a session variable to indicate admin role
            cursor.execute("SET LOCAL app.user_role = 'admin';")
        elif role == "regular":
            # For regular user, we need to set a specific regular user ID
            cursor.execute("RESET ROLE;")
            cursor.execute("SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000002';")
            cursor.execute("SET LOCAL app.user_role = 'user';")
        
        # Perform the specified operation
        if operation == "SELECT":
            cursor.execute(f"SELECT * FROM {table} LIMIT 10")
            result["data"] = cursor.fetchall()
            result["rows_affected"] = len(result["data"])
            result["success"] = True
        elif operation == "INSERT" and data:
            columns = ", ".join(data.keys())
            placeholders = ", ".join(["%s"] * len(data))
            cursor.execute(
                f"INSERT INTO {table} ({columns}) VALUES ({placeholders}) RETURNING id",
                list(data.values())
            )
            result["data"] = cursor.fetchone()
            result["rows_affected"] = cursor.rowcount
            result["success"] = True
        elif operation == "UPDATE" and data:
            set_clause = ", ".join([f"{k} = %s" for k in data.keys()])
            cursor.execute(
                f"UPDATE {table} SET {set_clause} WHERE id = %s RETURNING id",
                list(data.values()) + [data.get("id")]
            )
            result["data"] = cursor.fetchone()
            result["rows_affected"] = cursor.rowcount
            result["success"] = True
        elif operation == "DELETE" and data:
            cursor.execute(
                f"DELETE FROM {table} WHERE id = %s RETURNING id",
                [data.get("id")]
            )
            result["data"] = cursor.fetchone()
            result["rows_affected"] = cursor.rowcount
            result["success"] = True
        
    except Exception as e:
        result["error"] = str(e)
        logger.warning(f"Error during {operation} as {role} on {table}: {str(e)}")
    finally:
        # Reset role
        cursor.execute("RESET ROLE;")
        cursor.execute("RESET auth.uid;")
        cursor.close()
    
    return result

def verify_rls_effectiveness(conn, project_id=None):
    """Verify the effectiveness of RLS policies for different roles and operations."""
    logger.info("Verifying RLS policy effectiveness...")
    
    # Get all tables and views to test
    cursor = conn.cursor()
    
    # Get all tables
    cursor.execute(
        "SELECT tablename FROM pg_tables WHERE schemaname = 'public'"
    )
    tables = [row[0] for row in cursor.fetchall()]
    
    # Get all views
    cursor.execute(
        "SELECT viewname FROM pg_views WHERE schemaname = 'public'"
    )
    views = [row[0] for row in cursor.fetchall()]
    
    # Combine tables and views for testing
    tables_to_test = tables + views
    
    cursor.close()
    
    # Roles to test
    roles = ["service_role", "admin", "regular", "anonymous"]
    
    # Operations to test
    operations = ["SELECT", "INSERT", "UPDATE", "DELETE"]
    
    # Store results
    results = {
        "timestamp": datetime.now().isoformat(),
        "project_id": project_id,
        "tables": {},
        "performance": {}
    }
    
    # Test each table with each role and operation
    for table in tables_to_test:
        results["tables"][table] = {}
        
        for role in roles:
            results["tables"][table][role] = {}
            
            for operation in operations:
                # Skip operations that don't make sense for views
                if table.endswith("_with_results") or table.endswith("_with_components") or table.endswith("_with_properties"):
                    if operation != "SELECT":
                        results["tables"][table][role][operation] = {
                            "skipped": True,
                            "reason": "Operation not applicable to views"
                        }
                        continue
                
                # Prepare test data for write operations
                test_data = None
                if operation in ["INSERT", "UPDATE", "DELETE"]:
                    if table == "user_profile":
                        test_data = {
                            "id": "00000000-0000-0000-0000-000000000003",
                            "auth_user_id": "00000000-0000-0000-0000-000000000003",
                            "name": "Test User",
                            "email": "test@example.com",
                            "role": "user",
                            "project_id": project_id or "00000000-0000-0000-0000-000000000001"
                        }
                    elif table == "migrations":
                        test_data = {
                            "id": 9999,
                            "name": "test_migration",
                            "applied_at": "2023-01-01 00:00:00"
                        }
                    # Add more test data for other tables as needed
                
                # Simulate role access
                result = simulate_role_access(conn, role, operation, table, test_data)
                results["tables"][table][role][operation] = {
                    "success": result["success"],
                    "rows_affected": result["rows_affected"],
                    "error": result["error"]
                }
                
                # Log the result
                status = "✅" if result["success"] else "❌"
                logger.info(f"{status} {role} {operation} on {table}: {result['rows_affected']} rows affected")
    
    # Measure query performance for SELECT operations with different roles
    performance_queries = [
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
        },
        {
            "description": "Select all user profiles",
            "sql": "SELECT * FROM user_profile LIMIT 100"
        }
    ]
    
    for role in roles:
        results["performance"][role] = []
        
        # Set role for performance testing
        cursor = conn.cursor()
        try:
            if role == "anonymous":
                cursor.execute("RESET ROLE;")
                cursor.execute("SET LOCAL auth.uid = NULL;")
            elif role == "service_role":
                cursor.execute("SET ROLE service_role;")
            elif role == "admin":
                cursor.execute("RESET ROLE;")
                cursor.execute("SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';")
                cursor.execute("SET LOCAL app.user_role = 'admin';")
            elif role == "regular":
                cursor.execute("RESET ROLE;")
                cursor.execute("SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000002';")
                cursor.execute("SET LOCAL app.user_role = 'user';")
            
            # Run performance tests
            for query in performance_queries:
                try:
                    perf_result = measure_query_performance(
                        conn, 
                        query["sql"], 
                        f"{role}: {query['description']}"
                    )
                    results["performance"][role].append(perf_result)
                except Exception as e:
                    logger.error(f"Error measuring performance for {role} on {query['description']}: {str(e)}")
                    results["performance"][role].append({
                        "description": query["description"],
                        "error": str(e)
                    })
        finally:
            cursor.execute("RESET ROLE;")
            cursor.execute("RESET auth.uid;")
            cursor.close()
    
    # Generate verification report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/security/rls_effectiveness_report_{timestamp}.json"
    
    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"RLS effectiveness report saved to {report_path}")
    
    # Generate human-readable report
    generate_human_readable_report(results, timestamp)
    
    return results

def generate_human_readable_report(results, timestamp):
    """Generate a human-readable report from the verification results."""
    report_path = f"reports/security/rls_effectiveness_report_{timestamp}.md"
    
    with open(report_path, 'w') as f:
        f.write("# CryoProtect RLS Effectiveness Verification Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Project ID: {results['project_id'] or 'Default'}\n\n")
        
        f.write("## Access Pattern Verification\n\n")
        
        # Create a table for each operation
        for operation in ["SELECT", "INSERT", "UPDATE", "DELETE"]:
            f.write(f"### {operation} Operation\n\n")
            
            # Create table header
            f.write("| Table | Service Role | Admin | Regular User | Anonymous |\n")
            f.write("|-------|--------------|-------|--------------|----------|\n")
            
            # Add rows for each table
            for table in results["tables"]:
                row = f"| {table} |"
                
                for role in ["service_role", "admin", "regular", "anonymous"]:
                    if operation in results["tables"][table][role]:
                        if "skipped" in results["tables"][table][role][operation]:
                            row += " N/A |"
                        elif results["tables"][table][role][operation]["success"]:
                            row += f" ✅ ({results['tables'][table][role][operation]['rows_affected']} rows) |"
                        else:
                            error = results["tables"][table][role][operation]["error"]
                            short_error = error[:30] + "..." if error and len(error) > 30 else error
                            row += f" ❌ ({short_error}) |"
                    else:
                        row += " ? |"
                
                f.write(row + "\n")
            
            f.write("\n")
        
        f.write("## Performance Measurements\n\n")
        f.write("| Query | Role | Avg Time (s) | Median Time (s) | Min Time (s) | Max Time (s) | Rows |\n")
        f.write("|-------|------|--------------|-----------------|--------------|--------------|------|\n")
        
        for role in results["performance"]:
            for perf in results["performance"][role]:
                if "error" in perf:
                    f.write(f"| {perf['description']} | {role} | Error: {perf['error']} | | | | |\n")
                else:
                    f.write(f"| {perf['description']} | {role} | {perf['avg_execution_time']:.6f} | {perf['median_execution_time']:.6f} | {perf['min_execution_time']:.6f} | {perf['max_execution_time']:.6f} | {perf['row_count']} |\n")
        
        f.write("\n## Summary of Findings\n\n")
        
        # Generate summary based on results
        access_issues = []
        for table in results["tables"]:
            for role in ["service_role", "admin", "regular", "anonymous"]:
                for operation in ["SELECT", "INSERT", "UPDATE", "DELETE"]:
                    if operation in results["tables"][table][role] and not results["tables"][table][role][operation].get("skipped", False):
                        expected_success = False
                        
                        # Define expected behavior
                        if role == "service_role":
                            # Service role should have access to everything
                            expected_success = True
                        elif role == "admin":
                            # Admin should have SELECT access to everything, and full access to migrations
                            if operation == "SELECT" or table == "migrations":
                                expected_success = True
                        elif role == "regular":
                            # Regular users should have SELECT access to data views but not migrations
                            if operation == "SELECT" and table != "migrations":
                                expected_success = True
                        # Anonymous users should not have access to anything
                        
                        actual_success = results["tables"][table][role][operation]["success"]
                        if actual_success != expected_success:
                            access_issues.append({
                                "table": table,
                                "role": role,
                                "operation": operation,
                                "expected": expected_success,
                                "actual": actual_success
                            })
        
        if access_issues:
            f.write("### Access Control Issues\n\n")
            for issue in access_issues:
                expected = "have" if issue["expected"] else "not have"
                actual = "has" if issue["actual"] else "does not have"
                f.write(f"- **{issue['role']}** should {expected} {issue['operation']} access to **{issue['table']}**, but {actual} access\n")
        else:
            f.write("### Access Control\n\n")
            f.write("✅ All access controls are working as expected\n\n")
        
        # Performance analysis
        f.write("### Performance Impact\n\n")
        
        # Compare service_role performance to other roles
        if "service_role" in results["performance"] and results["performance"]["service_role"]:
            f.write("#### Performance Overhead by Role\n\n")
            f.write("| Query | Role | Avg Time (s) | Overhead vs Service Role |\n")
            f.write("|-------|------|--------------|-------------------------|\n")
            
            for i, query in enumerate(results["performance"]["service_role"]):
                if "error" not in query:
                    service_time = query["avg_execution_time"]
                    query_desc = query["description"]
                    
                    for role in ["admin", "regular", "anonymous"]:
                        if role in results["performance"] and i < len(results["performance"][role]):
                            role_query = results["performance"][role][i]
                            if "error" not in role_query:
                                role_time = role_query["avg_execution_time"]
                                if service_time > 0:
                                    overhead_pct = ((role_time - service_time) / service_time) * 100
                                    overhead_str = f"{overhead_pct:.2f}%"
                                else:
                                    overhead_str = "N/A"
                                
                                f.write(f"| {query_desc} | {role} | {role_time:.6f} | {overhead_str} |\n")
        
        f.write("\n## Recommendations\n\n")
        
        # Generate recommendations based on findings
        if access_issues:
            f.write("### Access Control Fixes\n\n")
            for issue in access_issues:
                if issue["expected"] and not issue["actual"]:
                    f.write(f"- Add policy to grant **{issue['role']}** {issue['operation']} access to **{issue['table']}**\n")
                elif not issue["expected"] and issue["actual"]:
                    f.write(f"- Modify policy to prevent **{issue['role']}** from {issue['operation']} access to **{issue['table']}**\n")
        
        # Performance recommendations
        high_overhead_found = False
        for role in ["admin", "regular", "anonymous"]:
            if role in results["performance"]:
                for query in results["performance"][role]:
                    if "error" not in query and "service_role" in results["performance"]:
                        # Find matching service role query
                        service_query = next((q for q in results["performance"]["service_role"] if q["description"] == query["description"]), None)
                        if service_query and "error" not in service_query:
                            service_time = service_query["avg_execution_time"]
                            role_time = query["avg_execution_time"]
                            if service_time > 0:
                                overhead_pct = ((role_time - service_time) / service_time) * 100
                                if overhead_pct > 50:  # Threshold for high overhead
                                    if not high_overhead_found:
                                        f.write("\n### Performance Optimizations\n\n")
                                        high_overhead_found = True
                                    f.write(f"- Optimize RLS policy for **{query['description']}** when accessed as **{role}** (current overhead: {overhead_pct:.2f}%)\n")
    
    logger.info(f"Human-readable report generated: {report_path}")
    return report_path

def main():
    """Main function to verify RLS policy effectiveness."""
    parser = argparse.ArgumentParser(description="Verify RLS policy effectiveness in CryoProtect Supabase")
    parser.add_argument("--project-id", help="Supabase project ID to verify", default="tsdlmynydfuypiugmkev")
    args = parser.parse_args()
    
    ensure_logs_dir()
    ensure_reports_dir()
    
    logger.info(f"Starting RLS policy verification for project {args.project_id}")
    
    try:
        # Get database connection
        conn = get_db_connection(args.project_id)
        logger.info("Successfully connected to the database")
        
        # Verify RLS policy effectiveness
        with transaction(conn):
            verify_rls_effectiveness(conn, args.project_id)
        
        logger.info("RLS policy verification completed successfully")
        
    except Exception as e:
        logger.error(f"Error during RLS policy verification: {e}")
        sys.exit(1)
    finally:
        if 'conn' in locals():
            conn.close()

if __name__ == "__main__":
    main()
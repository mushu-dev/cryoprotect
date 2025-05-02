#!/usr/bin/env python
"""
Test RLS Policies in CryoProtect Supabase Project

This script tests the effectiveness of RLS policies by:
1. Simulating different user roles (service role, admin, regular user, anonymous)
2. Testing access to tables and views
3. Measuring query performance with RLS enabled

Usage:
    python test_rls_policies.py --project-id <project_id>
"""

import os
import sys
import json
import time
import logging
import argparse
import statistics
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/test_rls_policies.log', 'a')
    ]
)
logger = logging.getLogger(__name__)

def ensure_directory(path):
    """Ensure a directory exists."""
    if not os.path.exists(path):
        os.makedirs(path)

def execute_sql_via_mcp(project_id, sql_query):
    """Execute SQL query via Supabase MCP server."""
    try:
        # Import the MCP client here to avoid dependency issues
        from supabase_mcp_tools import execute_sql_on_supabase
        
        # Execute the query
        result = execute_sql_on_supabase(project_id, sql_query)
        
        return result
    except ImportError:
        logger.error("supabase_mcp_tools module not found. Using fallback method.")
        return execute_sql_via_mcp_fallback(project_id, sql_query)

def execute_sql_via_mcp_fallback(project_id, sql_query):
    """Fallback method to execute SQL via MCP using the use_mcp_tool."""
    try:
        # Use the MCP tool to execute the SQL
        from supabase_mcp_client import SupabaseMCPClient
        
        # Create a client
        client = SupabaseMCPClient(project_id)
        
        # Execute the query
        result = client.execute_sql(sql_query)
        
        return result
    except Exception as e:
        logger.error(f"Error in fallback MCP execution: {str(e)}")
        
        # Try using the MCP server directly
        try:
            import subprocess
            import tempfile
            
            # Create a temporary file with the SQL query
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sql', delete=False) as temp_file:
                temp_file.write(sql_query)
                temp_file_path = temp_file.name
            
            # Use the MCP tool to execute the SQL
            cmd = [
                "python", "execute_rls_sql_via_mcp.py",
                "--project-id", project_id,
                "--sql-file", temp_file_path
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # Clean up the temporary file
            os.unlink(temp_file_path)
            
            if result.returncode != 0:
                logger.error(f"Error executing SQL via MCP: {result.stderr}")
                return None
            
            # Parse the output as JSON
            try:
                return json.loads(result.stdout)
            except json.JSONDecodeError:
                logger.error(f"Error parsing MCP result as JSON: {result.stdout}")
                return result.stdout
        except Exception as e2:
            logger.error(f"Error in subprocess MCP execution: {str(e2)}")
            return None

def test_role_access(project_id, role, operation, table):
    """Test access to a table as a specific role."""
    logger.info(f"Testing {role} {operation} access to {table}...")
    
    # Set up the SQL query based on the operation
    if operation == "SELECT":
        sql = f"""
        -- Set role
        {get_role_setup_sql(role)}
        
        -- Try to select from the table
        SELECT COUNT(*) FROM {table};
        """
    elif operation == "INSERT":
        sql = f"""
        -- Set role
        {get_role_setup_sql(role)}
        
        -- Try to insert into the table
        DO $$
        BEGIN
            BEGIN
                -- This is a test insert that will be rolled back
                INSERT INTO {table} (id, name)
                VALUES ('00000000-0000-0000-0000-000000000099', 'Test Record');
                RAISE NOTICE 'Insert succeeded';
            EXCEPTION WHEN OTHERS THEN
                RAISE NOTICE 'Insert failed: %', SQLERRM;
            END;
            -- Roll back the insert
            ROLLBACK;
        END
        $$;
        
        -- Return success/failure
        SELECT 'Insert test completed' AS result;
        """
    elif operation == "UPDATE":
        sql = f"""
        -- Set role
        {get_role_setup_sql(role)}
        
        -- Try to update the table
        DO $$
        BEGIN
            BEGIN
                -- This is a test update that will be rolled back
                UPDATE {table}
                SET name = 'Updated Name'
                WHERE id = '00000000-0000-0000-0000-000000000001';
                RAISE NOTICE 'Update succeeded';
            EXCEPTION WHEN OTHERS THEN
                RAISE NOTICE 'Update failed: %', SQLERRM;
            END;
            -- Roll back the update
            ROLLBACK;
        END
        $$;
        
        -- Return success/failure
        SELECT 'Update test completed' AS result;
        """
    elif operation == "DELETE":
        sql = f"""
        -- Set role
        {get_role_setup_sql(role)}
        
        -- Try to delete from the table
        DO $$
        BEGIN
            BEGIN
                -- This is a test delete that will be rolled back
                DELETE FROM {table}
                WHERE id = '00000000-0000-0000-0000-000000000001';
                RAISE NOTICE 'Delete succeeded';
            EXCEPTION WHEN OTHERS THEN
                RAISE NOTICE 'Delete failed: %', SQLERRM;
            END;
            -- Roll back the delete
            ROLLBACK;
        END
        $$;
        
        -- Return success/failure
        SELECT 'Delete test completed' AS result;
        """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        success = True
        
        # Check for error messages in the result
        if isinstance(result, list) and len(result) > 0:
            if "error" in str(result).lower() or "permission denied" in str(result).lower():
                success = False
        
        return {
            "role": role,
            "operation": operation,
            "table": table,
            "success": success,
            "result": result
        }
    except Exception as e:
        logger.error(f"Error testing {role} {operation} access to {table}: {str(e)}")
        return {
            "role": role,
            "operation": operation,
            "table": table,
            "success": False,
            "error": str(e)
        }

def get_role_setup_sql(role):
    """Get SQL to set up a specific role."""
    if role == "service_role":
        return "SET ROLE service_role;"
    elif role == "admin":
        return """
        RESET ROLE;
        SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
        SET LOCAL app.user_role = 'admin';
        """
    elif role == "regular":
        return """
        RESET ROLE;
        SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000002';
        SET LOCAL app.user_role = 'user';
        """
    elif role == "anonymous":
        return """
        RESET ROLE;
        SET LOCAL auth.uid = NULL;
        """
    else:
        return "RESET ROLE;"

def measure_query_performance(project_id, role, query, description, iterations=5):
    """Measure query performance with a specific role."""
    logger.info(f"Measuring performance of {description} as {role}...")
    
    sql = f"""
    -- Set role
    {get_role_setup_sql(role)}
    
    -- Execute the query
    EXPLAIN ANALYZE {query};
    """
    
    try:
        times = []
        
        for i in range(iterations):
            start_time = time.time()
            result = execute_sql_via_mcp(project_id, sql)
            end_time = time.time()
            
            times.append(end_time - start_time)
        
        # Calculate statistics
        avg_time = sum(times) / len(times)
        median_time = statistics.median(times)
        min_time = min(times)
        max_time = max(times)
        std_dev = statistics.stdev(times) if len(times) > 1 else 0
        
        logger.info(f"Performance - {description} as {role}: avg={avg_time:.6f}s, median={median_time:.6f}s")
        
        return {
            "role": role,
            "description": description,
            "query": query,
            "avg_execution_time": avg_time,
            "median_execution_time": median_time,
            "min_execution_time": min_time,
            "max_execution_time": max_time,
            "std_deviation": std_dev,
            "iterations": iterations
        }
    except Exception as e:
        logger.error(f"Error measuring performance of {description} as {role}: {str(e)}")
        return {
            "role": role,
            "description": description,
            "query": query,
            "error": str(e)
        }

def test_rls_policies(project_id):
    """Test RLS policies for different roles and operations."""
    logger.info("Testing RLS policies...")
    
    # Tables and views to test
    tables_to_test = [
        "public.molecules",
        "public.mixtures",
        "public.experiments",
        "public.migrations",
        "public.scientific_data_audit",
        "public.molecule_with_properties",
        "public.mixture_with_components",
        "public.experiment_with_results"
    ]
    
    # Roles to test
    roles = ["service_role", "admin", "regular", "anonymous"]
    
    # Operations to test
    operations = ["SELECT", "INSERT", "UPDATE", "DELETE"]
    
    # Store results
    results = {
        "timestamp": datetime.now().isoformat(),
        "project_id": project_id,
        "access_tests": [],
        "performance_tests": []
    }
    
    # Test access for each role, operation, and table
    for role in roles:
        for operation in operations:
            for table in tables_to_test:
                # Skip operations that don't make sense for views
                if table.endswith("_with_properties") or table.endswith("_with_components") or table.endswith("_with_results"):
                    if operation != "SELECT":
                        continue
                
                result = test_role_access(project_id, role, operation, table)
                results["access_tests"].append(result)
    
    # Test query performance
    performance_queries = [
        {
            "description": "Select all molecules with properties",
            "query": "SELECT * FROM public.molecule_with_properties LIMIT 100"
        },
        {
            "description": "Select all mixtures with components",
            "query": "SELECT * FROM public.mixture_with_components LIMIT 100"
        },
        {
            "description": "Select all experiments with results",
            "query": "SELECT * FROM public.experiment_with_results LIMIT 100"
        },
        {
            "description": "Select all user profiles",
            "query": "SELECT * FROM public.user_profile LIMIT 100"
        }
    ]
    
    for role in roles:
        for query_info in performance_queries:
            result = measure_query_performance(
                project_id,
                role,
                query_info["query"],
                query_info["description"]
            )
            results["performance_tests"].append(result)
    
    # Generate test report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/security/rls_test_report_{timestamp}.md"
    
    ensure_directory("reports/security")
    
    with open(report_path, 'w') as f:
        f.write("# RLS Policy Test Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Access Tests
        f.write("## Access Tests\n\n")
        
        # Group by operation
        for operation in operations:
            f.write(f"### {operation} Operation\n\n")
            
            # Create table header
            f.write("| Table | Service Role | Admin | Regular User | Anonymous |\n")
            f.write("|-------|--------------|-------|--------------|----------|\n")
            
            # Add rows for each table
            for table in tables_to_test:
                # Skip operations that don't make sense for views
                if table.endswith("_with_properties") or table.endswith("_with_components") or table.endswith("_with_results"):
                    if operation != "SELECT":
                        continue
                
                row = f"| {table} |"
                
                for role in roles:
                    # Find the test result for this combination
                    test_result = next(
                        (t for t in results["access_tests"] if t["role"] == role and t["operation"] == operation and t["table"] == table),
                        None
                    )
                    
                    if test_result:
                        if test_result["success"]:
                            row += " ✅ |"
                        else:
                            row += " ❌ |"
                    else:
                        row += " N/A |"
                
                f.write(row + "\n")
            
            f.write("\n")
        
        # Performance Tests
        f.write("## Performance Tests\n\n")
        f.write("| Query | Role | Avg Time (s) | Median Time (s) | Min Time (s) | Max Time (s) |\n")
        f.write("|-------|------|--------------|-----------------|--------------|--------------|n")
        
        for perf_test in results["performance_tests"]:
            if "error" in perf_test:
                f.write(f"| {perf_test['description']} | {perf_test['role']} | Error: {perf_test['error']} | | | |\n")
            else:
                f.write(f"| {perf_test['description']} | {perf_test['role']} | {perf_test['avg_execution_time']:.6f} | {perf_test['median_execution_time']:.6f} | {perf_test['min_execution_time']:.6f} | {perf_test['max_execution_time']:.6f} |\n")
        
        f.write("\n")
        
        # Performance Comparison
        f.write("## Performance Comparison\n\n")
        f.write("| Query | Role | Avg Time (s) | Overhead vs Service Role |\n")
        f.write("|-------|------|--------------|-------------------------|\n")
        
        # Group by query description
        query_descriptions = set(perf_test["description"] for perf_test in results["performance_tests"] if "description" in perf_test)
        
        for description in query_descriptions:
            # Get service role time for this query
            service_role_test = next(
                (t for t in results["performance_tests"] if "description" in t and t["description"] == description and t["role"] == "service_role"),
                None
            )
            
            if service_role_test and "avg_execution_time" in service_role_test:
                service_role_time = service_role_test["avg_execution_time"]
                
                # Compare other roles to service role
                for role in [r for r in roles if r != "service_role"]:
                    role_test = next(
                        (t for t in results["performance_tests"] if "description" in t and t["description"] == description and t["role"] == role),
                        None
                    )
                    
                    if role_test and "avg_execution_time" in role_test:
                        role_time = role_test["avg_execution_time"]
                        
                        if service_role_time > 0:
                            overhead_pct = ((role_time - service_role_time) / service_role_time) * 100
                            overhead_str = f"{overhead_pct:.2f}%"
                        else:
                            overhead_str = "N/A"
                        
                        f.write(f"| {description} | {role} | {role_time:.6f} | {overhead_str} |\n")
        
        f.write("\n")
        
        # Save raw results as JSON
        json_path = f"reports/security/rls_test_results_{timestamp}.json"
        with open(json_path, 'w') as json_file:
            json.dump(results, json_file, indent=2)
        
        f.write(f"Raw test results saved to: {json_path}\n")
    
    logger.info(f"Test report generated: {report_path}")
    return report_path

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Test RLS policies in CryoProtect Supabase project")
    parser.add_argument("--project-id", help="Supabase project ID", default="tsdlmynydfuypiugmkev")
    args = parser.parse_args()
    
    # Ensure directories exist
    ensure_directory("logs")
    ensure_directory("reports/security")
    
    logger.info(f"Starting RLS policy tests for project {args.project_id}")
    
    # Test RLS policies
    report_path = test_rls_policies(args.project_id)
    
    logger.info(f"RLS policy tests completed. Report saved to {report_path}")

if __name__ == "__main__":
    main()
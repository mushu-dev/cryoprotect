#!/usr/bin/env python3
"""
Demonstrate RLS Optimization Performance

This script demonstrates the performance improvements of the RLS optimizations
by running test queries before and after using the optimized functions.
"""

import os
import sys
import time
import json
import logging
import argparse
from datetime import datetime
from prettytable import PrettyTable

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f"rls_optimization_demo_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
    ]
)
logger = logging.getLogger(__name__)

# Test queries
ORIGINAL_QUERIES = {
    "property_range_query": """
        SELECT m.id, m.name, m.molecular_formula
        FROM molecules m
        JOIN molecular_properties mp ON m.id = mp.molecule_id
        WHERE mp.property_name = 'Molecular Weight'
        AND mp.property_value::numeric BETWEEN 100 AND 500
        AND (m.is_public = true OR m.created_by = auth.uid() OR 
             EXISTS (
                 SELECT 1 FROM project_molecules pm
                 JOIN team_projects tp ON pm.project_id = tp.project_id
                 JOIN user_profile up ON tp.team_id = up.team_id
                 WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
             ))
        LIMIT 10;
    """,
    
    "molecules_with_properties_query": """
        SELECT 
            m.id,
            m.name,
            m.smiles,
            m.molecular_formula,
            COUNT(mp.id) AS property_count
        FROM 
            molecules m
        LEFT JOIN
            molecular_properties mp ON m.id = mp.molecule_id
        WHERE 
            m.is_public = true 
            OR m.created_by = auth.uid()
            OR EXISTS (
                SELECT 1 FROM project_molecules pm
                JOIN team_projects tp ON pm.project_id = tp.project_id
                JOIN user_profile up ON tp.team_id = up.team_id
                WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
            )
        GROUP BY
            m.id, m.name, m.smiles, m.molecular_formula
        ORDER BY
            m.name
        LIMIT 10;
    """,
    
    "mixtures_with_components_query": """
        SELECT 
            m.id,
            m.name,
            m.description,
            COUNT(mc.id) AS component_count
        FROM 
            mixtures m
        LEFT JOIN
            mixture_components mc ON m.id = mc.mixture_id
        WHERE 
            m.is_public = true 
            OR m.created_by = auth.uid()
            OR EXISTS (
                SELECT 1 FROM project_mixtures pm
                JOIN team_projects tp ON pm.project_id = tp.project_id
                JOIN user_profile up ON tp.team_id = up.team_id
                WHERE pm.mixture_id = m.id AND up.auth_user_id = auth.uid()
            )
        GROUP BY
            m.id, m.name, m.description
        ORDER BY
            m.name
        LIMIT 10;
    """
}

# Optimized queries
OPTIMIZED_QUERIES = {
    "property_range_query": """
        SELECT m.id, m.name, m.molecular_formula
        FROM find_molecules_by_property_range('Molecular Weight', 100, 500) mol_ids
        JOIN molecules m ON m.id = mol_ids
        LIMIT 10;
    """,
    
    "molecules_with_properties_query": """
        SELECT * FROM get_molecules_with_properties(10, 0);
    """,
    
    "mixtures_with_components_query": """
        SELECT * FROM get_mixtures_with_components(10, 0);
    """
}

# Check if we have PrettyTable for nice output formatting
try:
    from prettytable import PrettyTable
    HAVE_PRETTYTABLE = True
except ImportError:
    HAVE_PRETTYTABLE = False
    logger.warning("PrettyTable not installed. Tables will be formatted plainly.")
    logger.warning("Install with: pip install prettytable")


def execute_query_mcp(project_id, query):
    """Execute a query using Supabase MCP."""
    try:
        from supabase_mcp_tools import execute_sql_on_supabase
        
        start_time = time.time()
        result = execute_sql_on_supabase(project_id, query)
        execution_time = time.time() - start_time
        
        return result, execution_time * 1000
    except Exception as e:
        logger.error(f"Error executing query with MCP: {e}")
        return None, 0


def execute_query_supabase(supabase, query):
    """Execute a query using Supabase client."""
    try:
        start_time = time.time()
        result = supabase.rpc('exec_sql', {'query': query}).execute()
        execution_time = time.time() - start_time
        
        if hasattr(result, 'data'):
            return result.data, execution_time * 1000
        else:
            return result, execution_time * 1000
    except Exception as e:
        logger.error(f"Error executing query with Supabase: {e}")
        return None, 0


def execute_query_direct(conn, query):
    """Execute a query using direct database connection."""
    try:
        cursor = conn.cursor()
        
        start_time = time.time()
        cursor.execute(query)
        result = cursor.fetchall()
        execution_time = time.time() - start_time
        
        cursor.close()
        return result, execution_time * 1000
    except Exception as e:
        logger.error(f"Error executing query with direct connection: {e}")
        return None, 0


def run_performance_test(connection, project_id=None, iterations=3):
    """Run performance tests on original and optimized queries."""
    results = {
        "original": {},
        "optimized": {},
        "summary": {}
    }
    
    # Helper function to execute query based on connection type
    def run_query(query):
        if project_id:  # MCP
            return execute_query_mcp(project_id, query)
        elif hasattr(connection, 'rpc'):  # Supabase
            return execute_query_supabase(connection, query)
        else:  # Direct
            return execute_query_direct(connection, query)
    
    # Test original queries
    logger.info("Testing original queries...")
    for name, query in ORIGINAL_QUERIES.items():
        logger.info(f"Testing {name}...")
        times = []
        row_count = 0
        
        for i in range(iterations):
            result, exec_time = run_query(query)
            if result is not None:
                times.append(exec_time)
                row_count = len(result) if isinstance(result, list) else 0
        
        if times:
            avg_time = sum(times) / len(times)
            results["original"][name] = {
                "avg_execution_time": avg_time,
                "row_count": row_count,
                "execution_times": times
            }
            logger.info(f"  Average execution time: {avg_time:.2f}ms")
            logger.info(f"  Rows returned: {row_count}")
        else:
            logger.warning(f"  No successful executions for {name}")
    
    # Test optimized queries
    logger.info("\nTesting optimized queries...")
    for name, query in OPTIMIZED_QUERIES.items():
        logger.info(f"Testing optimized {name}...")
        times = []
        row_count = 0
        
        for i in range(iterations):
            result, exec_time = run_query(query)
            if result is not None:
                times.append(exec_time)
                row_count = len(result) if isinstance(result, list) else 0
        
        if times:
            avg_time = sum(times) / len(times)
            results["optimized"][name] = {
                "avg_execution_time": avg_time,
                "row_count": row_count,
                "execution_times": times
            }
            logger.info(f"  Average execution time: {avg_time:.2f}ms")
            logger.info(f"  Rows returned: {row_count}")
        else:
            logger.warning(f"  No successful executions for optimized {name}")
    
    # Calculate improvement statistics
    logger.info("\nCalculating performance improvements...")
    for name in ORIGINAL_QUERIES.keys():
        if name in results["original"] and name in results["optimized"]:
            orig_time = results["original"][name]["avg_execution_time"]
            opt_time = results["optimized"][name]["avg_execution_time"]
            
            if orig_time > 0:
                improvement_pct = ((orig_time - opt_time) / orig_time) * 100
                speedup_factor = orig_time / opt_time if opt_time > 0 else float('inf')
                
                results["summary"][name] = {
                    "original_ms": orig_time,
                    "optimized_ms": opt_time,
                    "improvement_pct": improvement_pct,
                    "speedup_factor": speedup_factor
                }
                
                logger.info(f"Query {name} - Improvement: {improvement_pct:.2f}% "
                           f"({orig_time:.2f}ms → {opt_time:.2f}ms) - Speedup: {speedup_factor:.2f}x")
    
    return results


def connect_to_database(args):
    """Connect to the database using the specified method."""
    if args.mcp:
        logger.info(f"Using MCP with project ID: {args.project_id}")
        return None, args.project_id
    
    elif args.supabase:
        try:
            try:
                logger.info("Using service_role_helper to get Supabase client")
                from service_role_helper import get_supabase_client
                client = get_supabase_client()
            except ImportError:
                logger.info("Importing Supabase client directly")
                from supabase import create_client
                
                # Try to get URL and key from environment variables
                import os
                url = os.environ.get("SUPABASE_URL")
                key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY")
                
                if not url or not key:
                    logger.error("SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY environment variables are required")
                    return None, None
                
                client = create_client(url, key)
                
            logger.info("Connected to Supabase successfully")
            return client, None
        except Exception as e:
            logger.error(f"Error connecting to Supabase: {e}")
            return None, None
    
    else:  # Direct
        try:
            import psycopg2
            logger.info(f"Connecting to database {args.db_name} on {args.db_host}")
            conn = psycopg2.connect(
                host=args.db_host,
                port=args.db_port,
                dbname=args.db_name,
                user=args.db_user,
                password=args.db_password
            )
            logger.info("Connected to database successfully")
            return conn, None
        except Exception as e:
            logger.error(f"Error connecting to database: {e}")
            return None, None


def display_results(results):
    """Display performance test results in a nice table."""
    if HAVE_PRETTYTABLE:
        table = PrettyTable()
        table.field_names = ["Query", "Original (ms)", "Optimized (ms)", "Improvement", "Speedup"]
        
        for name, stats in results["summary"].items():
            table.add_row([
                name, 
                f"{stats['original_ms']:.2f}", 
                f"{stats['optimized_ms']:.2f}", 
                f"{stats['improvement_pct']:.2f}%", 
                f"{stats['speedup_factor']:.2f}x"
            ])
        
        print("\nPerformance Improvement Summary:")
        print(table)
    else:
        print("\nPerformance Improvement Summary:")
        print("-" * 80)
        print("| {:<30} | {:>12} | {:>12} | {:>12} | {:>12} |".format(
            "Query", "Original (ms)", "Optimized (ms)", "Improvement", "Speedup"))
        print("-" * 80)
        
        for name, stats in results["summary"].items():
            print("| {:<30} | {:>12.2f} | {:>12.2f} | {:>11.2f}% | {:>11.2f}x |".format(
                name, 
                stats["original_ms"], 
                stats["optimized_ms"], 
                stats["improvement_pct"], 
                stats["speedup_factor"]
            ))
        
        print("-" * 80)


def save_results(results, output_file):
    """Save performance test results to a file."""
    try:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Test results saved to {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error saving test results: {e}")
        return False


def main():
    """Main function to run the performance demonstration."""
    parser = argparse.ArgumentParser(description='Demonstrate RLS Query Optimization Performance')
    
    # Connection options
    connection_group = parser.add_mutually_exclusive_group(required=True)
    connection_group.add_argument('--direct', action='store_true', help='Use direct database connection')
    connection_group.add_argument('--supabase', action='store_true', help='Use Supabase client')
    connection_group.add_argument('--mcp', action='store_true', help='Use Supabase MCP')
    
    # MCP project ID
    parser.add_argument('--project-id', default="tsdlmynydfuypiugmkev", help='Supabase project ID (required for MCP)')
    
    # Direct database connection parameters
    parser.add_argument('--db-host', help='Database host (for direct connection)')
    parser.add_argument('--db-port', type=int, default=5432, help='Database port (for direct connection)')
    parser.add_argument('--db-name', help='Database name (for direct connection)')
    parser.add_argument('--db-user', help='Database user (for direct connection)')
    parser.add_argument('--db-password', help='Database password (for direct connection)')
    
    # Test parameters
    parser.add_argument('--iterations', type=int, default=3, help='Number of iterations per query (default: 3)')
    parser.add_argument('--output', help='Output JSON file path for test results')
    
    args = parser.parse_args()
    
    # Check if MCP with project ID
    if args.mcp and not args.project_id:
        logger.error("Project ID is required for MCP mode")
        return 1
    
    # Check if direct with connection params
    if args.direct and not all([args.db_host, args.db_name, args.db_user, args.db_password]):
        logger.error("Database connection parameters are required for direct connection")
        return 1
    
    # Connect to the database
    connection, project_id = connect_to_database(args)
    
    if not connection and not project_id:
        logger.error("Failed to establish database connection")
        return 1
    
    # Run performance tests
    print("\n" + "=" * 80)
    print("RUNNING RLS PERFORMANCE OPTIMIZATION TESTS")
    print("=" * 80)
    
    results = run_performance_test(connection, project_id, args.iterations)
    
    # Display performance summary
    display_results(results)
    
    # Save results if output file is specified
    if args.output:
        save_results(results, args.output)
    
    # Close the connection if needed
    if connection and not args.mcp and not args.supabase:
        connection.close()
    
    print("\n" + "=" * 80)
    print("PERFORMANCE TEST COMPLETED")
    print("=" * 80 + "\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
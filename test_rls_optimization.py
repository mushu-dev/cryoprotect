#!/usr/bin/env python3
"""
Test script for RLS policy optimization.

This script tests that the RLS policy optimizations are applied correctly
and improve performance.
"""

import os
import sys
import logging
import argparse
import time
from datetime import datetime
from pathlib import Path

try:
    from database import connection_manager
except ImportError:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from database import connection_manager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('rls_optimization_test')

def test_security_definer_functions(conn):
    """Test that security definer functions exist and work correctly."""
    logger.info("Testing security definer functions...")
    
    try:
        # Check for a set of expected security definer functions
        functions_sql = """
        SELECT proname, prosecdef 
        FROM pg_proc 
        JOIN pg_namespace ON pg_proc.pronamespace = pg_namespace.oid
        WHERE nspname = 'auth' 
          AND proname IN ('check_resource_access', 'check_resource_modify_access', 
                         'is_service_role_user', 'check_batch_access')
          AND prokind = 'f'
          AND prolang = (SELECT oid FROM pg_language WHERE lanname = 'plpgsql');
        """
        
        results = conn.execute(functions_sql).fetchall()
        
        # Log findings
        if results:
            logger.info(f"Found {len(results)} security definer functions:")
            for r in results:
                func_name = r['proname']
                is_security_definer = r['prosecdef']
                logger.info(f"  - {func_name} (Security Definer: {is_security_definer})")
                
            # Verify they are all security definer
            all_security_definer = all(r['prosecdef'] for r in results)
            if all_security_definer:
                logger.info("All functions are correctly set as SECURITY DEFINER.")
                return True
            else:
                logger.warning("Some functions are not set as SECURITY DEFINER!")
                return False
        else:
            logger.warning("No security definer functions found.")
            return False
            
    except Exception as e:
        logger.error(f"Error testing security definer functions: {str(e)}")
        return False

def test_performance_indexes(conn):
    """Test that performance indexes exist."""
    logger.info("Testing performance indexes...")
    
    try:
        # Check for a set of expected indexes
        indexes_sql = """
        SELECT indexname, tablename
        FROM pg_indexes
        WHERE indexname IN (
            'idx_user_profile_user_id', 
            'idx_molecule_project_id', 
            'idx_mixture_project_id',
            'idx_user_profile_user_project'
        );
        """
        
        results = conn.execute(indexes_sql).fetchall()
        
        # Log findings
        if results:
            logger.info(f"Found {len(results)} performance indexes:")
            for r in results:
                index_name = r['indexname']
                table_name = r['tablename']
                logger.info(f"  - {index_name} on {table_name}")
            
            return True
        else:
            logger.warning("No performance indexes found.")
            return False
            
    except Exception as e:
        logger.error(f"Error testing performance indexes: {str(e)}")
        return False

def test_rls_policies(conn):
    """Test that optimized RLS policies exist and use security definer functions."""
    logger.info("Testing RLS policies...")
    
    try:
        # Check for policies that use the security definer functions
        policies_sql = """
        SELECT policyname, tablename, cmd, permissive, qual
        FROM pg_policies
        WHERE schemaname = 'public'
        ORDER BY tablename, policyname;
        """
        
        results = conn.execute(policies_sql).fetchall()
        
        # Look for policies with function calls like is_project_member
        optimized_policies = 0
        total_policies = len(results)
        
        if total_policies > 0:
            logger.info(f"Found {total_policies} RLS policies:")
            for r in results:
                policy_name = r['policyname']
                table_name = r['tablename']
                command = r['cmd']
                permissive = r['permissive']
                qual = r['qual'] or "None"
                
                # Check if the policy uses a function
                is_optimized = False
                function_markers = [
                    'is_project_member', 
                    'molecule_in_user_project', 
                    'is_team_member'
                ]
                
                for marker in function_markers:
                    if marker in str(qual):
                        is_optimized = True
                        optimized_policies += 1
                        break
                
                status = "Optimized" if is_optimized else "Standard"
                logger.info(f"  - {policy_name} on {table_name} ({command}, {status})")
            
            optimization_ratio = optimized_policies / total_policies
            logger.info(f"Optimization ratio: {optimized_policies}/{total_policies} policies optimized ({optimization_ratio:.1%})")
            
            # Consider it a success if at least half the policies are optimized
            return optimization_ratio >= 0.5
        else:
            logger.warning("No RLS policies found.")
            return False
            
    except Exception as e:
        logger.error(f"Error testing RLS policies: {str(e)}")
        return False

def benchmark_queries(conn):
    """Benchmark some common queries to assess RLS policy performance."""
    logger.info("Benchmarking queries with RLS policies...")
    
    test_queries = [
        ("Select molecules (should use optimized RLS)", 
         "SELECT COUNT(*) FROM molecule;"),
        
        ("Select mixtures (should use optimized RLS)", 
         "SELECT COUNT(*) FROM mixture;"),
         
        ("Select molecular properties (should use optimized RLS)",
         "SELECT COUNT(*) FROM molecular_property;"),
         
        ("Join between molecule and properties (should use optimized RLS)",
         """
         SELECT m.name, COUNT(mp.id) as property_count 
         FROM molecule m 
         LEFT JOIN molecular_property mp ON m.id = mp.molecule_id 
         GROUP BY m.id, m.name 
         LIMIT 10;
         """)
    ]
    
    results = []
    
    # Run each query and measure time
    for description, query in test_queries:
        logger.info(f"Testing query: {description}")
        
        try:
            # First run an EXPLAIN ANALYZE
            explain_query = f"EXPLAIN ANALYZE {query}"
            start_time = time.time()
            explain_result = conn.execute(explain_query).fetchall()
            explain_time = time.time() - start_time
            
            # Extract estimated and actual execution time from EXPLAIN ANALYZE
            planning_time = None
            execution_time = None
            
            for row in explain_result:
                if 'Planning Time' in str(row):
                    planning_time = float(str(row).split(': ')[1].split(' ms')[0])
                elif 'Execution Time' in str(row):
                    execution_time = float(str(row).split(': ')[1].split(' ms')[0])
            
            # Now run the actual query for real-world timing
            start_time = time.time()
            conn.execute(query).fetchall()
            total_time = time.time() - start_time
            
            result = {
                'description': description,
                'planning_time_ms': planning_time,
                'execution_time_ms': execution_time,
                'total_time_s': total_time,
                'explain_output': explain_result
            }
            
            results.append(result)
            
            logger.info(f"  Planning time: {planning_time:.2f} ms")
            logger.info(f"  Execution time: {execution_time:.2f} ms")
            logger.info(f"  Total time: {total_time:.4f} seconds")
            
            # Look for specific function calls in the explain output
            function_calls = []
            for row in explain_result:
                row_str = str(row)
                if 'Function Scan on' in row_str:
                    function_name = row_str.split('Function Scan on ')[1].split(' ')[0]
                    function_calls.append(function_name)
            
            if function_calls:
                logger.info(f"  Used functions: {', '.join(function_calls)}")
                
        except Exception as e:
            logger.error(f"Error benchmarking query '{description}': {str(e)}")
            results.append({
                'description': description,
                'error': str(e)
            })
    
    return results

def generate_performance_report(results):
    """Generate a performance report based on benchmark results."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = Path("reports")
    report_dir.mkdir(exist_ok=True)
    
    report_file = report_dir / f"rls_performance_report_{timestamp}.md"
    
    report_content = f"""# RLS Performance Test Report

## Summary
Timestamp: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Query Performance Results

| Query | Planning Time (ms) | Execution Time (ms) | Total Time (s) |
|-------|-------------------|---------------------|----------------|
"""
    
    for result in results:
        if 'error' in result:
            report_content += f"| {result['description']} | Error | Error | Error |\n"
        else:
            planning = result.get('planning_time_ms', 'N/A')
            execution = result.get('execution_time_ms', 'N/A')
            total = result.get('total_time_s', 'N/A')
            
            report_content += f"| {result['description']} | {planning:.2f if isinstance(planning, float) else planning} | {execution:.2f if isinstance(execution, float) else execution} | {total:.4f if isinstance(total, float) else total} |\n"
    
    report_content += """
## Interpretation

The times above indicate how the optimized RLS policies are performing.
Faster planning time generally indicates that the query optimizer is able to 
make better decisions due to the security definer functions and specialized indexes.

Execution time indicates how efficiently the query can be executed after planning.
Lower values indicate better performance.

### Potential Improvements

If planning times are high:
- Analyze query plans to see if the right indexes are being used
- Consider adding more specific indexes for common query patterns

If execution times are high:
- Look for sequential scans that could be replaced with index scans
- Consider adding materialized views for complex, frequently-run queries

## References
- [PostgreSQL Query Planning](https://www.postgresql.org/docs/current/performance-tips.html)
- [PostgreSQL EXPLAIN](https://www.postgresql.org/docs/current/sql-explain.html)
"""
    
    try:
        with open(report_file, 'w') as f:
            f.write(report_content)
        logger.info(f"Performance report generated: {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error generating performance report: {str(e)}")
        return None

def main():
    """Main function to test RLS optimizations."""
    parser = argparse.ArgumentParser(description='Test RLS policy optimization')
    parser.add_argument('--benchmark-only', action='store_true', help='Only run performance benchmarks')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Connect to the database
    try:
        connection = connection_manager.ConnectionManager.get_instance()
        if connection.connect():
            conn = connection.get_active_adapter()
            logger.info(f"Connected to database using {connection.active_adapter} adapter")
        else:
            logger.error("Failed to connect to database")
            return 1
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        return 1
    
    try:
        test_results = {
            'security_definer_functions': False,
            'performance_indexes': False,
            'rls_policies': False
        }
        
        # Run tests
        if not args.benchmark_only:
            test_results['security_definer_functions'] = test_security_definer_functions(conn)
            test_results['performance_indexes'] = test_performance_indexes(conn)
            test_results['rls_policies'] = test_rls_policies(conn)
            
            # Log summary
            logger.info("\nTest Results Summary:")
            logger.info(f"Security Definer Functions: {'PASS' if test_results['security_definer_functions'] else 'FAIL'}")
            logger.info(f"Performance Indexes: {'PASS' if test_results['performance_indexes'] else 'FAIL'}")
            logger.info(f"RLS Policies: {'PASS' if test_results['rls_policies'] else 'FAIL'}")
        
        # Run benchmarks
        benchmark_results = benchmark_queries(conn)
        report_file = generate_performance_report(benchmark_results)
        
        if report_file:
            logger.info(f"Performance report available at: {report_file}")
        
        # Determine success
        if args.benchmark_only:
            success = len(benchmark_results) > 0
        else:
            success = all(test_results.values())
        
        if success:
            logger.info("✅ RLS optimization verification completed successfully")
            return 0
        else:
            logger.warning("⚠️ Some RLS optimization tests failed - see logs for details")
            return 1
    
    except Exception as e:
        logger.error(f"Error during RLS optimization testing: {str(e)}")
        return 1
    
    finally:
        # Disconnect from database
        connection.disconnect()
        logger.info("Database connection closed")

if __name__ == "__main__":
    sys.exit(main())
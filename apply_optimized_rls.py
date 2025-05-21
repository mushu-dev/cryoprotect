#!/usr/bin/env python3
"""
Apply Optimized RLS Policies Migration

This script applies the optimized RLS policies migration to improve
performance of Row Level Security policies in the database.
"""

import os
import sys
import logging
import argparse
import time
from pathlib import Path
from typing import Optional, List, Dict, Any

# Try to import the optimized connection pool
try:
    from db_pool import ConnectionPool, execute_query, execute_transaction
    USE_POOL = True
except ImportError:
    try:
        from optimized_connection_pool import (
            execute_query_with_retry as execute_query,
            transaction_context,
            ConnectionManager
        )
        USE_POOL = False
    except ImportError:
        print("Error: Cannot import database connection modules. Please ensure either db_pool.py or optimized_connection_pool.py is available.")
        sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/apply_optimized_rls.log', 'w')
    ]
)
logger = logging.getLogger(__name__)

def ensure_directory(path):
    """Ensure a directory exists."""
    if not os.path.exists(path):
        os.makedirs(path)

def read_migration_file(file_path: str) -> str:
    """
    Read the migration SQL file.
    
    Args:
        file_path: Path to the migration file
        
    Returns:
        SQL content as string
    """
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading migration file {file_path}: {str(e)}")
        raise

def execute_migration(sql: str) -> bool:
    """
    Execute the migration SQL.
    
    Args:
        sql: SQL migration to execute
        
    Returns:
        True if successful, False otherwise
    """
    try:
        if USE_POOL:
            result = execute_transaction(sql)
            logger.info("Migration executed successfully using connection pool")
            return True
        else:
            with transaction_context() as conn:
                with conn.cursor() as cursor:
                    cursor.execute(sql)
                
            logger.info("Migration executed successfully")
            return True
    except Exception as e:
        logger.error(f"Error executing migration: {str(e)}")
        return False

def verify_security_definer_functions() -> Dict[str, Any]:
    """
    Verify the security definer functions are present.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        'success': False,
        'functions_found': 0,
        'functions_expected': 8,
        'function_list': []
    }
    
    # List of functions to check for
    functions_to_check = [
        'auth.check_resource_access',
        'auth.check_resource_modify_access',
        'auth.is_service_role_user',
        'auth.check_verification_access',
        'auth.check_verification_modify_access',
        'auth.check_batch_access',
        'auth.get_user_permissions',
        'auth.user_has_permission'
    ]
    
    try:
        # Query to check for the functions
        sql = """
        SELECT 
            n.nspname as schema,
            p.proname as name,
            pg_get_functiondef(p.oid) as definition
        FROM 
            pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
        WHERE 
            n.nspname = 'auth' AND 
            p.proname IN (
                'check_resource_access',
                'check_resource_modify_access',
                'is_service_role_user',
                'check_verification_access',
                'check_verification_modify_access',
                'check_batch_access',
                'get_user_permissions',
                'user_has_permission'
            );
        """
        
        result = execute_query(sql)
        
        if result:
            # Process results
            for row in result:
                function_name = f"{row['schema']}.{row['name']}"
                results['function_list'].append(function_name)
            
            results['functions_found'] = len(results['function_list'])
            results['success'] = results['functions_found'] == results['functions_expected']
            
            logger.info(f"Found {results['functions_found']} of {results['functions_expected']} expected security definer functions")
            
            # Log specific results
            for function in functions_to_check:
                if function in results['function_list']:
                    logger.info(f"✓ Function {function} found")
                else:
                    logger.warning(f"✗ Function {function} not found!")
        else:
            logger.error("Failed to query database for security definer functions")
    
    except Exception as e:
        logger.error(f"Error verifying security definer functions: {str(e)}")
    
    return results

def verify_practical_functions() -> Dict[str, Any]:
    """
    Verify the practical RLS functions are present.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        'success': False,
        'functions_found': 0,
        'functions_expected': 6,
        'function_list': []
    }
    
    # List of functions to check for
    functions_to_check = [
        'public.is_project_member',
        'public.user_projects',
        'public.is_team_member',
        'public.is_project_owner',
        'public.molecule_in_user_project',
        'public.mixture_in_user_project'
    ]
    
    try:
        # Query to check for the functions
        sql = """
        SELECT 
            n.nspname as schema,
            p.proname as name,
            pg_get_functiondef(p.oid) as definition
        FROM 
            pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
        WHERE 
            n.nspname = 'public' AND 
            p.proname IN (
                'is_project_member',
                'user_projects',
                'is_team_member',
                'is_project_owner',
                'molecule_in_user_project',
                'mixture_in_user_project'
            );
        """
        
        result = execute_query(sql)
        
        if result:
            # Process results
            for row in result:
                function_name = f"{row['schema']}.{row['name']}"
                results['function_list'].append(function_name)
            
            results['functions_found'] = len(results['function_list'])
            results['success'] = results['functions_found'] == results['functions_expected']
            
            logger.info(f"Found {results['functions_found']} of {results['functions_expected']} expected practical RLS functions")
            
            # Log specific results
            for function in functions_to_check:
                if function in results['function_list']:
                    logger.info(f"✓ Function {function} found")
                else:
                    logger.warning(f"✗ Function {function} not found!")
        else:
            logger.error("Failed to query database for practical RLS functions")
    
    except Exception as e:
        logger.error(f"Error verifying practical RLS functions: {str(e)}")
    
    return results

def verify_performance_indexes() -> Dict[str, Any]:
    """
    Verify the performance indexes are present.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        'success': False,
        'indexes_found': 0,
        'indexes_expected': 12,
        'index_list': []
    }
    
    # List of indexes to check for
    indexes_to_check = [
        'idx_user_profile_user_id',
        'idx_user_profile_project_id',
        'idx_user_profile_team_id',
        'idx_user_profile_role',
        'idx_molecule_project_id',
        'idx_mixture_project_id',
        'idx_mixture_component_mixture_id',
        'idx_molecular_property_molecule_id',
        'idx_experiment_project_id',
        'idx_prediction_molecule_id',
        'idx_experiment_property_experiment_id',
        'idx_user_profile_user_project'
    ]
    
    try:
        # Query to check for the indexes
        sql = """
        SELECT 
            tablename,
            indexname
        FROM 
            pg_indexes
        WHERE 
            schemaname = 'public' AND
            indexname IN (
                'idx_user_profile_user_id',
                'idx_user_profile_project_id',
                'idx_user_profile_team_id',
                'idx_user_profile_role',
                'idx_molecule_project_id',
                'idx_mixture_project_id',
                'idx_mixture_component_mixture_id',
                'idx_molecular_property_molecule_id',
                'idx_experiment_project_id',
                'idx_prediction_molecule_id',
                'idx_experiment_property_experiment_id',
                'idx_user_profile_user_project'
            );
        """
        
        result = execute_query(sql)
        
        if result:
            # Process results
            for row in result:
                results['index_list'].append(row['indexname'])
            
            results['indexes_found'] = len(results['index_list'])
            results['success'] = results['indexes_found'] > 0
            
            logger.info(f"Found {results['indexes_found']} of {results['indexes_expected']} expected performance indexes")
            
            # Log specific results
            for index in indexes_to_check:
                if index in results['index_list']:
                    logger.info(f"✓ Index {index} found")
                else:
                    logger.info(f"✗ Index {index} not found!")
        else:
            logger.error("Failed to query database for performance indexes")
    
    except Exception as e:
        logger.error(f"Error verifying performance indexes: {str(e)}")
    
    return results

def verify_service_role_policies() -> Dict[str, Any]:
    """
    Verify the service role policies are present.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        'success': False,
        'tables_with_policy': 0,
        'tables_checked': 0,
        'tables_with_policy_list': []
    }
    
    try:
        # Query to get all tables
        sql_tables = """
        SELECT 
            tablename
        FROM 
            pg_tables
        WHERE 
            schemaname = 'public';
        """
        
        # Get list of tables
        tables_result = execute_query(sql_tables)
        
        if not tables_result:
            logger.error("Failed to query database for tables")
            return results
        
        # Process tables
        tables = [row['tablename'] for row in tables_result]
        results['tables_checked'] = len(tables)
        
        # Check each table for service role policy
        for table in tables:
            # Query to check for service role policy
            sql_policy = f"""
            SELECT 
                policyname
            FROM 
                pg_policies
            WHERE 
                schemaname = 'public' AND
                tablename = '{table}' AND
                policyname = 'service_role_all_access';
            """
            
            policy_result = execute_query(sql_policy)
            
            if policy_result and len(policy_result) > 0:
                results['tables_with_policy_list'].append(table)
        
        results['tables_with_policy'] = len(results['tables_with_policy_list'])
        results['success'] = results['tables_with_policy'] > 0
        
        # Calculate percentage of tables with service role policy
        if results['tables_checked'] > 0:
            coverage_pct = (results['tables_with_policy'] / results['tables_checked']) * 100
            logger.info(f"Service role policy coverage: {coverage_pct:.1f}% ({results['tables_with_policy']} of {results['tables_checked']} tables)")
        
    except Exception as e:
        logger.error(f"Error verifying service role policies: {str(e)}")
    
    return results

def run_benchmark() -> Dict[str, Any]:
    """
    Run a simple benchmark to measure the performance improvement.
    
    Returns:
        Dictionary with benchmark results
    """
    results = {
        'success': False,
        'queries': [],
        'summary': {}
    }
    
    # Define benchmark queries
    benchmark_queries = [
        {
            'name': 'Select all molecules for a user',
            'sql': """
            SET ROLE authenticated;
            SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
            SELECT COUNT(*) FROM molecule;
            """
        },
        {
            'name': 'Select all mixture components for a user',
            'sql': """
            SET ROLE authenticated;
            SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
            SELECT COUNT(*) FROM mixture_component;
            """
        },
        {
            'name': 'Select all molecular properties for a user',
            'sql': """
            SET ROLE authenticated;
            SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
            SELECT COUNT(*) FROM molecular_property;
            """
        }
    ]
    
    try:
        # Run each benchmark query
        for query in benchmark_queries:
            query_result = {
                'name': query['name'],
                'execution_time': None,
                'error': None
            }
            
            try:
                # Measure execution time
                start_time = time.time()
                
                # Execute query
                execute_query(query['sql'])
                
                # Calculate execution time
                execution_time = time.time() - start_time
                query_result['execution_time'] = execution_time
                
                logger.info(f"Benchmark: {query['name']} - {execution_time:.6f} seconds")
            except Exception as e:
                query_result['error'] = str(e)
                logger.error(f"Error in benchmark query {query['name']}: {str(e)}")
            
            results['queries'].append(query_result)
        
        # Calculate summary statistics
        successful_queries = [q for q in results['queries'] if q['execution_time'] is not None]
        
        if successful_queries:
            avg_time = sum(q['execution_time'] for q in successful_queries) / len(successful_queries)
            results['summary']['average_execution_time'] = avg_time
            results['summary']['successful_queries'] = len(successful_queries)
            results['summary']['total_queries'] = len(benchmark_queries)
            results['success'] = True
            
            logger.info(f"Benchmark summary: Average execution time {avg_time:.6f} seconds")
        else:
            logger.error("No benchmark queries completed successfully")
    
    except Exception as e:
        logger.error(f"Error in benchmark: {str(e)}")
    
    return results

def generate_report(verify_security_functions_results: Dict[str, Any], 
                    verify_practical_functions_results: Dict[str, Any],
                    verify_indexes_results: Dict[str, Any],
                    verify_service_role_results: Dict[str, Any],
                    benchmark_results: Dict[str, Any]) -> str:
    """
    Generate a verification report.
    
    Args:
        verify_security_functions_results: Results of security functions verification
        verify_practical_functions_results: Results of practical functions verification
        verify_indexes_results: Results of index verification
        verify_service_role_results: Results of service role policy verification
        benchmark_results: Results of benchmark
        
    Returns:
        Path to the report file
    """
    # Ensure directory exists
    ensure_directory("reports")
    
    # Generate timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/rls_optimization_report_{timestamp}.md"
    
    with open(report_path, 'w') as f:
        f.write("# RLS Optimization Verification Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Security Definer Functions
        f.write("## Security Definer Functions\n\n")
        f.write(f"Status: {'✅ All functions found' if verify_security_functions_results['success'] else '❌ Some functions missing'}\n\n")
        f.write(f"Found {verify_security_functions_results['functions_found']} of {verify_security_functions_results['functions_expected']} expected functions\n\n")
        
        if verify_security_functions_results['function_list']:
            f.write("| Function | Status |\n")
            f.write("|---------|--------|\n")
            
            # List of expected functions
            functions_to_check = [
                'auth.check_resource_access',
                'auth.check_resource_modify_access',
                'auth.is_service_role_user',
                'auth.check_verification_access',
                'auth.check_verification_modify_access',
                'auth.check_batch_access',
                'auth.get_user_permissions',
                'auth.user_has_permission'
            ]
            
            for function in functions_to_check:
                status = "✅ Found" if function in verify_security_functions_results['function_list'] else "❌ Missing"
                f.write(f"| {function} | {status} |\n")
        
        f.write("\n")
        
        # Practical Functions
        f.write("## Practical RLS Functions\n\n")
        f.write(f"Status: {'✅ All functions found' if verify_practical_functions_results['success'] else '❌ Some functions missing'}\n\n")
        f.write(f"Found {verify_practical_functions_results['functions_found']} of {verify_practical_functions_results['functions_expected']} expected functions\n\n")
        
        if verify_practical_functions_results['function_list']:
            f.write("| Function | Status |\n")
            f.write("|---------|--------|\n")
            
            # List of expected functions
            functions_to_check = [
                'public.is_project_member',
                'public.user_projects',
                'public.is_team_member',
                'public.is_project_owner',
                'public.molecule_in_user_project',
                'public.mixture_in_user_project'
            ]
            
            for function in functions_to_check:
                status = "✅ Found" if function in verify_practical_functions_results['function_list'] else "❌ Missing"
                f.write(f"| {function} | {status} |\n")
        
        f.write("\n")
        
        # Performance Indexes
        f.write("## Performance Indexes\n\n")
        f.write(f"Status: {'✅ Some indexes found' if verify_indexes_results['success'] else '❌ No indexes found'}\n\n")
        f.write(f"Found {verify_indexes_results['indexes_found']} of {verify_indexes_results['indexes_expected']} expected indexes\n\n")
        
        if verify_indexes_results['index_list']:
            f.write("| Index | Status |\n")
            f.write("|-------|--------|\n")
            
            # List of expected indexes
            indexes_to_check = [
                'idx_user_profile_user_id',
                'idx_user_profile_project_id',
                'idx_user_profile_team_id',
                'idx_user_profile_role',
                'idx_molecule_project_id',
                'idx_mixture_project_id',
                'idx_mixture_component_mixture_id',
                'idx_molecular_property_molecule_id',
                'idx_experiment_project_id',
                'idx_prediction_molecule_id',
                'idx_experiment_property_experiment_id',
                'idx_user_profile_user_project'
            ]
            
            for index in indexes_to_check:
                status = "✅ Found" if index in verify_indexes_results['index_list'] else "❌ Missing"
                f.write(f"| {index} | {status} |\n")
        
        f.write("\n")
        
        # Service Role Policies
        f.write("## Service Role Policies\n\n")
        f.write(f"Status: {'✅ Some policies found' if verify_service_role_results['success'] else '❌ No policies found'}\n\n")
        f.write(f"Found service role policies on {verify_service_role_results['tables_with_policy']} of {verify_service_role_results['tables_checked']} tables\n\n")
        
        if verify_service_role_results['tables_with_policy_list']:
            f.write("### Tables with service role policy:\n\n")
            for table in verify_service_role_results['tables_with_policy_list']:
                f.write(f"- {table}\n")
        
        f.write("\n")
        
        # Benchmark Results
        f.write("## Performance Benchmark\n\n")
        f.write(f"Status: {'✅ Completed' if benchmark_results['success'] else '❌ Failed'}\n\n")
        
        if benchmark_results['success']:
            f.write(f"Average execution time: {benchmark_results['summary']['average_execution_time']:.6f} seconds\n")
            f.write(f"Successful queries: {benchmark_results['summary']['successful_queries']} of {benchmark_results['summary']['total_queries']}\n\n")
            
            if benchmark_results['queries']:
                f.write("| Query | Execution Time (s) |\n")
                f.write("|-------|-------------------|\n")
                
                for query in benchmark_results['queries']:
                    if query['execution_time'] is not None:
                        f.write(f"| {query['name']} | {query['execution_time']:.6f} |\n")
                    else:
                        f.write(f"| {query['name']} | Error: {query['error']} |\n")
        
        f.write("\n")
        
        # Overall Status
        overall_success = (
            verify_security_functions_results['success'] or
            verify_practical_functions_results['success'] or
            verify_indexes_results['success'] or
            verify_service_role_results['success']
        )
        
        f.write("## Overall Status\n\n")
        f.write(f"Status: {'✅ Migration applied successfully' if overall_success else '❌ Migration not fully applied'}\n\n")
        
        # Summary Table
        f.write("| Component | Status |\n")
        f.write("|-----------|--------|\n")
        f.write(f"| Security Definer Functions | {'✅ Success' if verify_security_functions_results['success'] else '❌ Failed'} |\n")
        f.write(f"| Practical RLS Functions | {'✅ Success' if verify_practical_functions_results['success'] else '❌ Failed'} |\n")
        f.write(f"| Performance Indexes | {'✅ Success' if verify_indexes_results['success'] else '❌ Failed'} |\n")
        f.write(f"| Service Role Policies | {'✅ Success' if verify_service_role_results['success'] else '❌ Failed'} |\n")
        f.write(f"| Performance Benchmark | {'✅ Success' if benchmark_results['success'] else '❌ Failed'} |\n")
        
    logger.info(f"Report generated: {report_path}")
    return report_path

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Apply Optimized RLS Policies Migration")
    parser.add_argument("--security-definer-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/019_rls_security_definer_functions.sql",
                    help="Path to the security definer functions migration file")
    parser.add_argument("--optimized-policies-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/020_optimized_rls_policies.sql",
                    help="Path to the optimized RLS policies migration file")
    parser.add_argument("--verify-only", action="store_true", 
                    help="Only verify the migration (don't execute it)")
    parser.add_argument("--security-definer-only", action="store_true",
                    help="Only apply/verify the security definer functions")
    parser.add_argument("--optimized-policies-only", action="store_true",
                    help="Only apply/verify the optimized RLS policies")
    args = parser.parse_args()
    
    # Ensure the logs directory exists
    ensure_directory("logs")
    
    # Validate files exist
    if not args.verify_only:
        if not args.optimized_policies_only and not os.path.exists(args.security_definer_file):
            logger.error(f"Security definer functions migration file not found: {args.security_definer_file}")
            return 1
        
        if not args.security_definer_only and not os.path.exists(args.optimized_policies_file):
            logger.error(f"Optimized RLS policies migration file not found: {args.optimized_policies_file}")
            return 1
    
    # Verify only mode
    if args.verify_only:
        logger.info("Verify-only mode: Checking if migration has been applied")
        
        # Verify security definer functions
        security_functions_results = verify_security_definer_functions()
        
        # Verify practical functions
        practical_functions_results = verify_practical_functions()
        
        # Verify performance indexes
        indexes_results = verify_performance_indexes()
        
        # Verify service role policies
        service_role_results = verify_service_role_policies()
        
        # Run benchmark
        benchmark_results = run_benchmark()
        
        # Generate report
        report_path = generate_report(
            security_functions_results,
            practical_functions_results,
            indexes_results,
            service_role_results,
            benchmark_results
        )
        
        # Determine overall success
        overall_success = (
            security_functions_results['success'] or
            practical_functions_results['success'] or
            indexes_results['success'] or
            service_role_results['success']
        )
        
        logger.info(f"Verification report generated: {report_path}")
        
        return 0 if overall_success else 1
    
    # Execute migration mode
    try:
        # Security definer functions migration
        if not args.optimized_policies_only:
            logger.info(f"Applying security definer functions migration from {args.security_definer_file}")
            
            # Read migration file
            security_definer_sql = read_migration_file(args.security_definer_file)
            
            # Execute migration
            logger.info("Executing security definer functions migration...")
            security_definer_success = execute_migration(security_definer_sql)
            
            if not security_definer_success:
                logger.error("Security definer functions migration failed")
                return 1
            
            logger.info("Security definer functions migration applied successfully")
        
        # Optimized RLS policies migration
        if not args.security_definer_only:
            logger.info(f"Applying optimized RLS policies migration from {args.optimized_policies_file}")
            
            # Read migration file
            optimized_policies_sql = read_migration_file(args.optimized_policies_file)
            
            # Execute migration
            logger.info("Executing optimized RLS policies migration...")
            optimized_policies_success = execute_migration(optimized_policies_sql)
            
            if not optimized_policies_success:
                logger.error("Optimized RLS policies migration failed")
                return 1
            
            logger.info("Optimized RLS policies migration applied successfully")
        
        # Verify the migration
        logger.info("Verifying migration...")
        
        # Verify security definer functions
        security_functions_results = verify_security_definer_functions()
        
        # Verify practical functions
        practical_functions_results = verify_practical_functions()
        
        # Verify performance indexes
        indexes_results = verify_performance_indexes()
        
        # Verify service role policies
        service_role_results = verify_service_role_policies()
        
        # Run benchmark
        benchmark_results = run_benchmark()
        
        # Generate report
        report_path = generate_report(
            security_functions_results,
            practical_functions_results,
            indexes_results,
            service_role_results,
            benchmark_results
        )
        
        logger.info(f"Migration verification report generated: {report_path}")
        
        # Determine overall success
        if args.security_definer_only:
            overall_success = security_functions_results['success']
        elif args.optimized_policies_only:
            overall_success = (
                practical_functions_results['success'] and
                indexes_results['success'] and
                service_role_results['success']
            )
        else:
            overall_success = (
                security_functions_results['success'] and
                practical_functions_results['success'] and
                indexes_results['success'] and
                service_role_results['success']
            )
        
        if overall_success:
            logger.info("Migration completed and verified successfully")
            return 0
        else:
            logger.warning("Migration applied but verification found some issues")
            return 1
        
    except Exception as e:
        logger.error(f"Error in migration process: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
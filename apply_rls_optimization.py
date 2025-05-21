#!/usr/bin/env python3
"""
Apply RLS optimization to the database including:
- Security definer functions for common access patterns
- Performance indexes for RLS policies
- Materialized views for frequently accessed data
- Optimized RLS policies

This script improves database query performance while maintaining security
by applying optimized RLS (Row-Level Security) policies.
"""
import os
import sys
import argparse
import logging
import time
from pathlib import Path

# Add the project root to the Python module search path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import database modules
from database.adapter import get_db_adapter
from database.connection import ConnectionManager
import config

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define paths to SQL files
MIGRATIONS_DIR = Path(__file__).parent / "migrations" / "rls_helpers"

# SQL file paths
SECURITY_DEFINER_FUNCTIONS_FILE = MIGRATIONS_DIR / "security_definer_functions.sql"
PERFORMANCE_INDEXES_FILE = MIGRATIONS_DIR / "performance_indexes.sql"
MATERIALIZED_VIEWS_FILE = MIGRATIONS_DIR / "materialized_views.sql"
RLS_POLICIES_FILE = MIGRATIONS_DIR / "rls_policies.sql"

def read_sql_file(file_path):
    """Read SQL file content."""
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading SQL file {file_path}: {e}")
        raise

def apply_sql(conn, sql):
    """Apply SQL to database connection."""
    cursor = conn.cursor()
    try:
        cursor.execute(sql)
        conn.commit()
        return True
    except Exception as e:
        conn.rollback()
        logger.error(f"Error applying SQL: {e}")
        logger.error(f"SQL statement: {sql[:150]}...")  # Log just the beginning of the SQL
        return False
    finally:
        cursor.close()

def apply_optimization_step(conn, step_name, sql_file):
    """Apply a specific optimization step."""
    logger.info(f"Applying {step_name}...")
    start_time = time.time()
    
    try:
        sql = read_sql_file(sql_file)
        success = apply_sql(conn, sql)
        
        if success:
            elapsed = time.time() - start_time
            logger.info(f"Successfully applied {step_name} in {elapsed:.2f} seconds")
            return True
        else:
            logger.error(f"Failed to apply {step_name}")
            return False
    except Exception as e:
        logger.error(f"Error during {step_name}: {e}")
        return False

def verify_optimization(conn):
    """Verify that optimizations have been applied."""
    logger.info("Verifying optimizations...")
    
    # Verify security definer functions
    cursor = conn.cursor()
    try:
        # Check for a sample of functions
        functions_to_check = [
            "is_team_member", 
            "has_molecule_access", 
            "user_has_clearance"
        ]
        
        missing_functions = []
        for func in functions_to_check:
            cursor.execute(
                "SELECT COUNT(*) FROM pg_proc WHERE proname = %s",
                (func,)
            )
            count = cursor.fetchone()[0]
            if count == 0:
                missing_functions.append(func)
        
        if missing_functions:
            logger.warning(f"Missing security definer functions: {', '.join(missing_functions)}")
            return False
        
        # Check for a sample of indexes
        indexes_to_check = [
            "idx_user_profile_auth_user_id",
            "idx_molecules_created_by",
            "idx_molecules_public"
        ]
        
        missing_indexes = []
        for idx in indexes_to_check:
            cursor.execute(
                "SELECT COUNT(*) FROM pg_indexes WHERE indexname = %s",
                (idx,)
            )
            count = cursor.fetchone()[0]
            if count == 0:
                missing_indexes.append(idx)
        
        if missing_indexes:
            logger.warning(f"Missing performance indexes: {', '.join(missing_indexes)}")
            return False
        
        # Check for materialized views
        views_to_check = [
            "public_molecules_summary",
            "public_molecular_properties"
        ]
        
        missing_views = []
        for view in views_to_check:
            cursor.execute(
                "SELECT COUNT(*) FROM pg_matviews WHERE matviewname = %s",
                (view,)
            )
            count = cursor.fetchone()[0]
            if count == 0:
                missing_views.append(view)
        
        if missing_views:
            logger.warning(f"Missing materialized views: {', '.join(missing_views)}")
            return False
        
        logger.info("All optimizations were verified successfully")
        return True
    except Exception as e:
        logger.error(f"Error during verification: {e}")
        return False
    finally:
        cursor.close()

def report_query_performance(conn):
    """Run test queries and report performance improvements."""
    logger.info("Running performance tests...")
    
    test_queries = [
        ("Get molecules for authenticated user", 
         "SELECT * FROM molecules WHERE has_molecule_access(id) LIMIT 100"),
        
        ("Get public molecules using materialized view", 
         "SELECT * FROM public_molecules_summary LIMIT 100"),
        
        ("Get properties for a molecule", 
         "SELECT * FROM molecular_properties WHERE has_molecule_access(molecule_id) LIMIT 100")
    ]
    
    cursor = conn.cursor()
    for description, query in test_queries:
        try:
            start_time = time.time()
            cursor.execute(f"EXPLAIN ANALYZE {query}")
            result = cursor.fetchall()
            elapsed = time.time() - start_time
            
            # Extract query execution time from EXPLAIN ANALYZE output
            execution_time = None
            for row in result:
                if isinstance(row[0], str) and "execution time" in row[0].lower():
                    time_str = row[0].split(":")[1].strip()
                    if "ms" in time_str:
                        execution_time = float(time_str.split(" ")[0])
                        break
            
            if execution_time:
                logger.info(f"Query: {description}")
                logger.info(f"Execution time: {execution_time} ms")
            else:
                logger.info(f"Query: {description}")
                logger.info(f"Total time (including Python overhead): {elapsed:.2f} seconds")
        except Exception as e:
            logger.error(f"Error testing query '{description}': {e}")
    
    cursor.close()

def create_refresh_function_job(conn):
    """Creates or updates the scheduled job to refresh materialized views."""
    logger.info("Setting up materialized view refresh job...")
    
    try:
        cursor = conn.cursor()
        
        # Check if pg_cron extension is available
        cursor.execute("SELECT COUNT(*) FROM pg_extension WHERE extname = 'pg_cron'")
        has_cron = cursor.fetchone()[0] > 0
        
        if has_cron:
            # Check if the job already exists
            cursor.execute("SELECT COUNT(*) FROM cron.job WHERE jobname = 'refresh_materialized_views'")
            job_exists = cursor.fetchone()[0] > 0
            
            if job_exists:
                # Update existing job
                cursor.execute(
                    "UPDATE cron.job SET schedule = '0 * * * *', command = 'SELECT refresh_materialized_views()' "
                    "WHERE jobname = 'refresh_materialized_views'"
                )
                logger.info("Updated existing refresh job schedule")
            else:
                # Create new job
                cursor.execute(
                    "SELECT cron.schedule('refresh_materialized_views', '0 * * * *', 'SELECT refresh_materialized_views()')"
                )
                logger.info("Created new refresh job schedule (hourly)")
            
            conn.commit()
        else:
            logger.warning("pg_cron extension not available. Materialized views will need to be refreshed manually.")
            logger.info("You can refresh views by running: SELECT refresh_materialized_views()")
        
        cursor.close()
        return True
    except Exception as e:
        conn.rollback()
        logger.error(f"Error setting up refresh job: {e}")
        return False

def generate_report(steps_status, performance_results=None):
    """Generate a detailed report of the RLS optimization process."""
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = Path("reports")
    report_path.mkdir(exist_ok=True)
    report_file = report_path / f"rls_optimization_report_{timestamp}.md"

    report_content = f"""# RLS Policy Optimization Report

## Summary
Timestamp: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

| Optimization Step | Status |
|-------------------|---------|
| Security Definer Functions | {"✅ Applied" if steps_status.get('functions', False) else "❌ Failed"} |
| Performance Indexes | {"✅ Applied" if steps_status.get('indexes', False) else "❌ Failed"} |
| Materialized Views | {"✅ Applied" if steps_status.get('views', False) else "❌ Failed"} |
| Optimized RLS Policies | {"✅ Applied" if steps_status.get('policies', False) else "❌ Failed"} |
| Verification | {"✅ Passed" if steps_status.get('verify', False) else "❌ Failed"} |

## Overview
This report documents the application of optimized Row-Level Security (RLS) policies
to the CryoProtect database. The optimizations aim to improve query performance while
maintaining the security guarantees provided by RLS.

## Optimization Details

### Security Definer Functions
Security definer functions were implemented to optimize common access patterns.
These functions run with elevated privileges but enforce security checks internally,
reducing the overhead of repetitive RLS policy evaluations.

### Performance Indexes
Specialized indexes were added to support efficient execution of RLS policies.
These indexes target the columns commonly used in RLS policy conditions.

### Materialized Views
Materialized views were created for frequently accessed data combinations.
These views pre-compute complex joins and aggregations to improve query performance.

### Optimized RLS Policies
The RLS policies were updated to leverage the security definer functions and
performance indexes, reducing query planning and execution time.
"""

    # Add performance test results if available
    if performance_results:
        report_content += """\n## Performance Test Results\n
| Query | Execution Time |
|-------|----------------|
"""
        for query, time in performance_results:
            report_content += f"| {query} | {time} |\n"

    report_content += """\n## Next Steps

1. Monitor query performance with the optimized policies
2. Collect metrics on query execution times
3. Fine-tune indexes based on actual usage patterns
4. Consider additional optimizations for specific access patterns

## References
- [PostgreSQL RLS Documentation](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
- [PostgreSQL Security Definer Functions](https://www.postgresql.org/docs/current/sql-createfunction.html)
- [PostgreSQL Index Types](https://www.postgresql.org/docs/current/indexes-types.html)
- [PostgreSQL Materialized Views](https://www.postgresql.org/docs/current/rules-materializedviews.html)
"""

    try:
        with open(report_file, 'w') as f:
            f.write(report_content)
        logger.info(f"Generated optimization report: {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        return None

def main():
    """Main function to apply RLS optimizations."""
    parser = argparse.ArgumentParser(description='Apply RLS optimization to the database')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without making changes')
    parser.add_argument('--skip-functions', action='store_true', help='Skip creating security definer functions')
    parser.add_argument('--skip-indexes', action='store_true', help='Skip creating performance indexes')
    parser.add_argument('--skip-views', action='store_true', help='Skip creating materialized views')
    parser.add_argument('--skip-policies', action='store_true', help='Skip creating RLS policies')
    parser.add_argument('--verify', action='store_true', help='Only verify that optimizations have been applied')
    parser.add_argument('--performance-test', action='store_true', help='Run performance tests after applying optimizations')
    parser.add_argument('--no-security-definer', action='store_true', help='Skip applying security definer functions (alias for --skip-functions)')
    parser.add_argument('--no-indexes', action='store_true', help='Skip applying performance indexes (alias for --skip-indexes)')
    parser.add_argument('--no-rls-policies', action='store_true', help='Skip applying optimized RLS policies (alias for --skip-policies)')
    args = parser.parse_args()

    # Handle alias arguments
    if args.no_security_definer:
        args.skip_functions = True
    if args.no_indexes:
        args.skip_indexes = True
    if args.no_rls_policies:
        args.skip_policies = True
    
    # Get database adapter and connection
    try:
        db_adapter = get_db_adapter()
        
        if args.dry_run:
            logger.info("DRY RUN MODE - No changes will be made")
            logger.info(f"Would apply security definer functions from: {SECURITY_DEFINER_FUNCTIONS_FILE}")
            logger.info(f"Would apply performance indexes from: {PERFORMANCE_INDEXES_FILE}")
            logger.info(f"Would apply materialized views from: {MATERIALIZED_VIEWS_FILE}")
            logger.info(f"Would apply RLS policies from: {RLS_POLICIES_FILE}")
            return 0
        
        if args.verify:
            with ConnectionManager() as conn:
                success = verify_optimization(conn)
                return 0 if success else 1
        
        # Apply each optimization step
        with ConnectionManager() as conn:
            all_steps_success = True
            
            # Apply security definer functions
            if not args.skip_functions:
                success = apply_optimization_step(
                    conn, 
                    "security definer functions", 
                    SECURITY_DEFINER_FUNCTIONS_FILE
                )
                all_steps_success = all_steps_success and success
            
            # Apply performance indexes
            if not args.skip_indexes:
                success = apply_optimization_step(
                    conn, 
                    "performance indexes", 
                    PERFORMANCE_INDEXES_FILE
                )
                all_steps_success = all_steps_success and success
            
            # Apply materialized views
            if not args.skip_views:
                success = apply_optimization_step(
                    conn, 
                    "materialized views", 
                    MATERIALIZED_VIEWS_FILE
                )
                all_steps_success = all_steps_success and success
                
                # Set up refresh job for materialized views
                if success:
                    create_refresh_function_job(conn)
            
            # Apply RLS policies
            if not args.skip_policies:
                success = apply_optimization_step(
                    conn, 
                    "optimized RLS policies", 
                    RLS_POLICIES_FILE
                )
                all_steps_success = all_steps_success and success
            
            # Verify optimizations
            verify_success = verify_optimization(conn)
            
            # Run performance tests if requested
            if args.performance_test:
                report_query_performance(conn)
            
            # Collect steps status
            steps_status = {
                'functions': not args.skip_functions and all_steps_success,
                'indexes': not args.skip_indexes and all_steps_success,
                'views': not args.skip_views and all_steps_success,
                'policies': not args.skip_policies and all_steps_success,
                'verify': verify_success
            }

            # Generate performance results if test was run
            performance_results = None
            if args.performance_test:
                # Just a placeholder - actual results come from report_query_performance
                performance_results = []

            # Generate report
            report_file = generate_report(steps_status, performance_results)
            if report_file:
                logger.info(f"Full report available at: {report_file}")

            if all_steps_success:
                logger.info("Successfully applied all RLS optimizations")
                return 0
            else:
                logger.error("Some RLS optimization steps failed")
                return 1
            
    except Exception as e:
        logger.error(f"Error applying RLS optimizations: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
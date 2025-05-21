#!/usr/bin/env python3
"""
Script to consolidate duplicative RLS policies, with a focus on service_role policies.
This script applies the SQL migration and verifies the results.
"""

import os
import sys
import json
import logging
from datetime import datetime
import argparse
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to consolidate service_role policies."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '029_consolidate_service_role_policies.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to consolidate service_role policies")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def get_policy_stats_before_after():
    """Get statistics on RLS policies before and after the migration."""
    try:
        # Use the new auth.rls_policy_overview view
        query = """
        SELECT * FROM auth.rls_policy_overview;
        """
        
        policies = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Group policies by table
        policies_by_table = {}
        policy_types = {}
        
        for policy in policies:
            table = policy['tablename']
            policy_type = policy['policy_type']
            
            # Count policy types
            if policy_type not in policy_types:
                policy_types[policy_type] = 0
            policy_types[policy_type] += 1
            
            # Group by table
            if table not in policies_by_table:
                policies_by_table[table] = []
            policies_by_table[table].append(policy)
        
        # Get RLS coverage analysis
        coverage_query = """
        SELECT * FROM auth.analyze_rls_coverage();
        """
        
        coverage = db_utils.execute_query(coverage_query, cursor_factory=RealDictCursor)
        
        # Calculate stats
        stats = {
            'total_policies': len(policies),
            'tables_with_policies': len(policies_by_table),
            'policy_types': policy_types,
            'coverage': {
                'tables_analyzed': len(coverage),
                'tables_with_rls': sum(1 for c in coverage if c['has_rls']),
                'tables_with_service_role_policy': sum(1 for c in coverage if c['has_service_role_policy']),
                'tables_with_admin_policy': sum(1 for c in coverage if c['has_admin_policy']),
                'tables_with_owner_policy': sum(1 for c in coverage if c['has_owner_policy']),
                'tables_with_public_policy': sum(1 for c in coverage if c['has_public_policy']),
                'average_coverage_score': sum(c['coverage_score'] for c in coverage) / len(coverage) if coverage else 0,
                'tables_needing_attention': [c for c in coverage if c['coverage_score'] < 75]
            }
        }
        
        return {
            'policies': policies,
            'policies_by_table': {
                table: [{'name': p['policyname'], 'type': p['policy_type']} for p in table_policies]
                for table, table_policies in policies_by_table.items()
            },
            'coverage': coverage,
            'stats': stats
        }
    except Exception as e:
        logger.error(f"Error getting RLS policy stats: {e}")
        return {
            'policies': [],
            'policies_by_table': {},
            'coverage': [],
            'stats': {
                'total_policies': 0,
                'tables_with_policies': 0,
                'policy_types': {},
                'coverage': {
                    'tables_analyzed': 0,
                    'tables_with_rls': 0,
                    'tables_with_service_role_policy': 0,
                    'tables_with_admin_policy': 0,
                    'tables_with_owner_policy': 0,
                    'tables_with_public_policy': 0,
                    'average_coverage_score': 0,
                    'tables_needing_attention': []
                }
            }
        }

def verify_cached_function_performance():
    """Benchmark the performance of cached vs. uncached functions."""
    benchmark_query = """
    DO $$
    DECLARE
        start_time TIMESTAMPTZ;
        end_time TIMESTAMPTZ;
        cached_duration INTERVAL;
        uncached_duration INTERVAL;
        i INT;
        result BOOLEAN;
    BEGIN
        -- Benchmark uncached function
        start_time := clock_timestamp();
        FOR i IN 1..1000 LOOP
            SELECT auth.is_service_role() INTO result;
        END LOOP;
        end_time := clock_timestamp();
        uncached_duration := end_time - start_time;
        
        -- Benchmark cached function
        start_time := clock_timestamp();
        FOR i IN 1..1000 LOOP
            SELECT auth.is_service_role_cached() INTO result;
        END LOOP;
        end_time := clock_timestamp();
        cached_duration := end_time - start_time;
        
        -- Log results
        RAISE NOTICE 'Uncached duration: %', uncached_duration;
        RAISE NOTICE 'Cached duration: %', cached_duration;
        RAISE NOTICE 'Performance improvement: % times faster', 
            EXTRACT(EPOCH FROM uncached_duration) / NULLIF(EXTRACT(EPOCH FROM cached_duration), 0);
    END $$;
    """
    
    try:
        db_utils.execute_query(benchmark_query, fetch=False)
        
        # Now fetch actual performance metrics
        performance_query = """
        WITH benchmark AS (
            SELECT
                (SELECT auth.is_service_role()) AS uncached_result,
                (SELECT auth.is_service_role_cached()) AS cached_result,
                clock_timestamp() AS start_time
        ), uncached_test AS (
            SELECT auth.is_service_role() AS result
            FROM generate_series(1, 1000)
        ), uncached_done AS (
            SELECT clock_timestamp() AS end_time
        ), cached_test AS (
            SELECT auth.is_service_role_cached() AS result
            FROM generate_series(1, 1000)
        ), cached_done AS (
            SELECT clock_timestamp() AS end_time
        )
        SELECT
            b.start_time AS start_time,
            ud.end_time AS uncached_end_time,
            cd.end_time AS cached_end_time,
            ud.end_time - b.start_time AS uncached_duration,
            cd.end_time - ud.end_time AS cached_duration,
            EXTRACT(EPOCH FROM (ud.end_time - b.start_time)) / 
                NULLIF(EXTRACT(EPOCH FROM (cd.end_time - ud.end_time)), 0) AS improvement_factor
        FROM
            benchmark b,
            uncached_done ud,
            cached_done cd;
        """
        
        result = db_utils.execute_query(performance_query, cursor_factory=RealDictCursor)[0]
        
        return {
            'uncached_duration_ms': result['uncached_duration'].total_seconds() * 1000,
            'cached_duration_ms': result['cached_duration'].total_seconds() * 1000,
            'improvement_factor': result['improvement_factor'],
            'success': True
        }
    except Exception as e:
        logger.error(f"Error benchmarking function performance: {e}")
        return {
            'uncached_duration_ms': 0,
            'cached_duration_ms': 0,
            'improvement_factor': 0,
            'success': False,
            'error': str(e)
        }

def save_report(before_stats, after_stats, performance_test, filename_prefix="service_role_policy_consolidation_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'before': before_stats,
        'after': after_stats,
        'performance_test': performance_test,
        'summary': {
            'policy_reduction': before_stats['stats']['total_policies'] - after_stats['stats']['total_policies'],
            'policy_reduction_percentage': (before_stats['stats']['total_policies'] - after_stats['stats']['total_policies']) / before_stats['stats']['total_policies'] * 100 if before_stats['stats']['total_policies'] > 0 else 0,
            'coverage_improvement': after_stats['stats']['coverage']['average_coverage_score'] - before_stats['stats']['coverage']['average_coverage_score'],
            'performance_improvement': performance_test['improvement_factor'] if performance_test['success'] else 0
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Consolidate duplicative RLS policies, with a focus on service_role policies.")
    parser.add_argument("--test-only", action="store_true", help="Only test the existing policies without applying the migration")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Get RLS policy stats before migration
    logger.info("Getting RLS policy stats before migration...")
    before_stats = get_policy_stats_before_after()
    
    if not args.test_only:
        # Execute the migration
        logger.info("Executing migration to consolidate service_role policies...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Get RLS policy stats after migration
    logger.info("Getting RLS policy stats after migration...")
    after_stats = get_policy_stats_before_after()
    
    # Test cached function performance
    logger.info("Testing cached function performance...")
    performance_test = verify_cached_function_performance()
    
    if performance_test['success']:
        logger.info(f"Cached function is {performance_test['improvement_factor']:.2f}x faster than uncached function")
    else:
        logger.warning("Failed to test cached function performance")
    
    # Print summary
    before_policies = before_stats['stats']['total_policies']
    after_policies = after_stats['stats']['total_policies']
    policy_reduction = before_policies - after_policies
    
    if policy_reduction > 0:
        logger.info(f"Successfully consolidated {policy_reduction} policies ({policy_reduction / before_policies * 100:.1f}%)")
    else:
        logger.info("No policy consolidation was needed")
    
    # Save report
    report_file = save_report(before_stats, after_stats, performance_test)
    if report_file:
        logger.info(f"Report saved to {report_file}")
    
    # Check if we've met the goal of consolidating service_role policies
    service_role_before = sum(1 for p in before_stats['policies'] 
                              if p['policy_type'] in ('Service Role (uncached)', 'Service Role (by name)'))
    service_role_after = sum(1 for p in after_stats['policies'] 
                             if p['policy_type'] in ('Service Role (uncached)', 'Service Role (by name)'))
    service_role_cached_after = sum(1 for p in after_stats['policies'] 
                                    if p['policy_type'] == 'Service Role (cached)')
    
    if service_role_after < service_role_before:
        logger.info(f"Successfully consolidated service_role policies: {service_role_before} -> {service_role_after}")
        print(f"\n✅ Service role policy consolidation successful!")
        print(f"   - Reduced total policies: {before_policies} -> {after_policies} ({policy_reduction} fewer policies)")
        print(f"   - Reduced service_role policies: {service_role_before} -> {service_role_after}")
        print(f"   - Added cached service_role policies: {service_role_cached_after}")
        if performance_test['success']:
            print(f"   - Performance improvement: {performance_test['improvement_factor']:.2f}x faster")
        sys.exit(0)
    else:
        logger.warning(f"Failed to consolidate service_role policies: {service_role_before} -> {service_role_after}")
        print(f"\n⚠️ Service role policy consolidation did not reduce the number of policies")
        print(f"   - Total policies: {before_policies} -> {after_policies}")
        print(f"   - Service_role policies: {service_role_before} -> {service_role_after}")
        sys.exit(1)

if __name__ == "__main__":
    main()
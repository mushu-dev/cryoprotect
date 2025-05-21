#!/usr/bin/env python3
"""
Script to enhance RLS policies for better access control.
This script applies the SQL migration for RLS policy improvements
and tests the new functionality.
"""

import os
import sys
import argparse
import json
import logging
from datetime import datetime
import db_utils
from psycopg2.extras import RealDictCursor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to enhance RLS policies."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '028_enhance_rls_policies.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to enhance RLS policies")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def get_rls_policy_stats():
    """Get statistics on RLS policies before and after the migration."""
    query = """
    SELECT 
        schemaname, 
        tablename, 
        policyname, 
        roles, 
        cmd, 
        qual
    FROM 
        pg_policies
    WHERE 
        schemaname = 'public'
    ORDER BY 
        tablename, 
        policyname;
    """
    
    try:
        policies = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Analyze policies for potential issues
        analysis = {
            'total_policies': len(policies),
            'tables_with_policies': len(set(p['tablename'] for p in policies)),
            'potential_duplicate_policies': [],
            'tables_without_service_role': [],
            'tables_with_complex_policies': []
        }
        
        # Group policies by table
        policies_by_table = {}
        for policy in policies:
            table = policy['tablename']
            if table not in policies_by_table:
                policies_by_table[table] = []
            policies_by_table[table].append(policy)
        
        # Find potential duplicate policies (same table, same cmd, similar qual)
        for table, table_policies in policies_by_table.items():
            for i, policy1 in enumerate(table_policies):
                for j, policy2 in enumerate(table_policies):
                    if i >= j:
                        continue
                    
                    if policy1['cmd'] == policy2['cmd'] and ('service_role' in policy1['qual'] and 'service_role' in policy2['qual']):
                        analysis['potential_duplicate_policies'].append({
                            'table': table,
                            'policy1': policy1['policyname'],
                            'policy2': policy2['policyname'],
                            'cmd': policy1['cmd']
                        })
            
            # Check if this table has service_role policy
            has_service_role = any('service_role' in p['qual'] for p in table_policies)
            if not has_service_role:
                analysis['tables_without_service_role'].append(table)
            
            # Check for complex policies
            complex_policies = [p for p in table_policies if p['qual'] and len(p['qual']) > 200]
            if complex_policies:
                analysis['tables_with_complex_policies'].append({
                    'table': table,
                    'complex_policies': [p['policyname'] for p in complex_policies]
                })
        
        return {
            'policies': policies,
            'analysis': analysis
        }
    except Exception as e:
        logger.error(f"Error getting RLS policy stats: {e}")
        return {
            'policies': [],
            'analysis': {
                'total_policies': 0,
                'tables_with_policies': 0,
                'potential_duplicate_policies': [],
                'tables_without_service_role': [],
                'tables_with_complex_policies': []
            }
        }

def test_rls_functions():
    """Test the new RLS functions."""
    tests = []
    
    try:
        # Get a molecule ID to test with
        query = "SELECT id, name FROM molecules LIMIT 1"
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if not result:
            logger.error("No molecules found to test with")
            return []
        
        molecule_id = result[0]['id']
        molecule_name = result[0]['name']
        logger.info(f"Testing with molecule: {molecule_name} (ID: {molecule_id})")
        
        # Test has_access_to_molecule function
        query = f"""
        SELECT auth.has_access_to_molecule('{molecule_id}');
        """
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)[0]
        has_access = result['has_access_to_molecule']
        logger.info(f"Has access to molecule {molecule_name}: {has_access}")
        
        tests.append({
            'test': 'has_access_to_molecule',
            'molecule_id': str(molecule_id),
            'molecule_name': molecule_name,
            'has_access': has_access,
            'success': has_access is not None
        })
        
        # Get a mixture ID to test with
        query = "SELECT id, name FROM mixtures LIMIT 1"
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if result:
            mixture_id = result[0]['id']
            mixture_name = result[0]['name']
            logger.info(f"Testing with mixture: {mixture_name} (ID: {mixture_id})")
            
            # Test has_access_to_mixture function
            query = f"""
            SELECT auth.has_access_to_mixture('{mixture_id}');
            """
            result = db_utils.execute_query(query, cursor_factory=RealDictCursor)[0]
            has_access = result['has_access_to_mixture']
            logger.info(f"Has access to mixture {mixture_name}: {has_access}")
            
            tests.append({
                'test': 'has_access_to_mixture',
                'mixture_id': str(mixture_id),
                'mixture_name': mixture_name,
                'has_access': has_access,
                'success': has_access is not None
            })
        else:
            logger.warning("No mixtures found to test with")
        
        # Test test_rls_access function
        if molecule_id:
            query = f"""
            SELECT * FROM auth.test_rls_access('molecule', '{molecule_id}');
            """
            result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
            logger.info(f"RLS access test for molecule {molecule_name}: {len(result)} tests run")
            
            tests.append({
                'test': 'test_rls_access_molecule',
                'molecule_id': str(molecule_id),
                'results': [{
                    'access_type': r['access_type'],
                    'has_access': r['has_access'],
                    'reason': r['reason']
                } for r in result],
                'success': len(result) > 0
            })
        
        return tests
    except Exception as e:
        logger.error(f"Error testing RLS functions: {e}")
        return []

def save_report(stats, tests, filename_prefix="rls_policy_enhancement_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    all_tests_passed = all(test.get('success', False) for test in tests)
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'rls_policy_stats': stats,
        'tests': tests,
        'all_tests_passed': all_tests_passed,
        'total_tests': len(tests),
        'successful_tests': sum(1 for test in tests if test.get('success', False))
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
    parser = argparse.ArgumentParser(description="Enhance RLS policies for better access control.")
    parser.add_argument("--test-only", action="store_true", help="Only test the existing policies without applying the migration")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Get RLS policy stats before migration
    logger.info("Getting RLS policy stats before migration...")
    before_stats = get_rls_policy_stats()
    
    if not args.test_only:
        # Execute the migration
        logger.info("Executing migration to enhance RLS policies...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Get RLS policy stats after migration
    logger.info("Getting RLS policy stats after migration...")
    after_stats = get_rls_policy_stats()
    
    # Test the RLS functions
    logger.info("Testing RLS functions...")
    tests = test_rls_functions()
    
    if not tests:
        logger.error("Failed to run tests. Exiting.")
        sys.exit(1)
    
    # Print test results
    all_passed = all(test.get('success', False) for test in tests)
    logger.info(f"Test results: {sum(1 for test in tests if test.get('success', False))}/{len(tests)} tests passed")
    
    if all_passed:
        logger.info("All tests passed!")
    else:
        logger.warning("Some tests failed. Check the report for details.")
    
    # Prepare stats for report
    stats = {
        'before': before_stats,
        'after': after_stats
    }
    
    # Save report
    report_file = save_report(stats, tests)
    if report_file:
        logger.info(f"Report saved to {report_file}")
    
    if not all_passed and not args.test_only:
        logger.warning("Migration was applied but some tests failed.")
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()
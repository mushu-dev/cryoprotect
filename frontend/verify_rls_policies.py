#!/usr/bin/env python3
"""
Verification script for RLS policies.
This script verifies that the RLS policies are working as expected.
"""

import os
import sys
import json
import logging
from datetime import datetime
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def verify_rls_helper_functions():
    """Verify that the RLS helper functions are present and working."""
    functions_to_check = [
        'auth.has_access_to_molecule',
        'auth.has_access_to_mixture',
        'auth.get_user_clearance_level',
        'auth.clearance_level_value',
        'auth.has_clearance',
        'auth.test_rls_access'
    ]
    
    results = {}
    
    try:
        for function_name in functions_to_check:
            schema, name = function_name.split('.')
            query = """
            SELECT COUNT(*) 
            FROM pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
            WHERE n.nspname = %s AND p.proname = %s;
            """
            
            count = db_utils.execute_query(query, (schema, name))[0][0]
            results[function_name] = count > 0
            
            if count > 0:
                logger.info(f"Function {function_name} exists")
            else:
                logger.warning(f"Function {function_name} does not exist")
        
        return results
    except Exception as e:
        logger.error(f"Error verifying RLS helper functions: {e}")
        return {func: False for func in functions_to_check}

def verify_rls_policies():
    """Verify that the RLS policies are properly applied to all tables."""
    main_tables = [
        'molecules',
        'mixtures',
        'molecular_properties',
        'mixture_components'
    ]
    
    policy_results = {}
    
    try:
        for table in main_tables:
            query = """
            SELECT
                policyname,
                cmd,
                roles,
                qual
            FROM
                pg_policies
            WHERE
                schemaname = 'public'
                AND tablename = %s;
            """
            
            policies = db_utils.execute_query(query, (table,), cursor_factory=RealDictCursor)
            
            # Check if we have both access and modify policies
            has_access_policy = any('access_policy' in p['policyname'] for p in policies)
            has_modify_policy = any('modify_policy' in p['policyname'] for p in policies)
            
            policy_results[table] = {
                'policies': [{
                    'name': p['policyname'],
                    'cmd': p['cmd'],
                    'roles': p['roles'],
                    'qual_length': len(p['qual']) if p['qual'] else 0
                } for p in policies],
                'has_access_policy': has_access_policy,
                'has_modify_policy': has_modify_policy,
                'total_policies': len(policies)
            }
            
            if not has_access_policy:
                logger.warning(f"Table {table} does not have an access policy")
            
            if not has_modify_policy:
                logger.warning(f"Table {table} does not have a modify policy")
        
        return policy_results
    except Exception as e:
        logger.error(f"Error verifying RLS policies: {e}")
        return {table: {'error': str(e)} for table in main_tables}

def test_molecule_access():
    """Test access to molecules using the RLS helper functions."""
    try:
        # Get a few molecules to test with
        query = """
        SELECT id, name, is_public, created_by 
        FROM molecules 
        ORDER BY created_at DESC
        LIMIT 5;
        """
        
        molecules = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if not molecules:
            logger.warning("No molecules found to test with")
            return []
        
        test_results = []
        
        for molecule in molecules:
            molecule_id = molecule['id']
            molecule_name = molecule['name']
            is_public = molecule['is_public']
            
            # Test direct access
            query = f"""
            SELECT auth.has_access_to_molecule('{molecule_id}') as has_access;
            """
            access_result = db_utils.execute_query(query, cursor_factory=RealDictCursor)[0]
            has_access = access_result['has_access']
            
            # Test using the test_rls_access function
            query = f"""
            SELECT * FROM auth.test_rls_access('molecule', '{molecule_id}');
            """
            access_tests = db_utils.execute_query(query, cursor_factory=RealDictCursor)
            
            test_result = {
                'molecule_id': str(molecule_id),
                'molecule_name': molecule_name,
                'is_public': is_public,
                'has_access': has_access,
                'access_tests': [{
                    'access_type': t['access_type'],
                    'has_access': t['has_access'],
                    'reason': t['reason']
                } for t in access_tests],
                'success': has_access is not None
            }
            
            test_results.append(test_result)
            
            logger.info(f"Tested access to molecule {molecule_name}: {'✓' if has_access else '✗'}")
        
        return test_results
    except Exception as e:
        logger.error(f"Error testing molecule access: {e}")
        return []

def test_mixture_access():
    """Test access to mixtures using the RLS helper functions."""
    try:
        # Get a few mixtures to test with
        query = """
        SELECT id, name, is_public, created_by 
        FROM mixtures 
        ORDER BY created_at DESC
        LIMIT 5;
        """
        
        mixtures = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if not mixtures:
            logger.warning("No mixtures found to test with")
            return []
        
        test_results = []
        
        for mixture in mixtures:
            mixture_id = mixture['id']
            mixture_name = mixture['name']
            is_public = mixture['is_public']
            
            # Test direct access
            query = f"""
            SELECT auth.has_access_to_mixture('{mixture_id}') as has_access;
            """
            access_result = db_utils.execute_query(query, cursor_factory=RealDictCursor)[0]
            has_access = access_result['has_access']
            
            # Test using the test_rls_access function
            query = f"""
            SELECT * FROM auth.test_rls_access('mixture', '{mixture_id}');
            """
            access_tests = db_utils.execute_query(query, cursor_factory=RealDictCursor)
            
            test_result = {
                'mixture_id': str(mixture_id),
                'mixture_name': mixture_name,
                'is_public': is_public,
                'has_access': has_access,
                'access_tests': [{
                    'access_type': t['access_type'],
                    'has_access': t['has_access'],
                    'reason': t['reason']
                } for t in access_tests],
                'success': has_access is not None
            }
            
            test_results.append(test_result)
            
            logger.info(f"Tested access to mixture {mixture_name}: {'✓' if has_access else '✗'}")
        
        return test_results
    except Exception as e:
        logger.error(f"Error testing mixture access: {e}")
        return []

def save_verification_report(function_results, policy_results, 
                             molecule_tests, mixture_tests, 
                             filename_prefix="rls_verification_report"):
    """Save the verification report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    # Calculate success rates
    function_success_rate = sum(1 for _, exists in function_results.items() if exists) / len(function_results) if function_results else 0
    molecule_success_rate = sum(1 for test in molecule_tests if test.get('has_access', False)) / len(molecule_tests) if molecule_tests else 0
    mixture_success_rate = sum(1 for test in mixture_tests if test.get('has_access', False)) / len(mixture_tests) if mixture_tests else 0
    
    policy_coverage = sum(1 for _, result in policy_results.items() 
                         if result.get('has_access_policy', False) and result.get('has_modify_policy', False))
    policy_coverage_rate = policy_coverage / len(policy_results) if policy_results else 0
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'function_verification': {
            'results': function_results,
            'success_rate': function_success_rate
        },
        'policy_verification': {
            'results': policy_results,
            'policy_coverage_rate': policy_coverage_rate
        },
        'molecule_access_tests': molecule_tests,
        'mixture_access_tests': mixture_tests,
        'summary': {
            'function_success_rate': function_success_rate,
            'policy_coverage_rate': policy_coverage_rate,
            'molecule_success_rate': molecule_success_rate,
            'mixture_success_rate': mixture_success_rate,
            'overall_success': function_success_rate == 1.0 and policy_coverage_rate == 1.0
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Verification report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving verification report: {e}")
        return None

def main():
    """Main function to run the verification script."""
    parser = argparse.ArgumentParser(description="Verify RLS policies.")
    parser.add_argument("--quiet", action="store_true", help="Suppress INFO logging")
    args = parser.parse_args()
    
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Verify RLS helper functions
    logger.info("Verifying RLS helper functions...")
    function_results = verify_rls_helper_functions()
    
    # Verify RLS policies
    logger.info("Verifying RLS policies...")
    policy_results = verify_rls_policies()
    
    # Test molecule access
    logger.info("Testing molecule access...")
    molecule_tests = test_molecule_access()
    
    # Test mixture access
    logger.info("Testing mixture access...")
    mixture_tests = test_mixture_access()
    
    # Save verification report
    report_file = save_verification_report(
        function_results, 
        policy_results, 
        molecule_tests, 
        mixture_tests
    )
    
    if report_file:
        logger.info(f"Verification report saved to {report_file}")
        
        # Calculate success rates for final output
        function_success = all(function_results.values())
        policy_coverage = all(result.get('has_access_policy', False) and result.get('has_modify_policy', False) 
                             for _, result in policy_results.items())
        molecule_success = all(test.get('has_access', False) for test in molecule_tests) if molecule_tests else False
        mixture_success = all(test.get('has_access', False) for test in mixture_tests) if mixture_tests else False
        
        overall_success = function_success and policy_coverage and molecule_success and mixture_success
        
        if overall_success:
            logger.info("✓ All RLS policy verifications passed!")
            print("\n✓ RLS POLICIES VERIFICATION SUCCESSFUL")
            sys.exit(0)
        else:
            logger.warning("✗ Some RLS policy verifications failed. Check the report for details.")
            print("\n✗ RLS POLICIES VERIFICATION FAILED - See report for details")
            sys.exit(1)
    else:
        logger.error("Failed to save verification report.")
        sys.exit(1)

if __name__ == "__main__":
    main()
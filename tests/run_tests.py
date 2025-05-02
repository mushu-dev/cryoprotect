#!/usr/bin/env python
"""
CryoProtect Analyzer - Comprehensive Test Runner

This script runs all tests for the CryoProtect Analyzer project.
It executes Python unit tests, JavaScript database tests, API endpoint tests,
database schema validation, authentication tests, and performance benchmarks.
"""

import os
import sys
import json
import unittest
import subprocess
import argparse
import time
from datetime import datetime

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def run_python_tests(verbose=False, pattern=None):
    """Run Python unit tests."""
    print("\n===== Running Python Tests =====")
    
    # Discover and run tests
    loader = unittest.TestLoader()
    if pattern:
        test_suite = loader.discover('tests', pattern=pattern)
    else:
        test_suite = loader.discover('tests')
    
    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(test_suite)
    
    return result.wasSuccessful(), result

def run_js_tests(verbose=False):
    """Run JavaScript database tests."""
    print("\n===== Running JavaScript Database Tests =====")
    
    # Run the JavaScript test file using Node.js
    cmd = ['node', 'tests/test_database_schema.js']
    if verbose:
        cmd.append('--verbose')
    
    try:
        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if verbose or result.returncode != 0:
            print(result.stdout)
            if result.stderr:
                print(result.stderr)
        return result.returncode == 0, result.stdout
    except Exception as e:
        print(f"Error running JavaScript tests: {e}")
        return False, str(e)

def run_api_endpoint_tests(verbose=False):
    """Run API endpoint tests with different user roles."""
    print("\n===== Running API Endpoint Tests =====")
    
    # Run the API endpoint tests
    loader = unittest.TestLoader()
    test_suite = loader.discover('tests', pattern='test_api_endpoints.py')
    
    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(test_suite)
    
    return result.wasSuccessful(), result

def run_database_schema_validation(verbose=False):
    """Run database schema validation tests."""
    print("\n===== Running Database Schema Validation =====")
    
    # Run the database schema validation tests
    loader = unittest.TestLoader()
    test_suite = loader.discover('tests', pattern='test_database_schema_validation.py')
    
    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(test_suite)
    
    return result.wasSuccessful(), result

def run_authentication_tests(verbose=False):
    """Run authentication tests."""
    print("\n===== Running Authentication Tests =====")
    
    # Run the authentication tests
    loader = unittest.TestLoader()
    test_suite = loader.discover('tests', pattern='test_authentication.py')
    
    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(test_suite)
    
    return result.wasSuccessful(), result

def run_performance_tests(verbose=False):
    """Run performance benchmark tests."""
    print("\n===== Running Performance Benchmark Tests =====")
    
    # Run the performance benchmark tests
    loader = unittest.TestLoader()
    test_suite = loader.discover('tests', pattern='test_performance.py')
    
    runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
    result = runner.run(test_suite)
    
    return result.wasSuccessful(), result

def generate_report(test_results):
    """Generate a comprehensive test report."""
    print("\n===== Generating Test Report =====")
    
    report_dir = os.path.join(os.path.dirname(__file__), 'reports')
    os.makedirs(report_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_path = os.path.join(report_dir, f'test_report_{timestamp}.json')
    text_report_path = os.path.join(report_dir, f'test_report_{timestamp}.txt')
    
    # Create JSON report
    with open(report_path, 'w') as f:
        json.dump(test_results, f, indent=2)
    
    # Create text report
    with open(text_report_path, 'w') as f:
        f.write("CryoProtect Analyzer - Comprehensive Test Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Overall status
        overall_success = all(result['success'] for result in test_results.values())
        f.write(f"Overall Status: {'PASSED' if overall_success else 'FAILED'}\n\n")
        
        # Test results
        for test_type, result in test_results.items():
            f.write(f"{test_type}: {'PASSED' if result['success'] else 'FAILED'}\n")
            if 'details' in result:
                f.write(f"  Details: {result['details']}\n")
            if 'metrics' in result:
                f.write("  Metrics:\n")
                for metric, value in result['metrics'].items():
                    f.write(f"    {metric}: {value}\n")
            f.write("\n")
        
        # Recommendations
        f.write("Recommendations:\n")
        if not overall_success:
            f.write("  - Fix failing tests before proceeding to production\n")
        
        if 'performance_tests' in test_results and 'metrics' in test_results['performance_tests']:
            metrics = test_results['performance_tests']['metrics']
            if metrics.get('avg_response_time', 0) > 500:
                f.write("  - Optimize database queries to improve response time\n")
            if metrics.get('max_cpu_usage', 0) > 80:
                f.write("  - Consider scaling up database resources to handle the load\n")
    
    print(f"\nTest report generated: {report_path}")
    print(f"Text report generated: {text_report_path}")
    return report_path, text_report_path

def main():
    """Main function to run tests."""
    parser = argparse.ArgumentParser(description='Run CryoProtect Analyzer tests')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('-p', '--python-only', action='store_true', help='Run only Python tests')
    parser.add_argument('-j', '--js-only', action='store_true', help='Run only JavaScript tests')
    parser.add_argument('-a', '--api-only', action='store_true', help='Run only API endpoint tests')
    parser.add_argument('-d', '--db-only', action='store_true', help='Run only database schema validation tests')
    parser.add_argument('-u', '--auth-only', action='store_true', help='Run only authentication tests')
    parser.add_argument('-b', '--perf-only', action='store_true', help='Run only performance benchmark tests')
    parser.add_argument('-t', '--test', help='Run specific test pattern (e.g., test_rdkit*.py)')
    parser.add_argument('-r', '--report', action='store_true', help='Generate test report')
    
    args = parser.parse_args()
    
    # Initialize test results
    test_results = {}
    
    # Run Python tests if requested
    if args.test:
        success, result = run_python_tests(args.verbose, args.test)
        test_results['python_tests'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors"
        }
    elif not (args.js_only or args.api_only or args.db_only or args.auth_only or args.perf_only):
        success, result = run_python_tests(args.verbose)
        test_results['python_tests'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors"
        }
    
    # Run JavaScript tests if requested
    if not args.test and (args.js_only or not (args.python_only or args.api_only or args.db_only or args.auth_only or args.perf_only)):
        success, output = run_js_tests(args.verbose)
        test_results['js_tests'] = {
            'success': success,
            'details': output[:200] + "..." if len(output) > 200 else output
        }
    
    # Run API endpoint tests if requested
    if args.api_only or not (args.python_only or args.js_only or args.db_only or args.auth_only or args.perf_only or args.test):
        success, result = run_api_endpoint_tests(args.verbose)
        test_results['api_endpoint_tests'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors"
        }
    
    # Run database schema validation tests if requested
    if args.db_only or not (args.python_only or args.js_only or args.api_only or args.auth_only or args.perf_only or args.test):
        success, result = run_database_schema_validation(args.verbose)
        test_results['database_schema_validation'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors"
        }
    
    # Run authentication tests if requested
    if args.auth_only or not (args.python_only or args.js_only or args.api_only or args.db_only or args.perf_only or args.test):
        success, result = run_authentication_tests(args.verbose)
        test_results['authentication_tests'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors"
        }
    
    # Run performance benchmark tests if requested
    if args.perf_only or not (args.python_only or args.js_only or args.api_only or args.db_only or args.auth_only or args.test):
        start_time = time.time()
        success, result = run_performance_tests(args.verbose)
        end_time = time.time()
        
        # Extract performance metrics from the test results
        metrics = {
            'execution_time': round(end_time - start_time, 2),
            'avg_response_time': None,  # Will be populated by the performance tests
            'max_response_time': None,  # Will be populated by the performance tests
            'throughput': None,         # Will be populated by the performance tests
            'max_cpu_usage': None,      # Will be populated by the performance tests
            'max_memory_usage': None    # Will be populated by the performance tests
        }
        
        test_results['performance_tests'] = {
            'success': success,
            'details': f"Ran {result.testsRun} tests, {len(result.failures)} failures, {len(result.errors)} errors",
            'metrics': metrics
        }
    
    # Generate report if requested
    if args.report:
        report_path, text_report_path = generate_report(test_results)
    
    # Print summary
    print("\n===== Test Summary =====")
    for test_type, result in test_results.items():
        print(f"{test_type}: {'PASSED' if result['success'] else 'FAILED'}")
    
    # Return exit code based on test results
    return 0 if all(result['success'] for result in test_results.values()) else 1

if __name__ == '__main__':
    sys.exit(main())
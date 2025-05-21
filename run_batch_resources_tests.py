"""
Test runner script for batch_resources module.

This script runs the tests for the batch_resources module and measures the coverage
before and after running the tests. It generates a coverage report and displays
the coverage percentage for the module.
"""

import os
import sys
import coverage
import unittest
import importlib
from datetime import datetime

def measure_initial_coverage():
    """Measure initial coverage of the batch_resources module."""
    print("\n=== MEASURING INITIAL COVERAGE ===")
    
    # Create a coverage object with appropriate settings
    cov = coverage.Coverage(
        source=['api.batch_resources'],
        omit=['*/__pycache__/*', '*/.*/'],
        include=['*/api/batch_resources.py']
    )
    
    # Start coverage measurement
    cov.start()
    
    try:
        # Import the module to measure coverage
        import api.batch_resources
        
        # Force reload to ensure coverage is measured
        importlib.reload(api.batch_resources)
        
        # Create instances to verify imports
        resource = api.batch_resources.BatchOperationResource()
        service = api.batch_resources.BatchOperationService()
        
        # Print module info to verify import
        print(f"Successfully imported BatchOperationResource from {api.batch_resources.__file__}")
        print(f"Successfully imported BatchOperationService from {api.batch_resources.__file__}")
        
        # Print basic coverage info
        print("\nBasic module structure coverage:")
        print(f"- Module contains {len(api.batch_resources.__dict__)} items")
        print(f"- BatchOperationService has {len([m for m in dir(api.batch_resources.BatchOperationService) if not m.startswith('_')])} public methods")
        print(f"- BatchOperationResource has {len([m for m in dir(api.batch_resources.BatchOperationResource) if not m.startswith('_')])} public methods")
        
    finally:
        # Stop coverage measurement
        cov.stop()
        
        # Generate report
        print("\nInitial Coverage Report:")
        cov.report()
        
        # Save HTML report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_dir = f'reports/coverage_batch_resources_initial_{timestamp}'
        os.makedirs(report_dir, exist_ok=True)
        cov.html_report(directory=report_dir)
        
        print(f"\nInitial HTML coverage report saved to {report_dir}")
    
    return cov

def run_tests_and_measure_coverage():
    """Run tests and measure coverage after tests."""
    print("\n=== RUNNING TESTS AND MEASURING COVERAGE ===")
    
    # Create a coverage object with appropriate settings
    cov = coverage.Coverage(
        source=['api.batch_resources'],
        omit=['*/__pycache__/*', '*/.*/'],
        include=['*/api/batch_resources.py']
    )
    
    # Start coverage measurement
    cov.start()
    
    try:
        # Run the tests
        print("\nRunning tests for BatchOperationService...")
        from tests.test_batch_resources import TestBatchOperationService
        
        service_suite = unittest.TestLoader().loadTestsFromTestCase(TestBatchOperationService)
        service_result = unittest.TextTestRunner(verbosity=2).run(service_suite)
        
        print(f"\nBatchOperationService test results: {service_result.testsRun} tests run, "
              f"{len(service_result.errors)} errors, {len(service_result.failures)} failures")
        
        print("\nRunning Flask-specific tests for BatchOperationResource...")
        from tests.test_batch_resources_flask_only import TestBatchOperationResourceFlask
        
        flask_suite = unittest.TestLoader().loadTestsFromTestCase(TestBatchOperationResourceFlask)
        flask_result = unittest.TextTestRunner(verbosity=2).run(flask_suite)
        
        print(f"\nFlask-specific test results: {flask_result.testsRun} tests run, "
              f"{len(flask_result.errors)} errors, {len(flask_result.failures)} failures")
        
        # Calculate total test results
        total_tests = service_result.testsRun + flask_result.testsRun
        total_errors = len(service_result.errors) + len(flask_result.errors)
        total_failures = len(service_result.failures) + len(flask_result.failures)
        
        print(f"\nTotal test results: {total_tests} tests run, {total_errors} errors, {total_failures} failures")
        
    finally:
        # Stop coverage measurement
        cov.stop()
        
        # Generate report
        print("\nFinal Coverage Report:")
        cov.report()
        
        # Save HTML report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_dir = f'reports/coverage_batch_resources_final_{timestamp}'
        os.makedirs(report_dir, exist_ok=True)
        cov.html_report(directory=report_dir)
        
        print(f"\nFinal HTML coverage report saved to {report_dir}")
    
    return cov

def main():
    """Main function to run tests and measure coverage."""
    print("=== BATCH RESOURCES TEST RUNNER ===")
    print("This script measures coverage before and after running tests.")
    
    # Create reports directory if it doesn't exist
    os.makedirs('reports', exist_ok=True)
    
    # Measure initial coverage
    initial_cov = measure_initial_coverage()
    initial_coverage = initial_cov.get_data().line_percent('api/batch_resources.py')
    
    # Run tests and measure coverage
    final_cov = run_tests_and_measure_coverage()
    final_coverage = final_cov.get_data().line_percent('api/batch_resources.py')
    
    # Print summary
    print("\n=== COVERAGE SUMMARY ===")
    print(f"Initial coverage: {initial_coverage:.2f}%")
    print(f"Final coverage: {final_coverage:.2f}%")
    print(f"Coverage improvement: {final_coverage - initial_coverage:.2f}%")
    
    # Generate a summary report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_file = f'reports/batch_resources_coverage_summary_{timestamp}.txt'
    
    with open(summary_file, 'w') as f:
        f.write("=== BATCH RESOURCES COVERAGE SUMMARY ===\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Initial coverage: {initial_coverage:.2f}%\n")
        f.write(f"Final coverage: {final_coverage:.2f}%\n")
        f.write(f"Coverage improvement: {final_coverage - initial_coverage:.2f}%\n\n")
        f.write("Initial coverage report: reports/coverage_batch_resources_initial_" + timestamp + "\n")
        f.write("Final coverage report: reports/coverage_batch_resources_final_" + timestamp + "\n")
    
    print(f"\nSummary report saved to {summary_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
"""
Script to measure coverage of api/batch_resources.py before and after running tests.

This script demonstrates the improved testability of the refactored batch_resources.py
by showing coverage metrics before and after running the unit tests.
"""

import os
import sys
import coverage
import unittest
import importlib

def measure_module_coverage():
    """Measure coverage of the batch_resources module."""
    print("Measuring coverage of batch_resources.py...")
    
    # Create a coverage object with appropriate settings
    cov = coverage.Coverage(
        source=['api.batch_resources'],
        omit=['*/__pycache__/*', '*/\.*'],
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
        
        # Run some basic tests to improve coverage
        print("\nRunning basic functionality tests...")
        
        # Test validation
        validation_result = api.batch_resources.BatchOperationService.validate_request(
            operation="property_calculation",
            entity_type="molecule",
            ids=["test-id"]
        )
        print(f"- Validation test: {'Passed' if validation_result is None else 'Failed'}")
        
        # Test batch operation processing with invalid operation
        try:
            result = api.batch_resources.BatchOperationService.process_batch_operation(
                operation="invalid_operation",
                entity_type="molecule",
                ids=["test-id"]
            )
            print(f"- Invalid operation test: {'Passed' if result['status'] == 'ERROR' else 'Failed'}")
        except Exception as e:
            print(f"- Invalid operation test: Failed with exception: {e}")
        
        # Try to run the unit tests if they exist
        try:
            print("\nAttempting to run unit tests...")
            from tests.test_batch_resources import TestBatchOperationService
            
            # Run the tests
            test_suite = unittest.TestLoader().loadTestsFromTestCase(TestBatchOperationService)
            test_result = unittest.TextTestRunner(verbosity=2).run(test_suite)
            
            print(f"\nTest results: {test_result.testsRun} tests run, {len(test_result.errors)} errors, {len(test_result.failures)} failures")
            
        except ImportError as e:
            print(f"Could not run unit tests: {e}")
            print("This is expected if you haven't created the test file yet.")
    
    finally:
        # Stop coverage measurement
        cov.stop()
        
        # Generate report
        print("\nCoverage Report:")
        cov.report()
        
        # Save HTML report
        os.makedirs('reports', exist_ok=True)
        cov.html_report(directory='reports/coverage_batch_resources')
        
        print(f"\nHTML coverage report saved to reports/coverage_batch_resources")
    
    return cov

if __name__ == "__main__":
    # Measure module coverage
    cov = measure_module_coverage()
    
    # Print summary
    print("\nSummary:")
    print("--------")
    print("The coverage report shows how much of the code can be tested.")
    print("Higher coverage indicates better testability of the code.")
    print("The refactored code separates business logic from Flask-specific code,")
    print("making it easier to test without needing Flask's request context.")
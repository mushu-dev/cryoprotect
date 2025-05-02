"""
Test runner for predictive_models.py coverage

This script runs the tests in test_predictive_models.py
and generates a coverage report for api/predictive_models.py.
"""

import os
import sys
import unittest
import coverage

# Start coverage
cov = coverage.Coverage(source=['api.predictive_models'])
cov.start()

# Import our test module directly
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from tests.test_predictive_models import (
    TestPredictiveModel,
    TestModelManager,
    TestPredictiveFunctions
)

# Create a test suite
def create_test_suite():
    suite = unittest.TestSuite()
    
    # Add tests from test_predictive_models.py
    suite.addTest(unittest.makeSuite(TestPredictiveModel))
    suite.addTest(unittest.makeSuite(TestModelManager))
    suite.addTest(unittest.makeSuite(TestPredictiveFunctions))
    
    return suite

if __name__ == '__main__':
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(create_test_suite())
    
    # Stop coverage and generate report
    cov.stop()
    cov.save()
    
    print("\nCoverage Report:")
    cov.report()
    
    # Generate HTML report
    cov.html_report(directory='reports/predictive_models_coverage_html')
    print(f"\nHTML coverage report generated in reports/predictive_models_coverage_html")
    
    # Exit with appropriate code
    sys.exit(not result.wasSuccessful())
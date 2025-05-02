"""
Simple test runner for models.py coverage

This script runs the tests in test_models_focused.py and test_models_focused_2.py
and generates a coverage report for api/models.py.
"""

import os
import sys
import unittest
import coverage

# Start coverage
cov = coverage.Coverage(source=['api.models'])
cov.start()

# Import our test modules directly
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from tests.test_models_focused import (
    TestFlexibleDateTime,
    TestSchemas,
    TestPropertyType,
    TestMolecularProperty,
    TestUserProfile
)
from tests.test_models_focused_2 import (
    TestCalculationMethod,
    TestPrediction,
    TestExperiment,
    TestComparison,
    TestProject,
    TestProjectExperiment
)

# Create a test suite
def create_test_suite():
    suite = unittest.TestSuite()
    
    # Add tests from test_models_focused.py
    suite.addTest(unittest.makeSuite(TestFlexibleDateTime))
    suite.addTest(unittest.makeSuite(TestSchemas))
    suite.addTest(unittest.makeSuite(TestPropertyType))
    suite.addTest(unittest.makeSuite(TestMolecularProperty))
    suite.addTest(unittest.makeSuite(TestUserProfile))
    
    # Add tests from test_models_focused_2.py
    suite.addTest(unittest.makeSuite(TestCalculationMethod))
    suite.addTest(unittest.makeSuite(TestPrediction))
    suite.addTest(unittest.makeSuite(TestExperiment))
    suite.addTest(unittest.makeSuite(TestComparison))
    suite.addTest(unittest.makeSuite(TestProject))
    suite.addTest(unittest.makeSuite(TestProjectExperiment))
    
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
    cov.html_report(directory='reports/coverage_html')
    print(f"\nHTML coverage report generated in reports/coverage_html")
    
    # Exit with appropriate code
    sys.exit(not result.wasSuccessful())
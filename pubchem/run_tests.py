#!/usr/bin/env python
"""
Test runner for the PubChem client package.

This script runs all unit tests for the PubChem client package.
"""

import unittest
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Import test modules
from test_cache import TestPubChemCache
from test_rate_limiter import TestAdaptiveRateLimiter
from test_utils import TestCircuitBreaker, TestRetryWithBackoff, TestCreateRequestKey
from test_scheduler import TestWeekendJobScheduler
from test_client import TestResilientPubChemClient

if __name__ == "__main__":
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Add test cases
    test_suite.addTest(unittest.makeSuite(TestPubChemCache))
    test_suite.addTest(unittest.makeSuite(TestAdaptiveRateLimiter))
    test_suite.addTest(unittest.makeSuite(TestCircuitBreaker))
    test_suite.addTest(unittest.makeSuite(TestRetryWithBackoff))
    test_suite.addTest(unittest.makeSuite(TestCreateRequestKey))
    test_suite.addTest(unittest.makeSuite(TestWeekendJobScheduler))
    test_suite.addTest(unittest.makeSuite(TestResilientPubChemClient))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    print(f"\nTest Summary:")
    print(f"  Ran {result.testsRun} tests")
    print(f"  Failures: {len(result.failures)}")
    print(f"  Errors: {len(result.errors)}")
    print(f"  Skipped: {len(result.skipped)}")
    
    # Exit with appropriate code
    if result.wasSuccessful():
        print("\nAll tests passed successfully!")
        exit(0)
    else:
        print("\nSome tests failed!")
        exit(1)
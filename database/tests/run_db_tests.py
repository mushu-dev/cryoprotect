#!/usr/bin/env python
"""
Database Tests Runner

This script discovers and runs all database-related tests.
"""

import os
import sys
import unittest
import importlib
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))


def run_tests():
    """Discover and run all database tests."""
    print("=" * 70)
    print("CryoProtect Database Test Runner")
    print("=" * 70)
    
    # Get the directory containing this script
    test_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Create test suite
    test_suite = unittest.defaultTestLoader.discover(
        start_dir=test_dir,
        pattern='test_*.py'
    )
    
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    # Print summary
    print("\n" + "=" * 70)
    print(f"Tests run: {result.testsRun}, Errors: {len(result.errors)}, Failures: {len(result.failures)}")
    print("=" * 70)
    
    # Return success or failure
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    # Create tests directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(__file__)), exist_ok=True)
    
    # Run tests
    sys.exit(run_tests())
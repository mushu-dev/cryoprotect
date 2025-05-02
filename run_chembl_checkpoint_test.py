#!/usr/bin/env python
"""
Test runner for ChEMBL checkpoint tests.
"""

import os
import sys
import unittest

# Add the parent directory to the path so we can import the chembl package
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

def run_tests():
    """Run the checkpoint tests."""
    print("\n===== Running ChEMBL Checkpoint Tests =====")
    
    # Discover and run tests
    loader = unittest.TestLoader()
    test_suite = loader.discover('tests/unit/chembl', pattern='test_checkpoint.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)
    
    return result.wasSuccessful()

if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
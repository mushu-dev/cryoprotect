# Task 3.1: Fix ImportError Issues in Test Modules

## Objective
Fix the ImportError issues in test modules to ensure all tests can run successfully.

## Context
The test modules are currently showing import errors that prevent them from running correctly. This hinders our ability to verify that the codebase is functioning properly. We need to fix these import issues to enable comprehensive testing.

## Acceptance Criteria
- All ImportError issues in test modules are resolved
- Tests can be run without import-related errors
- A standardized import approach is implemented across test files
- Test coverage reporting works correctly

## Implementation Steps

1. First, identify the test files with import errors:
   ```bash
   cd tests
   python -m unittest discover . 2>&1 | grep "ImportError"
   ```

2. Create a common test utility module to standardize imports:
   ```python
   # tests/test_utils.py
   
   """
   Test Utilities for CryoProtect v2 tests
   
   This module provides common utilities and import helpers for test modules.
   """
   
   import os
   import sys
   import json
   from pathlib import Path
   
   # Add the parent directory to the path so we can import the api package
   project_root = Path(__file__).parent.parent
   sys.path.insert(0, str(project_root))
   
   # Import common testing libraries
   import unittest
   from unittest.mock import patch, MagicMock, Mock
   
   # Import application modules
   try:
       from config import Config
       from api.utils import get_supabase_client, _handle_json_serialization
       from api.models import BaseModel
   except ImportError as e:
       print(f"Error importing application modules: {e}")
       
   class MockSupabase:
       """Mock Supabase client for testing."""
       
       def __init__(self, data=None, error=None):
           self.data = data or {}
           self.error = error
           
       def from_(self, table):
           """Mock from_ method."""
           self._table = table
           return self
           
       def select(self, fields="*"):
           """Mock select method."""
           self._fields = fields
           return self
           
       def eq(self, field, value):
           """Mock eq method."""
           self._eq_field = field
           self._eq_value = value
           return self
           
       def execute(self):
           """Mock execute method."""
           class Response:
               def __init__(self, data, error):
                   self.data = data
                   self.error = error
                   
           return Response(self.data, self.error)
       
   def create_mock_supabase(data=None, error=None):
       """Create a mock Supabase client for testing."""
       return MockSupabase(data, error)
   
   def setup_test_environment():
       """Set up the test environment."""
       # Configure environment variables for testing
       os.environ['SUPABASE_URL'] = 'https://test.supabase.co'
       os.environ['SUPABASE_KEY'] = 'test-key'
       os.environ['SUPABASE_SERVICE_ROLE_KEY'] = 'test-service-role-key'
       os.environ['USE_SERVICE_ROLE'] = 'true'
       
       # Clear any existing cache or state
       if 'supabase' in globals():
           del globals()['supabase']
   
   def teardown_test_environment():
       """Tear down the test environment."""
       # Clean up environment variables
       for var in ['SUPABASE_URL', 'SUPABASE_KEY', 'SUPABASE_SERVICE_ROLE_KEY', 'USE_SERVICE_ROLE']:
           if var in os.environ:
               del os.environ[var]
   ```

3. Fix the test_models.py file:
   ```python
   # tests/test_models.py - updated version
   
   """
   Tests for the CryoProtect v2 data models
   """
   
   from tests.test_utils import (
       unittest, patch, MagicMock, setup_test_environment, 
       teardown_test_environment, create_mock_supabase
   )
   from api.models import (
       BaseModel, Molecule, Mixture, MixtureComponent,
       Prediction, Experiment, PropertyType
   )
   
   class TestBaseModel(unittest.TestCase):
       """Tests for the BaseModel class."""
       
       def setUp(self):
           """Set up test case."""
           setup_test_environment()
       
       def tearDown(self):
           """Tear down test case."""
           teardown_test_environment()
       
       @patch('api.models.get_supabase_client')
       def test_get_supabase(self, mock_get_supabase):
           """Test the get_supabase method."""
           mock_client = MagicMock()
           mock_get_supabase.return_value = mock_client
           
           client = BaseModel.get_supabase()
           self.assertEqual(client, mock_client)
           mock_get_supabase.assert_called_once()
       
       # Add more test methods...
   
   class TestMolecule(unittest.TestCase):
       """Tests for the Molecule class."""
       
       def setUp(self):
           """Set up test case."""
           setup_test_environment()
       
       def tearDown(self):
           """Tear down test case."""
           teardown_test_environment()
       
       @patch('api.models.BaseModel.get_supabase')
       def test_get(self, mock_get_supabase):
           """Test the get method."""
           # Mock data
           molecule_data = {
               'id': '123',
               'name': 'Test Molecule',
               'formula': 'C2H6O'
           }
           
           # Create mock response
           mock_client = create_mock_supabase([molecule_data])
           mock_get_supabase.return_value = mock_client
           
           # Call the method
           molecule = Molecule.get('123')
           
           # Assertions
           self.assertEqual(molecule['id'], '123')
           self.assertEqual(molecule['name'], 'Test Molecule')
           
       # Add more test methods...
   
   # Add more test classes for other models...
   
   if __name__ == '__main__':
       unittest.main()
   ```

4. Fix the base_test_case.py file to use the new test utilities:
   ```python
   # tests/base_test_case.py - updated version
   
   """
   Base Test Case for CryoProtect v2 tests
   """
   
   from tests.test_utils import (
       unittest, patch, MagicMock, setup_test_environment, 
       teardown_test_environment, create_mock_supabase
   )
   
   class BaseTestCase(unittest.TestCase):
       """Base class for all test cases."""
       
       def setUp(self):
           """Set up test case."""
           setup_test_environment()
       
       def tearDown(self):
           """Tear down test case."""
           teardown_test_environment()
       
       def create_mock_supabase(self, data=None, error=None):
           """Create a mock Supabase client for testing."""
           return create_mock_supabase(data, error)
   ```

5. Update the test runner to use the consistent approach:
   ```python
   # tests/run_tests.py - updated version
   
   #!/usr/bin/env python
   """
   Test Runner for CryoProtect v2 tests
   """
   
   import os
   import sys
   import unittest
   import argparse
   import json
   from datetime import datetime
   from pathlib import Path
   
   # Add the parent directory to the path
   project_root = Path(__file__).parent.parent
   sys.path.insert(0, str(project_root))
   
   # Import test utilities
   from tests.test_utils import setup_test_environment, teardown_test_environment
   
   def run_tests(pattern=None, verbose=False):
       """Run the tests matching the pattern."""
       setup_test_environment()
       
       try:
           # Discover and run tests
           loader = unittest.TestLoader()
           if pattern:
               suite = loader.discover('.', pattern=pattern)
           else:
               suite = loader.discover('.')
           
           runner = unittest.TextTestRunner(verbosity=2 if verbose else 1)
           result = runner.run(suite)
           
           return result
       finally:
           teardown_test_environment()
   
   def generate_report(result):
       """Generate a test report."""
       # Create report directory
       reports_dir = Path('reports')
       reports_dir.mkdir(exist_ok=True)
       
       # Generate report
       timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
       report_path = reports_dir / f'test_report_{timestamp}.json'
       
       report = {
           'timestamp': timestamp,
           'total_tests': result.testsRun,
           'failures': len(result.failures),
           'errors': len(result.errors),
           'skipped': len(result.skipped) if hasattr(result, 'skipped') else 0,
           'success': result.wasSuccessful()
       }
       
       with open(report_path, 'w') as f:
           json.dump(report, f, indent=2)
       
       return report_path
   
   def main():
       """Main function."""
       parser = argparse.ArgumentParser(description='Run CryoProtect v2 tests')
       parser.add_argument('-p', '--pattern', help='Test pattern (e.g., test_*.py)')
       parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
       parser.add_argument('-r', '--report', action='store_true', help='Generate report')
       
       args = parser.parse_args()
       
       result = run_tests(args.pattern, args.verbose)
       
       if args.report:
           report_path = generate_report(result)
           print(f"Report generated: {report_path}")
       
       return 0 if result.wasSuccessful() else 1
   
   if __name__ == '__main__':
       sys.exit(main())
   ```

6. Fix import errors in specific test files by adding the common import approach to each file:
   - Find test files with import errors
   - Update imports to use the test_utils module
   - Ensure consistent setUp and tearDown methods
   - Use the create_mock_supabase function for consistent mocking

## Files to Modify
- Create `tests/test_utils.py`
- Update `tests/base_test_case.py`
- Update `tests/run_tests.py`
- Fix import errors in failing test files

## Verification
1. Run the test discovery to see if import errors are resolved:
   ```bash
   cd tests
   python -m unittest discover .
   ```

2. Run specific tests that previously had import errors:
   ```bash
   python -m unittest tests.test_models
   ```

3. Run the updated test runner:
   ```bash
   python tests/run_tests.py
   ```

4. Generate a test report and verify it's created correctly:
   ```bash
   python tests/run_tests.py --report
   ```

## Notes for Roo Code Agent
- Focus on creating a consistent approach to imports across all test files
- The test_utils.py module acts as a central place for common imports and utilities
- Make sure to properly handle environment setup and teardown in each test
- The MockSupabase class provides a simple way to mock Supabase responses
- Update additional test files as needed to use the new approach
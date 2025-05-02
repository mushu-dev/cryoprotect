# Task: Create Unified Testing Framework

## Objective
Consolidate the scattered testing files into a structured, unified testing framework that reduces duplication and improves maintainability.

## Context
The codebase currently contains 65+ test files with significant duplication and overlapping functionality. A unified testing framework with shared fixtures, utilities, and organized structure would dramatically reduce complexity while improving test coverage.

## Acceptance Criteria
- A structured `tests` package is created with domain-specific submodules
- Common test utilities and fixtures are consolidated
- Test coverage is maintained or improved
- Tests run successfully with the new structure
- Documentation is updated to reflect the new testing approach

## Implementation Steps

1. Create the structured test directory organization:
   ```bash
   mkdir -p tests/api
   mkdir -p tests/database
   mkdir -p tests/rdkit
   mkdir -p tests/models
   mkdir -p tests/integration
   mkdir -p tests/fixtures
   mkdir -p tests/utils
   ```

2. Create the test utilities module:
   ```python
   # tests/utils/__init__.py
   
   """
   Test utilities for CryoProtect v2 tests.
   
   This module provides common utilities for setting up tests,
   creating mock objects, and other test-related functions.
   """
   
   import os
   import sys
   import json
   from pathlib import Path
   from unittest.mock import patch, MagicMock, Mock
   
   # Add the project root to the path
   project_root = Path(__file__).parent.parent.parent
   if str(project_root) not in sys.path:
       sys.path.insert(0, str(project_root))
   
   # Import common test dependencies
   import unittest
   import pytest
   
   # Import application modules that are commonly used in tests
   try:
       from config import Config
       from api.utils import get_supabase_client
   except ImportError as e:
       print(f"Warning: Could not import some application modules: {e}")
   
   
   def setup_test_environment():
       """Set up the test environment with required variables."""
       # Set environment variables for testing
       os.environ['SUPABASE_URL'] = 'https://test.supabase.co'
       os.environ['SUPABASE_KEY'] = 'test-key'
       os.environ['FLASK_ENV'] = 'testing'
       
       # Clear any existing cache
       if hasattr(get_supabase_client, '_instance'):
           delattr(get_supabase_client, '_instance')
   
   
   def teardown_test_environment():
       """Clean up the test environment."""
       # Remove environment variables
       for var in ['SUPABASE_URL', 'SUPABASE_KEY', 'FLASK_ENV']:
           if var in os.environ:
               del os.environ[var]
       
       # Clear any cache
       if hasattr(get_supabase_client, '_instance'):
           delattr(get_supabase_client, '_instance')
   
   
   class MockSupabaseResponse:
       """Mock Supabase response object."""
       
       def __init__(self, data=None, error=None):
           self.data = data if data is not None else []
           self.error = error
   
   
   class MockSupabase:
       """Mock Supabase client for testing."""
       
       def __init__(self, response_data=None, response_error=None):
           self.response_data = response_data
           self.response_error = response_error
           self._table = None
           self._select = "*"
           self._filters = []
       
       def from_(self, table):
           """Mock from_ method."""
           self._table = table
           return self
       
       def select(self, select="*"):
           """Mock select method."""
           self._select = select
           return self
       
       def eq(self, field, value):
           """Mock eq filter method."""
           self._filters.append(("eq", field, value))
           return self
       
       def execute(self):
           """Mock execute method."""
           return MockSupabaseResponse(self.response_data, self.response_error)
       
       def rpc(self, function, params=None):
           """Mock rpc method."""
           return self
   
   
   def create_mock_supabase_client(response_data=None, response_error=None):
       """
       Create a mock Supabase client for testing.
       
       Args:
           response_data: Data to return in the response
           response_error: Error to return in the response
           
       Returns:
           MockSupabase: A mock Supabase client
       """
       return MockSupabase(response_data, response_error)
   
   
   class BaseTestCase(unittest.TestCase):
       """Base test case with common setup and teardown."""
       
       def setUp(self):
           """Set up the test case."""
           setup_test_environment()
       
       def tearDown(self):
           """Tear down the test case."""
           teardown_test_environment()
       
       def create_mock_supabase(self, response_data=None, response_error=None):
           """Create a mock Supabase client."""
           return create_mock_supabase_client(response_data, response_error)
   ```

3. Create shared test fixtures:
   ```python
   # tests/fixtures/__init__.py
   
   """
   Test fixtures for CryoProtect v2 tests.
   
   This module provides common test fixtures that can be used
   across multiple test modules.
   """
   
   import os
   import json
   import pytest
   from pathlib import Path
   
   from api.utils import get_supabase_client
   from tests.utils import create_mock_supabase_client, MockSupabaseResponse
   
   
   @pytest.fixture
   def test_data_dir():
       """Return the path to the test data directory."""
       return Path(__file__).parent.parent / "data"
   
   
   @pytest.fixture
   def molecule_data(test_data_dir):
       """Load test molecule data."""
       data_file = test_data_dir / "molecules.json"
       if data_file.exists():
           with open(data_file, "r") as f:
               return json.load(f)
       return []
   
   
   @pytest.fixture
   def mixture_data(test_data_dir):
       """Load test mixture data."""
       data_file = test_data_dir / "mixtures.json"
       if data_file.exists():
           with open(data_file, "r") as f:
               return json.load(f)
       return []
   
   
   @pytest.fixture
   def mock_supabase(monkeypatch):
       """Provide a mock Supabase client for testing."""
       mock_client = create_mock_supabase_client()
       
       # Mock the get_supabase_client function
       def _mock_get_client():
           return mock_client
       
       monkeypatch.setattr("api.utils.get_supabase_client", _mock_get_client)
       return mock_client
   
   
   @pytest.fixture
   def flask_app():
       """Create a Flask app for testing."""
       from app import create_app
       app = create_app(testing=True)
       app.config['TESTING'] = True
       
       return app
   
   
   @pytest.fixture
   def flask_client(flask_app):
       """Create a Flask test client."""
       with flask_app.test_client() as client:
           yield client
   
   
   @pytest.fixture
   def auth_tokens():
       """Provide authentication tokens for different roles."""
       return {
           "admin": "admin-test-token",
           "user": "user-test-token",
           "service_role": "service-role-test-token",
       }
   ```

4. Create a common test data directory:
   ```bash
   mkdir -p tests/data
   ```

5. Create sample test data:
   ```python
   # tests/data/molecules.json
   
   [
     {
       "id": "test-molecule-1",
       "name": "Test Molecule 1",
       "formula": "C2H6O",
       "smiles": "CCO",
       "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
       "molecular_weight": 46.07
     },
     {
       "id": "test-molecule-2",
       "name": "Test Molecule 2",
       "formula": "C3H8O3",
       "smiles": "C(C(CO)O)O",
       "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
       "molecular_weight": 92.09
     }
   ]
   
   # tests/data/mixtures.json
   
   [
     {
       "id": "test-mixture-1",
       "name": "Test Mixture 1",
       "description": "A test mixture",
       "components": [
         {
           "molecule_id": "test-molecule-1",
           "concentration": 10.0,
           "concentration_unit": "mM"
         },
         {
           "molecule_id": "test-molecule-2",
           "concentration": 5.0,
           "concentration_unit": "mM"
         }
       ]
     }
   ]
   ```

6. Create a model test module using pytest:
   ```python
   # tests/models/test_molecule.py
   
   """
   Tests for the Molecule model.
   
   This module tests the Molecule model functionality including
   database operations and property calculations.
   """
   
   import pytest
   from tests.fixtures import mock_supabase, molecule_data
   
   from api.models import Molecule
   
   
   class TestMolecule:
       """Tests for the Molecule model."""
       
       def test_get_molecule_by_id(self, mock_supabase, molecule_data):
           """Test retrieving a molecule by ID."""
           # Set up the mock response
           mock_supabase.response_data = [molecule_data[0]]
           
           # Call the method
           molecule = Molecule.get(molecule_data[0]["id"])
           
           # Assertions
           assert molecule is not None
           assert molecule["id"] == molecule_data[0]["id"]
           assert molecule["name"] == molecule_data[0]["name"]
       
       def test_get_molecule_not_found(self, mock_supabase):
           """Test retrieving a molecule that doesn't exist."""
           # Set up the mock response
           mock_supabase.response_data = []
           
           # Call the method
           molecule = Molecule.get("nonexistent-id")
           
           # Assertions
           assert molecule is None
       
       def test_get_all_molecules(self, mock_supabase, molecule_data):
           """Test retrieving all molecules."""
           # Set up the mock response
           mock_supabase.response_data = molecule_data
           
           # Call the method
           molecules = Molecule.get_all()
           
           # Assertions
           assert molecules is not None
           assert len(molecules) == len(molecule_data)
           assert molecules[0]["id"] == molecule_data[0]["id"]
       
       def test_create_molecule(self, mock_supabase):
           """Test creating a molecule."""
           # Set up the mock response
           mock_supabase.response_data = [{"id": "new-molecule-id"}]
           
           # Create new molecule data
           new_molecule = {
               "name": "New Test Molecule",
               "formula": "C4H10O",
               "smiles": "CCCCO",
               "inchikey": "LBFFZMHNMNIGLB-UHFFFAOYSA-N",
               "molecular_weight": 74.12
           }
           
           # Call the method
           result = Molecule.create(new_molecule)
           
           # Assertions
           assert result is not None
           assert "id" in result
   ```

7. Create an API test module:
   ```python
   # tests/api/test_molecule_endpoints.py
   
   """
   Tests for the molecule-related API endpoints.
   
   This module tests the API endpoints related to molecules,
   including listing, retrieval, creation, and property calculation.
   """
   
   import pytest
   import json
   from tests.fixtures import flask_client, mock_supabase, molecule_data, auth_tokens
   
   
   class TestMoleculeEndpoints:
       """Tests for molecule-related API endpoints."""
       
       def test_get_molecules(self, flask_client, mock_supabase, molecule_data):
           """Test GET /api/v1/molecules endpoint."""
           # Set up the mock response
           mock_supabase.response_data = molecule_data
           
           # Make the request
           response = flask_client.get("/api/v1/molecules")
           
           # Assertions
           assert response.status_code == 200
           data = json.loads(response.data)
           assert len(data) == len(molecule_data)
       
       def test_get_molecule_by_id(self, flask_client, mock_supabase, molecule_data):
           """Test GET /api/v1/molecules/{id} endpoint."""
           # Set up the mock response
           mock_supabase.response_data = [molecule_data[0]]
           
           # Make the request
           response = flask_client.get(f"/api/v1/molecules/{molecule_data[0]['id']}")
           
           # Assertions
           assert response.status_code == 200
           data = json.loads(response.data)
           assert data["id"] == molecule_data[0]["id"]
           assert data["name"] == molecule_data[0]["name"]
       
       def test_get_molecule_not_found(self, flask_client, mock_supabase):
           """Test GET /api/v1/molecules/{id} for nonexistent molecule."""
           # Set up the mock response
           mock_supabase.response_data = []
           
           # Make the request
           response = flask_client.get("/api/v1/molecules/nonexistent-id")
           
           # Assertions
           assert response.status_code == 404
       
       def test_create_molecule(self, flask_client, mock_supabase, auth_tokens):
           """Test POST /api/v1/molecules endpoint."""
           # Set up the mock response
           mock_supabase.response_data = [{"id": "new-molecule-id"}]
           
           # Create new molecule data
           new_molecule = {
               "name": "New Test Molecule",
               "formula": "C4H10O",
               "smiles": "CCCCO",
               "inchikey": "LBFFZMHNMNIGLB-UHFFFAOYSA-N",
               "molecular_weight": 74.12
           }
           
           # Make the request
           response = flask_client.post(
               "/api/v1/molecules",
               data=json.dumps(new_molecule),
               content_type="application/json",
               headers={"Authorization": f"Bearer {auth_tokens['user']}"}
           )
           
           # Assertions
           assert response.status_code == 201
           data = json.loads(response.data)
           assert "id" in data
       
       def test_calculate_properties(self, flask_client, mock_supabase, molecule_data):
           """Test POST /api/v1/molecules/{id}/calculate-properties endpoint."""
           # Set up the mock response
           mock_supabase.response_data = [molecule_data[0]]
           
           # Make the request
           response = flask_client.post(
               f"/api/v1/molecules/{molecule_data[0]['id']}/calculate-properties"
           )
           
           # Assertions
           assert response.status_code == 200
           data = json.loads(response.data)
           assert "message" in data
   ```

8. Create a main test runner:
   ```python
   # tests/run_tests.py
   
   #!/usr/bin/env python3
   """
   Test runner for CryoProtect v2 tests.
   
   This script provides a unified interface for running tests
   using both unittest and pytest frameworks.
   """
   
   import os
   import sys
   import argparse
   import pytest
   import unittest
   import json
   from datetime import datetime
   from pathlib import Path
   
   # Add the project root to the path
   project_root = Path(__file__).parent.parent
   if str(project_root) not in sys.path:
       sys.path.insert(0, str(project_root))
   
   
   def run_tests(args):
       """
       Run tests based on command-line arguments.
       
       Args:
           args: Command-line arguments
           
       Returns:
           int: Exit code (0 for success, non-zero for failure)
       """
       # Set up environment variables for testing
       os.environ['FLASK_ENV'] = 'testing'
       os.environ['SUPABASE_URL'] = 'https://test.supabase.co'
       os.environ['SUPABASE_KEY'] = 'test-key'
       
       # Prepare report directory
       if args.report:
           report_dir = Path("reports/tests")
           report_dir.mkdir(parents=True, exist_ok=True)
       
       # Choose the test framework
       if args.framework == 'pytest' or (not args.framework and not args.unittest):
           # Run pytest tests
           pytest_args = []
           
           # Add verbosity
           if args.verbose:
               pytest_args.append('-v')
           
           # Add pattern if specified
           if args.pattern:
               pytest_args.append(args.pattern)
           else:
               pytest_args.append('tests/')
           
           # Add test coverage if requested
           if args.coverage:
               pytest_args.extend(['--cov=.', '--cov-report=term', '--cov-report=html:reports/coverage_html'])
           
           # Run the tests
           result = pytest.main(pytest_args)
           
           # Generate report if requested
           if args.report:
               # Create simple report - pytest has its own reporting mechanisms
               report = {
                   'timestamp': datetime.now().isoformat(),
                   'framework': 'pytest',
                   'exit_code': result,
                   'command': ' '.join(['pytest'] + pytest_args)
               }
               
               report_path = report_dir / f"pytest_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
               with open(report_path, 'w') as f:
                   json.dump(report, f, indent=2)
           
           return result
           
       elif args.framework == 'unittest' or args.unittest:
           # Run unittest tests
           loader = unittest.TestLoader()
           
           # Discover tests
           if args.pattern:
               suite = loader.discover('.', pattern=args.pattern)
           else:
               suite = loader.discover('tests')
           
           # Run the tests
           runner = unittest.TextTestRunner(verbosity=2 if args.verbose else 1)
           result = runner.run(suite)
           
           # Generate report if requested
           if args.report:
               report = {
                   'timestamp': datetime.now().isoformat(),
                   'framework': 'unittest',
                   'tests_run': result.testsRun,
                   'failures': len(result.failures),
                   'errors': len(result.errors),
                   'skipped': len(result.skipped) if hasattr(result, 'skipped') else 0,
                   'was_successful': result.wasSuccessful()
               }
               
               report_path = report_dir / f"unittest_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
               with open(report_path, 'w') as f:
                   json.dump(report, f, indent=2)
           
           return 0 if result.wasSuccessful() else 1
       
       return 0
   
   
   def main():
       """Main function for the test runner."""
       parser = argparse.ArgumentParser(description='Run CryoProtect v2 tests')
       parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
       parser.add_argument('-p', '--pattern', help='Test pattern (e.g., test_*.py)')
       parser.add_argument('-f', '--framework', choices=['pytest', 'unittest'], help='Test framework to use')
       parser.add_argument('-u', '--unittest', action='store_true', help='Use unittest framework')
       parser.add_argument('-c', '--coverage', action='store_true', help='Generate test coverage report')
       parser.add_argument('-r', '--report', action='store_true', help='Generate test report')
       parser.add_argument('-m', '--module', help='Run tests for specific module')
       
       args = parser.parse_args()
       
       # If module is specified, update pattern
       if args.module:
           if args.module == 'api':
               args.pattern = 'test_*api*.py'
           elif args.module == 'models':
               args.pattern = 'test_*model*.py'
           elif args.module == 'rdkit':
               args.pattern = 'test_*rdkit*.py'
           elif args.module == 'database':
               args.pattern = 'test_*database*.py'
           else:
               args.pattern = f'test_*{args.module}*.py'
       
       return run_tests(args)
   
   
   if __name__ == '__main__':
       sys.exit(main())
   ```

9. Create a conftest.py file for pytest:
   ```python
   # tests/conftest.py
   
   """
   Pytest configuration file for CryoProtect v2 tests.
   
   This file imports fixtures and sets up the pytest environment.
   """
   
   import os
   import sys
   from pathlib import Path
   
   # Make sure the project root is in the path
   project_root = Path(__file__).parent.parent
   if str(project_root) not in sys.path:
       sys.path.insert(0, str(project_root))
   
   # Import fixtures to make them available to all tests
   from tests.fixtures import (
       test_data_dir,
       molecule_data,
       mixture_data,
       mock_supabase,
       flask_app,
       flask_client,
       auth_tokens
   )
   
   
   def pytest_configure(config):
       """Configure pytest environment."""
       # Set environment variables for testing
       os.environ['FLASK_ENV'] = 'testing'
       os.environ['SUPABASE_URL'] = 'https://test.supabase.co'
       os.environ['SUPABASE_KEY'] = 'test-key'
   
   
   def pytest_unconfigure(config):
       """Clean up pytest environment."""
       # Remove test environment variables
       for var in ['FLASK_ENV', 'SUPABASE_URL', 'SUPABASE_KEY']:
           if var in os.environ:
               del os.environ[var]
   ```

10. Create a test_runner.py script at the project root:
    ```python
    # test_runner.py
    
    #!/usr/bin/env python3
    """
    Main test runner for CryoProtect v2 project.
    
    This script is a convenience wrapper around the tests/run_tests.py script.
    """
    
    import sys
    from pathlib import Path
    
    # Add the project root to the path
    project_root = Path(__file__).parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    # Import and run the test runner
    from tests.run_tests import main
    
    if __name__ == '__main__':
        sys.exit(main())
    ```

11. Make the test runner executable:
    ```bash
    chmod +x test_runner.py
    ```

12. Update documentation for the new testing framework:
    ```markdown
    # Testing Framework
    
    The CryoProtect v2 project uses a unified testing framework that supports both
    unittest and pytest. The framework provides shared utilities, fixtures, and
    a consistent approach to testing across the codebase.
    
    ## Running Tests
    
    Run all tests with the test runner:
    
    ```bash
    # Run all tests with pytest (default)
    ./test_runner.py
    
    # Run tests with verbose output
    ./test_runner.py -v
    
    # Run only API tests
    ./test_runner.py -m api
    
    # Run tests with a specific pattern
    ./test_runner.py -p "test_molecule*.py"
    
    # Generate test coverage report
    ./test_runner.py -c
    
    # Generate test report
    ./test_runner.py -r
    
    # Use unittest framework instead of pytest
    ./test_runner.py -u
    ```
    
    ## Test Organization
    
    Tests are organized by domain:
    
    - `tests/api/`: Tests for API endpoints
    - `tests/models/`: Tests for data models
    - `tests/database/`: Tests for database operations
    - `tests/rdkit/`: Tests for RDKit integration
    - `tests/integration/`: End-to-end integration tests
    
    ## Fixtures and Utilities
    
    Common fixtures and utilities are available in:
    
    - `tests/fixtures/__init__.py`: Shared test fixtures
    - `tests/utils/__init__.py`: Shared test utilities
    
    ## Test Data
    
    Test data is stored in the `tests/data/` directory.
    ```

## Files to Modify
- Create new directory structure and files as described above
- The existing test files can be gradually migrated to the new structure

## Verification
1. Run the new test runner:
   ```bash
   ./test_runner.py -v
   ```

2. Generate a test coverage report:
   ```bash
   ./test_runner.py -c
   ```

3. Try running specific test modules:
   ```bash
   ./test_runner.py -m api
   ```

4. Compare the test results with the original tests to ensure functionality is preserved.

## Notes for Roo Code Agent
- This implementation creates a flexible testing framework that supports both unittest and pytest
- The modular design allows for easy extension with new test types
- Shared fixtures and utilities reduce code duplication
- The approach makes it easy to gradually migrate existing tests to the new framework
- Full coverage should be maintained throughout the migration process
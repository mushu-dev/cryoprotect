# Task 2.1: Create Unified Test Directory Structure

## Objective
Reorganize the scattered test files into a unified, hierarchical structure that improves organization, reduces duplication, and makes tests easier to maintain and run.

## Context
The project currently has 65+ test files scattered throughout the codebase with inconsistent patterns and significant duplication. By creating a unified test directory structure, we'll improve test organization, make test discovery more consistent, and establish a foundation for better test coverage reporting.

## Acceptance Criteria
- A well-organized test directory structure that mirrors the project structure
- Consistent test naming patterns for easy discovery
- Common testing utilities and fixtures in a shared location
- Clear separation between different types of tests (unit, integration, etc.)
- Good documentation explaining the test organization
- No test functionality loss during restructuring
- All tests pass after reorganization

## Implementation Steps

1. Create the new test directory structure:
   ```
   tests/
   ├── __init__.py              # Package initialization
   ├── conftest.py              # Shared pytest fixtures
   ├── unit/                    # Unit tests
   │   ├── __init__.py
   │   ├── api/                 # Tests for API functionality
   │   │   ├── __init__.py
   │   │   ├── test_resources.py
   │   │   ├── test_models.py
   │   │   └── test_utils.py
   │   ├── database/            # Tests for database functionality
   │   │   ├── __init__.py
   │   │   ├── test_population.py
   │   │   ├── test_migrations.py
   │   │   └── test_verification.py
   │   └── utils/               # Tests for utility functions
   │       ├── __init__.py
   │       └── test_utils.py
   ├── integration/             # Integration tests
   │   ├── __init__.py
   │   ├── test_api_endpoints.py
   │   ├── test_database_operations.py
   │   └── test_authentication.py
   ├── functional/              # Functional tests
   │   ├── __init__.py
   │   ├── test_workflows.py
   │   └── test_user_scenarios.py
   ├── performance/             # Performance tests
   │   ├── __init__.py
   │   └── test_database_performance.py
   ├── fixtures/                # Test fixtures and test data
   │   ├── __init__.py
   │   ├── api_fixtures.py
   │   ├── database_fixtures.py
   │   └── data/                # Test data files
   │       ├── molecules.json
   │       └── mixtures.json
   └── utils/                   # Testing utilities
       ├── __init__.py
       ├── test_client.py       # Custom test client
       ├── mock_supabase.py     # Mocking utilities
       └── assertions.py        # Custom test assertions
   ```

2. Create a shared pytest configuration file (`conftest.py`):
   ```python
   """
   Shared pytest fixtures and configuration.
   
   This module provides pytest fixtures that can be used across all tests.
   """
   
   import os
   import pytest
   import tempfile
   import json
   from typing import Dict, Any, Generator
   
   from tests.utils.mock_supabase import MockSupabase
   
   @pytest.fixture
   def app_config() -> Dict[str, Any]:
       """
       Provide a test configuration for the application.
       
       Returns:
           Dictionary with test configuration
       """
       return {
           'TESTING': True,
           'DEBUG': True,
           'SUPABASE_URL': 'https://test-instance.supabase.co',
           'SUPABASE_KEY': 'test-key',
           'SUPABASE_SERVICE_ROLE_KEY': 'test-service-role-key',
       }
   
   @pytest.fixture
   def mock_supabase() -> Generator[MockSupabase, None, None]:
       """
       Provide a mock Supabase client for testing.
       
       Yields:
           MockSupabase instance
       """
       mock_client = MockSupabase()
       yield mock_client
       mock_client.reset()
   
   @pytest.fixture
   def temp_file() -> Generator[str, None, None]:
       """
       Create a temporary file for testing.
       
       Yields:
           Path to the temporary file
       """
       fd, path = tempfile.mkstemp()
       yield path
       os.close(fd)
       os.unlink(path)
       
   @pytest.fixture
   def sample_molecules_data() -> Dict:
       """
       Provide sample molecule data for testing.
       
       Returns:
           Dictionary with sample molecule data
       """
       data_path = os.path.join(
           os.path.dirname(__file__), 
           'fixtures', 
           'data', 
           'molecules.json'
       )
       with open(data_path, 'r') as f:
           return json.load(f)
       
   @pytest.fixture
   def sample_mixtures_data() -> Dict:
       """
       Provide sample mixture data for testing.
       
       Returns:
           Dictionary with sample mixture data
       """
       data_path = os.path.join(
           os.path.dirname(__file__), 
           'fixtures', 
           'data', 
           'mixtures.json'
       )
       with open(data_path, 'r') as f:
           return json.load(f)
   ```

3. Create a mock Supabase client for testing:
   ```python
   """
   Mock Supabase client for testing.
   
   This module provides a mock implementation of the Supabase client
   that can be used for testing without connecting to a real Supabase instance.
   """
   
   from typing import Dict, List, Any, Optional, Union, Callable
   
   class MockQuery:
       """Mock Supabase query builder."""
       
       def __init__(self, table: str, data: List[Dict]):
           self.table = table
           self.data = data
           self.filters = []
           self.selected_columns = None
           
       def select(self, columns: str = '*') -> 'MockQuery':
           """Mock select method."""
           self.selected_columns = columns
           return self
           
       def eq(self, column: str, value: Any) -> 'MockQuery':
           """Mock equality filter."""
           self.filters.append(('eq', column, value))
           return self
           
       def execute(self) -> Dict:
           """Mock query execution."""
           # Apply filters
           filtered_data = self.data
           for filter_type, column, value in self.filters:
               if filter_type == 'eq':
                   filtered_data = [
                       item for item in filtered_data
                       if item.get(column) == value
                   ]
                   
           return {'data': filtered_data, 'error': None}
           
       def insert(self, record: Dict) -> 'MockQuery':
           """Mock insert method."""
           # Add an ID if not present
           if 'id' not in record:
               record['id'] = f"mock-id-{len(self.data) + 1}"
               
           self.data.append(record)
           return self
           
       def update(self, record: Dict) -> 'MockQuery':
           """Mock update method."""
           for i, item in enumerate(self.data):
               match = True
               for filter_type, column, value in self.filters:
                   if filter_type == 'eq' and item.get(column) != value:
                       match = False
                       break
                       
               if match:
                   self.data[i] = {**item, **record}
                   
           return self
           
       def delete(self) -> 'MockQuery':
           """Mock delete method."""
           original_data = self.data.copy()
           for filter_type, column, value in self.filters:
               if filter_type == 'eq':
                   self.data = [
                       item for item in self.data
                       if item.get(column) != value
                   ]
                   
           return self
   
   class MockSupabase:
       """Mock Supabase client."""
       
       def __init__(self):
           self.tables = {}
           self.rpc_handlers = {}
           self.reset()
           
       def reset(self):
           """Reset mock state."""
           self.tables = {
               'molecules': [],
               'mixtures': [],
               'mixture_components': [],
               'migrations': []
           }
           
           # Add default RPC handlers
           self.rpc_handlers = {
               'has_table': lambda params: {'data': [params.get('table_name') in self.tables]}
           }
           
       def table(self, name: str) -> MockQuery:
           """Get a query builder for a table."""
           if name not in self.tables:
               self.tables[name] = []
               
           return MockQuery(name, self.tables[name])
           
       def rpc(self, function: str, params: Optional[Dict] = None) -> MockQuery:
           """Mock RPC call."""
           if function in self.rpc_handlers:
               result = self.rpc_handlers[function](params or {})
               return MockQuery('rpc', result.get('data', []))
           else:
               return MockQuery('rpc', [{'success': True}])
               
       def sql(self, query: str, params: Optional[Dict] = None) -> MockQuery:
           """Mock SQL query."""
           # Simple query parser for testing
           query = query.lower()
           
           if query.startswith('select'):
               # Extract table name (very basic parser)
               parts = query.split('from')
               if len(parts) > 1:
                   table_parts = parts[1].strip().split()
                   if table_parts:
                       table_name = table_parts[0].strip()
                       if table_name in self.tables:
                           return MockQuery(table_name, self.tables[table_name])
           
           # Default empty result for other queries
           return MockQuery('sql', [{'success': True}])
           
       def add_test_data(self, table: str, data: List[Dict]):
           """Add test data to a table."""
           if table not in self.tables:
               self.tables[table] = []
               
           self.tables[table].extend(data)
           
       def register_rpc_handler(self, function: str, handler: Callable):
           """Register a custom RPC handler."""
           self.rpc_handlers[function] = handler
   ```

4. Create a testing utilities module:
   ```python
   """
   Testing utilities.
   
   This module provides utility functions for writing tests.
   """
   
   import json
   from typing import Dict, List, Any, Optional, Union

   def load_test_data(file_path: str) -> Any:
       """
       Load test data from a JSON file.
       
       Args:
           file_path: Path to the JSON file
           
       Returns:
           Loaded data
       """
       with open(file_path, 'r') as f:
           return json.load(f)
           
   def generate_test_id(prefix: str = 'test') -> str:
       """
       Generate a unique test ID.
       
       Args:
           prefix: ID prefix
           
       Returns:
           Unique test ID
       """
       import uuid
       return f"{prefix}-{uuid.uuid4()}"
       
   def flatten_response_data(response: Dict) -> List[Dict]:
       """
       Flatten a nested response into a list of dictionaries.
       
       Args:
           response: Response dictionary
           
       Returns:
           Flattened list of dictionaries
       """
       if 'data' in response and isinstance(response['data'], list):
           return response['data']
       elif 'results' in response and isinstance(response['results'], list):
           return response['results']
       elif isinstance(response, list):
           return response
       else:
           return []
   ```

5. Create a README for the test directory:
   ```markdown
   # Test Structure

   This directory contains tests for the CryoProtect v2 project. The tests are organized
   into categories based on their scope and purpose.

   ## Test Categories

   - **Unit Tests** (`unit/`): Tests for individual functions and classes in isolation
   - **Integration Tests** (`integration/`): Tests for interactions between components
   - **Functional Tests** (`functional/`): Tests for end-to-end functionality
   - **Performance Tests** (`performance/`): Tests for performance and scalability

   ## Fixtures and Utilities

   - **Fixtures** (`fixtures/`): Shared test fixtures and test data
   - **Utilities** (`utils/`): Helper functions and classes for testing

   ## Running Tests

   ### Running All Tests

   ```bash
   pytest tests/
   ```

   ### Running Specific Test Categories

   ```bash
   # Run all unit tests
   pytest tests/unit/

   # Run all API unit tests
   pytest tests/unit/api/

   # Run all integration tests
   pytest tests/integration/
   ```

   ### Running with Coverage

   ```bash
   pytest --cov=database tests/unit/database/
   ```

   ## Writing Tests

   ### Naming Convention

   - Test modules should be named `test_*.py`
   - Test classes should be named `Test*`
   - Test methods should be named `test_*`

   ### Using Fixtures

   ```python
   def test_database_connection(mock_supabase):
       # Use the mock_supabase fixture
       result = mock_supabase.table('molecules').select().execute()
       assert 'data' in result
   ```

   ### Adding Test Data

   Add test data to the `fixtures/data/` directory in JSON format.
   ```

6. Create example test files for each category:

   a. Unit Test Example:
   ```python
   """
   Unit tests for database population module.
   """
   
   import pytest
   from database.population.molecules import populate_molecules
   
   def test_populate_molecules(mock_supabase, sample_molecules_data):
       """Test populating molecules table."""
       # Arrange
       # Add test data to mock
       
       # Act
       result = populate_molecules(mock_supabase)
       
       # Assert
       assert result > 0
       assert len(mock_supabase.tables['molecules']) == result
   ```

   b. Integration Test Example:
   ```python
   """
   Integration tests for API endpoints.
   """
   
   import pytest
   from app import create_app
   
   def test_molecules_endpoint(app_config):
       """Test molecules endpoint."""
       # Arrange
       app = create_app(app_config)
       client = app.test_client()
       
       # Act
       response = client.get('/api/molecules')
       data = response.get_json()
       
       # Assert
       assert response.status_code == 200
       assert 'data' in data
   ```

7. Create a pytest configuration file (`pytest.ini`):
   ```ini
   [pytest]
   testpaths = tests
   python_files = test_*.py
   python_classes = Test*
   python_functions = test_*
   markers =
       unit: Unit tests
       integration: Integration tests
       functional: Functional tests
       performance: Performance tests
       slow: Slow tests
   ```

8. Create a sample migration script for testing:
   ```python
   """
   Sample migration script for testing.
   
   This migration creates a test table.
   """
   
   def apply(conn, environment):
       """Apply the migration."""
       conn.sql("""
           CREATE TABLE IF NOT EXISTS test_table (
               id SERIAL PRIMARY KEY,
               name TEXT NOT NULL,
               value INTEGER
           )
       """).execute()
       
   def rollback(conn, environment):
       """Roll back the migration."""
       conn.sql("DROP TABLE IF EXISTS test_table").execute()
   ```

## Files to Modify

- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/conftest.py` 
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/README.md`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/pytest.ini`
- Create new subdirectory structure as outlined above
- Create example test files for each test category
- Create utility files for testing

## Verification
1. Verify that the directory structure is created correctly
2. Check that the shared fixtures are defined properly
3. Verify that the example tests can be discovered and run
4. Ensure test utilities work as expected
5. Check that test data can be loaded properly
6. Review the documentation for clarity and completeness
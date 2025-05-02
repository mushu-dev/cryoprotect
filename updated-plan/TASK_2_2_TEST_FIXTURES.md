# Task 2.2: Implement Shared Test Fixtures

## Objective
Create a comprehensive set of shared test fixtures that standardize test setup, reduce duplication, and make tests more maintainable across the entire project.

## Context
With our new test directory structure in place, we need reusable test fixtures to avoid duplicating setup code across test files. Currently, similar setup code appears in multiple test files, making maintenance difficult and reducing consistency. By implementing shared fixtures, we'll improve test reliability and reduce the effort needed to write and maintain tests.

## Acceptance Criteria
- Comprehensive fixtures for database, API, and authentication testing
- Fixtures for mocking external dependencies (Supabase, RDKit)
- Data fixtures for common test scenarios
- Clear documentation for all fixtures
- Integration with pytest's fixture system
- Tests demonstrating fixture usage
- No duplication in test setup code

## Implementation Steps

1. Create the fixtures directory structure:
   ```
   tests/fixtures/
   ├── __init__.py              # Package initialization
   ├── database/                # Database fixtures
   │   ├── __init__.py
   │   ├── connection.py        # Database connection fixtures
   │   ├── data.py              # Test data generation
   │   └── migrations.py        # Migration testing fixtures
   ├── api/                     # API fixtures
   │   ├── __init__.py
   │   ├── client.py            # API client fixtures
   │   └── auth.py              # Authentication fixtures
   ├── mocks/                   # Mock fixtures
   │   ├── __init__.py
   │   ├── supabase.py          # Supabase mocks
   │   └── rdkit.py             # RDKit mocks
   └── data/                    # Test data files
       ├── molecules.json
       ├── mixtures.json
       └── users.json
   ```

2. Create database connection fixtures:
   ```python
   """
   Database connection fixtures.
   
   This module provides fixtures for setting up database
   connections in tests.
   """
   
   import os
   import pytest
   from typing import Dict, Any, Generator, List
   
   from database.utils.connection import create_connection
   from tests.mocks.supabase import MockSupabase, patch_supabase_client
   
   @pytest.fixture
   def db_config() -> Dict[str, Any]:
       """
       Provide a database configuration for testing.
       
       Returns:
           Dictionary with test database configuration
       """
       return {
           'url': 'https://test-instance.supabase.co',
           'key': 'test-key',
           'service_role': 'test-service-role-key',
       }
   
   @pytest.fixture
   def mock_db() -> Generator[MockSupabase, None, None]:
       """
       Provide a mock database connection.
       
       Yields:
           MockSupabase instance
       """
       mock_client = MockSupabase()
       
       # Setup test data
       mock_client.add_test_data('molecules', [
           {'id': 'mol-1', 'name': 'Test Molecule 1', 'formula': 'C6H12O6'},
           {'id': 'mol-2', 'name': 'Test Molecule 2', 'formula': 'H2O'},
       ])
       
       mock_client.add_test_data('mixtures', [
           {'id': 'mix-1', 'name': 'Test Mixture 1'},
           {'id': 'mix-2', 'name': 'Test Mixture 2'},
       ])
       
       with patch_supabase_client(mock_client):
           yield mock_client
           
   @pytest.fixture
   def real_db_connection(db_config) -> Generator[Any, None, None]:
       """
       Provide a real database connection for integration tests.
       
       This fixture should only be used in integration tests.
       
       Yields:
           Database connection
       """
       # Only use this fixture when explicitly requested with a marker
       if not pytest.config.getoption('--use-real-db'):
           pytest.skip('Skipping test that requires a real database connection')
           
       conn = create_connection(db_config)
       yield conn
   
   @pytest.fixture
   def populated_mock_db(mock_db) -> Generator[MockSupabase, None, None]:
       """
       Provide a mock database populated with comprehensive test data.
       
       Yields:
           MockSupabase instance with test data
       """
       # Load and insert test data from JSON files
       from tests.fixtures.data import load_test_data
       
       molecules_data = load_test_data('molecules.json')
       mixtures_data = load_test_data('mixtures.json')
       users_data = load_test_data('users.json')
       
       mock_db.add_test_data('molecules', molecules_data)
       mock_db.add_test_data('mixtures', mixtures_data)
       mock_db.add_test_data('users', users_data)
       
       # Add relationships
       for i, mixture in enumerate(mixtures_data):
           mix_id = mixture['id']
           # Add components for each mixture
           for j in range(1, 3):  # 2 components per mixture
               mol_id = molecules_data[j % len(molecules_data)]['id']
               component = {
                   'id': f'comp-{i}-{j}',
                   'mixture_id': mix_id,
                   'molecule_id': mol_id,
                   'concentration': 10.0 * j,
                   'units': '%'
               }
               mock_db.add_test_data('mixture_components', [component])
       
       yield mock_db
   ```

3. Create data generation utilities:
   ```python
   """
   Test data generation utilities.
   
   This module provides functions for generating and loading
   test data for database tests.
   """
   
   import os
   import json
   import uuid
   from typing import Dict, List, Any, Optional
   
   def load_test_data(filename: str) -> List[Dict[str, Any]]:
       """
       Load test data from a JSON file in the data directory.
       
       Args:
           filename: Name of the JSON file
           
       Returns:
           List of dictionaries with test data
       """
       data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
       file_path = os.path.join(data_dir, filename)
       
       with open(file_path, 'r') as f:
           return json.load(f)
           
   def generate_molecule(
       name: Optional[str] = None,
       formula: Optional[str] = None,
       smiles: Optional[str] = None
   ) -> Dict[str, Any]:
       """
       Generate a test molecule.
       
       Args:
           name: Optional molecule name
           formula: Optional molecular formula
           smiles: Optional SMILES string
           
       Returns:
           Dictionary with molecule data
       """
       id_value = str(uuid.uuid4())
       return {
           'id': id_value,
           'name': name or f"Test Molecule {id_value[:6]}",
           'formula': formula or 'C6H12O6',
           'smiles': smiles or 'CCO',
       }
       
   def generate_mixture(
       name: Optional[str] = None,
       description: Optional[str] = None,
       component_count: int = 2
   ) -> Dict[str, Any]:
       """
       Generate a test mixture with components.
       
       Args:
           name: Optional mixture name
           description: Optional description
           component_count: Number of components to generate
           
       Returns:
           Dictionary with mixture data and components
       """
       id_value = str(uuid.uuid4())
       mixture = {
           'id': id_value,
           'name': name or f"Test Mixture {id_value[:6]}",
           'description': description or f"Description for test mixture",
           'components': []
       }
       
       # Generate components
       for i in range(component_count):
           molecule = generate_molecule()
           component = {
               'id': str(uuid.uuid4()),
               'mixture_id': mixture['id'],
               'molecule_id': molecule['id'],
               'concentration': 10.0 * (i + 1),
               'units': '%',
               'molecule': molecule
           }
           mixture['components'].append(component)
           
       return mixture
   ```

4. Create migration testing fixtures:
   ```python
   """
   Migration testing fixtures.
   
   This module provides fixtures for testing database migrations.
   """
   
   import os
   import pytest
   import tempfile
   from typing import Dict, Any, Generator, List
   
   from database.migrations import (
       initialize_migration_tracking,
       apply_migrations,
       rollback_migrations
   )
   
   @pytest.fixture
   def migration_db(mock_db) -> Generator[Any, None, None]:
       """
       Provide a mock database prepared for migration testing.
       
       Yields:
           MockSupabase instance set up for migrations
       """
       # Initialize migration tracking
       initialize_migration_tracking(mock_db)
       
       # Register additional RPC handlers for migration testing
       def has_column_handler(params):
           table = params.get('table_name', '')
           column = params.get('column_name', '')
           
           if table in mock_db.tables:
               # Simplified mock - assume columns exist in test tables
               return {'data': [True]}
           return {'data': [False]}
           
       mock_db.register_rpc_handler('has_column', has_column_handler)
       
       yield mock_db
       
   @pytest.fixture
   def temp_migration_script() -> Generator[str, None, None]:
       """
       Create a temporary migration script for testing.
       
       Yields:
           Path to the temporary migration script
       """
       # Create a temporary directory
       with tempfile.TemporaryDirectory() as temp_dir:
           # Create a migration script
           script_path = os.path.join(temp_dir, '999_test_migration.py')
           with open(script_path, 'w') as f:
               f.write('''
def apply(conn, environment):
    """Apply the migration."""
    conn.sql("""
        CREATE TABLE IF NOT EXISTS test_migration_table (
            id SERIAL PRIMARY KEY,
            name TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT NOW()
        )
    """).execute()
    
def rollback(conn, environment):
    """Roll back the migration."""
    conn.sql("DROP TABLE IF EXISTS test_migration_table").execute()
''')
           yield script_path
   ```

5. Create API client fixtures:
   ```python
   """
   API client fixtures.
   
   This module provides fixtures for testing API endpoints.
   """
   
   import os
   import pytest
   from typing import Dict, Any, Generator
   
   from flask import Flask
   from flask.testing import FlaskClient
   
   @pytest.fixture
   def app_config() -> Dict[str, Any]:
       """
       Provide a test configuration for the Flask application.
       
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
   def app(app_config, mock_db) -> Generator[Flask, None, None]:
       """
       Provide a Flask application instance for testing.
       
       Yields:
           Flask application
       """
       # Import here to avoid circular imports
       from app import create_app
       
       # Create the app with test configuration
       test_app = create_app(app_config)
       
       # Set up app context
       with test_app.app_context():
           yield test_app
           
   @pytest.fixture
   def client(app) -> FlaskClient:
       """
       Provide a Flask test client.
       
       Args:
           app: Flask application fixture
           
       Returns:
           Flask test client
       """
       return app.test_client()
       
   @pytest.fixture
   def api_headers() -> Dict[str, str]:
       """
       Provide headers for API requests.
       
       Returns:
           Dictionary with headers
       """
       return {
           'Content-Type': 'application/json',
           'Accept': 'application/json'
       }
   ```

6. Create authentication fixtures:
   ```python
   """
   Authentication fixtures.
   
   This module provides fixtures for testing authentication
   and authorization.
   """
   
   import os
   import pytest
   import jwt
   from datetime import datetime, timedelta
   from typing import Dict, Any, Generator
   
   @pytest.fixture
   def auth_config() -> Dict[str, Any]:
       """
       Provide authentication configuration for testing.
       
       Returns:
           Dictionary with auth configuration
       """
       return {
           'SECRET_KEY': 'test-secret-key',
           'JWT_ALGORITHM': 'HS256',
           'TOKEN_EXPIRY': 3600,  # 1 hour
       }
   
   @pytest.fixture
   def generate_token(auth_config) -> callable:
       """
       Provide a function to generate authentication tokens.
       
       Args:
           auth_config: Auth configuration fixture
           
       Returns:
           Function to generate tokens
       """
       def _generate_token(user_id: str, role: str = 'user', expires_in: int = None) -> str:
           """Generate a JWT token for testing."""
           now = datetime.utcnow()
           expires = now + timedelta(seconds=expires_in or auth_config['TOKEN_EXPIRY'])
           
           payload = {
               'sub': user_id,
               'role': role,
               'iat': now,
               'exp': expires
           }
           
           return jwt.encode(
               payload,
               auth_config['SECRET_KEY'],
               algorithm=auth_config['JWT_ALGORITHM']
           )
           
       return _generate_token
       
   @pytest.fixture
   def admin_token(generate_token) -> str:
       """
       Provide an admin authentication token.
       
       Args:
           generate_token: Token generator fixture
           
       Returns:
           Admin JWT token
       """
       return generate_token('admin-user-id', role='admin')
       
   @pytest.fixture
   def user_token(generate_token) -> str:
       """
       Provide a regular user authentication token.
       
       Args:
           generate_token: Token generator fixture
           
       Returns:
           User JWT token
       """
       return generate_token('regular-user-id', role='user')
       
   @pytest.fixture
   def auth_client(client, user_token, api_headers) -> Generator[Any, None, None]:
       """
       Provide an authenticated test client.
       
       Args:
           client: Flask test client fixture
           user_token: User token fixture
           api_headers: API headers fixture
           
       Yields:
           Flask test client with auth headers
       """
       # Update headers with auth token
       auth_headers = {**api_headers, 'Authorization': f'Bearer {user_token}'}
       
       # Create a wrapper around the client that adds auth headers
       class AuthClient:
           def __init__(self, client, headers):
               self.client = client
               self.headers = headers
               
           def get(self, path, **kwargs):
               kwargs.setdefault('headers', self.headers)
               return self.client.get(path, **kwargs)
               
           def post(self, path, **kwargs):
               kwargs.setdefault('headers', self.headers)
               return self.client.post(path, **kwargs)
               
           def put(self, path, **kwargs):
               kwargs.setdefault('headers', self.headers)
               return self.client.put(path, **kwargs)
               
           def delete(self, path, **kwargs):
               kwargs.setdefault('headers', self.headers)
               return self.client.delete(path, **kwargs)
               
       yield AuthClient(client, auth_headers)
   ```

7. Create Supabase mock utilities:
   ```python
   """
   Supabase mock utilities.
   
   This module provides tools for mocking Supabase in tests.
   """
   
   import pytest
   from unittest.mock import patch
   from contextlib import contextmanager
   from typing import Dict, List, Any, Optional, Union, Generator, Callable
   
   class MockQueryResult:
       """Mock Supabase query result."""
       
       def __init__(self, data: List[Dict] = None, error: Any = None):
           self.data = data or []
           self.error = error
           
   class MockQuery:
       """Mock Supabase query builder."""
       
       def __init__(self, table: str, data: List[Dict]):
           self.table = table
           self.data = data
           self.filters = []
           self.selected_columns = None
           self._order_by = None
           self._limit = None
           
       def select(self, columns: str = '*') -> 'MockQuery':
           """Mock select method."""
           self.selected_columns = columns
           return self
           
       def eq(self, column: str, value: Any) -> 'MockQuery':
           """Mock equality filter."""
           self.filters.append(('eq', column, value))
           return self
           
       def gt(self, column: str, value: Any) -> 'MockQuery':
           """Mock greater than filter."""
           self.filters.append(('gt', column, value))
           return self
           
       def lt(self, column: str, value: Any) -> 'MockQuery':
           """Mock less than filter."""
           self.filters.append(('lt', column, value))
           return self
           
       def like(self, column: str, pattern: str) -> 'MockQuery':
           """Mock LIKE filter."""
           self.filters.append(('like', column, pattern))
           return self
           
       def order(self, column: str, desc: bool = False) -> 'MockQuery':
           """Mock order by."""
           self._order_by = (column, desc)
           return self
           
       def limit(self, count: int) -> 'MockQuery':
           """Mock limit."""
           self._limit = count
           return self
           
       def execute(self) -> MockQueryResult:
           """Mock query execution."""
           # Apply filters
           filtered_data = self.data
           for filter_type, column, value in self.filters:
               if filter_type == 'eq':
                   filtered_data = [
                       item for item in filtered_data
                       if item.get(column) == value
                   ]
               elif filter_type == 'gt':
                   filtered_data = [
                       item for item in filtered_data
                       if item.get(column, 0) > value
                   ]
               elif filter_type == 'lt':
                   filtered_data = [
                       item for item in filtered_data
                       if item.get(column, 0) < value
                   ]
               elif filter_type == 'like':
                   # Simple like implementation - just check if value is in column
                   filtered_data = [
                       item for item in filtered_data
                       if value.replace('%', '') in str(item.get(column, ''))
                   ]
                   
           # Apply order by
           if self._order_by:
               column, desc = self._order_by
               filtered_data = sorted(
                   filtered_data,
                   key=lambda x: x.get(column, ''),
                   reverse=desc
               )
               
           # Apply limit
           if self._limit and self._limit < len(filtered_data):
               filtered_data = filtered_data[:self._limit]
                   
           return MockQueryResult(data=filtered_data)
           
       def insert(self, record: Dict) -> 'MockQuery':
           """Mock insert method."""
           # Add an ID if not present
           if 'id' not in record:
               import uuid
               record['id'] = str(uuid.uuid4())
               
           self.data.append(record)
           return self
           
       def update(self, record: Dict) -> 'MockQuery':
           """Mock update method."""
           updated_data = []
           
           for i, item in enumerate(self.data):
               match = True
               for filter_type, column, value in self.filters:
                   if filter_type == 'eq' and item.get(column) != value:
                       match = False
                       break
                       
               if match:
                   updated_item = {**item, **record}
                   updated_data.append(updated_item)
               else:
                   updated_data.append(item)
                   
           self.data = updated_data
           return self
           
       def delete(self) -> 'MockQuery':
           """Mock delete method."""
           original_data = self.data.copy()
           remaining_data = []
           
           for item in self.data:
               should_delete = True
               for filter_type, column, value in self.filters:
                   if filter_type == 'eq' and item.get(column) != value:
                       should_delete = False
                       break
                       
               if not should_delete:
                   remaining_data.append(item)
                   
           self.data = remaining_data
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
               'migrations': [],
               'users': [],
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
   
   @contextmanager
   def patch_supabase_client(mock_client):
       """
       Context manager to patch the Supabase client.
       
       Args:
           mock_client: Mock Supabase client
           
       Yields:
           The mock client
       """
       # Patch the create_connection function
       target = 'database.utils.connection.create_connection'
       with patch(target, return_value=mock_client):
           # Also patch any direct Supabase imports
           with patch('supabase.create_client', return_value=mock_client):
               yield mock_client
   ```

8. Create RDKit mock utilities:
   ```python
   """
   RDKit mock utilities.
   
   This module provides tools for mocking RDKit in tests.
   """
   
   import pytest
   from unittest.mock import MagicMock, patch
   from contextlib import contextmanager
   from typing import Dict, Any, Generator
   
   class MockMolecule:
       """Mock RDKit molecule."""
       
       def __init__(self, smiles: str = None, inchi: str = None):
           self.smiles = smiles or 'C'
           self.inchi = inchi or 'InChI=1S/CH4/h1H4'
           self.formula = 'C6H12O6'
           self.descriptors = {
               'MolWt': 180.16,
               'LogP': -0.5,
               'NumHDonors': 5,
               'NumHAcceptors': 6,
               'TPSA': 110.38,
           }
           
       def GetNumAtoms(self):
           """Get number of atoms."""
           return 12
           
       def GetNumBonds(self):
           """Get number of bonds."""
           return 12
           
       def GetNumHeavyAtoms(self):
           """Get number of heavy atoms."""
           return 6
   
   class MockRDKit:
       """Mock RDKit library."""
       
       def __init__(self):
           self.Chem = MagicMock()
           self.Descriptors = MagicMock()
           
           # Set up Chem module
           self.Chem.MolFromSmiles = self._mol_from_smiles
           self.Chem.MolFromInchi = self._mol_from_inchi
           self.Chem.MolToSmiles = self._mol_to_smiles
           self.Chem.MolToInchi = self._mol_to_inchi
           
           # Set up Descriptors module
           self.Descriptors.MolWt = self._mol_wt
           self.Descriptors.MolLogP = self._mol_logp
           
       def _mol_from_smiles(self, smiles):
           """Create molecule from SMILES."""
           return MockMolecule(smiles=smiles)
           
       def _mol_from_inchi(self, inchi):
           """Create molecule from InChI."""
           return MockMolecule(inchi=inchi)
           
       def _mol_to_smiles(self, mol):
           """Convert molecule to SMILES."""
           return mol.smiles
           
       def _mol_to_inchi(self, mol):
           """Convert molecule to InChI."""
           return mol.inchi
           
       def _mol_wt(self, mol):
           """Calculate molecular weight."""
           return mol.descriptors['MolWt']
           
       def _mol_logp(self, mol):
           """Calculate LogP."""
           return mol.descriptors['LogP']
           
   @contextmanager
   def patch_rdkit():
       """
       Context manager to patch RDKit.
       
       Yields:
           MockRDKit instance
       """
       mock_rdkit = MockRDKit()
       
       modules_to_patch = [
           'rdkit.Chem',
           'rdkit.Chem.Descriptors',
           'api.rdkit_utils.Chem',
           'api.rdkit_utils.Descriptors',
       ]
       
       # Start all patches
       patches = []
       for module in modules_to_patch:
           if 'Descriptors' in module:
               p = patch(module, mock_rdkit.Descriptors)
           else:
               p = patch(module, mock_rdkit.Chem)
           patches.append(p)
           p.start()
           
       try:
           yield mock_rdkit
       finally:
           # Stop all patches
           for p in patches:
               p.stop()
   ```

9. Create a README for fixtures:
   ```markdown
   # Test Fixtures and Utilities

   This directory contains fixtures and utilities for testing the CryoProtect v2 project.

   ## Fixtures Overview

   ### Database Fixtures
   - `db_config`: Configuration for test database connections
   - `mock_db`: Mock database connection using MockSupabase
   - `real_db_connection`: Real database connection for integration tests
   - `populated_mock_db`: Mock database pre-populated with test data
   - `migration_db`: Database prepared for migration testing
   - `temp_migration_script`: Temporary migration script for testing

   ### API Fixtures
   - `app_config`: Flask application configuration
   - `app`: Flask application instance
   - `client`: Flask test client
   - `api_headers`: Headers for API requests

   ### Authentication Fixtures
   - `auth_config`: Authentication configuration
   - `generate_token`: Function to generate JWT tokens
   - `admin_token`: Admin authentication token
   - `user_token`: Regular user authentication token
   - `auth_client`: Authenticated test client

   ### Mock Utilities
   - `patch_supabase_client`: Context manager to patch Supabase
   - `patch_rdkit`: Context manager to patch RDKit

   ## Usage Examples

   ### Mock Database

   ```python
   def test_molecule_query(mock_db):
       # Add test data
       mock_db.add_test_data('molecules', [
           {'id': '1', 'name': 'Glucose'},
           {'id': '2', 'name': 'Fructose'},
       ])
       
       # Test database operations
       result = mock_db.table('molecules').select().execute()
       assert len(result.data) == 2
   ```

   ### API Testing

   ```python
   def test_molecules_api(client, api_headers):
       response = client.get('/api/molecules', headers=api_headers)
       assert response.status_code == 200
       data = response.get_json()
       assert 'data' in data
   ```

   ### Authentication Testing

   ```python
   def test_protected_endpoint(auth_client):
       response = auth_client.get('/api/protected')
       assert response.status_code == 200
   ```

   ### RDKit Mocking

   ```python
   def test_molecule_processing():
       with patch_rdkit() as mock_rdkit:
           # RDKit functions are now mocked
           mol = mock_rdkit.Chem.MolFromSmiles('C')
           assert mol is not None
           weight = mock_rdkit.Descriptors.MolWt(mol)
           assert weight == 180.16
   ```
   ```

10. Create example data files:

    a. `molecules.json`:
    ```json
    [
      {
        "id": "mol-001",
        "name": "Glucose",
        "formula": "C6H12O6",
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"
      },
      {
        "id": "mol-002",
        "name": "Glycerol",
        "formula": "C3H8O3",
        "smiles": "C(C(CO)O)O"
      },
      {
        "id": "mol-003",
        "name": "Dimethyl Sulfoxide",
        "formula": "C2H6OS",
        "smiles": "CS(=O)C"
      },
      {
        "id": "mol-004",
        "name": "Ethylene Glycol",
        "formula": "C2H6O2",
        "smiles": "OCCO"
      },
      {
        "id": "mol-005",
        "name": "Propylene Glycol",
        "formula": "C3H8O2",
        "smiles": "CC(O)CO"
      }
    ]
    ```

    b. `mixtures.json`:
    ```json
    [
      {
        "id": "mix-001",
        "name": "Basic Cryoprotectant Mix 1",
        "description": "Standard mixture for general use"
      },
      {
        "id": "mix-002",
        "name": "High Concentration Mix",
        "description": "Higher concentration for improved protection"
      },
      {
        "id": "mix-003",
        "name": "Low Temperature Mix",
        "description": "Optimized for extremely low temperatures"
      },
      {
        "id": "mix-004",
        "name": "Cell Culture Mix",
        "description": "Designed for cell culture applications"
      }
    ]
    ```

    c. `users.json`:
    ```json
    [
      {
        "id": "user-001",
        "email": "admin@example.com",
        "role": "admin",
        "name": "Admin User"
      },
      {
        "id": "user-002",
        "email": "scientist@example.com",
        "role": "scientist",
        "name": "Scientist User"
      },
      {
        "id": "user-003",
        "email": "viewer@example.com",
        "role": "viewer",
        "name": "Viewer User"
      }
    ]
    ```

11. Update the main conftest.py to include all fixtures:
    ```python
    """
    Main pytest configuration and fixtures.
    
    This module imports and provides all test fixtures for the project.
    """
    
    # Import fixtures from fixture modules
    from tests.fixtures.database.connection import (
        db_config,
        mock_db,
        real_db_connection,
        populated_mock_db
    )
    
    from tests.fixtures.database.migrations import (
        migration_db,
        temp_migration_script
    )
    
    from tests.fixtures.api.client import (
        app_config,
        app,
        client,
        api_headers
    )
    
    from tests.fixtures.api.auth import (
        auth_config,
        generate_token,
        admin_token,
        user_token,
        auth_client
    )
    
    # Make them available to pytest
    __all__ = [
        # Database fixtures
        'db_config',
        'mock_db',
        'real_db_connection',
        'populated_mock_db',
        'migration_db',
        'temp_migration_script',
        
        # API fixtures
        'app_config',
        'app',
        'client',
        'api_headers',
        
        # Auth fixtures
        'auth_config',
        'generate_token',
        'admin_token',
        'user_token',
        'auth_client',
    ]
    
    # Add pytest command-line options
    def pytest_addoption(parser):
        parser.addoption(
            '--use-real-db',
            action='store_true',
            default=False,
            help='Run tests that require a real database connection'
        )
    ```

12. Create example tests using the fixtures:

    a. Example database test:
    ```python
    """
    Tests for database population module.
    """
    
    import pytest
    from database.population.molecules import populate_molecules
    
    def test_populate_molecules_empty(mock_db):
        """Test populating molecules table with empty database."""
        # Arrange
        assert len(mock_db.tables['molecules']) == 0
        
        # Act
        result = populate_molecules(mock_db)
        
        # Assert
        assert result == 0  # No data provided, so 0 molecules populated
        
    def test_populate_molecules_from_file(mock_db, monkeypatch, tmp_path):
        """Test populating molecules from a file."""
        # Arrange
        import json
        import os
        
        # Create a test data file
        test_data = [
            {'id': 'test-1', 'name': 'Test Molecule 1'},
            {'id': 'test-2', 'name': 'Test Molecule 2'},
        ]
        
        data_file = tmp_path / "test_molecules.json"
        with open(data_file, 'w') as f:
            json.dump(test_data, f)
            
        # Mock environment-specific file paths
        def mock_load_molecule_data(environment, file_path=None):
            if file_path:
                with open(file_path, 'r') as f:
                    return json.load(f)
            return []
            
        from database.population import molecules
        monkeypatch.setattr(molecules, '_load_molecule_data', mock_load_molecule_data)
        
        # Act
        result = populate_molecules(mock_db, file_path=str(data_file))
        
        # Assert
        assert result == 2
        assert len(mock_db.tables['molecules']) == 2
        assert mock_db.tables['molecules'][0]['name'] == 'Test Molecule 1'
        assert mock_db.tables['molecules'][1]['name'] == 'Test Molecule 2'
    
    def test_populate_molecules_with_pre_populated_data(populated_mock_db):
        """Test with pre-populated database."""
        # Arrange - populated_mock_db already has data
        
        # Act
        result = populate_molecules(populated_mock_db)
        
        # Assert - should not duplicate existing data
        assert result == 0
        
    def test_database_population_with_error_handling(mock_db, monkeypatch):
        """Test error handling during population."""
        # Arrange
        def mock_insert_that_fails(*args, **kwargs):
            raise Exception("Test exception")
            
        # Patch the mock to simulate a failure
        original_insert = mock_db.table('molecules').insert
        monkeypatch.setattr(mock_db.table('molecules'), 'insert', mock_insert_that_fails)
        
        # Act
        from database.population.molecules import populate_molecules
        result = populate_molecules(mock_db, environment='test')
        
        # Assert
        assert result == 0  # No molecules inserted due to error
    ```

    b. Example API test:
    ```python
    """
    Tests for molecule API endpoints.
    """
    
    import pytest
    import json
    
    def test_get_molecules(client, mock_db, api_headers):
        """Test getting all molecules."""
        # Arrange
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Test Molecule 1'},
            {'id': 'mol-2', 'name': 'Test Molecule 2'},
        ])
        
        # Act
        response = client.get('/api/molecules', headers=api_headers)
        
        # Assert
        assert response.status_code == 200
        data = json.loads(response.data)
        assert len(data['data']) == 2
        
    def test_get_molecule_by_id(client, mock_db, api_headers):
        """Test getting a specific molecule by ID."""
        # Arrange
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Test Molecule 1'},
        ])
        
        # Act
        response = client.get('/api/molecules/mol-1', headers=api_headers)
        
        # Assert
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['data']['name'] == 'Test Molecule 1'
        
    def test_create_molecule(client, mock_db, api_headers):
        """Test creating a new molecule."""
        # Arrange
        new_molecule = {
            'name': 'New Molecule',
            'formula': 'C6H12O6',
            'smiles': 'C(C1C(C(C(C(O1)O)O)O)O)O'
        }
        
        # Act
        response = client.post(
            '/api/molecules',
            data=json.dumps(new_molecule),
            headers=api_headers
        )
        
        # Assert
        assert response.status_code == 201
        data = json.loads(response.data)
        assert 'id' in data['data']
        assert data['data']['name'] == 'New Molecule'
        
        # Verify it was added to the database
        assert len(mock_db.tables['molecules']) == 1
        
    def test_protected_endpoint(auth_client):
        """Test an endpoint that requires authentication."""
        # Act
        response = auth_client.get('/api/user/profile')
        
        # Assert
        assert response.status_code == 200
        data = json.loads(response.data)
        assert 'user_id' in data
    ```

## Files to Modify

- Create new directories:
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/database/`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/api/`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/mocks/`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/data/`

- Create new files:
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/__init__.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/database/__init__.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/database/connection.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/database/data.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/database/migrations.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/api/__init__.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/api/client.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/api/auth.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/mocks/__init__.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/mocks/supabase.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/mocks/rdkit.py`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/data/molecules.json`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/data/mixtures.json`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/data/users.json`
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/fixtures/README.md`

- Update existing file:
  - `/mnt/c/Users/1edwa/documents/cryoprotect v2/tests/conftest.py`

## Verification
1. Verify that all fixture directories and files are created correctly
2. Check that the fixtures can be imported and used in tests
3. Run example tests to ensure fixtures work properly
4. Verify that mock objects behave as expected
5. Check that test data can be loaded and used
6. Review the documentation for clarity and completeness
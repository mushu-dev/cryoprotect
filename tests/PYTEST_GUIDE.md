# CryoProtect v2 Pytest Guide

This guide provides information on using pytest with the CryoProtect v2 testing framework, including available command-line options, fixtures, and best practices.

## Table of Contents

- [Command-Line Options](#command-line-options)
- [Fixture Categories](#fixture-categories)
  - [Database Fixtures](#database-fixtures)
  - [API Fixtures](#api-fixtures)
  - [Mock Object Fixtures](#mock-object-fixtures)
  - [Test Data Fixtures](#test-data-fixtures)
- [Test Categories and Markers](#test-categories-and-markers)
- [Best Practices](#best-practices)
- [Examples](#examples)

## Command-Line Options

The CryoProtect v2 testing framework provides several custom command-line options for pytest:

| Option | Description |
|--------|-------------|
| `--use-real-db` | Run tests that require a real database connection instead of mocks |
| `--include-slow` | Include tests marked as slow (which are skipped by default) |
| `--api-only` | Run only API tests (tests marked with `@pytest.mark.api`) |
| `--db-only` | Run only database tests (tests marked with `@pytest.mark.database`) |
| `--e2e` | Run end-to-end tests (tests marked with `@pytest.mark.e2e`) |
| `--skip-rdkit` | Skip tests that require RDKit functionality |

### Examples

```bash
# Run all tests
pytest

# Run only database tests
pytest --db-only

# Run tests with real database connection
pytest --use-real-db

# Run API tests and include slow tests
pytest --api-only --include-slow
```

## Fixture Categories

The CryoProtect v2 testing framework provides fixtures in several categories to support different testing needs.

### Database Fixtures

| Fixture | Description |
|---------|-------------|
| `db_config` | Provides test database configuration |
| `mock_db` | Provides a mock database connection |
| `real_db_connection` | Provides a real database connection (requires `--use-real-db`) |
| `populated_mock_db` | Provides a mock database populated with test data |
| `migration_db` | Provides a mock database prepared for migration testing |
| `temp_migration_script` | Creates a temporary migration script for testing |
| `init_test_data` | Initializes a mock database with standard test data |
| `clear_test_data` | Clears all test data from a mock database |

### API Fixtures

| Fixture | Description |
|---------|-------------|
| `app` | Provides a Flask app instance for testing |
| `client` | Provides a Flask test client |
| `api_client` | Provides an unauthenticated API client |
| `authenticated_client` | Provides an authenticated API client with a regular user token |
| `admin_client` | Provides an authenticated API client with an admin token |
| `scientist_client` | Provides an authenticated API client with a scientist token |
| `mock_api_client` | Provides a mock API client for testing without a real API |

#### Authentication Fixtures

| Fixture | Description |
|---------|-------------|
| `auth_token` | Generates a regular user authentication token |
| `admin_token` | Generates an admin user authentication token |
| `user_token` | Generates a regular user authentication token |
| `scientist_token` | Generates a scientist user authentication token |
| `expired_token` | Generates an expired authentication token |
| `invalid_token` | Generates an invalid authentication token |
| `mock_auth_middleware` | Mocks the authentication middleware |
| `mock_current_user` | Sets a mock current user for testing |
| `mock_admin_user` | Sets a mock admin user for testing |
| `mock_regular_user` | Sets a mock regular user for testing |
| `mock_scientist_user` | Sets a mock scientist user for testing |
| `mock_unauthenticated_user` | Sets no current user for testing unauthenticated scenarios |

### Mock Object Fixtures

| Fixture | Description |
|---------|-------------|
| `mock_rdkit` | Provides a mock RDKit implementation |
| `molecule_factory` | Provides a factory function for creating mock molecules |
| `mock_supabase` | Provides a mock Supabase client |

### Test Data Fixtures

| Fixture | Description |
|---------|-------------|
| `load_molecules` | Loads standard test molecules |
| `load_property_types` | Loads standard test property types |
| `load_mixtures` | Loads standard test mixtures |
| `load_experiments` | Loads standard test experiments |
| `load_calculation_methods` | Loads standard test calculation methods |
| `load_teams` | Loads standard test teams |
| `load_projects` | Loads standard test projects |
| `load_related_data` | Loads related data for a specific entity |

#### Data Generator Functions

| Function | Description |
|----------|-------------|
| `generate_molecule` | Generates a test molecule |
| `generate_mixture` | Generates a test mixture |
| `generate_experiment` | Generates a test experiment |
| `generate_prediction` | Generates a test prediction |
| `generate_user_profile` | Generates a test user profile |
| `generate_team` | Generates a test team |
| `generate_project` | Generates a test project |

## Test Categories and Markers

The CryoProtect v2 testing framework uses pytest markers to categorize tests:

| Marker | Description |
|--------|-------------|
| `@pytest.mark.slow` | Marks tests as slow (skipped by default) |
| `@pytest.mark.integration` | Marks tests that integrate multiple components |
| `@pytest.mark.database` | Marks tests that interact with the database |
| `@pytest.mark.api` | Marks tests that test API endpoints |
| `@pytest.mark.rdkit` | Marks tests that require RDKit functionality |
| `@pytest.mark.auth` | Marks tests related to authentication |
| `@pytest.mark.e2e` | Marks end-to-end tests |

Tests are automatically marked based on their location and naming patterns:
- Tests in the `integration/` directory are marked as `integration`
- Tests in the `unit/` directory are marked as `unit`
- Tests with `test_database` in the filename are marked as `database`
- Tests with `test_api` in the filename are marked as `api`
- Tests with `test_rdkit` in the filename are marked as `rdkit`
- Tests with `test_auth` in the filename are marked as `auth`

## Best Practices

### Fixture Isolation and Scope

Fixtures have different scopes that determine how often they are created:

- `scope="function"` (default): The fixture is created for each test function
- `scope="class"`: The fixture is created once per test class
- `scope="module"`: The fixture is created once per test module
- `scope="session"`: The fixture is created once per test session

Choose the appropriate scope for your fixtures to balance test isolation with performance:

```python
@pytest.fixture(scope="module")
def expensive_resource():
    # This fixture is created once per module
    resource = create_expensive_resource()
    yield resource
    cleanup_expensive_resource(resource)
```

### Using Fixtures in Tests

Fixtures are injected into tests by including them as function parameters:

```python
def test_molecule_creation(mock_db, molecule_factory):
    # Create a test molecule
    molecule = molecule_factory(name="Test Molecule")
    
    # Use the mock database to insert the molecule
    mock_db.table("molecules").insert(molecule).execute()
    
    # Query the molecule
    result = mock_db.table("molecules").select("*").eq("name", "Test Molecule").execute()
    
    # Assert the result
    assert len(result.data) == 1
    assert result.data[0]["name"] == "Test Molecule"
```

### Combining Fixtures

Fixtures can depend on other fixtures:

```python
@pytest.fixture
def authenticated_admin_api(api_client, admin_token):
    # Create an authenticated API client with admin token
    api_client.token = admin_token
    return api_client
```

## Examples

### Testing Database Operations

```python
@pytest.mark.database
def test_molecule_insertion(mock_db):
    # Create a test molecule
    molecule = {
        "id": "test-id",
        "name": "Test Molecule",
        "formula": "C6H12O6",
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"
    }
    
    # Insert the molecule
    mock_db.table("molecules").insert(molecule).execute()
    
    # Query the molecule
    result = mock_db.table("molecules").select("*").eq("id", "test-id").execute()
    
    # Assert the result
    assert len(result.data) == 1
    assert result.data[0]["name"] == "Test Molecule"
```

### Testing API Endpoints

```python
@pytest.mark.api
def test_get_molecule_api(authenticated_client):
    # Make a GET request to the API
    response = authenticated_client.get("/api/molecules/test-id")
    
    # Parse the response
    data = authenticated_client.parse_json(response)
    
    # Assert the response
    assert response.status_code == 200
    assert data["name"] == "Test Molecule"
```

### Testing with Mock RDKit

```python
@pytest.mark.rdkit
def test_molecule_properties(mock_rdkit, molecule_factory):
    # Create a mock molecule
    mock_mol = molecule_factory(smiles="CCO", mol_weight=46.07)
    
    # Use the mock RDKit to calculate properties
    mol = mock_rdkit.Chem.MolFromSmiles("CCO")
    
    # Assert the properties
    assert mock_rdkit.Chem.MolToSmiles(mol) == "CCO"
    assert mock_rdkit.Descriptors.MolWt(mol) == 46.07
```

### End-to-End Testing

```python
@pytest.mark.e2e
def test_molecule_workflow(authenticated_client, mock_db):
    # Create a test molecule
    molecule = {
        "name": "New Molecule",
        "formula": "C2H6O",
        "smiles": "CCO"
    }
    
    # Add the molecule via API
    response = authenticated_client.post("/api/molecules", json_data=molecule)
    data = authenticated_client.parse_json(response)
    
    # Assert the response
    assert response.status_code == 201
    assert "id" in data
    
    molecule_id = data["id"]
    
    # Get the molecule via API
    response = authenticated_client.get(f"/api/molecules/{molecule_id}")
    data = authenticated_client.parse_json(response)
    
    # Assert the response
    assert response.status_code == 200
    assert data["name"] == "New Molecule"
    
    # Verify the molecule was added to the database
    result = mock_db.table("molecules").select("*").eq("id", molecule_id).execute()
    assert len(result.data) == 1
    assert result.data[0]["name"] == "New Molecule"
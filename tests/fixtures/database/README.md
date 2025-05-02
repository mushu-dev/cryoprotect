# Database Test Fixtures

This directory contains fixtures for testing database-related functionality.

## Available Fixtures

### Connection Fixtures
- `db_config`: Configuration for test database connections
- `mock_db`: Mock database connection
- `real_db_connection`: Real database connection for integration tests
- `populated_mock_db`: Mock database pre-populated with test data

### Migration Fixtures
- `migration_db`: Database prepared for migration testing
- `temp_migration_script`: Temporary migration script for testing

## Usage Examples

### Basic Mock Database

```python
def test_database_operation(mock_db):
    # The mock_db fixture provides a mock Supabase client
    # with basic test data already loaded
    
    result = mock_db.table('molecules').select().execute()
    assert len(result.data) == 2
    
    # You can add more test data
    mock_db.add_test_data('molecules', [
        {'id': 'mol-3', 'name': 'New Molecule'}
    ])
    
    # And then use it in your test
    result = mock_db.table('molecules').select().execute()
    assert len(result.data) == 3
```

### Pre-populated Database

```python
def test_with_comprehensive_data(populated_mock_db):
    # The populated_mock_db fixture provides a more complete
    # dataset with relationships between tables
    
    # Query molecules
    molecules = populated_mock_db.table('molecules').select().execute()
    
    # Query mixtures
    mixtures = populated_mock_db.table('mixtures').select().execute()
    
    # Find components for a mixture
    mix_id = mixtures.data[0]['id']
    components = populated_mock_db.table('mixture_components')\
                  .select()\
                  .eq('mixture_id', mix_id)\
                  .execute()
```

### Migration Testing

```python
def test_migration(migration_db, temp_migration_script):
    # Test applying a migration
    from database.migrations import apply_migrations
    
    migrations_dir = os.path.dirname(temp_migration_script)
    result = apply_migrations(migration_db, migrations_dir=migrations_dir)
    
    # Verify the migration was applied
    assert 'test_migration_table' in migration_db.tables
```

### Data Generation

```python
def test_with_generated_data():
    from tests.fixtures.database.data import generate_molecule, generate_mixture
    
    # Generate a test molecule
    molecule = generate_molecule(name="Custom Molecule")
    
    # Generate a test mixture with components
    mixture = generate_mixture(component_count=3)
    assert len(mixture['components']) == 3
```

## Best Practices

1. **Use the simplest fixture** that meets your needs. For most tests, `mock_db` is sufficient.

2. **Avoid using `real_db_connection`** except in dedicated integration tests.

3. **Isolate your tests** by not relying on specific test data IDs or properties that might change.

4. **Add custom RPC handlers** to the mock database when testing functions that use RPC calls:
   ```python
   def test_with_custom_rpc(mock_db):
       # Add a custom RPC handler
       mock_db.register_rpc_handler('my_function', lambda params: {'data': [{'result': 'success'}]})
       
       # Now you can test code that calls this RPC
       result = your_function_that_calls_rpc(mock_db)
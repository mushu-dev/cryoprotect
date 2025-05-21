# Supabase Direct Connection

This module provides a direct PostgreSQL connection to Supabase using psycopg2 with connection pooling for improved performance over MCP.

## Implementation

The implementation includes:

1. `SupabaseDirectConnection` class with singleton pattern (`get_instance`)
2. Thread-safe connection pooling using `psycopg2.pool.ThreadedConnectionPool`
3. Reading credentials from environment variables
4. SQL execution functions with parameterized queries and transaction management
5. Connection management methods

## Files

- `supabase_direct.py`: Main implementation file
- `test_supabase_direct.py`: Manual test script for verifying the implementation
- `tests/unit/test_supabase_direct.py`: Unit tests using mocking

## Testing Requirements

To run the tests, you need to install the following dependencies:

```bash
# Install psycopg2 (requires Microsoft Visual C++ 14.0 or greater on Windows)
pip install psycopg2

# Or install the binary version (may be easier on some platforms)
pip install psycopg2-binary
```

### Running Unit Tests

```bash
# Run unit tests
python -m unittest tests/unit/test_supabase_direct.py
```

### Running Manual Tests

To run the manual tests, you need to set the following environment variables:

```bash
# Set environment variables
export SUPABASE_DB_HOST=your-supabase-host
export SUPABASE_DB_PORT=5432
export SUPABASE_DB_NAME=postgres
export SUPABASE_DB_USER=postgres
export SUPABASE_DB_PASSWORD=your-password
export SUPABASE_DB_MIN_CONNECTIONS=2
export SUPABASE_DB_MAX_CONNECTIONS=10

# Run manual tests
python test_supabase_direct.py
```

## Usage Example

```python
# Import the module
from supabase_direct import SupabaseDirectConnection

# Get connection instance
db = SupabaseDirectConnection.get_instance()

# Execute query
result = db.execute_query("SELECT 1 as test")
print(f"Query Result: {result}")

# Execute parameterized query
result = db.execute_query("SELECT %(value)s as param_test", {"value": "test"})
print(f"Parameterized Query Result: {result}")

# Execute batch queries
queries = [
    "CREATE TEMPORARY TABLE test_batch (id SERIAL PRIMARY KEY, name TEXT)",
    "INSERT INTO test_batch (name) VALUES ('test1')",
    "INSERT INTO test_batch (name) VALUES ('test2')"
]
db.execute_batch(queries)

# Close connections when done
db.close_all()
```

## Implementation Notes

1. The implementation uses a singleton pattern to ensure only one connection pool is created.
2. Connection pooling is thread-safe, allowing multiple threads to use the same pool.
3. Connections are automatically released back to the pool after use.
4. The implementation includes error handling and transaction management.
5. Query results are returned as a list of dictionaries for easy access.
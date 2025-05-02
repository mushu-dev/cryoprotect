# Direct Database Connection Guide for CryoProtect

This guide provides comprehensive instructions for implementing and using direct PostgreSQL connections to Supabase in the CryoProtect project, replacing the slower MCP-based approach.

## Why Direct Connections?

MCP-based database operations suffer from several limitations:
- **Performance**: MCP adds significant overhead, especially for bulk operations
- **Flexibility**: MCP doesn't support all PostgreSQL features
- **Reliability**: Connection failures are harder to diagnose and recover from
- **Efficiency**: Batch operations are inefficient through MCP

Direct PostgreSQL connections offer:
- **10-50x faster performance** for bulk database operations
- **Better error handling** and connection recovery
- **Full SQL feature support**, including advanced queries and transactions
- **Connection pooling** for optimal resource usage
- **Efficient batch operations** for importing thousands of records

## Prerequisites

1. **Environment Variables**
   ```
   SUPABASE_DB_HOST=db.tsdlmynydfuypiugmkev.supabase.co
   SUPABASE_DB_PORT=5432
   SUPABASE_DB_NAME=postgres
   SUPABASE_DB_USER=postgres
   SUPABASE_DB_PASSWORD=your-password
   DB_MIN_CONNECTIONS=1
   DB_MAX_CONNECTIONS=10
   ```

2. **Supabase Configuration**
   - Your project must be active (not paused)
   - Your IP address must be allowed in the Supabase dashboard

3. **Required Packages**
   ```
   pip install psycopg2-binary python-dotenv
   ```

## Implementation Components

The direct database connection implementation consists of two main components:

1. **`postgres_direct.py`** - Connection manager with pooling
2. **`sql_executor.py`** - SQL execution utilities

### PostgresDirectConnection Class

The `PostgresDirectConnection` class provides:
- Singleton pattern for connection management
- Thread-safe connection pooling
- Transaction management
- Error handling and retry mechanisms
- Performance statistics

### SQL Executor Module

The `sql_executor.py` module provides:
- Simple functions for common database operations
- Parameterized query execution
- Batch operation support
- Transaction management
- Bulk insert optimization

## Usage Examples

### Basic Query Execution

```python
from sql_executor import execute_query

# Simple query
result = execute_query("SELECT COUNT(*) as count FROM molecules")
molecule_count = result[0]['count']
print(f"Database contains {molecule_count} molecules")

# Parameterized query
result = execute_query(
    "SELECT * FROM molecules WHERE name LIKE %s LIMIT 10",
    ("%glycerol%",)
)
for molecule in result:
    print(f"Found molecule: {molecule['name']}")
```

### Bulk Insert Operations

```python
from sql_executor import bulk_insert

# Prepare data
molecules_to_insert = [
    {
        'id': '123e4567-e89b-12d3-a456-426614174000',
        'name': 'Glycerol',
        'formula': 'C3H8O3',
        'molecular_weight': 92.09,
        'chembl_id': 'CHEMBL388978'
    },
    {
        'id': '223e4567-e89b-12d3-a456-426614174001',
        'name': 'DMSO',
        'formula': 'C2H6OS',
        'molecular_weight': 78.13,
        'chembl_id': 'CHEMBL1098659'
    }
]

# Perform bulk insert
bulk_insert('molecules', molecules_to_insert)
```

### Transaction Management

```python
from sql_executor import begin_transaction, commit_transaction, rollback_transaction

# Begin a transaction
transaction = begin_transaction()

try:
    # Execute queries within the transaction
    # These will be automatically prepared with the transaction's connection
    cursor = transaction.cursor()
    cursor.execute("INSERT INTO molecules (id, name) VALUES (%s, %s)", 
                  ('123e4567-e89b-12d3-a456-426614174002', 'Trehalose'))
    cursor.execute("INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value) VALUES (%s, %s, %s)",
                  ('123e4567-e89b-12d3-a456-426614174002', '123', 342.3))
    
    # Commit the transaction
    commit_transaction(transaction)
except Exception as e:
    # Rollback on error
    rollback_transaction(transaction)
    raise
```

## Batch Operations for Database Population

For large-scale database operations, use the batch processing capabilities:

```python
from postgres_direct import PostgresDirectConnection
from sql_executor import bulk_insert

# Get database connection
db = PostgresDirectConnection.get_instance()

# Prepare property data
property_data = []
for molecule_id, properties in molecules_with_properties.items():
    for prop_name, value in properties.items():
        property_type_id = get_property_type_id(prop_name)
        if property_type_id:
            property_data.append({
                'molecule_id': molecule_id,
                'property_type_id': property_type_id,
                'numeric_value': float(value) if isinstance(value, (int, float)) else None,
                'text_value': str(value) if not isinstance(value, (int, float)) else None,
                'created_at': 'NOW()',
                'updated_at': 'NOW()'
            })

# Process in batches of 1000
batch_size = 1000
for i in range(0, len(property_data), batch_size):
    batch = property_data[i:i+batch_size]
    bulk_insert('molecular_properties', batch)
    print(f"Inserted batch {i//batch_size + 1}/{(len(property_data) + batch_size - 1)//batch_size}")
```

## Performance Optimization Tips

1. **Batch Operations**
   - Group related operations into batches
   - Use `bulk_insert` for inserting multiple rows
   - Use transactions for atomicity

2. **Connection Pooling**
   - Configure appropriate pool size (`DB_MIN_CONNECTIONS` and `DB_MAX_CONNECTIONS`)
   - Always release connections back to the pool
   - Close all connections when done

3. **Transaction Management**
   - Use transactions for related operations
   - Keep transactions as short as possible
   - Ensure proper error handling and rollback

4. **Database Indexes**
   - Create indexes on commonly queried columns
   - Monitor query performance with `EXPLAIN ANALYZE`
   - Balance indexing with write performance

## Troubleshooting

### Connection Issues

**Issue**: `Error: Could not resolve hostname db.tsdlmynydfuypiugmkev.supabase.co`

**Solution**: Use IP address directly instead of hostname
```
# Run supabase_postgres_connection.py to get the IP
python supabase_postgres_connection.py

# Update your .env file with the IP
SUPABASE_DB_HOST=172.64.xx.xx
```

**Issue**: `Error: Password authentication failed for user postgres`

**Solution**: Verify password in .env file and Supabase dashboard

**Issue**: `Error: Connection refused`

**Solution**: Check if your IP is allowed in Supabase dashboard network settings

### Performance Issues

**Issue**: Slow batch operations

**Solution**:
- Increase batch size for bulk operations
- Use prepared statements for repeated operations
- Reduce transaction scope
- Add appropriate indexes

**Issue**: Connection pool exhaustion

**Solution**:
- Increase `DB_MAX_CONNECTIONS`
- Ensure connections are properly released
- Use connection timeouts to prevent deadlocks

## Migrating from MCP

When migrating existing code from MCP to direct connections:

1. Replace MCP imports:
   ```python
   # Old MCP-based code
   from use_mcp_tool import execute_sql, get_project_id
   
   # New direct connection code
   from sql_executor import execute_query, bulk_insert
   ```

2. Update execution calls:
   ```python
   # Old MCP-based execution
   result = execute_sql(query, get_project_id())
   
   # New direct execution
   result = execute_query(query)
   ```

3. Replace formatted SQL with parameterized queries:
   ```python
   # Old MCP-based approach (susceptible to SQL injection)
   query = f"SELECT * FROM molecules WHERE name = '{name}'"
   result = execute_sql(query, get_project_id())
   
   # New secure parameterized approach
   query = "SELECT * FROM molecules WHERE name = %s"
   result = execute_query(query, (name,))
   ```

## Best Practices

1. **Security**
   - Always use parameterized queries, never string formatting
   - Store credentials in environment variables, not code
   - Use least-privilege database roles

2. **Error Handling**
   - Implement proper try/except blocks
   - Log database errors with context
   - Always release connections in finally blocks

3. **Resource Management**
   - Close all connections when your application exits
   - Monitor connection pool usage
   - Release connections as soon as possible

4. **Performance Monitoring**
   - Use the built-in statistics in `PostgresDirectConnection.get_stats()`
   - Add query timing for critical operations
   - Periodically analyze and optimize slow queries

## Reference

### PostgresDirectConnection Methods

| Method | Description |
|--------|-------------|
| `get_instance()` | Get singleton connection manager instance |
| `execute_query(query, params, fetch_one)` | Execute a SQL query with parameters |
| `execute_batch(queries, transaction)` | Execute multiple queries, optionally in a transaction |
| `bulk_insert(table, data, columns, return_ids)` | Insert multiple rows efficiently |
| `begin_transaction()` | Begin a new transaction |
| `commit_transaction(transaction)` | Commit a transaction |
| `rollback_transaction(transaction)` | Rollback a transaction |
| `get_stats()` | Get connection and query statistics |
| `test_connection()` | Test database connectivity |
| `close_all()` | Close all connections in the pool |

### SQL Executor Functions

| Function | Description |
|----------|-------------|
| `execute_query(query, params, fetch_one)` | Execute a SQL query with parameters |
| `execute_transaction(queries)` | Execute multiple queries in a transaction |
| `bulk_insert(table, data, columns, return_ids)` | Insert multiple rows efficiently |
| `begin_transaction()` | Begin a new transaction |
| `commit_transaction(transaction)` | Commit a transaction |
| `rollback_transaction(transaction)` | Rollback a transaction |

## Conclusion

Direct PostgreSQL connections provide substantial performance and reliability improvements over MCP-based database operations, especially for large-scale data imports. By implementing connection pooling, parameterized queries, and efficient batch operations, you can achieve 10-50x faster database operations while maintaining security and reliability.
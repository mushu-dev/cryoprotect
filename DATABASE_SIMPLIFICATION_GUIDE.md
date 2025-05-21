# Database Simplification Guide

## Overview

We've simplified the database connection architecture in CryoProtect to make it easier to maintain, debug, and extend. This guide explains how to migrate from the old complex adapter-based system to the new streamlined approach.

## Why We Simplified

The previous database architecture had several issues:
- Overly complex with multiple layers of abstraction
- Difficult to debug and maintain
- Excessive fallback mechanisms that were rarely used
- Mixing of concerns between connection management and business logic
- Overengineered for the actual use cases

## New Architecture

The new architecture focuses on simplicity and directness:

1. **Core Database Module (`db.py`)**: 
   - Direct PostgreSQL connection pooling
   - Basic query execution functions
   - Transaction management
   - Connection pool management

2. **Utility Functions (`utils.py`)**:
   - Higher-level database operations
   - Business logic for molecule and property management
   - Import/export utilities
   - Database statistics and health checks

## Migration Steps

### Step 1: Update Imports

Replace old imports with new ones:

```python
# OLD
from database.adapter import get_adapter
from database.connection_manager import get_connection_manager

# NEW
from database import db, utils
```

### Step 2: Replace Connection Management

```python
# OLD
db = get_connection_manager()
conn = db.get_connection()
# ... use conn ...

# NEW
db.init_connection_pool()  # Only call once at application start
result = db.execute_query("SELECT * FROM molecules")
```

### Step 3: Update Transaction Handling

```python
# OLD
db = get_connection_manager()
transaction = db.begin_transaction()
try:
    # ... execute queries ...
    db.commit_transaction(transaction)
except Exception as e:
    db.rollback_transaction(transaction)
    raise

# NEW
with db.transaction() as cursor:
    cursor.execute("INSERT INTO molecules ...")
    cursor.execute("INSERT INTO molecular_properties ...")
    # Auto-commits on success, auto-rollbacks on exception
```

### Step 4: Update Batch Operations

```python
# OLD
db = get_connection_manager()
results = db.execute_batch([
    "SELECT * FROM molecules LIMIT 10",
    "SELECT * FROM property_types"
])

# NEW
results = db.execute_batch([
    "SELECT * FROM molecules LIMIT 10",
    "SELECT * FROM property_types"
])
```

### Step 5: Use Utility Functions

```python
# OLD - Custom query logic
db = get_connection_manager()
results = db.execute_query(
    "SELECT * FROM molecules WHERE id = %s",
    (molecule_id,)
)

# NEW - Use utility functions
molecule = utils.get_molecule_by_id(molecule_id)
properties = utils.get_molecule_properties(molecule_id)
```

## ChEMBL and PubChem Import Scripts

For import scripts that need to insert large amounts of data:

1. Use `batch_save_properties()` for saving multiple properties in a single transaction
2. Use `db.transaction()` for custom transaction logic
3. Use `utils.update_import_status()` to track import progress

Example:

```python
# Initialize database
from database import db, utils
db.init_connection_pool()

# Update status
utils.update_import_status('pubchem', 'in_progress', progress=0)

# Process data in batches
for batch in data_batches:
    with db.transaction() as cursor:
        for item in batch:
            # Process each item
            cursor.execute("INSERT INTO ...")
    
    # Update progress
    progress = (current_count / total_count) * 100
    utils.update_import_status('pubchem', 'in_progress', progress=progress)

# Mark as complete
utils.update_import_status('pubchem', 'completed', progress=100)
```

## Configuration

Database connection is configured through environment variables:

```
DB_HOST=localhost
DB_PORT=5432
DB_NAME=cryoprotect
DB_USER=postgres
DB_PASSWORD=yourpassword
```

For Supabase connections, you can use the existing Supabase-specific variables:

```
SUPABASE_DB_HOST=aws-0-us-east-1.pooler.supabase.com
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres.projectid
SUPABASE_DB_PASSWORD=password
```

## Testing

A test script is provided to verify the database module works correctly:

```bash
python test_simplified_db.py
```

## Best Practices

1. **Initialize once**: Call `db.init_connection_pool()` only once at application start
2. **Use transactions**: Always use transactions for related operations
3. **Close connections**: Call `db.close_all_connections()` when shutting down
4. **Batch operations**: Use batch operations for better performance
5. **Error handling**: Use try/except blocks and provide specific error messages
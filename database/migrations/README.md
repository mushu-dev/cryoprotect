# Migrations

This directory manages database schema migrations and versioning. Migration scripts and tools for upgrading or downgrading the database schema are stored here.

## Overview

The migrations module provides a command line interface and programmatic API for managing database schema migrations for the CryoProtect v2 project.

## Command Line Usage

The migrations module provides a command line interface for common migration tasks:

```bash
# Show migration status
python -m database.migrations status

# Apply all pending migrations
python -m database.migrations apply

# Apply migrations up to a specific version
python -m database.migrations apply --target 005

# Roll back the most recent migration
python -m database.migrations rollback --target 004

# Roll back all migrations
python -m database.migrations rollback

# Dry run to see what would happen
python -m database.migrations apply --dry-run
```

## Programmatic Usage

The migrations module can also be used programmatically:

```python
from database.migrations import apply_migrations, rollback_migrations, get_migration_status

# Get migration status
status = get_migration_status()

# Apply migrations
apply_migrations(target_version="005", environment="development")

# Roll back migrations
rollback_migrations(target_version="004", environment="development")
```

## Creating New Migrations

To create a new migration:

1. Create a new Python file in the `database/migrations/scripts` directory
2. Follow the naming convention: `NNN_descriptive_name.py` (e.g., `006_add_user_preferences.py`)
3. Implement the required functions:

```python
def apply(conn, environment):
    # Code to apply the migration
    conn.sql("CREATE TABLE new_table (...)").execute()

def rollback(conn, environment):
    # Code to roll back the migration
    conn.sql("DROP TABLE new_table").execute()
```

## Best Practices

1. Always include both apply and rollback functions
2. Make migrations idempotent when possible
3. Test migrations thoroughly before applying to production
4. Keep migrations small and focused
5. Include appropriate logging

## Directory Structure

- `__init__.py`: Package exports
- `runner.py`: Main migration runner module
- `tracker.py`: Migration tracking utilities
- `scripts/`: Directory containing migration scripts
  - `__init__.py`: Package initialization
  - `001_initial_schema.py`: Example migration script
  - `...`: Additional migration scripts

## Migration Script Structure

Each migration script should:

1. Have a clear, descriptive name following the pattern `NNN_descriptive_name.py`
2. Include docstrings explaining the purpose of the migration
3. Implement both `apply` and `rollback` functions
4. Handle errors gracefully
5. Include appropriate logging

Example:

```python
"""
Add user preferences table.

This migration adds a table for storing user preferences.
"""

import logging

# Configure logging
logger = logging.getLogger(__name__)

def apply(conn, environment):
    """Apply the migration."""
    logger.info(f"Applying migration in {environment} environment")
    
    conn.sql("""
        CREATE TABLE user_preferences (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
            theme TEXT NOT NULL DEFAULT 'light',
            notifications BOOLEAN NOT NULL DEFAULT true,
            created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
            updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        CREATE UNIQUE INDEX user_preferences_user_id_idx ON user_preferences(user_id);
    """).execute()

def rollback(conn, environment):
    """Roll back the migration."""
    logger.info(f"Rolling back migration in {environment} environment")
    
    conn.sql("""
        DROP TABLE IF EXISTS user_preferences;
    """).execute()
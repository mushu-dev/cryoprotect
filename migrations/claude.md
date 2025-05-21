# Migrations Module Reference

## Overview

The migrations module manages database schema evolution for CryoProtect, providing tools for versioned schema changes, database upgrades/downgrades, and reliable schema management.

## Key Components

### Migration System
- **Migration Scripts**: Versioned SQL/Python scripts in `migrations/scripts/`
- **Runner**: Execution management in `migrations/runner.py`
- **Tracker**: Applied migration tracking in `migrations/tracker.py`
- **CLI Interface**: Command-line interface for migration management

### Migration Design
- **Sequential Versioning**: Migrations are numbered sequentially (001, 002, etc.)
- **Idempotent Execution**: Safe to apply multiple times
- **Atomic Operations**: All-or-nothing execution
- **Reversible Changes**: Support for rollback operations

## Migration Workflow

1. **Creation**: Create a new migration script with up/down operations
2. **Validation**: Verify migration script correctness
3. **Execution**: Apply migration to target database
4. **Verification**: Confirm successful application
5. **Tracking**: Record applied migration in tracking table

## Command-Line Usage

```bash
# Show migration status
python -m database.migrations status

# Apply all pending migrations
python -m database.migrations apply

# Apply migrations up to a specific version
python -m database.migrations apply --target 005

# Roll back the most recent migration
python -m database.migrations rollback

# Roll back to a specific version
python -m database.migrations rollback --target 004

# Dry run to see what would happen
python -m database.migrations apply --dry-run
```

## Creating New Migrations

### Migration Script Structure

Each migration script includes:
- **Metadata**: Version number and description
- **Apply Function**: Code to apply the migration
- **Rollback Function**: Code to undo the migration
- **Documentation**: Purpose and changes

### Example Migration Script

```python
"""
006_add_user_preferences.py - Add user preferences table

This migration adds a table for storing user preferences with theme settings and notifications.
"""

import logging

logger = logging.getLogger(__name__)

def apply(conn, environment):
    """Apply the migration."""
    logger.info(f"Applying migration 006 in {environment} environment")
    
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
    logger.info(f"Rolling back migration 006 in {environment} environment")
    
    conn.sql("""
        DROP TABLE IF EXISTS user_preferences;
    """).execute()
```

## Migration Types

### Schema Changes
- Table creation/modification/deletion
- Index creation/modification/deletion
- Constraint management
- View management

### Data Migrations
- Data transformation
- Data seeding
- Reference data setup
- Data cleanup

### Function/Procedure Changes
- Database function creation/updates
- Stored procedure management
- Trigger management

## Best Practices

1. **Keep Migrations Small**: Focus on related changes
2. **Include Both Apply/Rollback**: Always implement rollback functions
3. **Test Migrations Thoroughly**: Verify on dev before production
4. **Make Them Idempotent**: Safe to run multiple times
5. **Version Control**: Keep migrations in source control
6. **Document Purpose**: Clear descriptions of what and why
7. **Avoid Large Data Operations**: Split data migrations if needed

## Common Pitfalls

1. **Missing Dependencies**: Not respecting migration ordering
2. **Incomplete Rollbacks**: Not properly reversing all changes
3. **Transaction Issues**: Long-running migrations that timeout
4. **Performance Problems**: Schema changes that lock tables
5. **Missing Indexes**: Schema changes that affect performance
6. **Inadequate Testing**: Migrations that work in dev but fail in production

## Integration with CI/CD

- **Automated Testing**: Run migrations in test environment
- **Migration Verification**: Validate applied migrations
- **Blue/Green Deployment**: Apply migrations before switching
- **Rollback Planning**: Have rollback strategy for failed deployments
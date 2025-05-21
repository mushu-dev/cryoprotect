# CryoProtect Convex Migration Framework

This document describes the migration framework for Convex integration in the CryoProtect application. The framework allows for seamless management of database migrations between Supabase and Convex databases.

## Overview

The Convex Migration Framework provides:

1. A TypeScript-based migration system for Convex
2. A Python bridge for integrating with the existing Supabase migration system
3. CLI tools for managing migrations
4. Data synchronization utilities between Supabase and Convex

The framework is designed to enable the gradual transition from Supabase to Convex while maintaining backward compatibility.

## Architecture

The migration framework consists of the following components:

### 1. Convex Migration System

The core Convex migration components are located in the `convex/migrations` directory:

- `tracker.ts`: Functions for tracking applied migrations
- `runner.ts`: Functions for applying and rolling back migrations
- `schema.ts`: Schema definitions for the migrations table
- `cli.ts`: Command-line interface for managing migrations
- `examples/`: Directory containing example migration files

### 2. Python Bridge

The Python bridge integrates the Convex migration system with the existing Supabase migration system:

- `database/migrations/convex_bridge.py`: Bridge module with functions for interacting with the Convex migration system

### 3. Utility Scripts

Several utility scripts are provided for common tasks:

- `run_convex_migrations.sh`: Bash script for running Convex migrations
- `sync_supabase_to_convex.py`: Python script for syncing data from Supabase to Convex

## Migration Format

Convex migrations are defined as JSON files with the following structure:

```json
{
  "version": "001",
  "name": "initial_schema",
  "description": "Initial database schema",
  "schema": {
    "users": {
      "name": "string",
      "email": "string",
      "password_hash": "string?",
      "role": "string?",
      "created_at": "number",
      "updated_at": "number"
    }
  },
  "data": [
    {
      "type": "insert",
      "table": "users",
      "data": {
        "name": "System Administrator",
        "email": "admin@cryoprotect.com",
        "password_hash": "5f4dcc3b5aa765d61d8327deb882cf99",
        "role": "admin",
        "created_at": 1620000000000,
        "updated_at": 1620000000000
      }
    }
  ],
  "indexes": [
    {
      "table": "users",
      "name": "by_email",
      "fields": ["email"]
    }
  ]
}
```

Migration files are named with a sequential 3-digit prefix followed by a snake_case name (e.g., `001_initial_schema.json`).

## Usage

### Setting Up

Before using the migration framework, make sure Convex is enabled in your environment:

```bash
export USE_CONVEX=true
export CONVEX_URL=https://your-deployment.convex.cloud
```

You can add these to your `.env` file for persistence.

### Running Migrations

To check the status of migrations:

```bash
./run_convex_migrations.sh status
```

To apply pending migrations:

```bash
./run_convex_migrations.sh apply
```

To apply migrations up to a specific version:

```bash
./run_convex_migrations.sh apply --target=003
```

To roll back applied migrations:

```bash
./run_convex_migrations.sh rollback
```

To roll back migrations to a specific version:

```bash
./run_convex_migrations.sh rollback --target=002
```

To create a new migration:

```bash
./run_convex_migrations.sh create --name="add_new_feature"
```

### Running Universal Migrations

To run migrations on both Supabase and Convex databases:

```python
from database.migrations.convex_bridge import run_universal_migrations

results = run_universal_migrations(
    target_version="003",
    dry_run=True,
    environment="development"
)

print(results)
```

### Syncing Data

To sync data from Supabase to Convex:

```bash
# Sync all tables
./sync_supabase_to_convex.py

# Sync specific tables
./sync_supabase_to_convex.py users teams projects

# Dry run to see what would be synced
./sync_supabase_to_convex.py --dry-run

# Limit the number of records per table
./sync_supabase_to_convex.py --limit=500

# Output in JSON format
./sync_supabase_to_convex.py --format=json
```

## Migration Best Practices

1. **Version Numbers**: Use sequential 3-digit version numbers (001, 002, 003, etc.)
2. **Migration Names**: Use descriptive snake_case names
3. **Schema Changes**: Define schema changes in both the migration file and the main `schema.ts` file
4. **Data Migrations**: Include original data in updates and deletes for rollback support
5. **Testing**: Always test migrations with `--dry-run` first
6. **Sync After Migrations**: After applying schema migrations, sync any affected data

## Integration with Existing Migrations

The Convex Migration Framework integrates with the existing Supabase migration system through the Python bridge. This allows for:

1. Running migrations on both databases with a single command
2. Tracking the applied state of migrations in both databases
3. Syncing data between databases as needed

When running `run_universal_migrations()`, the system will apply migrations to both Supabase and Convex (if enabled) and provide a unified result.

## Fallback Mechanism

If Convex is not enabled (`USE_CONVEX=false`), the migration bridge will fall back to using only the Supabase migration system. This allows for gradual adoption of Convex without disrupting existing workflows.

## Troubleshooting

If you encounter issues with the migration framework, check the following:

1. **Environment Variables**: Make sure `USE_CONVEX` and `CONVEX_URL` are set correctly
2. **Node.js and ts-node**: Ensure Node.js and ts-node are installed for the CLI tools
3. **Migration Logs**: Check the logs in the `migration_logs` directory
4. **Network Connectivity**: Verify connectivity to the Convex deployment
5. **Permission Issues**: Make sure the scripts have execute permissions

For more detailed troubleshooting, run the commands with increased verbosity:

```bash
NODE_DEBUG=* ./run_convex_migrations.sh status
```

## Conclusion

The Convex Migration Framework provides a robust solution for managing database migrations in a hybrid Supabase/Convex environment. By following the best practices and using the provided tools, you can ensure a smooth transition to Convex while maintaining backward compatibility with existing code.
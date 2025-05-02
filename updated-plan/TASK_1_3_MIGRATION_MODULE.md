# Task 1.3: Create Migration Management Module

## Objective
Implement a centralized migration management module that provides tools for managing database schema migrations in a consistent, reliable manner.

## Context
The project currently uses multiple separate migration scripts and approaches without a unified system. This leads to inconsistent migration patterns, difficulties tracking migration history, and potential issues when deploying to different environments. A consolidated migration management module will provide a reliable way to apply, track, and roll back database schema changes.

## Acceptance Criteria
- Create a modular migration management package with clean API
- Support for applying, rolling back, and tracking migrations
- CLI interface for migration operations
- Version tracking for applied migrations
- Support for environment-specific migrations
- Migration verification and validation
- Clear documentation and examples

## Implementation Steps

1. Create the migration management module structure:
   ```
   database/migrations/
   ├── __init__.py          # Package exports
   ├── runner.py            # Main entry point
   ├── tracker.py           # Migration tracking utilities
   ├── scripts/             # Migration scripts directory
   │   └── __init__.py
   └── utils.py             # Migration utilities
   ```

2. Implement the migrations module `__init__.py`:
   ```python
   """
   Database migration module for CryoProtect v2.
   
   This module provides tools for managing database schema migrations
   including applying, rolling back, and tracking migration states.
   """
   
   from database.migrations.runner import (
       apply_migrations,
       rollback_migrations,
       get_migration_status,
       initialize_migration_tracking
   )
   
   __all__ = [
       'apply_migrations',
       'rollback_migrations',
       'get_migration_status',
       'initialize_migration_tracking'
   ]
   ```

3. Implement the main migration runner module:
   ```python
   """
   Main entry point for database migration operations.
   
   This module provides functions for applying and rolling back migrations,
   as well as tracking migration status.
   """
   
   import logging
   import os
   import importlib.util
   from datetime import datetime
   from typing import Dict, List, Optional, Tuple, Union, Any
   
   from database.utils.connection import create_connection
   from database.migrations.tracker import (
       get_applied_migrations,
       record_migration,
       remove_migration_record
   )
   
   logger = logging.getLogger(__name__)
   
   def _load_migration_script(script_path: str) -> Any:
       """
       Dynamically load a migration script module.
       
       Args:
           script_path: Path to the migration script
           
       Returns:
           Loaded module
       """
       script_name = os.path.basename(script_path)
       module_name = os.path.splitext(script_name)[0]
       
       spec = importlib.util.spec_from_file_location(module_name, script_path)
       module = importlib.util.module_from_spec(spec)
       spec.loader.exec_module(module)
       
       return module
   
   def _get_migration_scripts(
       migrations_dir: Optional[str] = None
   ) -> List[Tuple[str, str]]:
       """
       Get available migration scripts ordered by version.
       
       Args:
           migrations_dir: Directory containing migration scripts
           
       Returns:
           List of tuples (version, script_path)
       """
       if migrations_dir is None:
           # Default to scripts directory within this package
           migrations_dir = os.path.join(
               os.path.dirname(os.path.abspath(__file__)),
               'scripts'
           )
           
       if not os.path.exists(migrations_dir):
           logger.warning(f"Migrations directory not found: {migrations_dir}")
           return []
           
       scripts = []
       for filename in os.listdir(migrations_dir):
           if filename.endswith('.py') and not filename.startswith('__'):
               # Extract version from filename (e.g., 001_initial_schema.py -> 001)
               parts = filename.split('_', 1)
               if len(parts) >= 2 and parts[0].isdigit():
                   version = parts[0]
                   script_path = os.path.join(migrations_dir, filename)
                   scripts.append((version, script_path))
               
       # Sort by version
       return sorted(scripts, key=lambda x: x[0])
       
   def initialize_migration_tracking(
       conn: Any,
       config: Optional[Dict] = None
   ) -> bool:
       """
       Initialize migration tracking.
       
       Creates the migrations table if it doesn't exist.
       
       Args:
           conn: Database connection
           config: Optional configuration dictionary
           
       Returns:
           True if initialization was successful, False otherwise
       """
       logger.info("Initializing migration tracking")
       
       if conn is None:
           conn = create_connection(config)
       
       try:
           # Check if the migrations table exists
           result = conn.rpc('has_table', {'table_name': 'migrations'}).execute()
           table_exists = result.data[0] if result.data else False
           
           if not table_exists:
               # Create migrations table
               conn.sql("""
                   CREATE TABLE migrations (
                       id SERIAL PRIMARY KEY,
                       version VARCHAR(50) NOT NULL,
                       name VARCHAR(255) NOT NULL,
                       applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                       status VARCHAR(50) DEFAULT 'applied'
                   );
                   CREATE UNIQUE INDEX migrations_version_idx ON migrations(version);
               """).execute()
               logger.info("Created migrations table")
           
           return True
       except Exception as e:
           logger.error(f"Error initializing migration tracking: {str(e)}")
           return False
   
   def get_migration_status(
       conn: Optional[Any] = None,
       config: Optional[Dict] = None,
       migrations_dir: Optional[str] = None
   ) -> Dict[str, Dict]:
       """
       Get status of all migrations.
       
       Args:
           conn: Database connection
           config: Optional configuration dictionary
           migrations_dir: Directory containing migration scripts
           
       Returns:
           Dictionary mapping versions to status dictionaries
       """
       logger.info("Getting migration status")
       
       if conn is None:
           conn = create_connection(config)
           
       # Initialize migration tracking if needed
       initialize_migration_tracking(conn)
       
       # Get available migrations
       available_migrations = _get_migration_scripts(migrations_dir)
       
       # Get applied migrations
       applied_migrations = get_applied_migrations(conn)
       applied_versions = {m['version']: m for m in applied_migrations}
       
       # Build status dictionary
       status = {}
       for version, script_path in available_migrations:
           script_name = os.path.basename(script_path)
           name = script_name.split('_', 1)[1].rsplit('.', 1)[0] if '_' in script_name else script_name
           
           if version in applied_versions:
               applied_info = applied_versions[version]
               status[version] = {
                   'version': version,
                   'name': name,
                   'script': script_path,
                   'status': applied_info['status'],
                   'applied_at': applied_info['applied_at']
               }
           else:
               status[version] = {
                   'version': version,
                   'name': name,
                   'script': script_path,
                   'status': 'pending',
                   'applied_at': None
               }
               
       return status
   
   def apply_migrations(
       conn: Optional[Any] = None,
       config: Optional[Dict] = None,
       migrations_dir: Optional[str] = None,
       target_version: Optional[str] = None,
       environment: str = 'development',
       dry_run: bool = False
   ) -> List[Dict]:
       """
       Apply pending migrations up to target version.
       
       Args:
           conn: Database connection
           config: Optional configuration dictionary
           migrations_dir: Directory containing migration scripts
           target_version: Target version to migrate to, or None for all
           environment: Target environment
           dry_run: If True, don't actually apply migrations
           
       Returns:
           List of applied migration info dictionaries
       """
       logger.info(f"Applying migrations to {target_version or 'latest'} in {environment} environment")
       
       if conn is None:
           conn = create_connection(config)
           
       # Initialize migration tracking if needed
       initialize_migration_tracking(conn)
       
       # Get migration status
       status = get_migration_status(conn, config, migrations_dir)
       
       # Determine which migrations to apply
       pending_migrations = [
           (version, info['script'])
           for version, info in status.items()
           if info['status'] == 'pending' and (target_version is None or version <= target_version)
       ]
       pending_migrations.sort(key=lambda x: x[0])  # Sort by version
       
       if not pending_migrations:
           logger.info("No pending migrations to apply")
           return []
           
       # Apply migrations
       applied = []
       for version, script_path in pending_migrations:
           script_name = os.path.basename(script_path)
           name = script_name.split('_', 1)[1].rsplit('.', 1)[0] if '_' in script_name else script_name
           
           logger.info(f"Applying migration {version}: {name}")
           
           # Load migration script
           try:
               module = _load_migration_script(script_path)
               
               if not hasattr(module, 'apply') or not callable(getattr(module, 'apply')):
                   logger.error(f"Migration {version} does not have an 'apply' function")
                   continue
                   
               # Apply migration
               if not dry_run:
                   try:
                       module.apply(conn, environment)
                       record_migration(conn, version, name)
                       
                       migration_info = {
                           'version': version,
                           'name': name,
                           'status': 'applied',
                           'applied_at': datetime.now().isoformat()
                       }
                       applied.append(migration_info)
                       logger.info(f"Applied migration {version}: {name}")
                   except Exception as e:
                       logger.error(f"Error applying migration {version}: {str(e)}")
                       break
               else:
                   logger.info(f"Would apply migration {version}: {name} (dry run)")
                   migration_info = {
                       'version': version,
                       'name': name,
                       'status': 'would_apply',
                       'applied_at': None
                   }
                   applied.append(migration_info)
           except Exception as e:
               logger.error(f"Error loading migration {version}: {str(e)}")
               break
               
       return applied
   
   def rollback_migrations(
       conn: Optional[Any] = None,
       config: Optional[Dict] = None,
       migrations_dir: Optional[str] = None,
       target_version: Optional[str] = None,
       environment: str = 'development',
       dry_run: bool = False
   ) -> List[Dict]:
       """
       Roll back applied migrations down to target version.
       
       Args:
           conn: Database connection
           config: Optional configuration dictionary
           migrations_dir: Directory containing migration scripts
           target_version: Target version to rollback to, or None for all
           environment: Target environment
           dry_run: If True, don't actually roll back migrations
           
       Returns:
           List of rolled back migration info dictionaries
       """
       logger.info(f"Rolling back migrations to {target_version or 'beginning'} in {environment} environment")
       
       if conn is None:
           conn = create_connection(config)
           
       # Initialize migration tracking if needed
       initialize_migration_tracking(conn)
       
       # Get migration status
       status = get_migration_status(conn, config, migrations_dir)
       
       # Determine which migrations to roll back
       applied_migrations = [
           (version, info['script'])
           for version, info in status.items()
           if info['status'] == 'applied' and (target_version is None or version > target_version)
       ]
       applied_migrations.sort(key=lambda x: x[0], reverse=True)  # Sort by version, descending
       
       if not applied_migrations:
           logger.info("No migrations to roll back")
           return []
           
       # Roll back migrations
       rolled_back = []
       for version, script_path in applied_migrations:
           script_name = os.path.basename(script_path)
           name = script_name.split('_', 1)[1].rsplit('.', 1)[0] if '_' in script_name else script_name
           
           logger.info(f"Rolling back migration {version}: {name}")
           
           # Load migration script
           try:
               module = _load_migration_script(script_path)
               
               if not hasattr(module, 'rollback') or not callable(getattr(module, 'rollback')):
                   logger.error(f"Migration {version} does not have a 'rollback' function")
                   continue
                   
               # Roll back migration
               if not dry_run:
                   try:
                       module.rollback(conn, environment)
                       remove_migration_record(conn, version)
                       
                       migration_info = {
                           'version': version,
                           'name': name,
                           'status': 'rolled_back',
                           'applied_at': None
                       }
                       rolled_back.append(migration_info)
                       logger.info(f"Rolled back migration {version}: {name}")
                   except Exception as e:
                       logger.error(f"Error rolling back migration {version}: {str(e)}")
                       break
               else:
                   logger.info(f"Would roll back migration {version}: {name} (dry run)")
                   migration_info = {
                       'version': version,
                       'name': name,
                       'status': 'would_rollback',
                       'applied_at': None
                   }
                   rolled_back.append(migration_info)
           except Exception as e:
               logger.error(f"Error loading migration {version}: {str(e)}")
               break
               
       return rolled_back
   
   def main():
       """
       CLI entry point for database migrations.
       """
       import argparse
       import json
       
       parser = argparse.ArgumentParser(description='Manage CryoProtect database migrations')
       parser.add_argument(
           'action',
           choices=['apply', 'rollback', 'status'],
           help='Migration action to perform'
       )
       parser.add_argument(
           '--env', 
           choices=['development', 'staging', 'production'],
           default='development',
           help='Target environment'
       )
       parser.add_argument(
           '--target',
           help='Target version for apply/rollback'
       )
       parser.add_argument(
           '--dir',
           help='Directory containing migration scripts'
       )
       parser.add_argument(
           '--dry-run',
           action='store_true',
           help='Print what would be done without actually doing it'
       )
       parser.add_argument(
           '--format',
           choices=['text', 'json'],
           default='text',
           help='Output format'
       )
       
       args = parser.parse_args()
       
       conn = create_connection()
       
       if args.action == 'status':
           status = get_migration_status(conn, migrations_dir=args.dir)
           if args.format == 'json':
               print(json.dumps(status, indent=2, default=str))
           else:
               print("Migration Status:")
               print("-" * 60)
               for version, info in sorted(status.items()):
                   applied_at = info['applied_at'] or 'N/A'
                   print(f"{version}: {info['name']} - {info['status']} ({applied_at})")
       elif args.action == 'apply':
           results = apply_migrations(
               conn, 
               migrations_dir=args.dir,
               target_version=args.target,
               environment=args.env,
               dry_run=args.dry_run
           )
           if args.format == 'json':
               print(json.dumps(results, indent=2, default=str))
           else:
               print("Applied Migrations:")
               print("-" * 60)
               for info in results:
                   print(f"{info['version']}: {info['name']} - {info['status']}")
       elif args.action == 'rollback':
           results = rollback_migrations(
               conn, 
               migrations_dir=args.dir,
               target_version=args.target,
               environment=args.env,
               dry_run=args.dry_run
           )
           if args.format == 'json':
               print(json.dumps(results, indent=2, default=str))
           else:
               print("Rolled Back Migrations:")
               print("-" * 60)
               for info in results:
                   print(f"{info['version']}: {info['name']} - {info['status']}")
       
   if __name__ == '__main__':
       main()
   ```

4. Implement the migration tracker module:
   ```python
   """
   Migration tracking utilities.
   
   This module provides functions for tracking applied migrations.
   """
   
   import logging
   from typing import Any, Dict, List, Optional
   
   logger = logging.getLogger(__name__)
   
   def get_applied_migrations(conn: Any) -> List[Dict]:
       """
       Get list of applied migrations.
       
       Args:
           conn: Database connection
           
       Returns:
           List of applied migration dictionaries
       """
       try:
           result = conn.table('migrations').select('*').order('version').execute()
           return result.data if result.data else []
       except Exception as e:
           logger.error(f"Error getting applied migrations: {str(e)}")
           return []
   
   def record_migration(
       conn: Any,
       version: str,
       name: str,
       status: str = 'applied'
   ) -> bool:
       """
       Record a migration as applied.
       
       Args:
           conn: Database connection
           version: Migration version
           name: Migration name
           status: Migration status
           
       Returns:
           True if recording was successful, False otherwise
       """
       try:
           conn.table('migrations').insert({
               'version': version,
               'name': name,
               'status': status
           }).execute()
           return True
       except Exception as e:
           logger.error(f"Error recording migration {version}: {str(e)}")
           return False
   
   def remove_migration_record(conn: Any, version: str) -> bool:
       """
       Remove a migration record.
       
       Args:
           conn: Database connection
           version: Migration version
           
       Returns:
           True if removal was successful, False otherwise
       """
       try:
           conn.table('migrations').delete().eq('version', version).execute()
           return True
       except Exception as e:
           logger.error(f"Error removing migration record {version}: {str(e)}")
           return False
   ```

5. Create a sample migration script template:
   ```python
   """
   Migration: 001_initial_schema
   
   Creates the initial database schema.
   """
   
   import logging
   from typing import Any
   
   logger = logging.getLogger(__name__)
   
   def apply(conn: Any, environment: str):
       """
       Apply the migration.
       
       Args:
           conn: Database connection
           environment: Target environment
       """
       logger.info(f"Applying migration 001_initial_schema in {environment} environment")
       
       # Example: create tables
       conn.sql("""
           CREATE TABLE IF NOT EXISTS molecules (
               id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
               name VARCHAR(255) NOT NULL,
               formula VARCHAR(255),
               smiles TEXT,
               created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
               updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
           );
           
           CREATE TABLE IF NOT EXISTS mixtures (
               id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
               name VARCHAR(255) NOT NULL,
               description TEXT,
               created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
               updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
           );
           
           CREATE TABLE IF NOT EXISTS mixture_components (
               id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
               mixture_id UUID REFERENCES mixtures(id) ON DELETE CASCADE,
               molecule_id UUID REFERENCES molecules(id) ON DELETE RESTRICT,
               concentration NUMERIC(10, 5),
               units VARCHAR(50),
               created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
               updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
           );
       """).execute()
   
   def rollback(conn: Any, environment: str):
       """
       Roll back the migration.
       
       Args:
           conn: Database connection
           environment: Target environment
       """
       logger.info(f"Rolling back migration 001_initial_schema in {environment} environment")
       
       # Drop tables in reverse order
       conn.sql("""
           DROP TABLE IF EXISTS mixture_components;
           DROP TABLE IF EXISTS mixtures;
           DROP TABLE IF EXISTS molecules;
       """).execute()
   ```

6. Create a README for the migrations module:
   ```markdown
   # Database Migrations Module

   This module provides tools for managing database schema migrations for the CryoProtect v2 project.

   ## Usage

   ### Command Line Interface

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

   ### Programmatic Usage

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

   ## Creating Migrations

   To create a new migration:

   1. Create a new Python file in the `database/migrations/scripts` directory
   2. Follow the naming convention: `NNN_descriptive_name.py` (e.g., `006_add_user_preferences.py`)
   3. Implement both `apply` and `rollback` functions:

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
   ```

7. Update the main database `__init__.py` to include the migrations module:
   ```python
   """
   Database package for CryoProtect v2.
   
   This package centralizes all database-related functionality including
   models, operations, migrations, and verification tools.
   """
   
   from database.main import (
       initialize_database,
       populate_database,
       verify_database,
       backup_database,
       restore_database
   )
   
   # Include the population module exports
   from database.population import (
       populate_all,
       populate_specific,
       populate_from_file
   )
   
   # Include the migrations module exports
   from database.migrations import (
       apply_migrations,
       rollback_migrations,
       get_migration_status,
       initialize_migration_tracking
   )
   
   __all__ = [
       'initialize_database',
       'populate_database',
       'verify_database',
       'backup_database',
       'restore_database',
       'populate_all',
       'populate_specific',
       'populate_from_file',
       'apply_migrations',
       'rollback_migrations',
       'get_migration_status',
       'initialize_migration_tracking'
   ]
   ```

## Files to Modify

- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/`
- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/scripts/`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/runner.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/tracker.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/scripts/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/scripts/001_initial_schema.py` (example)
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/migrations/README.md`
- Update existing file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/__init__.py`

## Verification
1. Verify that all modules import properly without errors
2. Check that migration tracking table can be created
3. Verify CLI entry point works with all parameters
4. Test apply and rollback functionality with a sample migration
5. Ensure migration status reporting works correctly
6. Verify error handling and logging functionality
7. Check compatibility with existing migration scripts
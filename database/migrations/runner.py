"""
Main entry point for database migration operations.

This module provides functions for applying and rolling back migrations,
as well as tracking migration status.
"""

import os
import sys
import argparse
import importlib.util
import logging
from typing import Any, Dict, List, Optional, Tuple

from database.utils.connection import create_connection
from database.migrations.tracker import (
    get_applied_migrations,
    record_migration,
    remove_migration_record
)

# Configure logging
logger = logging.getLogger(__name__)

def _load_migration_script(script_path: str) -> Any:
    """
    Dynamically load a migration script module.

    Args:
        script_path: Path to the migration script

    Returns:
        Loaded module object
    """
    try:
        module_name = os.path.basename(script_path).replace('.py', '')
        spec = importlib.util.spec_from_file_location(module_name, script_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load spec for {script_path}")
        
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        logger.error(f"Error loading migration script {script_path}: {str(e)}")
        raise

def _get_migration_scripts(
    migrations_dir: Optional[str] = None
) -> List[Tuple[str, str]]:
    """
    Get available migration scripts ordered by version.

    Args:
        migrations_dir: Directory containing migration scripts

    Returns:
        List of (version, script_path) tuples
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
            try:
                parts = filename.split('_', 1)
                if len(parts) >= 2:
                    version = parts[0]
                    script_path = os.path.join(migrations_dir, filename)
                    scripts.append((version, script_path))
            except Exception as e:
                logger.warning(f"Invalid migration script filename: {filename}")
    
    return sorted(scripts, key=lambda x: x[0])

def initialize_migration_tracking(
    conn: Any,
) -> bool:
    """
    Initialize migration tracking.

    Creates the migrations table if it doesn't exist.

    Args:
        conn: Database connection object

    Returns:
        True if successful, False otherwise
    """
    logger.info("Initializing migration tracking")

    if conn is None:
        conn = create_connection()

    try:
        # Check if the migrations table exists
        result = conn.rpc('has_table', {'table_name': 'migrations'}).execute()
        table_exists = result.data[0] if result.data else False

        if not table_exists:
            # Create migrations table
            conn.sql("""
                CREATE TABLE migrations (
                    id SERIAL PRIMARY KEY,
                    version TEXT NOT NULL,
                    name TEXT NOT NULL,
                    applied_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
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
        conn: Database connection object
        config: Optional configuration dictionary
        migrations_dir: Directory containing migration scripts

    Returns:
        Dictionary of migration status information
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

    status = {}
    for version, script_path in available_migrations:
        script_name = os.path.basename(script_path)
        name = script_name.replace('.py', '')
        
        is_applied = version in applied_versions
        applied_at = applied_versions[version]['applied_at'] if is_applied else None
        
        status[version] = {
            'version': version,
            'name': name,
            'script': script_path,
            'applied': is_applied,
            'applied_at': applied_at
        }
    
    # Add any applied migrations that don't have scripts
    for version, migration in applied_versions.items():
        if version not in status:
            status[version] = {
                'version': version,
                'name': migration['name'],
                'script': None,
                'applied': True,
                'applied_at': migration['applied_at']
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
        conn: Database connection object
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
        if not info['applied'] and info['script'] is not None and
           (target_version is None or version <= target_version)
    ]
    pending_migrations.sort(key=lambda x: x[0])  # Sort by version

    if not pending_migrations:
        logger.info("No pending migrations to apply")
        return []

    # Apply migrations
    applied = []
    for version, script_path in pending_migrations:
        script_name = os.path.basename(script_path)
        name = script_name.replace('.py', '')

        logger.info(f"Applying migration {version}: {name}")

        # Load migration script
        try:
            module = _load_migration_script(script_path)

            # Check if module has apply function
            if not hasattr(module, 'apply'):
                logger.error(f"Migration {version} does not have an apply function")
                break

            # Apply migration
            if not dry_run:
                try:
                    module.apply(conn, environment)
                    record_migration(conn, version, name)

                    migration_info = {
                        'version': version,
                        'name': name,
                        'script': script_path,
                        'status': 'applied'
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
                    'script': script_path,
                    'status': 'would_apply'
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
        conn: Database connection object
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
        if info['applied'] and info['script'] is not None and
           (target_version is None or version > target_version)
    ]
    applied_migrations.sort(key=lambda x: x[0], reverse=True)  # Sort by version, descending

    if not applied_migrations:
        logger.info("No migrations to roll back")
        return []

    # Roll back migrations
    rolled_back = []
    for version, script_path in applied_migrations:
        script_name = os.path.basename(script_path)
        name = script_name.replace('.py', '')

        logger.info(f"Rolling back migration {version}: {name}")

        # Load migration script
        try:
            module = _load_migration_script(script_path)

            # Check if module has rollback function
            if not hasattr(module, 'rollback'):
                logger.error(f"Migration {version} does not have a rollback function")
                break

            # Roll back migration
            if not dry_run:
                try:
                    module.rollback(conn, environment)
                    remove_migration_record(conn, version)

                    migration_info = {
                        'version': version,
                        'name': name,
                        'script': script_path,
                        'status': 'rolled_back'
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
                    'script': script_path,
                    'status': 'would_roll_back'
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
    parser = argparse.ArgumentParser(description='Manage CryoProtect database migrations')
    parser.add_argument(
        'action',
        choices=['apply', 'rollback', 'status'],
        help='Migration action to perform'
    )
    parser.add_argument(
        '--target',
        help='Target migration version'
    )
    parser.add_argument(
        '--env',
        default='development',
        help='Target environment (development, staging, production)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Perform a dry run without making changes'
    )
    parser.add_argument(
        '--format',
        choices=['text', 'json'],
        default='text',
        help='Output format'
    )
    parser.add_argument(
        '--dir',
        help='Directory containing migration scripts'
    )

    args = parser.parse_args()

    # Create database connection
    conn = create_connection()

    if args.action == 'status':
        status = get_migration_status(conn, migrations_dir=args.dir)
        if args.format == 'json':
            import json
            print(json.dumps(status, indent=2, default=str))
        else:
            print("\nMigration Status:\n")
            for version, info in sorted(status.items()):
                applied = "Applied" if info['applied'] else "Pending"
                applied_at = f" at {info['applied_at']}" if info['applied'] else ""
                print(f"{version}: {info['name']} - {applied}{applied_at}")
            print()

    elif args.action == 'apply':
        results = apply_migrations(
            conn,
            migrations_dir=args.dir,
            target_version=args.target,
            environment=args.env,
            dry_run=args.dry_run
        )
        if args.format == 'json':
            import json
            print(json.dumps(results, indent=2, default=str))
        else:
            if args.dry_run:
                print("\nDry run - would apply the following migrations:\n")
            else:
                print("\nApplied migrations:\n")
            
            for info in results:
                print(f"{info['version']}: {info['name']} - {info['status']}")
            print()

    elif args.action == 'rollback':
        results = rollback_migrations(
            conn,
            migrations_dir=args.dir,
            target_version=args.target,
            environment=args.env,
            dry_run=args.dry_run
        )
        if args.format == 'json':
            import json
            print(json.dumps(results, indent=2, default=str))
        else:
            if args.dry_run:
                print("\nDry run - would roll back the following migrations:\n")
            else:
                print("\nRolled back migrations:\n")
            
            for info in results:
                print(f"{info['version']}: {info['name']} - {info['status']}")
            print()

if __name__ == '__main__':
    main()
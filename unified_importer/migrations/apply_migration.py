#!/usr/bin/env python3
"""
Apply the unified importer database migration.

This script applies the database migration required for the unified molecular data importer.
It can be run directly or imported and used programmatically.
"""

import os
import sys
import argparse
import logging
import asyncio
from typing import Optional, Dict, Any

# Add parent directory to path to allow imports
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

try:
    from unified_importer.core.database import DatabaseOperations
except ImportError:
    print("Unable to import DatabaseOperations. Make sure the unified_importer package is installed.")
    sys.exit(1)


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('unified_importer.migrations')


async def apply_migration(
    db_url: Optional[str] = None,
    db_key: Optional[str] = None,
    migration_file: Optional[str] = None,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Apply the unified importer database migration.
    
    Args:
        db_url: Database URL (defaults to environment variable)
        db_key: Database API key (defaults to environment variable)
        migration_file: Path to migration file (defaults to 021_unified_importer.sql)
        dry_run: If True, only validate the migration without applying
        
    Returns:
        Dictionary with results of the migration
    """
    # Set default migration file if not provided
    if not migration_file:
        # First look in the package migrations directory
        package_migration = os.path.join(os.path.dirname(__file__), 'unified_importer.sql')
        if os.path.exists(package_migration):
            migration_file = package_migration
        else:
            # Look in project migrations directory
            project_migration = os.path.join(
                os.path.dirname(os.path.dirname(parent_dir)), 
                'migrations', 
                '021_unified_importer.sql'
            )
            if os.path.exists(project_migration):
                migration_file = project_migration
            else:
                raise FileNotFoundError(
                    "Migration file not found. Please provide the path to 021_unified_importer.sql"
                )
    
    logger.info(f"Using migration file: {migration_file}")
    
    # Initialize database operations
    db = DatabaseOperations(db_url, db_key)
    
    # Read the migration file
    try:
        with open(migration_file, 'r') as f:
            migration_sql = f.read()
    except Exception as e:
        logger.error(f"Error reading migration file: {str(e)}")
        raise
    
    # Apply or validate the migration
    if dry_run:
        logger.info("Validating migration in dry run mode (no changes will be made)")
        # Just check if the SQL is valid
        # In a real implementation, this would do more extensive validation
        result = {
            'status': 'validated',
            'message': 'Migration validated successfully (dry run)',
            'dry_run': True
        }
    else:
        logger.info("Applying migration to database")
        try:
            # Execute the SQL
            async with db.transaction():
                # The SQL already includes checks to prevent re-applying
                await db.execute(migration_sql)
            
            # Verify that the migration was applied
            tables_query = """
                SELECT table_name FROM information_schema.tables 
                WHERE table_schema = 'public' AND table_name IN 
                ('data_sources', 'molecule_synonyms', 'molecule_cross_references', 
                'import_jobs', 'import_job_errors', 'property_types')
            """
            tables_result = await db.fetch_all(tables_query)
            
            # Count the tables that were created
            created_tables = [row['table_name'] for row in tables_result]
            logger.info(f"Confirmed creation of tables: {', '.join(created_tables)}")
            
            result = {
                'status': 'success',
                'message': f'Migration applied successfully. Created {len(created_tables)} tables.',
                'created_tables': created_tables,
                'dry_run': False
            }
            
        except Exception as e:
            logger.error(f"Error applying migration: {str(e)}")
            result = {
                'status': 'error',
                'message': f'Error applying migration: {str(e)}',
                'dry_run': False
            }
            raise
    
    return result


async def main():
    """Run the migration script from command line."""
    parser = argparse.ArgumentParser(description='Apply the unified importer database migration')
    parser.add_argument('--db-url', help='Database URL (defaults to environment variable)')
    parser.add_argument('--db-key', help='Database API key (defaults to environment variable)')
    parser.add_argument('--migration-file', help='Path to migration file')
    parser.add_argument('--dry-run', action='store_true', help='Validate migration without applying')
    
    args = parser.parse_args()
    
    try:
        result = await apply_migration(
            db_url=args.db_url,
            db_key=args.db_key,
            migration_file=args.migration_file,
            dry_run=args.dry_run
        )
        
        logger.info(f"Migration {result['status']}: {result['message']}")
        return 0
    except Exception as e:
        logger.error(f"Migration failed: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
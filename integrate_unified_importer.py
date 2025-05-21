#!/usr/bin/env python3
"""
Integration script for the unified molecular data importer.

This script applies the database migration and creates a basic integration
of the unified importer with the main CryoProtect application.
"""

import os
import sys
import argparse
import logging
import asyncio
from typing import Dict, Any, List, Optional

from unified_importer.core.config import ImporterConfig
from unified_importer.core.logging import setup_logging
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker, ConsoleProgressReporter
from unified_importer.core.database import DatabaseOperations
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.sources.chembl_source import ChEMBLDataSource


async def apply_migration(db_url: str, db_key: str) -> None:
    """
    Apply the unified importer database migration.
    
    Args:
        db_url: Supabase URL
        db_key: Supabase key
    """
    logger = setup_logging(
        'cryoprotect.integration.migration',
        log_level='INFO',
        log_file=f"logs/migration_{os.path.basename(__file__).split('.')[0]}.log"
    )
    
    logger.info("Applying database migration for unified importer")
    
    try:
        # Set up database connection
        db_operations = DatabaseOperations(
            connection_type="supabase",
            connection_params={
                'url': db_url,
                'key': db_key
            },
            logger=logger
        )
        
        # Read migration SQL
        migration_path = os.path.join(
            os.path.dirname(__file__),
            'migrations',
            '021_unified_importer.sql'
        )
        
        with open(migration_path, 'r') as f:
            migration_sql = f.read()
        
        # Execute migration
        logger.info("Executing migration SQL")
        await db_operations.execute_query(migration_sql)
        
        logger.info("Migration applied successfully")
    except Exception as e:
        logger.error(f"Error applying migration: {str(e)}")
        raise


async def integrate_with_app(db_url: str, db_key: str) -> None:
    """
    Integrate the unified importer with the main application.
    
    Args:
        db_url: Supabase URL
        db_key: Supabase key
    """
    logger = setup_logging(
        'cryoprotect.integration.app',
        log_level='INFO',
        log_file=f"logs/app_integration_{os.path.basename(__file__).split('.')[0]}.log"
    )
    
    logger.info("Integrating unified importer with the main application")
    
    try:
        # Create integration directory
        os.makedirs('integration', exist_ok=True)
        
        # Create integration module
        integration_path = os.path.join('integration', 'unified_importer.py')
        
        with open(integration_path, 'w') as f:
            f.write("""'''
Unified importer integration module.

This module provides integration between the unified importer
and the main CryoProtect application.
'''

import os
import logging
from typing import Dict, Any, List, Optional, Union

from unified_importer.core.config import ImporterConfig
from unified_importer.core.logging import setup_logging
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker
from unified_importer.core.database import DatabaseOperations
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.sources.chembl_source import ChEMBLDataSource


class MoleculeImporter:
    '''
    Integration class for the unified molecular importer.
    
    This class provides a simplified interface for importing
    molecules from various sources into the CryoProtect database.
    '''
    
    def __init__(
        self,
        db_url: str,
        db_key: str,
        log_level: str = 'INFO',
        checkpoint_dir: str = 'checkpoints',
        log_dir: str = 'logs'
    ):
        '''
        Initialize the molecule importer.
        
        Args:
            db_url: Supabase URL
            db_key: Supabase key
            log_level: Logging level
            checkpoint_dir: Directory for checkpoint files
            log_dir: Directory for log files
        '''
        # Create directories
        os.makedirs(checkpoint_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)
        
        # Set up logger
        self.logger = setup_logging(
            'cryoprotect.importer',
            log_level=log_level,
            log_file=os.path.join(log_dir, 'molecule_importer.log')
        )
        
        # Set up database connection
        self.db_operations = DatabaseOperations(
            connection_type="supabase",
            connection_params={
                'url': db_url,
                'key': db_key
            },
            logger=self.logger
        )
        
        # Store configuration
        self.checkpoint_dir = checkpoint_dir
        self.log_dir = log_dir
    
    async def import_pubchem_compounds(
        self,
        compound_ids: List[str],
        checkpoint_file: Optional[str] = None,
        dry_run: bool = False
    ) -> Dict[str, Any]:
        '''
        Import compounds from PubChem.
        
        Args:
            compound_ids: List of PubChem CIDs
            checkpoint_file: Path to checkpoint file (optional)
            dry_run: Whether to run in dry run mode
            
        Returns:
            Dictionary with import results
        '''
        # Set up checkpoint manager
        if not checkpoint_file:
            import time
            checkpoint_file = os.path.join(
                self.checkpoint_dir,
                f"pubchem_import_{time.strftime('%Y%m%d_%H%M%S')}.json"
            )
        
        checkpoint_manager = CheckpointManager(
            checkpoint_file=checkpoint_file,
            logger=self.logger
        )
        
        # Set up progress tracker
        progress_tracker = ProgressTracker(
            total_items=len(compound_ids),
            logger=self.logger
        )
        
        # Create PubChem data source
        source = PubChemDataSource(
            db_operations=self.db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            config={'dry_run': dry_run},
            logger=self.logger
        )
        
        # Import compounds
        success, failure, failures = await source.import_compounds(compound_ids)
        
        # Return results
        return {
            'success_count': success,
            'failure_count': failure,
            'failures': [
                {'id': cid, 'reason': reason}
                for cid, reason in failures
            ],
            'checkpoint_file': checkpoint_file
        }
    
    async def import_chembl_compounds(
        self,
        compound_ids: List[str],
        checkpoint_file: Optional[str] = None,
        dry_run: bool = False
    ) -> Dict[str, Any]:
        '''
        Import compounds from ChEMBL.
        
        Args:
            compound_ids: List of ChEMBL IDs
            checkpoint_file: Path to checkpoint file (optional)
            dry_run: Whether to run in dry run mode
            
        Returns:
            Dictionary with import results
        '''
        # Set up checkpoint manager
        if not checkpoint_file:
            import time
            checkpoint_file = os.path.join(
                self.checkpoint_dir,
                f"chembl_import_{time.strftime('%Y%m%d_%H%M%S')}.json"
            )
        
        checkpoint_manager = CheckpointManager(
            checkpoint_file=checkpoint_file,
            logger=self.logger
        )
        
        # Set up progress tracker
        progress_tracker = ProgressTracker(
            total_items=len(compound_ids),
            logger=self.logger
        )
        
        # Create ChEMBL data source
        source = ChEMBLDataSource(
            db_operations=self.db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            config={'dry_run': dry_run},
            logger=self.logger
        )
        
        # Import compounds
        success, failure, failures = await source.import_compounds(compound_ids)
        
        # Return results
        return {
            'success_count': success,
            'failure_count': failure,
            'failures': [
                {'id': cid, 'reason': reason}
                for cid, reason in failures
            ],
            'checkpoint_file': checkpoint_file
        }
    
    async def search_and_import(
        self,
        source: str,
        query: str,
        max_results: Optional[int] = None,
        checkpoint_file: Optional[str] = None,
        dry_run: bool = False
    ) -> Dict[str, Any]:
        '''
        Search for compounds and import them.
        
        Args:
            source: Data source ('pubchem' or 'chembl')
            query: Search query
            max_results: Maximum number of results to import
            checkpoint_file: Path to checkpoint file (optional)
            dry_run: Whether to run in dry run mode
            
        Returns:
            Dictionary with import results
        '''
        # Set up checkpoint manager
        if not checkpoint_file:
            import time
            checkpoint_file = os.path.join(
                self.checkpoint_dir,
                f"{source}_search_{time.strftime('%Y%m%d_%H%M%S')}.json"
            )
        
        checkpoint_manager = CheckpointManager(
            checkpoint_file=checkpoint_file,
            logger=self.logger
        )
        
        # Set up progress tracker
        progress_tracker = ProgressTracker(logger=self.logger)
        
        # Create data source
        if source.lower() == 'pubchem':
            data_source = PubChemDataSource(
                db_operations=self.db_operations,
                checkpoint_manager=checkpoint_manager,
                progress_tracker=progress_tracker,
                config={'dry_run': dry_run},
                logger=self.logger
            )
        elif source.lower() == 'chembl':
            data_source = ChEMBLDataSource(
                db_operations=self.db_operations,
                checkpoint_manager=checkpoint_manager,
                progress_tracker=progress_tracker,
                config={'dry_run': dry_run},
                logger=self.logger
            )
        else:
            raise ValueError(f"Unsupported data source: {source}")
        
        # Get estimated count
        count = await data_source.get_compound_count(query)
        if count > 0:
            self.logger.info(f"Found approximately {count} matching compounds")
            limit = min(count, max_results) if max_results else count
            progress_tracker.total_items = limit
        
        # Search and import
        success, failure, failures = await data_source.search_and_import(
            query=query,
            max_results=max_results
        )
        
        # Return results
        return {
            'success_count': success,
            'failure_count': failure,
            'failures': [
                {'id': cid, 'reason': reason}
                for cid, reason in failures
            ],
            'checkpoint_file': checkpoint_file,
            'total_found': count
        }


# Decorator for database operations
def import_molecules(func):
    '''
    Decorator for molecule import operations.
    
    This decorator wraps import functions to handle common
    tasks like connecting to the database and logging.
    '''
    async def wrapper(*args, **kwargs):
        # Get database configuration from environment
        db_url = os.environ.get('SUPABASE_URL')
        db_key = os.environ.get('SUPABASE_KEY')
        
        if not db_url or not db_key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY environment variables must be set")
        
        # Create importer
        importer = MoleculeImporter(db_url, db_key)
        
        # Call the wrapped function
        return await func(importer, *args, **kwargs)
    
    return wrapper


# Example usage as a decorator
@import_molecules
async def import_aspirin(importer):
    '''Import aspirin from PubChem.'''
    return await importer.import_pubchem_compounds(['2244'])


# Example usage as a function
async def import_cryoprotectants(identifiers):
    '''Import multiple cryoprotectants.'''
    # Get database configuration from environment
    db_url = os.environ.get('SUPABASE_URL')
    db_key = os.environ.get('SUPABASE_KEY')
    
    if not db_url or not db_key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY environment variables must be set")
    
    # Create importer
    importer = MoleculeImporter(db_url, db_key)
    
    # Import compounds
    return await importer.import_pubchem_compounds(identifiers)
""")
        
        logger.info(f"Created integration module at {integration_path}")
        
        # Create __init__.py file
        init_path = os.path.join('integration', '__init__.py')
        
        with open(init_path, 'w') as f:
            f.write("""'''
Integration modules for CryoProtect.

This package contains integration modules for various
components of the CryoProtect application.
'''

from .unified_importer import MoleculeImporter, import_molecules

__all__ = ['MoleculeImporter', 'import_molecules']
""")
        
        logger.info(f"Created integration __init__.py at {init_path}")
        
        # Create example usage script
        example_path = os.path.join('integration', 'import_example.py')
        
        with open(example_path, 'w') as f:
            f.write("""#!/usr/bin/env python3
'''
Example usage of the unified importer integration.

This script demonstrates how to use the unified importer
integration with the main CryoProtect application.
'''

import os
import sys
import asyncio
import argparse
from typing import List, Optional

from unified_importer.core.logging import setup_logging
from integration.unified_importer import MoleculeImporter


async def main():
    '''Main entry point.'''
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Example usage of unified importer")
    
    parser.add_argument(
        "--mode",
        choices=["pubchem", "chembl", "search"],
        required=True,
        help="Import mode"
    )
    
    parser.add_argument(
        "--ids",
        help="Comma-separated list of compound IDs"
    )
    
    parser.add_argument(
        "--query",
        help="Search query (for search mode)"
    )
    
    parser.add_argument(
        "--source",
        choices=["pubchem", "chembl"],
        help="Data source for search mode"
    )
    
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Maximum number of compounds to import"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate without inserting data"
    )
    
    args = parser.parse_args()
    
    # Set up logger
    logger = setup_logging('example', log_level='INFO')
    
    # Get database configuration from environment
    db_url = os.environ.get('SUPABASE_URL')
    db_key = os.environ.get('SUPABASE_KEY')
    
    if not db_url or not db_key:
        logger.error("SUPABASE_URL and SUPABASE_KEY environment variables must be set")
        return 1
    
    # Create importer
    importer = MoleculeImporter(db_url, db_key)
    
    try:
        if args.mode == "pubchem":
            if not args.ids:
                logger.error("--ids is required for pubchem mode")
                return 1
                
            # Import from PubChem
            ids = args.ids.split(',')
            logger.info(f"Importing {len(ids)} compounds from PubChem")
            
            result = await importer.import_pubchem_compounds(
                ids,
                dry_run=args.dry_run
            )
            
            logger.info(f"Import completed: {result['success_count']} successful, {result['failure_count']} failed")
            
        elif args.mode == "chembl":
            if not args.ids:
                logger.error("--ids is required for chembl mode")
                return 1
                
            # Import from ChEMBL
            ids = args.ids.split(',')
            logger.info(f"Importing {len(ids)} compounds from ChEMBL")
            
            result = await importer.import_chembl_compounds(
                ids,
                dry_run=args.dry_run
            )
            
            logger.info(f"Import completed: {result['success_count']} successful, {result['failure_count']} failed")
            
        elif args.mode == "search":
            if not args.query:
                logger.error("--query is required for search mode")
                return 1
                
            if not args.source:
                logger.error("--source is required for search mode")
                return 1
                
            # Search and import
            logger.info(f"Searching for compounds matching '{args.query}' in {args.source}")
            
            result = await importer.search_and_import(
                args.source,
                args.query,
                max_results=args.limit,
                dry_run=args.dry_run
            )
            
            logger.info(f"Found {result.get('total_found', 0)} matching compounds")
            logger.info(f"Import completed: {result['success_count']} successful, {result['failure_count']} failed")
    
    except Exception as e:
        logger.error(f"Error in import process: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
""")
        
        logger.info(f"Created example usage script at {example_path}")
        os.chmod(example_path, 0o755)  # Make executable
        
        logger.info("Integration completed successfully")
    except Exception as e:
        logger.error(f"Error integrating with application: {str(e)}")
        raise


async def main():
    """Main entry point."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Integrate unified molecular importer")
    
    parser.add_argument(
        "--apply-migration",
        action="store_true",
        help="Apply database migration"
    )
    
    parser.add_argument(
        "--integrate",
        action="store_true",
        help="Integrate with main application"
    )
    
    parser.add_argument(
        "--db-url",
        help="Supabase URL"
    )
    
    parser.add_argument(
        "--db-key",
        help="Supabase key"
    )
    
    args = parser.parse_args()
    
    # Get database configuration from environment if not provided
    db_url = args.db_url or os.environ.get('SUPABASE_URL')
    db_key = args.db_key or os.environ.get('SUPABASE_KEY')
    
    if (args.apply_migration or args.integrate) and (not db_url or not db_key):
        print("Error: Database URL and key are required for migration and integration")
        print("  Provide them as arguments or set SUPABASE_URL and SUPABASE_KEY environment variables")
        return 1
    
    # Create logs directory
    os.makedirs('logs', exist_ok=True)
    
    # Apply migration if requested
    if args.apply_migration:
        try:
            await apply_migration(db_url, db_key)
        except Exception as e:
            print(f"Error applying migration: {str(e)}")
            return 1
    
    # Integrate with application if requested
    if args.integrate:
        try:
            await integrate_with_app(db_url, db_key)
        except Exception as e:
            print(f"Error integrating with application: {str(e)}")
            return 1
    
    # If no action specified, show help
    if not args.apply_migration and not args.integrate:
        parser.print_help()
        return 0
    
    print("Integration completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
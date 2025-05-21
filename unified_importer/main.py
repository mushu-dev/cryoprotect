#!/usr/bin/env python3
"""
Unified Molecular Data Importer

Main entry point for the unified importer system that integrates
multiple chemical data sources into a standardized database.
"""

import os
import sys
import json
import argparse
import logging
import asyncio
import time
from typing import Dict, Any, List, Optional, Union, Set, Tuple

from .core.config import load_config
from .core.database import DatabaseOperations
from .core.enhanced_database import EnhancedDatabaseOperations  # Enhanced version with advanced connection pooling
from .core.enhanced_batch_processing import AdaptiveBatchProcessor, ChemicalDataBatchStrategy
from .core.batch_integration import BatchProcessingManager
from .core.checkpoint import CheckpointManager
from .core.progress import ProgressTracker, ConsoleProgressReporter
from .core.validation import MoleculeValidator

from .sources.source_base import MolecularDataSource
from .sources.chembl_source import ChEMBLDataSource
from .sources.pubchem_source import PubChemDataSource

from .transforms.molecule_transform import MoleculeTransformer
from .transforms.property_transform import PropertyTransformer


class MolecularImporter:
    """
    Main class for the unified molecular data importer.
    
    This class coordinates the import process, manages data sources,
    and provides a high-level interface for the import functionality.
    """
    
    def __init__(
        self,
        config: Optional[Dict[str, Any]] = None,
        config_file: Optional[str] = None,
        db_url: Optional[str] = None,
        db_key: Optional[str] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the molecular importer.
        
        Args:
            config: Configuration dictionary
            config_file: Path to configuration file
            db_url: Database URL (overrides config)
            db_key: Database API key (overrides config)
            logger: Logger instance
        """
        # Load configuration
        self.config = load_config(config, config_file)
        
        # Set up logging
        self.logger = logger or self._setup_logging()
        
        # Initialize database connection
        db_url = db_url or self.config.get('database', {}).get('url')
        db_key = db_key or self.config.get('database', {}).get('key')
        
        # Configure database with pool settings from config
        db_config = self.config.get('database', {})

        # Map configuration parameters to what DatabaseOperations expects
        connection_type = "supabase" if db_url else "direct"
        connection_params = {
            "url": db_url,
            "key": db_key,
            "host": db_config.get('host', 'localhost'),
            "port": db_config.get('port', 5432),
            "database": db_config.get('database', 'postgres'),
            "user": db_config.get('user', 'postgres'),
            "password": db_config.get('password', '')
        }

        # Use enhanced connection pool if configured
        use_enhanced_pool = self.config.get('use_enhanced_connection_pool', True)

        if use_enhanced_pool:
            # Use the enhanced database operations with advanced connection pooling
            self.db = EnhancedDatabaseOperations(
                connection_type=connection_type,
                connection_params=connection_params,
                batch_size=db_config.get('batch_size', 100),
                max_retries=db_config.get('max_retries', 3),
                retry_delay=db_config.get('retry_delay', 2.0),
                pool_min_size=db_config.get('pool_min_size', db_config.get('pool_size', 5)),
                pool_max_size=db_config.get('pool_max_size', db_config.get('pool_size', 10) * 2),
                pool_target_utilization=db_config.get('pool_target_utilization', 0.7),
                pool_timeout=db_config.get('pool_timeout', 30.0),
                pool_max_lifetime=db_config.get('pool_max_lifetime', 3600.0),
                pool_validation_interval=db_config.get('pool_validation_interval', 300.0),
                pool_health_check_interval=db_config.get('pool_health_check_interval', 60.0),
                logger=self.logger
            )
            self.logger.info("Using enhanced database operations with advanced connection pooling")
        else:
            # Use the standard database operations
            self.db = DatabaseOperations(
                connection_type=connection_type,
                connection_params=connection_params,
                batch_size=db_config.get('batch_size', 100),
                max_retries=db_config.get('max_retries', 3),
                retry_delay=db_config.get('retry_delay', 2.0),
                pool_min_size=db_config.get('pool_min_size', db_config.get('pool_size', 2)),
                pool_max_size=db_config.get('pool_max_size', db_config.get('pool_size', 10) * 2),
                pool_timeout=db_config.get('pool_timeout', 30.0),
                max_pool_connections=db_config.get('max_pool_connections', 20),
                logger=self.logger
            )
            self.logger.info("Using standard database operations")
        
        # Initialize transformer instances
        transformer_config = self.config.get('transforms', {})
        self.molecule_transformer = MoleculeTransformer(
            logger=self.logger,
            config=transformer_config.get('molecule_transform', {})
        )

        self.property_transformer = PropertyTransformer(
            logger=self.logger,
            config=transformer_config.get('property_transform', {})
        )

        # Set up validator
        self.validator = MoleculeValidator(logger=self.logger)

        # Configure batch processing
        batch_config = self.config.get('batch_processing', {})
        self.batch_manager = BatchProcessingManager(
            db=self.db,
            max_workers=batch_config.get('max_workers'),
            adaptive_batch_sizing=batch_config.get('adaptive_batch_sizing', True),
            parallel_processing=batch_config.get('parallel_processing', True),
            logger=self.logger
        )
        self.logger.info(
            f"Initialized BatchProcessingManager with "
            f"adaptive_batch_sizing={batch_config.get('adaptive_batch_sizing', True)}, "
            f"parallel_processing={batch_config.get('parallel_processing', True)}"
        )

        # Initialize data sources
        self.sources = {}
        self._initialize_sources()
    
    def _setup_logging(self) -> logging.Logger:
        """
        Set up logging configuration.
        
        Returns:
            Configured logger
        """
        log_config = self.config.get('logging', {})
        log_level = getattr(logging, log_config.get('level', 'INFO'))
        log_format = log_config.get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        log_file = log_config.get('file')
        
        # Configure logger
        logger = logging.getLogger('unified_importer')
        logger.setLevel(log_level)
        
        # Remove existing handlers
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
        
        # Create console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(console_handler)
        
        # Create file handler if specified
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(log_level)
            file_handler.setFormatter(logging.Formatter(log_format))
            logger.addHandler(file_handler)
        
        return logger
    
    def _initialize_sources(self) -> None:
        """Initialize all configured data sources."""
        sources_config = self.config.get('sources', {})

        # Initialize ChEMBL source if configured
        if 'chembl' in sources_config:
            self.sources['chembl'] = ChEMBLDataSource(
                db_operations=self.db,
                config=sources_config['chembl'],
                validator=self.validator,
                batch_manager=self.batch_manager,  # Add batch manager to source
                logger=self.logger
            )
            self.logger.info("Initialized ChEMBL data source")

        # Initialize PubChem source if configured
        if 'pubchem' in sources_config:
            self.sources['pubchem'] = PubChemDataSource(
                db_operations=self.db,
                config=sources_config['pubchem'],
                validator=self.validator,
                batch_manager=self.batch_manager,  # Add batch manager to source
                logger=self.logger
            )
            self.logger.info("Initialized PubChem data source")
    
    def get_source(self, source_name: str) -> Optional[MolecularDataSource]:
        """
        Get a data source by name.
        
        Args:
            source_name: Name of the data source
            
        Returns:
            Data source instance or None if not found
        """
        source = self.sources.get(source_name.lower())
        if not source:
            self.logger.warning(f"Data source '{source_name}' not found")
        return source
    
    async def import_from_source(
        self,
        source_name: str,
        query: Optional[str] = None,
        identifiers: Optional[List[str]] = None,
        max_results: Optional[int] = None,
        checkpoint_file: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Import data from a specific source.
        
        Args:
            source_name: Name of the data source
            query: Search query
            identifiers: List of specific identifiers to import
            max_results: Maximum number of results to import
            checkpoint_file: Path to checkpoint file
            
        Returns:
            Import results dictionary
        """
        # Get the data source
        source = self.get_source(source_name)
        if not source:
            return {
                'status': 'error',
                'message': f"Data source '{source_name}' not found",
                'success_count': 0,
                'failure_count': 0
            }
        
        # Set up checkpoint manager if specified
        if checkpoint_file:
            source.checkpoint_manager = CheckpointManager(
                checkpoint_file,
                backup_interval=self.config.get('checkpoints', {}).get('backup_interval', 5),
                logger=self.logger
            )
        
        # Set up progress tracker
        source.progress_tracker = ProgressTracker(
            window_size=20,
            logger=self.logger
        )
        
        # Add console reporter
        reporter = ConsoleProgressReporter(
            source.progress_tracker,
            update_interval=5.0,
            logger=self.logger
        )
        
        start_time = time.time()
        
        try:
            # Import based on provided parameters
            if identifiers:
                self.logger.info(f"Importing {len(identifiers)} specific identifiers from {source_name}")
                success_count, failure_count, failures = await source.import_compounds(identifiers)
            elif query:
                self.logger.info(f"Searching for '{query}' in {source_name}")
                success_count, failure_count, failures = await source.search_and_import(query, max_results)
            else:
                self.logger.info(f"Importing all compounds from {source_name} (limit: {max_results})")
                success_count, failure_count, failures = await source.import_all(max_results)
            
            # Generate results
            elapsed_time = time.time() - start_time
            
            # Save progress report if progress tracker is available
            report_path = None
            if source.progress_tracker:
                report_path = f"{source_name}_import_report_{time.strftime('%Y%m%d_%H%M%S')}.json"
                source.progress_tracker.save_report(report_path)
            
            results = {
                'status': 'success',
                'message': f"Imported {success_count} compounds from {source_name}",
                'success_count': success_count,
                'failure_count': failure_count,
                'elapsed_seconds': elapsed_time,
                'report_file': report_path
            }
            
            if failures:
                # Limit the number of failures in the results
                max_failures = 10
                failures_to_report = failures[:max_failures]
                
                results['failures'] = [
                    {'id': failure[0], 'reason': failure[1]}
                    for failure in failures_to_report
                ]
                
                if len(failures) > max_failures:
                    results['additional_failures'] = len(failures) - max_failures
            
            return results
            
        except Exception as e:
            self.logger.exception(f"Error importing from {source_name}: {str(e)}")
            return {
                'status': 'error',
                'message': f"Error importing from {source_name}: {str(e)}",
                'success_count': 0,
                'failure_count': 0
            }
        finally:
            # Clean up
            if reporter:
                reporter.stop()

            # Get batch processing stats if available
            if source and hasattr(source, 'batch_manager') and source.batch_manager:
                try:
                    batch_stats = source.batch_manager.get_stats()
                    self.logger.info(f"Batch processing stats for {source_name}: {json.dumps(batch_stats, indent=2)}")
                except Exception as e:
                    self.logger.warning(f"Error getting batch stats: {str(e)}")
    
    async def import_from_multiple_sources(
        self,
        sources: List[str],
        query: Optional[str] = None,
        identifiers: Optional[Dict[str, List[str]]] = None,
        max_results: Optional[int] = None,
        checkpoint_dir: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Import data from multiple sources.
        
        Args:
            sources: List of source names
            query: Search query (applied to all sources)
            identifiers: Dictionary mapping source names to identifier lists
            max_results: Maximum results per source
            checkpoint_dir: Directory for checkpoint files
            
        Returns:
            Combined import results
        """
        if not sources:
            self.logger.warning("No sources specified for import")
            return {
                'status': 'error',
                'message': 'No sources specified',
                'results': {}
            }
        
        # Create checkpoint directory if it doesn't exist
        if checkpoint_dir and not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
        
        # Import from each source
        all_results = {}
        for source_name in sources:
            # Get identifiers for this source
            source_identifiers = None
            if identifiers and source_name in identifiers:
                source_identifiers = identifiers[source_name]
            
            # Create checkpoint file
            checkpoint_file = None
            if checkpoint_dir:
                checkpoint_file = os.path.join(
                    checkpoint_dir,
                    f"{source_name}_checkpoint_{time.strftime('%Y%m%d_%H%M%S')}.json"
                )
            
            # Import from this source
            results = await self.import_from_source(
                source_name,
                query=query,
                identifiers=source_identifiers,
                max_results=max_results,
                checkpoint_file=checkpoint_file
            )
            
            all_results[source_name] = results
        
        # Combine results
        total_success = sum(results.get('success_count', 0) for results in all_results.values())
        total_failure = sum(results.get('failure_count', 0) for results in all_results.values())
        
        return {
            'status': 'success',
            'message': f"Imported {total_success} compounds from {len(sources)} sources",
            'total_success_count': total_success,
            'total_failure_count': total_failure,
            'results': all_results
        }
    
    async def standardize_molecule(
        self,
        molecule_data: Dict[str, Any],
        resolve_ids: bool = True
    ) -> Dict[str, Any]:
        """
        Standardize a molecule using the configured transformers.
        
        Args:
            molecule_data: Molecule data dictionary
            resolve_ids: Whether to resolve cross-database IDs
            
        Returns:
            Standardized molecule data
        """
        return await self.molecule_transformer.standardize_molecule(
            molecule_data,
            resolve_ids=resolve_ids
        )
    
    def run_import(
        self, 
        sources: Optional[List[str]] = None,
        query: Optional[str] = None,
        identifiers: Optional[Dict[str, List[str]]] = None,
        limit: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Run the import process synchronously.
        
        This is a convenience method that runs the async import
        in a new event loop. For more control, use the async methods directly.
        
        Args:
            sources: List of source names
            query: Search query
            identifiers: Dictionary mapping source names to identifier lists
            limit: Maximum results per source
            
        Returns:
            Import results
        """
        # Use all sources if not specified
        if not sources:
            sources = list(self.sources.keys())
        
        # Get checkpoint directory from config
        checkpoint_dir = self.config.get('checkpoints', {}).get('directory', 'checkpoints')
        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
        
        # Run the import
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
        
        results = loop.run_until_complete(
            self.import_from_multiple_sources(
                sources=sources,
                query=query,
                identifiers=identifiers,
                max_results=limit,
                checkpoint_dir=checkpoint_dir
            )
        )
        
        return results


async def main():
    """Main entry point when run as a script."""
    parser = argparse.ArgumentParser(description='Unified Molecular Data Importer')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--sources', help='Comma-separated list of sources (pubchem,chembl)')
    parser.add_argument('--query', help='Search query')
    parser.add_argument('--identifiers', help='Comma-separated list of identifiers to import')
    parser.add_argument('--limit', type=int, help='Maximum number of results per source')
    parser.add_argument('--checkpoint-dir', help='Directory for checkpoint files')
    parser.add_argument('--db-url', help='Database URL')
    parser.add_argument('--db-key', help='Database API key')
    parser.add_argument('--log-level', default='INFO', help='Logging level (DEBUG, INFO, WARNING, ERROR)')
    parser.add_argument('--log-file', help='Log file path')
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = getattr(logging, args.log_level)
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=args.log_file
    )
    logger = logging.getLogger('unified_importer')
    
    # Create importer
    config_file = args.config or os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'config',
        'config_example.json'
    )
    
    try:
        importer = MolecularImporter(
            config_file=config_file if os.path.exists(config_file) else None,
            db_url=args.db_url,
            db_key=args.db_key,
            logger=logger
        )
        
        # Parse sources
        sources = args.sources.split(',') if args.sources else list(importer.sources.keys())
        
        # Parse identifiers
        identifiers = None
        if args.identifiers:
            # If identifiers are provided but no specific source, use the first source
            if len(sources) == 1:
                identifiers = {sources[0]: args.identifiers.split(',')}
        
        # Run import
        results = await importer.import_from_multiple_sources(
            sources=sources,
            query=args.query,
            identifiers=identifiers,
            max_results=args.limit,
            checkpoint_dir=args.checkpoint_dir
        )
        
        # Print results
        print(json.dumps(results, indent=2))
        return 0
        
    except Exception as e:
        logger.exception(f"Error running importer: {str(e)}")
        print(f"Error: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
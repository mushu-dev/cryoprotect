#!/usr/bin/env python3
"""
Unified Molecular Data Importer Example

This script demonstrates how to use the unified importer framework
to import data from PubChem and ChEMBL.
"""

import os
import sys
import asyncio
import argparse
import logging
import time
from typing import Dict, Any, List, Optional

# Add the project directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from unified_importer.core.config import ImporterConfig
from unified_importer.core.logging import setup_logging
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker, ConsoleProgressReporter
from unified_importer.core.database import DatabaseOperations
from unified_importer.core.validation import MoleculeValidator
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.sources.chembl_source import ChEMBLDataSource


async def import_from_pubchem(chemids: List[str], db_config: Dict[str, Any]) -> None:
    """
    Import compounds from PubChem.
    
    Args:
        chemids: List of PubChem CIDs to import
        db_config: Database configuration
    """
    # Set up logger
    logger = setup_logging(
        'example.pubchem',
        log_level='INFO',
        log_file=f"logs/pubchem_example_{time.strftime('%Y%m%d_%H%M%S')}.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(
        checkpoint_file=f"checkpoints/pubchem_example_{time.strftime('%Y%m%d_%H%M%S')}.json",
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(
        total_items=len(chemids),
        logger=logger
    )
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=2.0,
        logger=logger
    )
    
    # Create database operations
    db_operations = DatabaseOperations(
        connection_type=db_config.get('type', 'direct'),
        connection_params=db_config.get('connection', {}),
        logger=logger
    )
    
    # Create validator
    validator = MoleculeValidator(logger=logger)
    
    # Create PubChem data source
    source = PubChemDataSource(
        db_operations=db_operations,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker,
        validator=validator,
        config={
            'batch_size': 10,
            'api_delay': 1.0,
            'max_retries': 3,
            'retry_delay': 2.0,
            'dry_run': True  # Set to False to actually insert data
        },
        logger=logger
    )
    
    logger.info(f"Starting import of {len(chemids)} compounds from PubChem")
    
    try:
        # Import compounds
        success, failure, failures = await source.import_compounds(chemids)
        
        # Print results
        logger.info(f"Import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for chemid, reason in failures:
                logger.info(f"  - {chemid}: {reason}")
    finally:
        # Stop the reporter
        reporter.stop()


async def import_from_chembl(chemblids: List[str], db_config: Dict[str, Any]) -> None:
    """
    Import compounds from ChEMBL.
    
    Args:
        chemblids: List of ChEMBL IDs to import
        db_config: Database configuration
    """
    # Set up logger
    logger = setup_logging(
        'example.chembl',
        log_level='INFO',
        log_file=f"logs/chembl_example_{time.strftime('%Y%m%d_%H%M%S')}.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(
        checkpoint_file=f"checkpoints/chembl_example_{time.strftime('%Y%m%d_%H%M%S')}.json",
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(
        total_items=len(chemblids),
        logger=logger
    )
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=2.0,
        logger=logger
    )
    
    # Create database operations
    db_operations = DatabaseOperations(
        connection_type=db_config.get('type', 'direct'),
        connection_params=db_config.get('connection', {}),
        logger=logger
    )
    
    # Create validator
    validator = MoleculeValidator(logger=logger)
    
    # Create ChEMBL data source
    source = ChEMBLDataSource(
        db_operations=db_operations,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker,
        validator=validator,
        config={
            'batch_size': 10,
            'api_delay': 1.0,
            'max_retries': 3,
            'retry_delay': 2.0,
            'dry_run': True  # Set to False to actually insert data
        },
        logger=logger
    )
    
    logger.info(f"Starting import of {len(chemblids)} compounds from ChEMBL")
    
    try:
        # Import compounds
        success, failure, failures = await source.import_compounds(chemblids)
        
        # Print results
        logger.info(f"Import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for chemblid, reason in failures:
                logger.info(f"  - {chemblid}: {reason}")
    finally:
        # Stop the reporter
        reporter.stop()


async def search_and_import(
    source_name: str,
    query: str,
    max_results: int,
    db_config: Dict[str, Any]
) -> None:
    """
    Search for compounds and import the results.
    
    Args:
        source_name: Data source name ('pubchem' or 'chembl')
        query: Search query
        max_results: Maximum number of results to import
        db_config: Database configuration
    """
    # Set up logger
    logger = setup_logging(
        f'example.search.{source_name}',
        log_level='INFO',
        log_file=f"logs/{source_name}_search_{time.strftime('%Y%m%d_%H%M%S')}.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(
        checkpoint_file=f"checkpoints/{source_name}_search_{time.strftime('%Y%m%d_%H%M%S')}.json",
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(logger=logger)
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=2.0,
        logger=logger
    )
    
    # Create database operations
    db_operations = DatabaseOperations(
        connection_type=db_config.get('type', 'direct'),
        connection_params=db_config.get('connection', {}),
        logger=logger
    )
    
    # Create validator
    validator = MoleculeValidator(logger=logger)
    
    # Create appropriate data source
    if source_name == 'pubchem':
        source = PubChemDataSource(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            validator=validator,
            config={
                'batch_size': 10,
                'api_delay': 1.0,
                'max_retries': 3,
                'retry_delay': 2.0,
                'dry_run': True  # Set to False to actually insert data
            },
            logger=logger
        )
    elif source_name == 'chembl':
        source = ChEMBLDataSource(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            validator=validator,
            config={
                'batch_size': 10,
                'api_delay': 1.0,
                'max_retries': 3,
                'retry_delay': 2.0,
                'dry_run': True  # Set to False to actually insert data
            },
            logger=logger
        )
    else:
        logger.error(f"Unsupported data source: {source_name}")
        return
    
    logger.info(f"Starting search and import from {source_name} for query: {query}")
    
    try:
        # Get estimated count
        count = await source.get_compound_count(query)
        if count > 0:
            logger.info(f"Found approximately {count} matching compounds")
            progress_tracker.total_items = min(count, max_results)
        
        # Search and import
        success, failure, failures = await source.search_and_import(
            query=query,
            max_results=max_results
        )
        
        # Print results
        logger.info(f"Import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for compound_id, reason in failures[:10]:  # Show only first 10
                logger.info(f"  - {compound_id}: {reason}")
                
            if len(failures) > 10:
                logger.info(f"  ... and {len(failures) - 10} more")
    finally:
        # Stop the reporter
        reporter.stop()


def main() -> None:
    """Main entry point for the example script."""
    parser = argparse.ArgumentParser(description="Unified Molecular Data Importer Example")
    
    parser.add_argument(
        "--mode",
        choices=["pubchem", "chembl", "search"],
        required=True,
        help="Import mode"
    )
    
    parser.add_argument(
        "--ids",
        help="Comma-separated list of compound IDs to import"
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
        default=100,
        help="Maximum number of compounds to import"
    )
    
    args = parser.parse_args()
    
    # Use dry run mode to test without database connection
    db_config = {
        'type': 'mock',  # Use mock database type to avoid actual DB connection
        'connection': {},
        'dry_run': True  # Enable dry run mode for testing
    }
    
    if args.mode == "pubchem":
        if not args.ids:
            print("Error: --ids is required for pubchem mode")
            return
            
        chemids = args.ids.split(',')
        asyncio.run(import_from_pubchem(chemids, db_config))
    elif args.mode == "chembl":
        if not args.ids:
            print("Error: --ids is required for chembl mode")
            return
            
        chemblids = args.ids.split(',')
        asyncio.run(import_from_chembl(chemblids, db_config))
    elif args.mode == "search":
        if not args.query:
            print("Error: --query is required for search mode")
            return
            
        if not args.source:
            print("Error: --source is required for search mode")
            return
            
        asyncio.run(search_and_import(args.source, args.query, args.limit, db_config))


if __name__ == "__main__":
    main()
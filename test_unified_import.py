#!/usr/bin/env python3
"""
Test the unified molecular data importer with real data.

This script tests importing a small set of compounds from PubChem and ChEMBL
using the unified importer to verify its functionality with real data.
"""

import os
import sys
import asyncio
import logging
from typing import Dict, List, Any, Optional

from unified_importer.core.config import ImporterConfig
from unified_importer.core.logging import setup_logging
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker, ConsoleProgressReporter
from unified_importer.core.database import DatabaseOperations
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.sources.chembl_source import ChEMBLDataSource


# Common cryoprotectants to test with
PUBCHEM_TEST_COMPOUNDS = [
    '702', # Ethylene glycol (basic cryoprotectant)
    '753', # Glycerol (common cryoprotectant)
    '6059', # DMSO (dimethyl sulfoxide, widely used cryoprotectant)
    '1030', # Propylene glycol (used in some cryopreservation protocols)
    '6568', # Trehalose (natural cryoprotectant used in some organisms)
]

CHEMBL_TEST_COMPOUNDS = [
    'CHEMBL422679', # Glycerol
    'CHEMBL1512', # DMSO
    'CHEMBL91403', # Ethylene glycol
    'CHEMBL550', # Propylene glycol
    'CHEMBL2106876' # Trehalose
]


async def test_pubchem_import(dry_run: bool = True) -> None:
    """
    Test importing compounds from PubChem.
    
    Args:
        dry_run: Whether to run in dry run mode without database writes
    """
    logger = setup_logging(
        'test.pubchem',
        log_level='INFO',
        log_file="logs/pubchem_test.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    logger.info(f"Testing PubChem import with compounds: {PUBCHEM_TEST_COMPOUNDS}")
    
    # Create checkpoint manager
    checkpoint_file = "checkpoints/pubchem_test.json"
    checkpoint_manager = CheckpointManager(
        checkpoint_file=checkpoint_file,
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(
        total_items=len(PUBCHEM_TEST_COMPOUNDS),
        logger=logger
    )
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=1.0,
        logger=logger
    )
    
    try:
        # Create database connection (mock)
        logger.info("Setting up mock database connection for testing")
        db_operations = DatabaseOperations(
            connection_type="mock",
            logger=logger
        )
        
        # Create PubChem data source
        source = PubChemDataSource(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            config={
                'batch_size': 2,
                'api_delay': 1.0,
                'max_retries': 2,
                'retry_delay': 1.0,
                'dry_run': dry_run
            },
            logger=logger
        )
        
        # Import compounds
        success, failure, failures = await source.import_compounds(PUBCHEM_TEST_COMPOUNDS)
        
        # Log results
        logger.info(f"PubChem import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for compound_id, reason in failures:
                logger.info(f"  - {compound_id}: {reason}")
        
        # Return success status
        return success > 0 and failure == 0
        
    except Exception as e:
        logger.error(f"Error in PubChem import test: {str(e)}")
        return False
    finally:
        # Stop the reporter
        reporter.stop()
        

async def test_chembl_import(dry_run: bool = True) -> None:
    """
    Test importing compounds from ChEMBL.
    
    Args:
        dry_run: Whether to run in dry run mode without database writes
    """
    logger = setup_logging(
        'test.chembl',
        log_level='INFO',
        log_file="logs/chembl_test.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    logger.info(f"Testing ChEMBL import with compounds: {CHEMBL_TEST_COMPOUNDS}")
    
    # Create checkpoint manager
    checkpoint_file = "checkpoints/chembl_test.json"
    checkpoint_manager = CheckpointManager(
        checkpoint_file=checkpoint_file,
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(
        total_items=len(CHEMBL_TEST_COMPOUNDS),
        logger=logger
    )
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=1.0,
        logger=logger
    )
    
    try:
        # Create database connection (mock)
        logger.info("Setting up mock database connection for testing")
        db_operations = DatabaseOperations(
            connection_type="mock",
            logger=logger
        )
        
        # Create ChEMBL data source
        source = ChEMBLDataSource(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            config={
                'batch_size': 2,
                'api_delay': 1.0,
                'max_retries': 2,
                'retry_delay': 1.0,
                'dry_run': dry_run
            },
            logger=logger
        )
        
        # Import compounds
        success, failure, failures = await source.import_compounds(CHEMBL_TEST_COMPOUNDS)
        
        # Log results
        logger.info(f"ChEMBL import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for compound_id, reason in failures:
                logger.info(f"  - {compound_id}: {reason}")
        
        # Return success status
        return success > 0 and failure == 0
        
    except Exception as e:
        logger.error(f"Error in ChEMBL import test: {str(e)}")
        return False
    finally:
        # Stop the reporter
        reporter.stop()


async def test_search_import(dry_run: bool = True) -> None:
    """
    Test search and import functionality.
    
    Args:
        dry_run: Whether to run in dry run mode without database writes
    """
    logger = setup_logging(
        'test.search',
        log_level='INFO',
        log_file="logs/search_test.log"
    )
    
    # Ensure directories exist
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)
    
    # Search query for cryoprotectants
    query = "cryoprotectant"
    max_results = 5
    
    logger.info(f"Testing search and import with query: '{query}' (max {max_results} results)")
    
    # Create checkpoint manager
    checkpoint_file = "checkpoints/search_test.json"
    checkpoint_manager = CheckpointManager(
        checkpoint_file=checkpoint_file,
        logger=logger
    )
    
    # Create progress tracker
    progress_tracker = ProgressTracker(logger=logger)
    
    # Set up console reporter
    reporter = ConsoleProgressReporter(
        tracker=progress_tracker,
        update_interval=1.0,
        logger=logger
    )
    
    try:
        # Create database connection (mock)
        logger.info("Setting up mock database connection for testing")
        db_operations = DatabaseOperations(
            connection_type="mock",
            logger=logger
        )
        
        # Create PubChem data source for searching
        source = PubChemDataSource(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            config={
                'batch_size': 2,
                'api_delay': 1.0,
                'max_retries': 2,
                'retry_delay': 1.0,
                'dry_run': dry_run
            },
            logger=logger
        )
        
        # Get estimated count
        try:
            count = await source.get_compound_count(query)
            if count > 0:
                logger.info(f"Found approximately {count} matching compounds")
                progress_tracker.total_items = min(count, max_results)
            else:
                logger.warning(f"No compounds found matching query: {query}")
        except Exception as e:
            logger.warning(f"Could not get compound count: {str(e)}")
        
        # Search and import
        success, failure, failures = await source.search_and_import(
            query=query,
            max_results=max_results
        )
        
        # Log results
        logger.info(f"Search and import completed: {success} successful, {failure} failed")
        
        if failures:
            logger.info("Failed compounds:")
            for compound_id, reason in failures:
                logger.info(f"  - {compound_id}: {reason}")
        
        # Return success status
        return success > 0
        
    except Exception as e:
        logger.error(f"Error in search and import test: {str(e)}")
        return False
    finally:
        # Stop the reporter
        reporter.stop()


async def run_tests() -> bool:
    """
    Run all tests.
    
    Returns:
        True if all tests passed, False otherwise
    """
    logger = setup_logging('test.main', log_level='INFO')
    logger.info("Starting unified importer tests")
    
    try:
        # Run PubChem import test
        logger.info("=== Testing PubChem Import ===")
        pubchem_success = await test_pubchem_import(dry_run=True)
        
        # Run ChEMBL import test
        logger.info("=== Testing ChEMBL Import ===")
        chembl_success = await test_chembl_import(dry_run=True)
        
        # Run search import test
        logger.info("=== Testing Search Import ===")
        search_success = await test_search_import(dry_run=True)
        
        # Check results
        all_passed = pubchem_success and chembl_success and search_success
        
        if all_passed:
            logger.info("✅ All tests passed successfully")
        else:
            logger.error("❌ Some tests failed")
            if not pubchem_success:
                logger.error("  - PubChem import test failed")
            if not chembl_success:
                logger.error("  - ChEMBL import test failed")
            if not search_success:
                logger.error("  - Search import test failed")
        
        return all_passed
        
    except Exception as e:
        logger.error(f"Error running tests: {str(e)}")
        return False


if __name__ == "__main__":
    success = asyncio.run(run_tests())
    sys.exit(0 if success else 1)
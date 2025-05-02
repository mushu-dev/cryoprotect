#!/usr/bin/env python3
"""
Example script demonstrating how to use the enhanced ChEMBL import script
with direct PostgreSQL connection.

This example shows how to:
1. Configure the direct PostgreSQL connection
2. Import a small set of compounds from ChEMBL
3. Use batch processing for efficiency
4. Handle checkpointing for resumable operations
5. Verify the imported data
"""

import os
import sys
import logging
import argparse
import time
from datetime import datetime

# Add parent directory to path to import modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import required modules
from postgres_direct import PostgresDirectConnection
from sql_executor import execute_query, get_db, get_connection_metrics
from property_utils import PropertyManager
from chembl.client import ChEMBLClient
from chembl.checkpoint import CheckpointManager
from chembl.worker import ChEMBLWorker

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_example.log')
    ]
)
logger = logging.getLogger(__name__)

def setup_environment():
    """
    Set up the environment for the example.
    
    This function checks if the required environment variables are set
    and creates necessary directories.
    """
    # Check if environment variables are set
    required_vars = [
        'SUPABASE_DB_HOST',
        'SUPABASE_DB_PORT',
        'SUPABASE_DB_NAME',
        'SUPABASE_DB_USER',
        'SUPABASE_DB_PASSWORD'
    ]
    
    missing_vars = [var for var in required_vars if not os.getenv(var)]
    if missing_vars:
        logger.error(f"Missing required environment variables: {', '.join(missing_vars)}")
        logger.info("Please set these variables in a .env file or in your environment")
        sys.exit(1)
    
    # Create necessary directories
    os.makedirs('logs', exist_ok=True)
    os.makedirs('checkpoints', exist_ok=True)

def verify_database_connection():
    """
    Verify the direct PostgreSQL connection.
    
    Returns:
        bool: True if connection is successful, False otherwise
    """
    logger.info("Verifying direct PostgreSQL connection...")
    try:
        # Get database connection
        db = get_db()
        
        # Test connection with a simple query
        result = execute_query("SELECT 1 as connection_test")
        if result and result[0]['connection_test'] == 1:
            logger.info("Successfully connected to PostgreSQL database")
            
            # Get connection metrics
            metrics = get_connection_metrics()
            logger.info(f"Connection metrics: {metrics}")
            
            return True
        else:
            logger.warning("Database connection test failed")
            return False
    except Exception as e:
        logger.error(f"Database connection verification failed: {e}")
        return False

def import_example_compounds(search_terms, batch_size=10, num_workers=2):
    """
    Import a small set of example compounds from ChEMBL.
    
    Args:
        search_terms: List of search terms to use
        batch_size: Number of compounds to process in each batch
        num_workers: Number of worker threads to use
        
    Returns:
        dict: Statistics about the import process
    """
    logger.info(f"Starting example import with {len(search_terms)} search terms")
    logger.info(f"Using batch size {batch_size} and {num_workers} workers")
    
    # Initialize ChEMBL client
    client = ChEMBLClient()
    
    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager('./checkpoints', prefix='chembl_example')
    
    # Check for existing checkpoint
    if checkpoint_manager.load_checkpoint():
        logger.info("Resuming from existing checkpoint")
    else:
        logger.info("Starting fresh import")
    
    # Get compound IDs to process
    compound_ids = []
    for term in search_terms:
        logger.info(f"Searching for compounds with term: {term}")
        results = client.search_compounds(term, limit=5)  # Limit to 5 compounds per term for example
        logger.info(f"Found {len(results)} compounds for term '{term}'")
        compound_ids.extend(results)
    
    # Deduplicate
    compound_ids = list(set(compound_ids))
    logger.info(f"Total unique compounds after deduplication: {len(compound_ids)}")
    
    # Skip already processed compounds if resuming
    if checkpoint_manager.state.get("processed_compounds"):
        processed = checkpoint_manager.state.get("processed_compounds", [])
        skipped = 0
        new_compound_ids = []
        for cid in compound_ids:
            if cid in processed:
                skipped += 1
            else:
                new_compound_ids.append(cid)
        compound_ids = new_compound_ids
        logger.info(f"Skipping {skipped} already processed compounds")
    
    # Store configuration in checkpoint
    checkpoint_manager.state["config"] = {
        "search_terms": search_terms,
        "batch_size": batch_size,
        "num_workers": num_workers,
        "total_compounds": len(compound_ids),
        "start_time": datetime.now().isoformat()
    }
    
    # Save initial checkpoint
    checkpoint_manager.save_checkpoint()
    
    # Set up worker pool
    from queue import Queue
    task_queue = Queue()
    result_queue = Queue()
    
    # Create workers
    workers = []
    for i in range(num_workers):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        workers.append(worker)
    
    # Start workers
    for worker in workers:
        worker.start()
    
    logger.info(f"Started {len(workers)} workers")
    
    # Add tasks to queue
    for compound_id in compound_ids:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": False,  # Set to True to skip database operations
            "checkpoint_manager": checkpoint_manager,
            "batch_size": batch_size
        })
    
    logger.info(f"Added {len(compound_ids)} compounds to task queue")
    
    # Statistics
    stats = {
        "total": len(compound_ids),
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": datetime.now(),
        "errors": {},
        "batch_times": [],
        "search_terms": search_terms
    }
    
    # Process results
    try:
        # Enhanced batch processing
        total_processed = 0
        max_batches = (len(compound_ids) + batch_size - 1) // batch_size
        
        logger.info(f"Processing {len(compound_ids)} compounds in up to {max_batches} batches")
        
        while total_processed < len(compound_ids):
            batch_results = []
            batch_start_time = time.time()
            
            # Process a batch of results
            for _ in range(min(batch_size, len(compound_ids) - total_processed)):
                try:
                    result = result_queue.get(timeout=60)  # 1 minute timeout per compound
                    batch_results.append(result)
                    result_queue.task_done()
                except Exception as e:
                    logger.warning(f"Error getting result: {e}")
                    break
            
            # Process all results in this batch
            for result in batch_results:
                # Update statistics
                stats["processed"] += 1
                
                if result.get("status") == "success":
                    stats["success"] += 1
                else:
                    stats["error"] += 1
                    error_category = result.get("error_category", "UNKNOWN")
                    
                    # Track error by category
                    if error_category not in stats["errors"]:
                        stats["errors"][error_category] = 0
                    stats["errors"][error_category] += 1
            
            # Calculate batch statistics
            batch_duration = time.time() - batch_start_time
            stats["batch_times"].append(batch_duration)
            
            compounds_per_second = len(batch_results) / batch_duration if batch_duration > 0 else 0
            logger.info(f"Batch completed: {len(batch_results)} compounds in {batch_duration:.2f}s " +
                      f"({compounds_per_second:.2f} compounds/s)")
            
            # Save checkpoint after each batch
            checkpoint_manager.save_checkpoint()
            
            # Update total processed
            total_processed += len(batch_results)
            
            if len(batch_results) == 0:
                logger.warning("No compounds processed in last batch, breaking")
                break
    
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
    finally:
        # Stop workers
        for worker in workers:
            task_queue.put(None)  # Send shutdown signal
        
        for worker in workers:
            worker.stop()
    
    # Calculate final statistics
    elapsed_seconds = (datetime.now() - stats["start_time"]).total_seconds()
    compounds_per_second = stats["processed"] / elapsed_seconds if elapsed_seconds > 0 else 0
    success_rate = stats["success"] / stats["processed"] * 100 if stats["processed"] > 0 else 0
    
    # Log summary
    logger.info(f"Import Summary:")
    logger.info(f"  Total compounds: {stats['total']}")
    logger.info(f"  Processed: {stats['processed']}")
    logger.info(f"  Success: {stats['success']} ({success_rate:.1f}%)")
    logger.info(f"  Errors: {stats['error']}")
    logger.info(f"  Total time: {elapsed_seconds:.2f} seconds")
    logger.info(f"  Performance: {compounds_per_second:.2f} compounds/second")
    
    # Log error categories if any
    if stats["errors"]:
        logger.info("  Error categories:")
        for category, count in stats["errors"].items():
            logger.info(f"    {category}: {count}")
    
    return stats

def verify_imported_data(compound_ids):
    """
    Verify that the compounds were imported correctly.
    
    Args:
        compound_ids: List of ChEMBL IDs to verify
        
    Returns:
        bool: True if verification is successful, False otherwise
    """
    logger.info(f"Verifying {len(compound_ids)} imported compounds")
    
    try:
        # Query the database to check if compounds exist
        placeholders = ', '.join(['%s'] * len(compound_ids))
        query = f"""
            SELECT chembl_id, name, formula, molecular_weight
            FROM molecules
            WHERE chembl_id IN ({placeholders})
        """
        
        result = execute_query(query, compound_ids)
        
        if not result:
            logger.warning("No compounds found in database")
            return False
        
        logger.info(f"Found {len(result)} compounds in database")
        
        # Check if all compounds were imported
        found_ids = [row['chembl_id'] for row in result]
        missing_ids = [cid for cid in compound_ids if cid not in found_ids]
        
        if missing_ids:
            logger.warning(f"Missing {len(missing_ids)} compounds in database: {missing_ids}")
        else:
            logger.info("All compounds were imported successfully")
        
        # Check if properties were imported
        property_manager = PropertyManager()
        
        # Get molecule IDs
        molecule_ids = [row['id'] for row in result if 'id' in row]
        
        if molecule_ids:
            # Get properties for a sample of molecules
            sample_ids = molecule_ids[:min(5, len(molecule_ids))]
            properties = property_manager.batch_get_properties(sample_ids)
            
            if properties:
                logger.info(f"Found properties for {len(properties)} molecules")
                
                # Log a sample of properties
                for molecule_id, props in list(properties.items())[:2]:
                    logger.info(f"Properties for molecule {molecule_id}: {len(props)} properties")
                    if props:
                        logger.info(f"Sample properties: {list(props.keys())[:5]}")
            else:
                logger.warning("No properties found for imported molecules")
        
        return True
    
    except Exception as e:
        logger.error(f"Error verifying imported data: {e}")
        return False

def main():
    """Main function for the example script."""
    parser = argparse.ArgumentParser(description="Example script for ChEMBL import with direct PostgreSQL connection")
    parser.add_argument("--search-terms", nargs="+", default=["cryoprotectant", "antifreeze"],
                        help="Search terms for finding compounds")
    parser.add_argument("--batch-size", type=int, default=5,
                        help="Number of compounds per batch")
    parser.add_argument("--workers", type=int, default=2,
                        help="Number of worker threads")
    args = parser.parse_args()
    
    # Set up environment
    setup_environment()
    
    # Verify database connection
    if not verify_database_connection():
        logger.error("Failed to connect to database. Exiting.")
        return 1
    
    # Import example compounds
    stats = import_example_compounds(
        search_terms=args.search_terms,
        batch_size=args.batch_size,
        num_workers=args.workers
    )
    
    # Verify imported data if any compounds were successfully imported
    if stats["success"] > 0:
        # Get the list of successfully imported compounds from the checkpoint
        checkpoint_manager = CheckpointManager('./checkpoints', prefix='chembl_example')
        if checkpoint_manager.load_checkpoint():
            successful_compounds = []
            for compound_id, status in checkpoint_manager.state.get("compound_status", {}).items():
                if status.get("success"):
                    successful_compounds.append(compound_id)
            
            if successful_compounds:
                verify_imported_data(successful_compounds)
    
    logger.info("Example completed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())
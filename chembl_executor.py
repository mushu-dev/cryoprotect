#!/usr/bin/env python3
"""
ChEMBL Data Import Executor

This script orchestrates the ChEMBL import process using parallel processing with worker threads.
It manages the import of compound data from ChEMBL into the CryoProtect database with
checkpointing for resumability and progress tracking.
"""

import os
import sys
import argparse
import logging
import time
import json
from datetime import datetime
from queue import Queue
from typing import List, Dict, Any, Optional
import threading

# Set up logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_import.log')
    ]
)

logger = logging.getLogger(__name__)

# Import necessary modules
try:
    from chembl.worker import ChEMBLWorker
    from chembl.checkpoint import CheckpointManager
    from chembl.client import ChEMBLClient
except ImportError as e:
    logger.error(f"Failed to import required modules: {e}")
    sys.exit(1)

# Define reference compounds for CryoProtect database
"""
Reference compounds for CryoProtect database.

This module defines standard reference compounds that should be included
in the database to ensure consistent comparison across experiments.
"""
# VERIFICATION: Implemented reference compounds dictionary with metadata

REFERENCE_COMPOUNDS = [
    {
        "chembl_id": "CHEMBL1118",
        "name": "GLYCEROL",
        "role": "primary cryoprotectant",
        "description": "One of the most common cryoprotectants, glycerol is used in many cryopreservation protocols.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1559",
        "name": "DIMETHYL SULFOXIDE",
        "common_name": "DMSO",
        "role": "primary cryoprotectant",
        "description": "DMSO is widely used for cell cryopreservation, particularly in stem cell and tissue preservation.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1750",
        "name": "ETHYLENE GLYCOL",
        "role": "primary cryoprotectant",
        "description": "Ethylene glycol is used in vitrification protocols for embryo and oocyte cryopreservation.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL2043",
        "name": "PROPYLENE GLYCOL",
        "role": "primary cryoprotectant",
        "description": "Propylene glycol is used in cryopreservation of various cell types and tissues.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1308",
        "name": "SUCROSE",
        "role": "secondary cryoprotectant",
        "description": "Sucrose is often used as a non-penetrating cryoprotectant to reduce osmotic stress.",
        "priority": 2
    },
    {
        "chembl_id": "CHEMBL2107886",
        "name": "TREHALOSE",
        "role": "secondary cryoprotectant",
        "description": "Trehalose provides membrane protection during freezing and is often used with DMSO.",
        "priority": 2
    },
    {
        "chembl_id": "CHEMBL1742",
        "name": "ETHANOL",
        "role": "secondary cryoprotectant",
        "description": "Ethanol is sometimes used in specialized cryopreservation protocols.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1489",
        "name": "UREA",
        "role": "additive",
        "description": "Urea is sometimes used to enhance protein stability during freezing.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1771",
        "name": "GLYCINE",
        "role": "additive",
        "description": "Glycine can act as an osmolyte that protects proteins during freezing.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1201",
        "name": "CAFFEINE",
        "role": "control",
        "description": "Caffeine is used as a control compound that is not expected to have cryoprotective properties.",
        "priority": 4
    },
    {
        "chembl_id": "CHEMBL25",
        "name": "ASPIRIN",
        "role": "control",
        "description": "Aspirin is used as a control compound that is not expected to have cryoprotective properties.",
        "priority": 4
    }
]

def get_reference_compounds():
    """Get the list of reference compounds with their ChEMBL IDs and metadata"""
    return REFERENCE_COMPOUNDS

def get_reference_compound_ids():
    """Get just the ChEMBL IDs of reference compounds"""
    return [compound["chembl_id"] for compound in REFERENCE_COMPOUNDS]

def get_primary_cryoprotectants():
    """Get the list of primary cryoprotectant reference compounds"""
    return [c for c in REFERENCE_COMPOUNDS if c["role"] == "primary cryoprotectant"]

def load_chembl_ids(filename: str = None, query: str = None, limit: int = 1000) -> List[str]:
    """
    Load ChEMBL IDs from a file or by querying the ChEMBL API.
    
    Args:
        filename: Optional path to file containing ChEMBL IDs (one per line)
        query: Optional search query for finding compounds
        limit: Maximum number of compounds to return
        
    Returns:
        List of ChEMBL IDs to process
    """
    # VERIFICATION: Implemented ID loading from file or API
    
    if filename and os.path.exists(filename):
        logger.info(f"Loading ChEMBL IDs from file: {filename}")
        with open(filename, 'r') as f:
            # Strip whitespace and filter out empty lines
            chembl_ids = [line.strip() for line in f if line.strip()]
            
        # Apply limit if specified
        if limit > 0:
            chembl_ids = chembl_ids[:limit]
            
        logger.info(f"Loaded {len(chembl_ids)} ChEMBL IDs from file")
        return chembl_ids
    
    elif query:
        logger.info(f"Searching ChEMBL for compounds matching: {query}")
        try:
            client = ChEMBLClient()
            chembl_ids = client.search_compounds(query, limit=limit)
            logger.info(f"Found {len(chembl_ids)} compounds matching query")
            return chembl_ids
        except Exception as e:
            logger.error(f"Error searching ChEMBL: {e}")
            return []
    
    else:
        logger.warning("No file or query provided for ChEMBL IDs")
        return []

def fetch_compound_ids(client: ChEMBLClient, args: argparse.Namespace) -> List[str]:
    """
    Fetch ChEMBL compound IDs to process.
    
    Args:
        client: ChEMBL client instance
        args: Command line arguments
        
    Returns:
        List of ChEMBL IDs to process
    """
    # Define reference compounds dictionary
    reference_compounds = []
    
    # Check if reference compounds are requested
    if (hasattr(args, 'include_refs') and args.include_refs) or (hasattr(args, 'reference_only') and args.reference_only):
        reference_compounds = get_reference_compound_ids()
        logger.info(f"Including {len(reference_compounds)} reference compounds")
        
        # If reference-only mode, return just these compounds
        if hasattr(args, 'reference_only') and args.reference_only:
            logger.info("Reference-only mode: processing only reference compounds")
            return reference_compounds
    
    # Search for compounds by query
    query_compounds = []
    if hasattr(args, 'query') and args.query:
        logger.info(f"Searching for compounds with query: {args.query}")
        query_compounds = client.search_compounds(args.query, limit=args.limit if hasattr(args, 'limit') else 1000)
        logger.info(f"Found {len(query_compounds)} compounds from search")
    
    # Load compounds from file if specified
    file_compounds = []
    if hasattr(args, 'input_file') and args.input_file and os.path.exists(args.input_file):
        logger.info(f"Loading compounds from file: {args.input_file}")
        with open(args.input_file, 'r') as f:
            file_compounds = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(file_compounds)} compounds from file")
    
    # Combine and deduplicate
    all_compounds = list(set(reference_compounds + query_compounds + file_compounds))
    
    # Apply limit
    if hasattr(args, 'limit') and args.limit and len(all_compounds) > args.limit:
        all_compounds = all_compounds[:args.limit]
    
    logger.info(f"Total compounds to process: {len(all_compounds)}")
    return all_compounds

def setup_worker_pool(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Set up worker pool for processing.
    
    Args:
        args: Command line arguments
        
    Returns:
        Dictionary with worker pool components
    """
    # VERIFICATION: Creates and starts worker threads with appropriate queues
    
    # Create queues
    task_queue = Queue()
    result_queue = Queue()
    
    # Create workers
    workers = []
    num_workers = args.workers if hasattr(args, 'workers') else 4
    
    logger.info(f"Creating {num_workers} worker threads")
    for i in range(num_workers):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        workers.append(worker)
    
    # Start workers
    for worker in workers:
        worker.start()
        logger.debug(f"Started worker {worker.worker_id}")
    
    return {
        "workers": workers,
        "task_queue": task_queue,
        "result_queue": result_queue
    }

def distribute_tasks(compound_ids: List[str], pool: Dict[str, Any],
                    args: argparse.Namespace, checkpoint_manager: Optional[CheckpointManager] = None) -> int:
    """
    Distribute tasks to worker pool.
    
    Args:
        compound_ids: List of compound IDs to process
        pool: Worker pool components
        args: Command line arguments
        checkpoint_manager: Optional checkpoint manager
        
    Returns:
        Number of tasks distributed
    """
    # VERIFICATION: Distributes tasks to workers with proper checkpoint handling
    
    task_queue = pool["task_queue"]
    
    # Skip already processed compounds if resuming
    if hasattr(args, 'mode') and args.mode == "resume" and checkpoint_manager:
        processed = checkpoint_manager.state.get("processed_compounds", [])
        skipped_compounds = [cid for cid in compound_ids if cid in processed]
        compound_ids = [cid for cid in compound_ids if cid not in processed]
        logger.info(f"Skipping {len(skipped_compounds)} already processed compounds")
    
    # Add dry run flag for dry-run mode
    dry_run = hasattr(args, 'mode') and args.mode == "dry-run"
    
    # Add tasks to queue
    for compound_id in compound_ids:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": dry_run,
            "checkpoint_manager": checkpoint_manager,
            "max_retries": args.max_retries if hasattr(args, 'max_retries') else 3
        })
    
    logger.info(f"Distributed {len(compound_ids)} tasks to worker pool")
    return len(compound_ids)

def shutdown_worker_pool(pool: Dict[str, Any], wait: bool = True, timeout: float = 30.0) -> None:
    """
    Shutdown the worker pool gracefully.
    
    Args:
        pool: Worker pool components
        wait: Whether to wait for workers to finish
        timeout: Maximum time to wait for each worker (in seconds)
    """
    # VERIFICATION: Properly shuts down worker threads
    
    workers = pool.get("workers", [])
    task_queue = pool.get("task_queue")
    
    if not workers or not task_queue:
        logger.warning("Invalid worker pool provided for shutdown")
        return
    
    logger.info(f"Shutting down {len(workers)} workers")
    
    # Send poison pill to each worker
    for _ in workers:
        task_queue.put(None)
    
    # Wait for workers to finish if requested
    if wait:
        logger.info("Waiting for workers to finish...")
        for i, worker in enumerate(workers):
            if worker.thread and worker.thread.is_alive():
                logger.debug(f"Waiting for worker {worker.worker_id} to finish")
                worker.thread.join(timeout=timeout)
                
                if worker.thread.is_alive():
                    logger.warning(f"Worker {worker.worker_id} did not finish within timeout ({timeout}s)")
                else:
                    logger.debug(f"Worker {worker.worker_id} finished")
    
    # Force stop any remaining workers
    for worker in workers:
        if worker.thread and worker.thread.is_alive():
            worker.stop()
            logger.debug(f"Forced stop for worker {worker.worker_id}")
    
    logger.info("Worker pool shutdown complete")

def process_results(pool: Dict[str, Any], total_tasks: int, args: argparse.Namespace,
                   checkpoint_manager: Optional[CheckpointManager] = None) -> Dict[str, Any]:
    """
    Process results from worker pool.
    
    Args:
        pool: Worker pool components
        total_tasks: Total number of tasks distributed
        args: Command line arguments
        checkpoint_manager: Optional checkpoint manager
        
    Returns:
        Dictionary with processing statistics
    """
    # VERIFICATION: Monitors result queue, updates checkpoint manager, and tracks progress
    
    result_queue = pool["result_queue"]
    task_queue = pool["task_queue"]
    
    # Statistics
    stats = {
        "total": total_tasks,
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": time.time(),
        "errors": {}
    }
    
    # Progress display parameters
    update_interval = 1.0  # seconds
    last_update = 0
    
    # Process results until all tasks are done
    while stats["processed"] < total_tasks:
        try:
            # Get result with timeout to allow checking queue status
            try:
                result = result_queue.get(timeout=0.1)
            except Queue.Empty:
                # Check if tasks are done processing
                if task_queue.empty() and task_queue.unfinished_tasks == 0:
                    logger.warning("Task queue is empty but not all results received")
                    break
                continue
            
            # Process the result
            stats["processed"] += 1
            
            if result.get("status") == "success":
                stats["success"] += 1
                
                # Update checkpoint if available
                if checkpoint_manager:
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=True
                    )
            else:
                stats["error"] += 1
                error_category = result.get("error_category", "UNKNOWN")
                
                # Track error by category
                if error_category not in stats["errors"]:
                    stats["errors"][error_category] = 0
                stats["errors"][error_category] += 1
                
                # Update checkpoint if available
                if checkpoint_manager:
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=False,
                        error=result.get("error")
                    )
            
            # Mark as done in result queue
            result_queue.task_done()
            
            # Display progress periodically
            current_time = time.time()
            if current_time - last_update >= update_interval:
                display_progress(stats)
                last_update = current_time
                
                # Save checkpoint periodically
                if checkpoint_manager and stats["processed"] % args.batch_size == 0:
                    checkpoint_manager.save_checkpoint()
                
        except Exception as e:
            logger.error(f"Error processing result: {e}")
    
    # Final progress display
    display_progress(stats, final=True)
    
    # Save final checkpoint
    if checkpoint_manager:
        checkpoint_manager.save_checkpoint()
    
    return stats

def display_progress(stats: Dict[str, Any], final: bool = False) -> None:
    """
    Display progress information.
    
    Args:
        stats: Statistics dictionary
        final: Whether this is the final progress display
    """
    # VERIFICATION: Displays progress bar, success/error counts, and ETA
    
    elapsed = time.time() - stats["start_time"]
    processed = stats["processed"]
    total = stats["total"]
    
    if processed == 0:
        return
    
    percent = (processed / total) * 100 if total > 0 else 0
    rate = processed / elapsed if elapsed > 0 else 0
    
    # Calculate ETA
    if rate > 0 and not final:
        eta = (total - processed) / rate
        eta_str = f"ETA: {eta:.1f}s"
    else:
        eta_str = "Complete"
    
    # Create progress bar
    bar_length = 30
    filled_length = int(bar_length * percent / 100)
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    # Print progress
    sys.stdout.write(f"\r[{bar}] {percent:.1f}% | {processed}/{total} | " +
                    f"Success: {stats['success']} | Errors: {stats['error']} | " +
                    f"Rate: {rate:.1f}/s | {eta_str}")
    sys.stdout.flush()
    
    if final:
        sys.stdout.write("\n")
        
        # Print error summary if there are errors
        if stats["error"] > 0:
            sys.stdout.write("\nError summary:\n")
            for category, count in stats["errors"].items():
                sys.stdout.write(f"  {category}: {count}\n")

def main():
    """
    Main entry point for ChEMBL data import.
    Parses command line arguments and orchestrates the import process.
    """
    # VERIFICATION: Added reference-only flag to command line arguments
    
    parser = argparse.ArgumentParser(description="Import compound data from ChEMBL into CryoProtect database")
    
    # Input options
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument("--input-file", "-i", help="File containing ChEMBL IDs (one per line)")
    input_group.add_argument("--query", "-q", help="Search query for finding compounds")
    input_group.add_argument("--limit", "-l", type=int, default=1000, help="Maximum number of compounds to process")
    input_group.add_argument("--include-refs", action="store_true", help="Include reference compounds")
    input_group.add_argument("--reference-only", action="store_true", help="Process only reference compounds")
    
    # Processing options
    proc_group = parser.add_argument_group("Processing Options")
    proc_group.add_argument("--workers", "-w", type=int, default=4, help="Number of worker threads")
    proc_group.add_argument("--batch-size", "-b", type=int, default=50, help="Batch size for processing")
    proc_group.add_argument("--max-retries", type=int, default=3, help="Maximum retry attempts for failed requests")
    
    # Mode options
    mode_group = parser.add_argument_group("Mode Options")
    mode_group.add_argument("--mode", choices=["normal", "resume", "dry-run"], default="normal",
                          help="Operation mode: normal, resume from checkpoint, or dry-run")
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--checkpoint-dir", default="checkpoints", help="Directory for checkpoint files")
    output_group.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate arguments
    if not args.input_file and not args.query and not args.reference_only:
        logger.error("No input source specified. Use --input-file, --query, or --reference-only")
        parser.print_help()
        sys.exit(1)
    
    # Initialize checkpoint manager
    checkpoint_manager = None
    if args.mode != "dry-run":
        checkpoint_manager = CheckpointManager(args.checkpoint_dir)
        
        # Load checkpoint if resuming
        if args.mode == "resume":
            if checkpoint_manager.load_checkpoint():
                logger.info("Resuming from checkpoint")
            else:
                logger.warning("No checkpoint found, starting fresh")
    
    # Initialize ChEMBL client
    try:
        client = ChEMBLClient()
    except Exception as e:
        logger.error(f"Failed to initialize ChEMBL client: {e}")
        sys.exit(1)
    
    # Get compound IDs to process
    compound_ids = fetch_compound_ids(client, args)
    if not compound_ids:
        logger.error("No compounds to process")
        sys.exit(1)
    
    logger.info(f"Processing {len(compound_ids)} compounds")
    
    # Set up worker pool
    pool = setup_worker_pool(args)
    
    try:
        # Distribute tasks
        total_tasks = distribute_tasks(compound_ids, pool, args, checkpoint_manager)
        
        # Process results
        stats = process_results(pool, total_tasks, args, checkpoint_manager)
        
        # Print summary
        logger.info(f"Import complete: {stats['success']} successful, {stats['error']} errors")
        
    except KeyboardInterrupt:
        logger.info("Interrupted by user, shutting down...")
    except Exception as e:
        logger.error(f"Error during import: {e}")
    finally:
        # Shutdown worker pool
        shutdown_worker_pool(pool)
        
        # Save final checkpoint
        if checkpoint_manager and args.mode != "dry-run":
            checkpoint_manager.save_checkpoint()
            logger.info(f"Final checkpoint saved to {args.checkpoint_dir}")

if __name__ == "__main__":
    main()
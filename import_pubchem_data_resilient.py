#!/usr/bin/env python3
"""
CryoProtect PubChem Data Importer with Resilient PubChem Client

This script imports cryoprotectant data from PubChem using the CID-Synonym-curated file
and stores it in the Supabase database using direct Supabase connection. It uses the
ResilientPubChemClient from the pubchem directory, which provides:

- Adaptive rate limiting (slower on weekdays, faster on weekends)
- Multi-level caching (in-memory and disk-based)
- Exponential backoff retry logic
- Circuit breaker to prevent repeated failures
- Fallback to pre-cached data when API is unavailable

Usage:
    python import_pubchem_data_resilient.py [--batch-size BATCH_SIZE] [--target TARGET]
                                           [--workers WORKERS] [--db-batch-size DB_BATCH_SIZE]
                                           [--resume] [--cache-dir CACHE_DIR]

Options:
    --batch-size: Number of compounds to process in each batch (default: 50)
    --target: Maximum number of compounds to import (default: 5000)
    --workers: Number of worker threads (default: 3)
    --db-batch-size: Database batch size (default: 25)
    --resume: Resume from last checkpoint
    --cache-dir: Directory to store cache files (default: cache/pubchem)
"""

import os
import sys
import time
import json
import logging
import argparse
import threading
import traceback
import concurrent.futures
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional, Set

import requests
from postgrest.exceptions import APIError

# Import the ResilientPubChemClient
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pubchem.client import ResilientPubChemClient

# Set up logging
LOG_FILE = "logs/pubchem_import_resilient.log"
SKIPPED_CID_LOG = "logs/skipped_cids_resilient.log"
os.makedirs("logs", exist_ok=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# Default parameters
DEFAULT_BATCH_SIZE = 50
DEFAULT_TARGET = 5000
DEFAULT_WORKERS = 3
DEFAULT_DB_BATCH_SIZE = 25
DEFAULT_CACHE_DIR = "cache/pubchem"

# Input file
CID_FILE = "CID-Synonym-curated"

# Checkpoint file path
CHECKPOINT_FILE = "checkpoints/pubchem_import_resilient.json"

# Core filtering criteria for cryoprotectants
CORE_CRITERIA = {
    "mw_range": (32, 500),  # Molecular weight range
    "logP_range": (-3, 3),  # LogP range
    "TPSA_range": (20, 140)  # TPSA range
}

# Queues for tracking progress
progress_queue = []
error_queue = []
rate_limit_queue = []

# Supabase client
def get_supabase_client():
    """Get a Supabase client."""
    from supabase import create_client, Client
    
    # Get Supabase URL and key from environment variables
    url = os.environ.get("SUPABASE_URL")
    key = os.environ.get("SUPABASE_KEY")
    
    if not url or not key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY environment variables must be set")
    
    # Create and return Supabase client
    return create_client(url, key)

def get_cid_list():
    """Retrieve all CIDs from the curated PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.error(f"CID file not found: {CID_FILE}")
        sys.exit(1)
    
    with open(CID_FILE, "r") as f:
        lines = f.readlines()
    
    cids = [line.strip().split("\t")[0] for line in lines if line.strip()]
    
    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

def meets_criteria(molecule):
    """Check if a molecule meets the criteria for a cryoprotectant."""
    if "Error" in molecule:
        logger.warning(f"Skipped CID {molecule['CID']}: {molecule['Error']}")
        return False
    
    if not molecule.get("Molecular Formula"):
        logger.warning(f"Skipped CID {molecule['CID']}: Missing Molecular Formula")
        return False
    
    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError):
        logger.warning(f"Skipped CID {molecule['CID']}: Invalid numerical values")
        return False
    
    if mw is None or not (CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1]):
        logger.warning(f"Skipped CID {molecule['CID']}: Molecular weight {mw} outside range {CORE_CRITERIA['mw_range']}")
        return False
    if logp is None or not (CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1]):
        logger.warning(f"Skipped CID {molecule['CID']}: LogP {logp} outside range {CORE_CRITERIA['logP_range']}")
        return False
    if tpsa is None or not (CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
        logger.warning(f"Skipped CID {molecule['CID']}: TPSA {tpsa} outside range {CORE_CRITERIA['TPSA_range']}")
        return False
    
    return True

def prepare_molecule_for_insert(molecule):
    """Prepare a molecule for insertion into the database."""
    # Extract and format the molecule data
    pubchem_cid = molecule["CID"]
    name = molecule.get("Title", "")
    smiles = molecule.get("SMILES", "")
    molecular_weight = molecule.get("Molecular Weight")
    formula = molecule.get("Molecular Formula", "")
    inchi = molecule.get("InChI", "")
    inchikey = molecule.get("InChIKey", "")
    data_source = "PubChem"
    created_at = datetime.utcnow().isoformat() + "Z"
    
    return {
        "pubchem_cid": pubchem_cid,
        "name": name,
        "smiles": smiles,
        "molecular_weight": molecular_weight,
        "formula": formula,
        "inchi": inchi,
        "inchikey": inchikey,
        "data_source": data_source,
        "created_at": created_at
    }

def batch_insert_molecules(molecules):
    """Insert a batch of molecules into the database using direct Supabase connection."""
    if not molecules:
        return 0
    
    try:
        # Get Supabase client
        supabase = get_supabase_client()
        
        # Prepare data for insertion
        data = []
        for molecule in molecules:
            data.append({
                "pubchem_cid": str(molecule["pubchem_cid"]),
                "name": molecule["name"],
                "smiles": molecule["smiles"],
                "molecular_weight": molecule["molecular_weight"],
                "formula": molecule["formula"],
                "inchi": molecule["inchi"],
                "inchikey": molecule["inchikey"],
                "data_source": molecule["data_source"],
                "created_at": molecule["created_at"]
            })
        
        # Insert data using upsert (insert or update)
        response = supabase.table("molecules").upsert(
            data, 
            on_conflict="pubchem_cid"
        ).execute()
        
        # Check for errors
        if hasattr(response, 'error') and response.error:
            logger.error(f"Database error for batch insert: {response.error}")
            return 0
        
        # Count successful inserts
        successful_inserts = len(response.data) if hasattr(response, 'data') else 0
        logger.info(f"SUCCESS: Inserted {successful_inserts} molecules in batch")
        return successful_inserts
    
    except Exception as e:
        logger.error(f"Error in batch insert: {str(e)}")
        logger.error(traceback.format_exc())
        return 0

def process_batch(batch_cids, pubchem_client, db_batch_size):
    """Process a batch of CIDs using the ResilientPubChemClient."""
    processed = 0
    imported = 0
    skipped = 0
    errors = 0
    
    # Fetch properties for all CIDs in the batch
    molecules = []
    for cid in batch_cids:
        try:
            # Get molecule properties using the resilient client
            molecule = pubchem_client.get_molecule_properties(cid)
            molecules.append(molecule)
            processed += 1
        except Exception as e:
            logger.error(f"Error fetching properties for CID {cid}: {str(e)}")
            errors += 1
            error_queue.append(1)
    
    # Filter molecules that meet criteria
    valid_molecules = []
    for molecule in molecules:
        if meets_criteria(molecule):
            valid_molecules.append(molecule)
        else:
            skipped += 1
            
            # Log skipped CIDs
            with open(SKIPPED_CID_LOG, "a") as f:
                f.write(f"{molecule['CID']}\t{datetime.now().isoformat()}\n")
    
    # Prepare molecules for insertion
    molecules_to_insert = [prepare_molecule_for_insert(m) for m in valid_molecules]
    
    # Insert molecules in batches
    if molecules_to_insert:
        for i in range(0, len(molecules_to_insert), db_batch_size):
            batch = molecules_to_insert[i:i+db_batch_size]
            batch_imported = batch_insert_molecules(batch)
            imported += batch_imported
            
            # Count errors
            batch_errors = len(batch) - batch_imported
            if batch_errors > 0:
                errors += batch_errors
                error_queue.append(batch_errors)
    
    return processed, imported, skipped, errors

def process_batch_worker(batch_cids, pubchem_client, db_batch_size):
    """Worker function to process a batch of CIDs."""
    try:
        return process_batch(batch_cids, pubchem_client, db_batch_size)
    except Exception as e:
        logger.error(f"Error in worker process: {str(e)}")
        logger.error(traceback.format_exc())
        return 0, 0, len(batch_cids), 0

def save_checkpoint(last_batch, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, status="Running"):
    """Save the current import state to a checkpoint file."""
    # Ensure checkpoints directory exists
    Path("checkpoints").mkdir(exist_ok=True)
    
    checkpoint_data = {
        "last_completed_batch": last_batch,
        "total_processed": total_processed,
        "total_imported": total_imported,
        "total_skipped": total_skipped,
        "total_errors": total_errors,
        "batch_times": batch_times,
        "start_time": start_time,
        "status": status,
        "timestamp": datetime.now().isoformat()
    }
    
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(checkpoint_data, f, indent=2)
    
    logger.info(f"Checkpoint saved: Batch {last_batch}, {total_processed} processed, {total_imported} imported")

def load_checkpoint():
    """Load the checkpoint file if it exists."""
    if os.path.exists(CHECKPOINT_FILE):
        try:
            with open(CHECKPOINT_FILE, "r") as f:
                checkpoint = json.load(f)
            
            logger.info(f"Checkpoint loaded: Batch {checkpoint['last_completed_batch']}, {checkpoint['total_processed']} processed")
            return checkpoint
        except Exception as e:
            logger.error(f"Error loading checkpoint: {str(e)}")
    
    logger.info("No checkpoint found. Starting from the beginning.")
    return None

def progress_monitor(start_time, total_cids, num_batches, stop_event):
    """Monitor and log progress."""
    last_processed = 0
    
    while not stop_event.is_set():
        # Sleep for a bit
        time.sleep(5)
        
        # Get current progress
        processed = sum(progress_queue)
        
        if processed > 0:
            # Calculate progress percentage
            progress = processed / total_cids * 100
            
            # Calculate ETA
            elapsed = time.time() - start_time
            cids_per_second = processed / elapsed if elapsed > 0 else 0
            remaining_cids = total_cids - processed
            eta_seconds = remaining_cids / cids_per_second if cids_per_second > 0 else 0
            
            # Format ETA
            if eta_seconds < 60:
                eta = f"{eta_seconds:.0f}s"
            elif eta_seconds < 3600:
                eta = f"{eta_seconds/60:.0f}m {eta_seconds%60:.0f}s"
            else:
                eta = f"{eta_seconds/3600:.0f}h {(eta_seconds%3600)/60:.0f}m"
            
            # Calculate processing rate
            new_processed = processed - last_processed
            rate = new_processed / 5 if new_processed > 0 else 0
            
            # Log progress
            logger.info(f"Progress: {progress:.1f}% ({processed}/{total_cids}), ETA: {eta}, Rate: {rate:.1f} CIDs/s")
            
            # Update last processed
            last_processed = processed
        
        # Log stats
        errors = sum(error_queue)
        rate_limits = sum(rate_limit_queue)
        logger.info(f"Stats: {processed} processed, {errors} errors, {rate_limits} rate limits")

def generate_import_report(total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, end_time):
    """Generate a report of the import process."""
    # Calculate statistics
    elapsed_time = end_time - start_time
    elapsed_minutes = elapsed_time / 60
    
    cids_per_minute = total_processed / elapsed_minutes if elapsed_minutes > 0 else 0
    success_rate = total_imported / total_processed * 100 if total_processed > 0 else 0
    
    # Calculate average batch time
    avg_batch_time = sum(batch_times) / len(batch_times) if batch_times else 0
    
    # Create report
    report = {
        "timestamp": datetime.now().isoformat(),
        "elapsed_time_seconds": elapsed_time,
        "elapsed_time_formatted": f"{int(elapsed_time // 3600)}h {int((elapsed_time % 3600) // 60)}m {int(elapsed_time % 60)}s",
        "total_processed": total_processed,
        "total_imported": total_imported,
        "total_skipped": total_skipped,
        "total_errors": total_errors,
        "cids_per_minute": cids_per_minute,
        "success_rate_percent": success_rate,
        "average_batch_time_seconds": avg_batch_time
    }
    
    # Save the report to a file
    report_file = f"reports/pubchem_import_report_resilient_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    Path("reports").mkdir(exist_ok=True)
    
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Import report generated: {report_file}")
    
    return report

def main():
    """Main function to run the PubChem data import."""
    parser = argparse.ArgumentParser(description="Import PubChem data into Supabase database using ResilientPubChemClient.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help=f"Batch size for processing (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--target", type=int, default=DEFAULT_TARGET, help=f"Target number of compounds to import (default: {DEFAULT_TARGET})")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help=f"Number of worker threads (default: {DEFAULT_WORKERS})")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    parser.add_argument("--db-batch-size", type=int, default=DEFAULT_DB_BATCH_SIZE, help=f"Database batch size (default: {DEFAULT_DB_BATCH_SIZE})")
    parser.add_argument("--cache-dir", type=str, default=DEFAULT_CACHE_DIR, help=f"Cache directory (default: {DEFAULT_CACHE_DIR})")
    args = parser.parse_args()
    
    # Initialize the ResilientPubChemClient
    pubchem_client = ResilientPubChemClient(
        cache_dir=args.cache_dir,
        weekday_requests_per_second=2.0,  # Conservative rate limit for weekdays
        weekend_requests_per_second=5.0,  # Higher rate limit for weekends
        max_retries=5,
        failure_threshold=3,
        recovery_timeout=60,
        cache_ttl=86400 * 30,  # 30 days
        memory_cache_size=1000
    )
    
    # Get list of CIDs
    cids = get_cid_list()
    
    # Limit to target number
    if args.target and args.target < len(cids):
        cids = cids[:args.target]
    
    # Initialize counters
    total_processed = 0
    total_imported = 0
    total_skipped = 0
    total_errors = 0
    batch_times = []
    
    # Calculate batches
    num_batches = (len(cids) + args.batch_size - 1) // args.batch_size
    
    # Load checkpoint if resuming
    start_batch = 0
    checkpoint = None
    if args.resume:
        checkpoint = load_checkpoint()
        if checkpoint:
            start_batch = checkpoint["last_completed_batch"] + 1
            total_processed = checkpoint["total_processed"]
            total_imported = checkpoint["total_imported"]
            total_skipped = checkpoint["total_skipped"]
            total_errors = checkpoint["total_errors"]
            batch_times = checkpoint["batch_times"]
            logger.info(f"Resuming from batch {start_batch} of {num_batches}")
    
    start_time = time.time()
    
    # Create a stop event for the progress monitor
    stop_event = threading.Event()
    
    # Start the progress monitor in a separate thread
    progress_thread = threading.Thread(target=progress_monitor, args=(start_time, len(cids), num_batches, stop_event))
    progress_thread.daemon = True
    progress_thread.start()
    
    try:
        # Process batches using thread pool
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = []
            
            # Submit all batches to the thread pool
            for i in range(start_batch, num_batches):
                # Get batch of CIDs
                start_idx = i * args.batch_size
                end_idx = min(start_idx + args.batch_size, len(cids))
                batch_cids = cids[start_idx:end_idx]
                
                logger.info(f"Processing batch {i+1}/{num_batches} with {len(batch_cids)} compounds")
                
                # Submit batch to thread pool
                future = executor.submit(process_batch_worker, batch_cids, pubchem_client, args.db_batch_size)
                futures.append((i, future))
            
            # Process results as they complete
            for i, future in futures:
                batch_start_time = time.time()
                
                # Get results from the future
                processed, imported, skipped, errors = future.result()
                
                # Update totals
                total_processed += processed
                total_imported += imported
                total_skipped += skipped
                total_errors += errors
                
                # Track batch time
                batch_time = time.time() - batch_start_time
                batch_times.append(batch_time)
                
                # Update progress queue
                progress_queue.append(processed)
                
                # Save checkpoint
                save_checkpoint(i, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time)
                
                logger.info(f"Batch {i+1}/{num_batches} completed in {batch_time:.2f}s: {processed} processed, {imported} imported, {skipped} skipped, {errors} errors")
        
        # Generate final report
        end_time = time.time()
        report = generate_import_report(total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, end_time)
        
        # Log completion
        logger.info(f"Import completed: {total_processed} processed, {total_imported} imported, {total_skipped} skipped, {total_errors} errors")
        logger.info(f"Elapsed time: {report['elapsed_time_formatted']}")
        logger.info(f"Success rate: {report['success_rate_percent']:.2f}%")
        
        # Save final checkpoint
        save_checkpoint(num_batches - 1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Completed")
        
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
        save_checkpoint(i if 'i' in locals() else -1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Interrupted")
    
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        logger.error(traceback.format_exc())
        save_checkpoint(i if 'i' in locals() else -1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Error")
    
    finally:
        # Stop the progress monitor
        stop_event.set()
        progress_thread.join(timeout=1)
        
        # Log client stats
        logger.info("PubChem client stats:")
        logger.info(f"Cache stats: {pubchem_client.get_cache_stats()}")
        logger.info(f"Rate limiter stats: {pubchem_client.get_rate_limiter_stats()}")
        logger.info(f"Circuit breaker stats: {pubchem_client.get_circuit_breaker_stats()}")

if __name__ == "__main__":
    main()
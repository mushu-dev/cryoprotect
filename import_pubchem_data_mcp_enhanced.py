#!/usr/bin/env python3
"""
CryoProtect PubChem Data Importer with MCP Integration - Enhanced Version

This script imports cryoprotectant data from PubChem using the CID-Synonym-curated file
and stores it in the Supabase database using MCP tools. It's optimized for Sunday's
higher PubChem API rate limits and includes the following enhancements:

1. Parallel processing for concurrent requests
2. Batch database operations
3. Enhanced checkpoint system
4. Improved error handling and reporting

Usage:
    python import_pubchem_data_mcp_enhanced.py [--batch-size BATCH_SIZE] [--api-delay API_DELAY] 
                                              [--target TARGET] [--workers WORKERS]
                                              [--resume] [--db-batch-size DB_BATCH_SIZE]

Parameters:
    --batch-size: Number of compounds to process in each batch (default: 100)
    --api-delay: Delay between PubChem API calls in seconds (default: 0.15)
    --target: Maximum number of compounds to import (default: 5000)
    --workers: Number of worker processes for parallel processing (default: 4)
    --resume: Resume from last checkpoint
    --db-batch-size: Number of compounds to insert in a single database operation (default: 25)
"""

import os
import sys
import time
import json
import argparse
import requests
import asyncio
import aiohttp
import concurrent.futures
from datetime import datetime, timedelta
import logging
from pathlib import Path
import queue
import threading
import traceback

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Set up logging
LOG_FILE = "logs/pubchem_import_enhanced.log"
SKIPPED_CID_LOG = "logs/skipped_cids.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Sunday-optimized parameters
DEFAULT_API_DELAY = 0.15  # Sunday-optimized delay
DEFAULT_BATCH_SIZE = 100  # Sunday-optimized batch size
DEFAULT_TARGET = 5000     # Target number of compounds
DEFAULT_WORKERS = 4       # Default number of worker processes
DEFAULT_DB_BATCH_SIZE = 25  # Default database batch size

# CID File
CID_FILE = "CID-Synonym-curated"

# Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-5, 5),
    "mw_range": (0, 1000),
    "TPSA_range": (0, 200)
}

# Checkpoint file path
CHECKPOINT_FILE = "checkpoints/pubchem_import_enhanced.json"

# Global variables for tracking progress
processed_queue = queue.Queue()
imported_queue = queue.Queue()
skipped_queue = queue.Queue()
error_queue = queue.Queue()
rate_limit_queue = queue.Queue()

def get_cid_list():
    """Retrieve all CIDs from the curated PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.warning(f"WARNING: CID file '{CID_FILE}' not found.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

async def get_molecule_properties_async(session, cid, api_delay, semaphore):
    """Fetch molecular properties and names from PubChem asynchronously."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
        "IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
    )
    
    async with semaphore:  # Limit concurrent requests
        try:
            async with session.get(url) as response:
                # Check for rate limiting
                if response.status == 429:
                    logger.warning(f"Rate limit hit for CID {cid}. Implementing exponential backoff.")
                    with open("logs/rate_limit_errors.log", "a") as f:
                        f.write(f"{datetime.now().isoformat()}: Rate limit hit for CID {cid}\n")
                    
                    # Implement exponential backoff
                    backoff_time = api_delay * 2
                    max_backoff = 10  # Maximum backoff in seconds
                    retries = 0
                    max_retries = 3
                    
                    while retries < max_retries:
                        logger.info(f"Retrying after {backoff_time:.2f}s (retry {retries+1}/{max_retries})")
                        await asyncio.sleep(backoff_time)
                        
                        async with session.get(url) as retry_response:
                            if retry_response.status == 200:
                                data = await retry_response.json()
                                break
                            
                            retries += 1
                            backoff_time = min(backoff_time * 2, max_backoff)
                    
                    if retries == max_retries:
                        logger.error(f"Failed to fetch CID {cid} after {max_retries} retries")
                        rate_limit_queue.put(1)  # Count rate limit errors
                        return {"CID": cid, "Error": "Rate limit exceeded", "RateLimit": True}
                
                # Standard delay for API rate limiting
                await asyncio.sleep(api_delay)
                
                if response.status == 200:
                    data = await response.json()
                    properties = data["PropertyTable"]["Properties"][0]
                    return {
                        "CID": cid,
                        "Molecular Formula": properties.get("MolecularFormula"),
                        "Molecular Weight": properties.get("MolecularWeight"),
                        "LogP": properties.get("XLogP"),
                        "TPSA": properties.get("TPSA"),
                        "H-Bond Donors": properties.get("HBondDonorCount"),
                        "H-Bond Acceptors": properties.get("HBondAcceptorCount"),
                        "SMILES": properties.get("IsomericSMILES"),
                        "InChI": properties.get("InChI"),
                        "InChIKey": properties.get("InChIKey"),
                        "IUPACName": properties.get("IUPACName"),
                        "Title": properties.get("Title"),
                        "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
                    }
                logger.warning(f"No molecular properties found for CID {cid}. Status code: {response.status}")
                return {"CID": cid, "Error": f"No data found (Status: {response.status})"}
        except Exception as e:
            logger.error(f"Error fetching properties for CID {cid}: {str(e)}")
            return {"CID": cid, "Error": str(e)}

def filter_molecule(molecule):
    """Initial filtering based on core cryoprotectant properties."""
    if "Error" in molecule:
        logger.warning(f"Skipped CID {molecule['CID']}: Error in molecule data")
        return False

    # Check for required fields for Supabase schema
    if not molecule.get("SMILES"):
        logger.warning(f"Skipped CID {molecule['CID']}: Missing SMILES")
        return False
    if not molecule.get("InChI"):
        logger.warning(f"Skipped CID {molecule['CID']}: Missing InChI")
        return False
    if not molecule.get("InChIKey"):
        logger.warning(f"Skipped CID {molecule['CID']}: Missing InChIKey")
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
    
    # Escape single quotes in string values
    if name:
        name = name.replace("'", "''")
    if smiles:
        smiles = smiles.replace("'", "''")
    if formula:
        formula = formula.replace("'", "''")
    if inchi:
        inchi = inchi.replace("'", "''")
    if inchikey:
        inchikey = inchikey.replace("'", "''")
    
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

def batch_insert_molecules_mcp(molecules, project_id):
    """Insert a batch of molecules into the database using MCP."""
    from use_mcp_tool import use_mcp_tool
    
    if not molecules:
        return 0
    
    try:
        # Prepare the values part of the SQL statement
        values_list = []
        for molecule in molecules:
            values = (
                f"('{molecule['pubchem_cid']}', "
                f"'{molecule['name']}', "
                f"'{molecule['smiles']}', "
                f"{molecule['molecular_weight'] if molecule['molecular_weight'] is not None else 'NULL'}, "
                f"'{molecule['formula']}', "
                f"'{molecule['inchi']}', "
                f"'{molecule['inchikey']}', "
                f"'{molecule['data_source']}', "
                f"'{molecule['created_at']}')"
            )
            values_list.append(values)
        
        values_str = ", ".join(values_list)
        
        # Build SQL statement with proper value formatting
        sql = f"""
        WITH inserted AS (
            INSERT INTO molecules (
                pubchem_cid, name, smiles, molecular_weight, formula,
                inchi, inchikey, data_source, created_at
            ) VALUES {values_str}
            ON CONFLICT (pubchem_cid) DO UPDATE SET
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                molecular_weight = EXCLUDED.molecular_weight,
                formula = EXCLUDED.formula,
                inchi = EXCLUDED.inchi,
                inchikey = EXCLUDED.inchikey,
                data_source = EXCLUDED.data_source,
                updated_at = NOW()
            RETURNING id, pubchem_cid
        )
        SELECT json_agg(
            json_build_object(
                'id', id,
                'pubchem_cid', pubchem_cid
            )
        ) as result
        FROM inserted;
        """
        
        # Execute the SQL using MCP
        result = use_mcp_tool("supabase", "execute_sql", {
            "project_id": project_id,
            "query": sql
        })
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Database error for batch insert: {result['error']}")
            return 0
        
        # Count successful inserts
        successful_inserts = len(result) if isinstance(result, list) else 0
        logger.info(f"SUCCESS: Inserted {successful_inserts} molecules in batch")
        return successful_inserts
    
    except Exception as e:
        logger.error(f"Error in batch insert: {str(e)}")
        logger.error(traceback.format_exc())
        return 0

async def process_batch_async(batch_cids, project_id, api_delay, db_batch_size):
    """Process a batch of CIDs asynchronously."""
    processed = 0
    imported = 0
    skipped = 0
    errors = 0
    
    # Create a semaphore to limit concurrent API requests
    semaphore = asyncio.Semaphore(10)  # Limit to 10 concurrent requests
    
    # Create an aiohttp session for async requests
    async with aiohttp.ClientSession() as session:
        # Fetch properties for all CIDs in the batch concurrently
        tasks = [get_molecule_properties_async(session, cid, api_delay, semaphore) for cid in batch_cids]
        molecules = await asyncio.gather(*tasks)
        
        # Process the results
        processed = len(molecules)
        processed_queue.put(processed)
        
        # Filter molecules
        valid_molecules = []
        for molecule in molecules:
            if filter_molecule(molecule):
                valid_molecules.append(prepare_molecule_for_insert(molecule))
            else:
                skipped += 1
                skipped_queue.put(1)
        
        # Insert molecules in smaller batches
        for i in range(0, len(valid_molecules), db_batch_size):
            batch = valid_molecules[i:i+db_batch_size]
            batch_imported = batch_insert_molecules_mcp(batch, project_id)
            imported += batch_imported
            imported_queue.put(batch_imported)
            
            # Count errors
            batch_errors = len(batch) - batch_imported
            if batch_errors > 0:
                errors += batch_errors
                error_queue.put(batch_errors)
    
    return processed, imported, skipped, errors

def process_batch_worker(batch_cids, project_id, api_delay, db_batch_size):
    """Worker function to process a batch of CIDs using asyncio."""
    try:
        # Run the async function in a new event loop
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        result = loop.run_until_complete(process_batch_async(batch_cids, project_id, api_delay, db_batch_size))
        loop.close()
        return result
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
        "elapsed_seconds": time.time() - start_time,
        "timestamp": datetime.now().isoformat(),
        "status": status
    }
    
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(checkpoint_data, f)
    
    logger.info(f"Checkpoint saved: Batch {last_batch}, {total_processed} processed, {total_imported} imported")

def load_checkpoint():
    """Load the last checkpoint if it exists."""
    if not os.path.exists(CHECKPOINT_FILE):
        logger.info("No checkpoint found. Starting from the beginning.")
        return None
    
    try:
        with open(CHECKPOINT_FILE, "r") as f:
            checkpoint = json.load(f)
        
        logger.info(f"Checkpoint loaded: Batch {checkpoint['last_completed_batch']}, {checkpoint['total_processed']} processed")
        return checkpoint
    
    except Exception as e:
        logger.error(f"Error loading checkpoint: {str(e)}")
        return None

def generate_import_report(total_processed, total_imported, total_skipped, total_errors, elapsed_time, rate_limit_errors=0):
    """Generate a final import report with statistics."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "statistics": {
            "total_processed": total_processed,
            "total_imported": total_imported,
            "total_skipped": total_skipped,
            "total_errors": total_errors,
            "rate_limit_errors": rate_limit_errors,
            "success_rate": round((total_imported / total_processed) * 100, 2) if total_processed > 0 else 0,
            "error_rate": round((total_errors / total_processed) * 100, 2) if total_processed > 0 else 0,
            "rate_limit_error_rate": round((rate_limit_errors / total_processed) * 100, 2) if total_processed > 0 else 0
        },
        "performance": {
            "elapsed_time_seconds": elapsed_time,
            "elapsed_time_formatted": str(timedelta(seconds=int(elapsed_time))),
            "compounds_per_second": round(total_processed / elapsed_time, 2) if elapsed_time > 0 else 0
        },
        "status": "Completed" if total_imported >= DEFAULT_TARGET else "Incomplete"
    }
    
    # Save the report to a file
    report_file = f"reports/pubchem_import_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    Path("reports").mkdir(exist_ok=True)
    
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Import report generated: {report_file}")
    return report

def progress_monitor(start_time, total_cids, num_batches, stop_event):
    """Monitor and report progress periodically."""
    last_processed = 0
    last_imported = 0
    last_skipped = 0
    last_errors = 0
    
    while not stop_event.is_set():
        try:
            # Get current counts from queues without blocking
            current_processed = last_processed
            current_imported = last_imported
            current_skipped = last_skipped
            current_errors = last_errors
            
            # Update counts from queues
            while not processed_queue.empty():
                current_processed += processed_queue.get_nowait()
            while not imported_queue.empty():
                current_imported += imported_queue.get_nowait()
            while not skipped_queue.empty():
                current_skipped += skipped_queue.get_nowait()
            while not error_queue.empty():
                current_errors += error_queue.get_nowait()
            
            # Calculate progress metrics
            elapsed_time = time.time() - start_time
            progress_percent = round((current_processed / total_cids) * 100, 1)
            
            # Calculate ETA
            if current_processed > 0:
                time_per_compound = elapsed_time / current_processed
                remaining_compounds = total_cids - current_processed
                eta_seconds = time_per_compound * remaining_compounds
                eta = str(timedelta(seconds=int(eta_seconds)))
            else:
                eta = "Unknown"
            
            # Log progress
            logger.info(f"Progress: {progress_percent}% ({current_processed}/{total_cids}), ETA: {eta}")
            logger.info(f"Stats: {current_imported} imported, {current_skipped} skipped, {current_errors} errors")
            
            # Update last counts
            last_processed = current_processed
            last_imported = current_imported
            last_skipped = current_skipped
            last_errors = current_errors
            
            # Sleep for 15 minutes (as per requirements)
            time.sleep(900)  # 15 minutes
        
        except Exception as e:
            logger.error(f"Error in progress monitor: {str(e)}")
            time.sleep(60)  # Sleep for 1 minute on error

def main():
    """Main function to run the PubChem data import."""
    parser = argparse.ArgumentParser(description="Import PubChem data into Supabase database using MCP with enhanced parallel processing.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help=f"Batch size for processing (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--api-delay", type=float, default=DEFAULT_API_DELAY, help=f"Delay between PubChem API calls in seconds (default: {DEFAULT_API_DELAY})")
    parser.add_argument("--target", type=int, default=DEFAULT_TARGET, help=f"Target number of compounds to import (default: {DEFAULT_TARGET})")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help=f"Number of worker processes (default: {DEFAULT_WORKERS})")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    parser.add_argument("--db-batch-size", type=int, default=DEFAULT_DB_BATCH_SIZE, help=f"Database batch size (default: {DEFAULT_DB_BATCH_SIZE})")
    args = parser.parse_args()
    
    # Get Supabase project ID
    from use_mcp_tool import use_mcp_tool
    
    # Get project ID from environment variable
    project_id = os.getenv("SUPABASE_PROJECT_ID")
    if not project_id:
        # Extract from URL if not found in environment
        supabase_url = os.getenv("SUPABASE_URL", "")
        if supabase_url:
            # Extract project ID from URL (format: https://[project_id].supabase.co)
            parts = supabase_url.split(".")
            if len(parts) >= 3 and "supabase" in parts[1]:
                project_id = parts[0].replace("https://", "")
    
    if not project_id:
        logger.error("Failed to get Supabase project ID. Make sure SUPABASE_PROJECT_ID or SUPABASE_URL is set in .env file.")
        sys.exit(1)
    
    logger.info(f"Using Supabase project ID: {project_id}")
    
    # Get CID list
    cids = get_cid_list()
    if not cids:
        logger.error("No CIDs found. Make sure the CID file exists.")
        sys.exit(1)
    
    # Limit to target number if specified
    if args.target and args.target < len(cids):
        logger.info(f"Limiting to {args.target} compounds (out of {len(cids)} available)")
        cids = cids[:args.target]
    
    # Initialize counters and tracking variables
    total_processed = 0
    total_imported = 0
    total_skipped = 0
    total_errors = 0
    rate_limit_errors = 0
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
        # Create a thread pool for parallel processing
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
            for i in range(start_batch, num_batches):
                # Get batch of CIDs
                start_idx = i * args.batch_size
                end_idx = min(start_idx + args.batch_size, len(cids))
                batch_cids = cids[start_idx:end_idx]
                
                logger.info(f"Processing batch {i+1}/{num_batches} with {len(batch_cids)} compounds")
                
                batch_start_time = time.time()
                
                # Process batch in parallel
                future = executor.submit(process_batch_worker, batch_cids, project_id, args.api_delay, args.db_batch_size)
                processed, imported, skipped, errors = future.result()
                
                # Update totals
                total_processed += processed
                total_imported += imported
                total_skipped += skipped
                total_errors += errors
                
                # Track batch time
                batch_time = time.time() - batch_start_time
                batch_times.append(batch_time)
                
                # Calculate progress metrics
                progress_percent = round((total_processed / len(cids)) * 100, 1)
                elapsed_time = time.time() - start_time
                avg_time_per_batch = sum(batch_times) / len(batch_times)
                remaining_batches = num_batches - (i + 1)
                eta_seconds = avg_time_per_batch * remaining_batches
                eta = str(timedelta(seconds=int(eta_seconds)))
                
                # Log progress
                logger.info(f"Batch {i+1}/{num_batches} complete. Progress: {progress_percent}%, ETA: {eta}")
                logger.info(f"Total: {total_processed} processed, {total_imported} imported, {total_skipped} skipped, {total_errors} errors")
                
                # Save checkpoint after each batch
                save_checkpoint(i, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time)
                
                # Check if we've reached the target
                if total_imported >= args.target:
                    logger.info(f"Target of {args.target} compounds reached. Stopping import.")
                    break
        
        # Get rate limit errors count
        while not rate_limit_queue.empty():
            rate_limit_errors += rate_limit_queue.get_nowait()
        
        # Generate final report
        elapsed_time = time.time() - start_time
        report = generate_import_report(total_processed, total_imported, total_skipped, total_errors, elapsed_time, rate_limit_errors)
        
        # Update checkpoint with completed status
        save_checkpoint(num_batches-1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Completed")
        
        logger.info(f"Import completed. Total: {total_processed} processed, {total_imported} imported, {total_skipped} skipped, {total_errors} errors")
        logger.info(f"Success rate: {report['statistics']['success_rate']}%, Error rate: {report['statistics']['error_rate']}%")
    
    except KeyboardInterrupt:
        logger.info("Import interrupted by user.")
        # Save checkpoint with paused status
        elapsed_time = time.time() - start_time
        save_checkpoint(i if 'i' in locals() else -1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Paused")
        generate_import_report(total_processed, total_imported, total_skipped, total_errors, elapsed_time, rate_limit_errors)
    
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        logger.error(traceback.format_exc())
        # Save checkpoint with error status
        elapsed_time = time.time() - start_time
        save_checkpoint(i if 'i' in locals() else -1, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, "Error")
        generate_import_report(total_processed, total_imported, total_skipped, total_errors, elapsed_time, rate_limit_errors)
        sys.exit(1)
    
    finally:
        # Stop the progress monitor
        stop_event.set()
        if progress_thread.is_alive():
            progress_thread.join(timeout=1)

if __name__ == "__main__":
    main()
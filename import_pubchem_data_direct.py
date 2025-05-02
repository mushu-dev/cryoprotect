#!/usr/bin/env python3
"""
CryoProtect PubChem Data Importer with Direct Supabase Connection

This script imports cryoprotectant data from PubChem using the CID-Synonym-curated file
and stores it in the Supabase database using direct Supabase connection. It's optimized for Sunday's
higher PubChem API rate limits and includes the following enhancements:

1. Direct Supabase connection for faster database operations
2. Parallel processing for concurrent requests
3. Batch database operations
4. Enhanced checkpoint system

Usage:
    python import_pubchem_data_direct.py [--batch-size BATCH_SIZE] [--api-delay API_DELAY] 
                                        [--target TARGET] [--workers WORKERS]
                                        [--resume] [--db-batch-size DB_BATCH_SIZE]

Parameters:
    --batch-size: Number of compounds to process in each batch (default: 100)
    --api-delay: Delay between PubChem API calls in seconds (default: 0.15)
    --target: Maximum number of compounds to import (default: 5000)
    --workers: Number of worker threads for parallel processing (default: 10)
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
import threading
import concurrent.futures
from datetime import datetime, timedelta
import logging
from pathlib import Path
import queue
import traceback
from typing import List, Dict, Any, Optional
from dotenv import load_dotenv

# Import new utility modules for improved resilience and performance
from db_connection_utils import get_db_connection, safe_transaction
from transaction_utils import with_transaction_retry, execute_in_transaction, is_transaction_active
from batch_utils import bulk_insert_properties, resumable_batch_import, batch_delete_properties

# Import legacy Supabase client for backward compatibility
from supabase import create_client, Client

# Load environment variables
load_dotenv()

# Get Supabase credentials
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Set up logging
LOG_FILE = "logs/pubchem_import_direct.log"
SKIPPED_CID_LOG = "logs/skipped_cids_direct.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Parameters adjusted to avoid PubChem API rate limiting
DEFAULT_API_DELAY = 0.5   # Increased delay to avoid 503 errors
DEFAULT_BATCH_SIZE = 50   # Reduced batch size
DEFAULT_TARGET = 5000     # Target number of compounds
DEFAULT_WORKERS = 3       # Reduced number of worker threads
DEFAULT_DB_BATCH_SIZE = 25  # Default database batch size
DEFAULT_MAX_RETRIES = 5   # Maximum number of retries for API calls
DEFAULT_RETRY_DELAY = 2   # Initial retry delay in seconds

# CID File
CID_FILE = "CID-Synonym-curated"

# Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-5, 5),
    "mw_range": (0, 1000),
    "TPSA_range": (0, 200)
}

# Checkpoint file path
CHECKPOINT_FILE = "checkpoints/pubchem_import_direct.json"

# Global variables for tracking progress
processed_queue = queue.Queue()
imported_queue = queue.Queue()
skipped_queue = queue.Queue()
error_queue = queue.Queue()
rate_limit_queue = queue.Queue()

# Global Supabase client
supabase_client = None

def get_supabase_client() -> Client:
    """
    Get a Supabase client with improved connection resilience.
    Uses the new db_connection_utils module for better error handling and retry logic.
    """
    global supabase_client
    
    if supabase_client is not None:
        return supabase_client
    
    if not SUPABASE_URL or not SUPABASE_KEY:
        logger.error("Supabase URL or key not found in environment variables.")
        sys.exit(1)
    
    # Try to get a connection with retry logic
    for attempt in range(3):  # Try up to 3 times
        try:
            # First try to use the new db_connection_utils
            try:
                # Check if we can get a direct database connection
                with get_db_connection() as conn:
                    logger.info("Successfully connected to database using db_connection_utils")
                    # If successful, still create a Supabase client for backward compatibility
            except Exception as db_err:
                logger.warning(f"Could not connect using db_connection_utils: {str(db_err)}")
            
            # Create the Supabase client
            supabase_client = create_client(SUPABASE_URL, SUPABASE_KEY)
            logger.info(f"Connected to Supabase at {SUPABASE_URL}")
            return supabase_client
            
        except Exception as e:
            logger.warning(f"Connection attempt {attempt+1} failed: {str(e)}")
            if attempt < 2:  # If not the last attempt
                wait_time = (attempt + 1) * 2  # Exponential backoff: 2s, 4s
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                logger.error(f"Failed to connect to Supabase after 3 attempts: {str(e)}")
                sys.exit(1)

def get_cid_list():
    """Retrieve all CIDs from the curated PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.warning(f"WARNING: CID file '{CID_FILE}' not found.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

async def get_molecule_properties_async(session, cid, api_delay, semaphore, max_retries=DEFAULT_MAX_RETRIES, retry_delay=DEFAULT_RETRY_DELAY):
    """Fetch molecular properties and names from PubChem asynchronously with improved error handling."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
        "IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
    )
    
    async with semaphore:  # Limit concurrent requests
        # Standard delay before making the request to avoid overwhelming the API
        await asyncio.sleep(api_delay)
        
        retries = 0
        backoff_time = retry_delay
        max_backoff = 30  # Maximum backoff in seconds
        
        while retries <= max_retries:
            try:
                async with session.get(url, timeout=30) as response:
                    # Handle different response status codes
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
                    
                    # Handle rate limiting (429) or server errors (5xx)
                    elif response.status == 429 or response.status >= 500:
                        if retries < max_retries:
                            error_type = "Rate limit" if response.status == 429 else "Server error"
                            logger.warning(f"{error_type} for CID {cid}. Status: {response.status}. Retry {retries+1}/{max_retries} after {backoff_time:.2f}s")
                            
                            # Log rate limit errors
                            if response.status == 429:
                                with open("logs/rate_limit_errors.log", "a") as f:
                                    f.write(f"{datetime.now().isoformat()}: Rate limit hit for CID {cid}\n")
                                rate_limit_queue.put(1)  # Count rate limit errors
                            
                            # Wait before retrying with exponential backoff
                            await asyncio.sleep(backoff_time)
                            retries += 1
                            backoff_time = min(backoff_time * 2, max_backoff)
                            continue
                        else:
                            logger.error(f"Failed to fetch CID {cid} after {max_retries} retries. Status: {response.status}")
                            return {"CID": cid, "Error": f"Failed after {max_retries} retries (Status: {response.status})"}
                    
                    # Handle other errors
                    else:
                        logger.warning(f"No molecular properties found for CID {cid}. Status code: {response.status}")
                        return {"CID": cid, "Error": f"No data found (Status: {response.status})"}
            
            except asyncio.TimeoutError:
                if retries < max_retries:
                    logger.warning(f"Timeout for CID {cid}. Retry {retries+1}/{max_retries} after {backoff_time:.2f}s")
                    await asyncio.sleep(backoff_time)
                    retries += 1
                    backoff_time = min(backoff_time * 2, max_backoff)
                    continue
                else:
                    logger.error(f"Timeout fetching CID {cid} after {max_retries} retries")
                    return {"CID": cid, "Error": f"Timeout after {max_retries} retries"}
            
            except Exception as e:
                if retries < max_retries:
                    logger.warning(f"Error fetching CID {cid}: {str(e)}. Retry {retries+1}/{max_retries} after {backoff_time:.2f}s")
                    await asyncio.sleep(backoff_time)
                    retries += 1
                    backoff_time = min(backoff_time * 2, max_backoff)
                    continue
                else:
                    logger.error(f"Error fetching CID {cid} after {max_retries} retries: {str(e)}")
                    return {"CID": cid, "Error": str(e)}
        
        # This should not be reached, but just in case
        return {"CID": cid, "Error": "Unknown error in retry loop"}

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
    """
    Insert a batch of molecules into the database using enhanced batch utilities.
    This function now uses the new batch_utils module for improved resilience and performance.
    """
    if not molecules:
        return 0
    
    try:
        # First try to use the new batch_utils module with direct database connection
        try:
            # Convert molecules to property records format for batch_utils
            property_records = []
            for molecule in molecules:
                # Create a base record with molecule ID and properties
                molecule_id = str(molecule["pubchem_cid"])
                
                # Add each property as a separate record
                for prop_name, prop_value in molecule.items():
                    if prop_name == "pubchem_cid":
                        continue  # Skip the ID field
                        
                    # Determine property type
                    if prop_name in ["molecular_weight", "logP", "TPSA"]:
                        property_records.append({
                            "molecule_id": molecule_id,
                            "property_name": prop_name,
                            "numeric_value": prop_value,
                            "text_value": None,
                            "created_at": molecule.get("created_at")
                        })
                    else:
                        property_records.append({
                            "molecule_id": molecule_id,
                            "property_name": prop_name,
                            "numeric_value": None,
                            "text_value": str(prop_value) if prop_value is not None else None,
                            "created_at": molecule.get("created_at")
                        })
            
            # Use resumable batch import for better resilience
            checkpoint_file = f"checkpoints/pubchem_batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            
            # Define a process function for the batch
            def process_batch(batch):
                with safe_transaction():
                    return bulk_insert_properties(batch, batch_size=25)
            
            # Execute the batch import with checkpointing
            successful_inserts = resumable_batch_import(
                property_records,
                process_func=process_batch,
                checkpoint_file=checkpoint_file,
                batch_size=25
            )
            
            logger.info(f"SUCCESS: Inserted {successful_inserts} molecule properties using batch_utils")
            return successful_inserts
            
        except Exception as batch_error:
            logger.warning(f"Error using batch_utils: {str(batch_error)}. Falling back to Supabase client.")
            
            # Fall back to the original Supabase client method
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
            logger.info(f"SUCCESS: Inserted {successful_inserts} molecules in batch using Supabase fallback")
            return successful_inserts
    
    except Exception as e:
        logger.error(f"Error in batch insert: {str(e)}")
        logger.error(traceback.format_exc())
        return 0

async def process_batch_async(batch_cids, api_delay, db_batch_size, max_retries=DEFAULT_MAX_RETRIES, retry_delay=DEFAULT_RETRY_DELAY):
    """Process a batch of CIDs asynchronously."""
    processed = 0
    imported = 0
    skipped = 0
    errors = 0
    
    # Create a semaphore to limit concurrent API requests
    semaphore = asyncio.Semaphore(5)  # Reduced to 5 concurrent requests to avoid overwhelming the API
    
    # Create an aiohttp session for async requests
    async with aiohttp.ClientSession() as session:
        # Fetch properties for all CIDs in the batch concurrently
        tasks = [get_molecule_properties_async(session, cid, api_delay, semaphore, max_retries, retry_delay) for cid in batch_cids]
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
            batch_imported = batch_insert_molecules(batch)
            imported += batch_imported
            imported_queue.put(batch_imported)
            
            # Count errors
            batch_errors = len(batch) - batch_imported
            if batch_errors > 0:
                errors += batch_errors
                error_queue.put(batch_errors)
    
    return processed, imported, skipped, errors

def process_batch_worker(batch_cids, api_delay, db_batch_size, max_retries=DEFAULT_MAX_RETRIES, retry_delay=DEFAULT_RETRY_DELAY):
    """Worker function to process a batch of CIDs using asyncio."""
    try:
        # Run the async function in a new event loop
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        result = loop.run_until_complete(process_batch_async(batch_cids, api_delay, db_batch_size, max_retries, retry_delay))
        loop.close()
        return result
    except Exception as e:
        logger.error(f"Error in worker process: {str(e)}")
        logger.error(traceback.format_exc())
        return 0, 0, len(batch_cids), 0

def save_checkpoint(last_batch, total_processed, total_imported, total_skipped, total_errors, batch_times, start_time, status="Running"):
    """
    Save the current import state to a checkpoint file.
    Enhanced with better error handling and using the batch_utils module's checkpoint format.
    """
    # Ensure checkpoints directory exists
    Path("checkpoints").mkdir(exist_ok=True)
    
    try:
        # Create checkpoint data in the format compatible with batch_utils
        checkpoint_data = {
            "position": last_batch * DEFAULT_BATCH_SIZE,  # Convert batch number to item position
            "processed": total_processed,
            "last_batch": last_batch,
            "total_imported": total_imported,
            "total_skipped": total_skipped,
            "total_errors": total_errors,
            "batch_times": batch_times,
            "elapsed_seconds": time.time() - start_time,
            "last_updated": datetime.now().isoformat(),
            "status": status
        }
        
        # Save checkpoint using safe transaction pattern
        with safe_transaction():
            with open(CHECKPOINT_FILE, "w") as f:
                json.dump(checkpoint_data, f, indent=2)
        
        logger.info(f"Checkpoint saved: Batch {last_batch}, {total_processed} processed, {total_imported} imported")
    except Exception as e:
        logger.error(f"Error saving checkpoint: {str(e)}")
        # Try to save to an alternative location as a backup
        try:
            backup_file = f"checkpoints/backup_checkpoint_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(backup_file, "w") as f:
                json.dump(checkpoint_data, f, indent=2)
            logger.info(f"Backup checkpoint saved to {backup_file}")
        except Exception as backup_error:
            logger.error(f"Failed to save backup checkpoint: {str(backup_error)}")

def load_checkpoint():
    """
    Load the last checkpoint if it exists.
    Enhanced to handle both old and new checkpoint formats for backward compatibility.
    """
    if not os.path.exists(CHECKPOINT_FILE):
        logger.info("No checkpoint found. Starting from the beginning.")
        return None
    
    try:
        # Try to load the checkpoint using safe_transaction for better error handling
        with safe_transaction():
            with open(CHECKPOINT_FILE, "r") as f:
                checkpoint = json.load(f)
        
        # Handle both old and new checkpoint formats
        if "position" in checkpoint:
            # New format (compatible with batch_utils)
            last_batch = checkpoint.get("last_batch", 0)
            processed = checkpoint.get("processed", 0)
            logger.info(f"Checkpoint loaded (new format): Batch {last_batch}, {processed} processed")
            
            # Convert to old format for backward compatibility
            if "last_completed_batch" not in checkpoint:
                checkpoint["last_completed_batch"] = last_batch
            if "total_processed" not in checkpoint and "processed" in checkpoint:
                checkpoint["total_processed"] = checkpoint["processed"]
        else:
            # Old format
            logger.info(f"Checkpoint loaded (old format): Batch {checkpoint.get('last_completed_batch', 0)}, {checkpoint.get('total_processed', 0)} processed")
        
        return checkpoint
    
    except Exception as e:
        logger.error(f"Error loading checkpoint: {str(e)}")
        
        # Try to find backup checkpoints
        try:
            backup_files = [f for f in os.listdir("checkpoints") if f.startswith("backup_checkpoint_")]
            if backup_files:
                # Sort by timestamp (newest first)
                backup_files.sort(reverse=True)
                latest_backup = os.path.join("checkpoints", backup_files[0])
                logger.info(f"Attempting to load backup checkpoint: {latest_backup}")
                
                with open(latest_backup, "r") as f:
                    backup_checkpoint = json.load(f)
                
                logger.info(f"Successfully loaded backup checkpoint: {latest_backup}")
                return backup_checkpoint
        except Exception as backup_error:
            logger.error(f"Failed to load backup checkpoint: {str(backup_error)}")
        
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
    report_file = f"reports/pubchem_import_report_direct_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
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
    parser = argparse.ArgumentParser(description="Import PubChem data into Supabase database using direct connection.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help=f"Batch size for processing (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--api-delay", type=float, default=DEFAULT_API_DELAY, help=f"Delay between PubChem API calls in seconds (default: {DEFAULT_API_DELAY})")
    parser.add_argument("--target", type=int, default=DEFAULT_TARGET, help=f"Target number of compounds to import (default: {DEFAULT_TARGET})")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help=f"Number of worker threads (default: {DEFAULT_WORKERS})")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    parser.add_argument("--db-batch-size", type=int, default=DEFAULT_DB_BATCH_SIZE, help=f"Database batch size (default: {DEFAULT_DB_BATCH_SIZE})")
    parser.add_argument("--max-retries", type=int, default=DEFAULT_MAX_RETRIES, help=f"Maximum number of retries for API calls (default: {DEFAULT_MAX_RETRIES})")
    parser.add_argument("--retry-delay", type=float, default=DEFAULT_RETRY_DELAY, help=f"Initial retry delay in seconds (default: {DEFAULT_RETRY_DELAY})")
    args = parser.parse_args()
    
    # Initialize Supabase client
    supabase = get_supabase_client()
    logger.info(f"Using Supabase URL: {SUPABASE_URL}")
    
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
                future = executor.submit(process_batch_worker, batch_cids, args.api_delay, args.db_batch_size, args.max_retries, args.retry_delay)
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
#!/usr/bin/env python3
"""
Full ChEMBL Import Script

This script imports cryoprotectant-related compounds from ChEMBL
into the CryoProtect database using the worker pool architecture.
Uses the connection factory for improved connection resilience.
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime, timedelta
from queue import Queue, Empty
import requests
from typing import List, Dict, Any, Optional, Tuple, Set
from dotenv import load_dotenv

# Import RDKit for molecular property filtering and SMARTS pattern matching
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, MolSurf
    from rdkit.Chem.Draw import MolToImage
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Molecular property filtering and SMARTS pattern matching will be disabled.")

# Import new utility modules for improved resilience and performance
from db_connection_utils import get_db_connection, safe_transaction
from transaction_utils import with_transaction_retry, execute_in_transaction, is_transaction_active
from batch_utils import bulk_insert_properties, resumable_batch_import, batch_delete_properties
from chembl_search_utils import (
    find_potential_cryoprotectants,
    find_similar_compounds,
    identify_compounds_by_chemical_class,
    store_cryoprotectant_candidates,
    store_similar_compounds,
    store_chemical_class_compounds
)

# Import legacy utilities for backward compatibility
from database.connection import get_db_connection_info
from sql_executor import (
    get_db,
    execute_query,
    bulk_insert,
    execute_batch,
    with_retry,
    process_in_batches,
    get_connection_metrics
)

# Import enhanced property manager
from property_utils import PropertyManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_import.log')
    ]
)

# Load environment variables
load_dotenv()

logger = logging.getLogger(__name__)

class ChEMBLClient:
    """
    Client for interacting with the ChEMBL API to fetch compound IDs.
    Uses retry mechanism for improved resilience.
    """
    
    def __init__(self, base_url: str = "https://www.ebi.ac.uk/chembl/api/data"):
        """
        Initialize the ChEMBL client.
        
        Args:
            base_url: Base URL for the ChEMBL API
        """
        self.base_url = base_url
        logger.info(f"ChEMBL client initialized with base URL: {base_url}")
    
    @with_retry(max_retries=3, retry_backoff=2.0)
    def search_compounds(self, term: str, limit: int = 1000) -> List[str]:
        """
        Search for compounds by term and return their ChEMBL IDs.
        
        Args:
            term: Search term
            limit: Maximum number of results to return
            
        Returns:
            List of ChEMBL IDs
        """
        logger.info(f"Searching for compounds with term: {term}")
        
        # Construct the search URL
        url = f"{self.base_url}/molecule"
        params = {
            "pref_name__icontains": term,
            "limit": limit
        }
        
        try:
            # Make the request
            headers = {"Accept": "application/json"}
            response = requests.get(url, params=params, headers=headers, timeout=30)
            response.raise_for_status()
            
            # Parse the response
            data = response.json()
            molecules = data.get("molecules", [])
            
            # Extract ChEMBL IDs
            chembl_ids = [mol.get("molecule_chembl_id") for mol in molecules if mol.get("molecule_chembl_id")]
            
            logger.info(f"Found {len(chembl_ids)} compounds for term '{term}'")
            return chembl_ids
            
        except Exception as e:
            logger.error(f"Error searching for compounds with term '{term}': {e}")
            raise  # Re-raise for retry decorator
            
    @with_retry(max_retries=3, retry_backoff=2.0)
    def get_similar_compounds(self, chembl_id: str, similarity: int = 70, limit: int = 20) -> List[str]:
        """
        Get similar compounds to a given ChEMBL ID based on structural similarity.
        
        Args:
            chembl_id: ChEMBL ID of the reference compound
            similarity: Minimum similarity percentage (0-100)
            limit: Maximum number of similar compounds to return
            
        Returns:
            List of ChEMBL IDs of similar compounds
        """
        logger.info(f"Searching for compounds similar to {chembl_id} with similarity >= {similarity}%")
        
        # Construct the search URL
        url = f"{self.base_url}/similarity/{chembl_id}/{similarity}"
        params = {
            "limit": limit
        }
        
        try:
            # Make the request
            headers = {"Accept": "application/json"}
            response = requests.get(url, params=params, headers=headers, timeout=30)
            response.raise_for_status()
            
            # Parse the response
            data = response.json()
            molecules = data.get("molecules", [])
            
            # Extract ChEMBL IDs
            similar_ids = []
            for mol in molecules:
                mol_id = mol.get("molecule_chembl_id")
                if mol_id and mol_id != chembl_id:  # Exclude the reference compound itself
                    similar_ids.append(mol_id)
            
            logger.info(f"Found {len(similar_ids)} compounds similar to {chembl_id}")
            return similar_ids
            
        except Exception as e:
            logger.error(f"Error searching for compounds similar to {chembl_id}: {e}")
            raise  # Re-raise for retry decorator
    
    @with_retry(max_retries=3, retry_backoff=2.0)
    def get_all_compound_ids(self, limit: Optional[int] = None) -> List[str]:
        """
        Fetch all compound IDs from ChEMBL.
        
        This method retrieves all compound IDs from ChEMBL using pagination.
        
        Args:
            limit: Optional maximum number of IDs to return
            
        Returns:
            List of ChEMBL IDs
        """
        logger.info("Fetching all compound IDs from ChEMBL")
        
        all_ids = []
        page_size = 1000
        offset = 0
        
        while True:
            # Construct the URL for this page
            url = f"{self.base_url}/molecule"
            params = {
                "limit": page_size,
                "offset": offset
            }
            
            try:
                # Make the request
                headers = {"Accept": "application/json"}
                response = requests.get(url, params=params, headers=headers, timeout=30)
                response.raise_for_status()
                
                # Parse the response
                data = response.json()
                molecules = data.get("molecules", [])
                
                # If no molecules returned, we've reached the end
                if not molecules:
                    break
                
                # Extract ChEMBL IDs
                page_ids = [mol.get("molecule_chembl_id") for mol in molecules if mol.get("molecule_chembl_id")]
                all_ids.extend(page_ids)
                
                logger.info(f"Fetched {len(page_ids)} compound IDs (total: {len(all_ids)})")
                
                # If we've reached the limit, stop
                if limit and len(all_ids) >= limit:
                    all_ids = all_ids[:limit]
                    break
                
                # Move to the next page
                offset += page_size
                
                # Add a small delay to avoid overwhelming the API
                time.sleep(0.5)
                
            except Exception as e:
                logger.error(f"Error fetching compound IDs at offset {offset}: {e}")
                # For pagination, we'll continue with the next page despite errors
                # rather than retrying the same page
                offset += page_size
        
        logger.info(f"Fetched a total of {len(all_ids)} compound IDs")
        return all_ids
        
    @with_retry(max_retries=3, retry_backoff=2.0)
    def get_compound(self, chembl_id: str) -> Dict[str, Any]:
        """
        Get compound data by ChEMBL ID.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Dictionary with compound data in the format expected by property filtering functions:
            {
                'molecule_structures': {
                    'canonical_smiles': '...'
                },
                'molecule_chembl_id': 'CHEMBLXXXX',
                ...
            }
        """
        logger.debug(f"Getting compound data for {chembl_id}")
        
        try:
            # Construct the URL
            url = f"{self.base_url}/molecule/{chembl_id}"
            
            # Make the request
            headers = {"Accept": "application/json"}
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            # Parse the response
            data = response.json()
            
            # Ensure the data has the expected structure
            if "molecule_structures" not in data:
                data["molecule_structures"] = {}
            
            # Make sure canonical_smiles is present
            if "canonical_smiles" not in data.get("molecule_structures", {}):
                # Try to get SMILES from other fields if available
                smiles = None
                
                # Check if we have a SMILES string in another location
                if "molecule_structures" in data and "smiles" in data["molecule_structures"]:
                    smiles = data["molecule_structures"]["smiles"]
                elif "smiles" in data:
                    smiles = data["smiles"]
                
                # Add it to the expected location
                if smiles:
                    data["molecule_structures"]["canonical_smiles"] = smiles
            
            return data
            
        except Exception as e:
            logger.warning(f"Error getting compound data for {chembl_id}: {e}")
            # Return a minimal structure that won't cause errors in property filtering
            return {
                "molecule_chembl_id": chembl_id,
                "error": str(e),
                "molecule_structures": {}
            }

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Import ChEMBL data into CryoProtect database")
    
    # Basic execution parameters
    parser.add_argument("--limit", type=int, default=1000,
                        help="Maximum number of compounds to import")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Number of compounds per batch")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of worker threads")
    parser.add_argument("--checkpoint-dir", default="./checkpoints",
                        help="Directory for checkpoint files")
    parser.add_argument("--mode", choices=["full", "resume", "dry-run"], default="resume",
                        help="Execution mode")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose logging")
    
    # Compound discovery strategy options
    strategy_group = parser.add_argument_group('Compound Discovery Strategy')
    strategy_group.add_argument("--search-terms", nargs="+",
                        default=["cryoprotectant", "antifreeze", "cryopreservation", "cell preservation"],
                        help="Search terms for finding compounds")
    strategy_group.add_argument("--search-terms-all", action="store_true",
                        help="Use all predefined search terms from categories")
    strategy_group.add_argument("--use-reference-seeds", action="store_true",
                        help="Use reference compounds as seeds for similarity search")
    strategy_group.add_argument("--similarity", type=int, default=70,
                        help="Minimum similarity percentage for structural similarity search (0-100)")
    strategy_group.add_argument("--expand-similar", action="store_true",
                        help="Expand search results with similar compounds")
    strategy_group.add_argument("--include-refs", action="store_true",
                        help="Include reference compounds")
    strategy_group.add_argument("--fetch-all", action="store_true",
                        help="Fetch all ChEMBL compounds instead of searching by terms")
    
    # Molecular property filtering options
    property_group = parser.add_argument_group('Molecular Property Filtering')
    property_group.add_argument("--filter-properties", action="store_true",
                        help="Filter compounds based on molecular properties")
    property_group.add_argument("--min-mw", type=float, default=30.0,
                        help="Minimum molecular weight")
    property_group.add_argument("--max-mw", type=float, default=500.0,
                        help="Maximum molecular weight")
    property_group.add_argument("--min-hba", type=int, default=2,
                        help="Minimum number of hydrogen bond acceptors")
    property_group.add_argument("--min-hbd", type=int, default=1,
                        help="Minimum number of hydrogen bond donors")
    property_group.add_argument("--max-alogp", type=float, default=3.0,
                        help="Maximum ALogP value")
    
    # SMARTS pattern matching options
    smarts_group = parser.add_argument_group('SMARTS Pattern Matching')
    smarts_group.add_argument("--use-smarts", action="store_true",
                        help="Use SMARTS pattern matching for specific chemical classes")
    smarts_group.add_argument("--smarts-classes", nargs="+",
                        choices=["polyol", "sugar", "amide", "amine", "alcohol", "ether", "all"],
                        default=["all"],
                        help="Chemical classes to match with SMARTS patterns")
    
    return parser.parse_args()

def main():
    """Main execution function"""
    args = parse_arguments()
    
    # Create logs directory if it doesn't exist
    os.makedirs('logs', exist_ok=True)
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting full ChEMBL import in {args.mode} mode")
    logger.info(f"Settings: limit={args.limit}, batch_size={args.batch_size}, workers={args.workers}")
    
    if args.fetch_all:
        logger.info("Fetching all ChEMBL compounds")
    else:
        logger.info(f"Search terms: {args.search_terms}")
    
    # Import necessary modules
    try:
        from chembl.worker import ChEMBLWorker
        from chembl.checkpoint import CheckpointManager
        from chembl.reference_compounds import get_reference_compound_ids
        from chembl.search_terms import get_all_search_terms
    except ImportError as e:
        logger.error(f"Failed to import required modules: {e}")
        return 1
    
    # Verify database connection using enhanced connection utilities
    logger.info("Verifying database connection using enhanced connection utilities...")
    try:
        with get_db_connection() as connection:
            # Test connection with a simple query
            result = connection.execute_query("SELECT 1 as connection_test")
            if result and result[0]['connection_test'] == 1:
                logger.info("Successfully connected to database")
                
                # Get connection info for backward compatibility
                conn_info = get_db_connection_info()
                logger.info(f"Connection info: {conn_info}")
                
                # Get connection metrics for backward compatibility
                metrics = get_connection_metrics()
                logger.info(f"Connection metrics: {metrics}")
            else:
                logger.warning("Database connection test returned unexpected result. Will attempt to reconnect when needed.")
    except Exception as e:
        logger.warning(f"Database connection verification failed: {e}")
        logger.info("Will attempt to connect when needed.")
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(args.checkpoint_dir, prefix="full_chembl_import")
    
    # Load checkpoint if in resume mode
    if args.mode == "resume":
        checkpoint_loaded = checkpoint_manager.load_checkpoint()
        if checkpoint_loaded:
            logger.info("Resuming from checkpoint")
        else:
            logger.info("No checkpoint found, starting fresh")
    
    # Initialize client for fetching compound IDs
    client = ChEMBLClient()
    
    # Get compound IDs to process using the new strategy
    compound_ids = []
    
    # Add reference compounds if requested
    if args.include_refs:
        reference_compounds = get_reference_compound_ids()
        logger.info(f"Including {len(reference_compounds)} reference compounds")
        compound_ids.extend(reference_compounds)
    
    # Strategy 1: Use reference compounds as seeds for similarity search
    if args.use_reference_seeds:
        logger.info(f"Using reference compounds as seeds for similarity search (similarity >= {args.similarity}%)")
        similar_to_refs = get_compounds_from_reference_seeds(client, similarity=args.similarity, limit=args.limit)
        compound_ids.extend(similar_to_refs)
        logger.info(f"Found {len(similar_to_refs)} compounds similar to reference compounds")
    
    # Strategy 2: Traditional search terms approach (as fallback or if specifically requested)
    elif not args.fetch_all:
        # Search for compounds by terms
        if args.search_terms_all:
            # Use all predefined search terms from categories
            search_terms = get_all_search_terms()
            logger.info(f"Using all predefined search terms: {len(search_terms)} terms")
            args.search_terms = search_terms
        
        for term in args.search_terms:
            logger.info(f"Searching for compounds with term: {term}")
            results = client.search_compounds(term)
            logger.info(f"Found {len(results)} compounds for term '{term}'")
            compound_ids.extend(results)
            
        # Expand search results with similar compounds if requested
        if args.expand_similar:
            logger.info(f"Expanding search results with similar compounds (similarity >= {args.similarity}%)...")
            compound_ids = expand_search_results(client, compound_ids, similarity=args.similarity, limit=args.limit)
            logger.info(f"Total compounds after expansion: {len(compound_ids)}")
    
    # Strategy 3: Fetch all compounds (if specifically requested)
    elif args.fetch_all:
        logger.info("Fetching all ChEMBL compound IDs")
        all_ids = client.get_all_compound_ids(limit=args.limit)
        compound_ids.extend(all_ids)
    
    # Deduplicate
    compound_ids = list(set(compound_ids))
    logger.info(f"Total unique compounds after deduplication: {len(compound_ids)}")
    
    # Strategy 4: Filter by molecular properties common to cryoprotectants
    if args.filter_properties and RDKIT_AVAILABLE:
        logger.info("Filtering compounds based on molecular properties...")
        filtered_ids = filter_by_molecular_properties(
            client,
            compound_ids,
            min_mw=args.min_mw,
            max_mw=args.max_mw,
            min_hba=args.min_hba,
            min_hbd=args.min_hbd,
            max_alogp=args.max_alogp
        )
        compound_ids = filtered_ids
        logger.info(f"Total compounds after property filtering: {len(compound_ids)}")
    
    # Strategy 5: Use SMARTS pattern matching for specific chemical classes
    if args.use_smarts and RDKIT_AVAILABLE:
        logger.info("Applying SMARTS pattern matching for specific chemical classes...")
        smarts_matches = match_smarts_patterns(client, compound_ids)
        
        # Filter by selected chemical classes
        if "all" not in args.smarts_classes:
            filtered_by_class = set()
            for class_name in args.smarts_classes:
                if class_name in smarts_matches:
                    filtered_by_class.update(smarts_matches[class_name])
            
            compound_ids = list(filtered_by_class)
            logger.info(f"Total compounds after SMARTS filtering: {len(compound_ids)}")
    
    # Apply limit
    if args.limit and len(compound_ids) > args.limit:
        logger.info(f"Limiting to {args.limit} compounds")
        compound_ids = compound_ids[:args.limit]
    
    # Skip already processed compounds if resuming
    if args.mode == "resume" and checkpoint_manager.state.get("processed_compounds"):
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
    
    # Store configuration in checkpoint with enhanced strategy parameters
    checkpoint_manager.state["config"] = {
        "mode": args.mode,
        "workers": args.workers,
        "batch_size": args.batch_size,
        "limit": args.limit,
        "search_terms": args.search_terms if not args.fetch_all else ["all"],
        "fetch_all": args.fetch_all,
        "include_refs": args.include_refs,
        "expand_similar": args.expand_similar,
        "search_terms_all": args.search_terms_all,
        "use_reference_seeds": args.use_reference_seeds,
        "similarity": args.similarity,
        "filter_properties": args.filter_properties,
        "property_filters": {
            "min_mw": args.min_mw,
            "max_mw": args.max_mw,
            "min_hba": args.min_hba,
            "min_hbd": args.min_hbd,
            "max_alogp": args.max_alogp
        },
        "use_smarts": args.use_smarts,
        "smarts_classes": args.smarts_classes,
        "rdkit_available": RDKIT_AVAILABLE,
        "total_compounds": len(compound_ids)
    }
    
    # Save initial checkpoint
    checkpoint_manager.save_checkpoint()
    
    # If dry-run mode, just log what would be done
    if args.mode == "dry-run":
        logger.info(f"DRY RUN: Would import {len(compound_ids)} compounds")
        for i, compound_id in enumerate(compound_ids[:10]):  # Show first 10 as example
            logger.info(f"DRY RUN: Would process compound {i+1}/{len(compound_ids)}: {compound_id}")
        if len(compound_ids) > 10:
            logger.info(f"DRY RUN: ... and {len(compound_ids) - 10} more compounds")
        logger.info("DRY RUN: Import simulation complete")
        return 0
    
    # Set up worker pool
    task_queue = Queue()
    result_queue = Queue()
    
    # Create workers
    workers = []
    for i in range(args.workers):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        workers.append(worker)
    
    # Start workers
    for worker in workers:
        worker.start()
        
    logger.info(f"Started {len(workers)} workers")
    
    # Add dry run flag for dry-run mode
    dry_run = (args.mode == "dry-run")
    
    # Add tasks to queue
    for compound_id in compound_ids:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": dry_run,
            "reference": False,  # Not a reference compound unless explicitly included
            "checkpoint_manager": checkpoint_manager
        })
    
    logger.info(f"Added {len(compound_ids)} compounds to task queue")
    
    # Statistics with enhanced strategy information
    stats = {
        "total": len(compound_ids),
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": datetime.now(),
        "errors": {},
        "batch_times": [],
        "search_terms": args.search_terms if not args.fetch_all else ["all"],
        "fetch_all": args.fetch_all,
        "use_reference_seeds": args.use_reference_seeds,
        "similarity": args.similarity,
        "expand_similar": args.expand_similar,
        "filter_properties": args.filter_properties,
        "property_filters": {
            "min_mw": args.min_mw,
            "max_mw": args.max_mw,
            "min_hba": args.min_hba,
            "min_hbd": args.min_hbd,
            "max_alogp": args.max_alogp
        },
        "use_smarts": args.use_smarts,
        "smarts_classes": args.smarts_classes,
        "rdkit_available": RDKIT_AVAILABLE
    }
    
    # Enhanced batch processing with optimized checkpointing
    batch_size = args.batch_size
    
    # Create a function to process results in batches
    def process_results_batch():
        batch_start_time = time.time()
        batch_results = []
        batch_timeout = 300 * batch_size  # Timeout proportional to batch size
        
        try:
            for _ in range(batch_size):
                if len(batch_results) >= len(compound_ids):
                    break  # All compounds processed
                    
                try:
                    result = result_queue.get(timeout=300)  # 5 minute timeout per compound
                    batch_results.append(result)
                    result_queue.task_done()
                except Empty:
                    logger.warning("Timeout waiting for result")
                    break
            
            # Process all results in this batch using resumable batch processing
            # This ensures we can recover from failures and continue where we left off
            def process_result(result):
                # Update statistics
                if result.get("status") == "success":
                    # Update checkpoint with success
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=True
                    )
                    return True
                else:
                    error_category = result.get("error_category", "UNKNOWN")
                    
                    # Track error by category
                    if error_category not in stats["errors"]:
                        stats["errors"][error_category] = 0
                    stats["errors"][error_category] += 1
                    
                    # Update checkpoint with error
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=False,
                        error=result.get("error")
                    )
                    return False
            
            # Use safe transaction context for batch processing
            with safe_transaction():
                # Process each result in the batch
                for result in batch_results:
                    stats["processed"] += 1
                    if process_result(result):
                        stats["success"] += 1
                    else:
                        stats["error"] += 1
            
            # Calculate batch statistics
            batch_duration = time.time() - batch_start_time
            stats["batch_times"].append(batch_duration)
            
            compounds_per_second = len(batch_results) / batch_duration if batch_duration > 0 else 0
            logger.info(f"Batch completed: {len(batch_results)} compounds in {batch_duration:.2f}s " +
                      f"({compounds_per_second:.2f} compounds/s)")
            
            # Save checkpoint after each batch
            checkpoint_manager.save_checkpoint()
            
            # Display progress
            display_progress(stats)
            
            return len(batch_results)
        except Exception as e:
            logger.error(f"Error processing batch: {e}")
            # Use transaction_utils for error handling
            if is_transaction_active(None):
                logger.warning("Rolling back active transaction due to batch processing error")
            return 0
    
    try:
        # Process compounds in batches until all are done
        total_processed = 0
        max_batches = (len(compound_ids) + batch_size - 1) // batch_size
        
        logger.info(f"Processing {len(compound_ids)} compounds in up to {max_batches} batches")
        
        while total_processed < len(compound_ids):
            processed_in_batch = process_results_batch()
            total_processed += processed_in_batch
            
            if processed_in_batch == 0:
                logger.warning("No compounds processed in last batch, breaking")
                break
                
            # Calculate and log estimated time remaining
            if stats["processed"] > 0:
                elapsed = (datetime.now() - stats["start_time"]).total_seconds()
                items_per_second = stats["processed"] / elapsed
                remaining_items = len(compound_ids) - stats["processed"]
                eta_seconds = remaining_items / items_per_second if items_per_second > 0 else 0
                eta = str(timedelta(seconds=int(eta_seconds)))
                logger.info(f"Estimated time remaining: {eta}")
                
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
    finally:
        # Stop workers
        for worker in workers:
            task_queue.put(None)  # Send shutdown signal
        
        for worker in workers:
            worker.stop()
        
        # Final progress display
        display_progress(stats, final=True)
        
        # Save final checkpoint
        checkpoint_manager.save_checkpoint()
        
        # Generate summary report
        report_file = generate_summary_report(stats, args.mode)
        
        # Run verification if not in dry-run mode
        if args.mode != "dry-run" and stats["success"] > 0:
            try:
                from verify_imported_data import verify_chembl_import
                logger.info("Running verification...")
                
                # Ensure direct PostgreSQL connection for verification
                db = get_db()
                try:
                    # Test connection with a simple query
                    test_result = execute_query("SELECT 1 as connection_test")
                    if test_result and test_result[0]['connection_test'] == 1:
                        logger.info("Database connection verified for verification")
                        verify_chembl_import()
                    else:
                        logger.warning("Database connection test failed for verification")
                        logger.info("Verification will be skipped")
                except Exception as conn_error:
                    logger.warning(f"Failed to connect to database for verification: {conn_error}")
                    logger.info("Verification will be skipped")
            except ImportError:
                logger.warning("Verification module not available")
            except Exception as e:
                logger.error(f"Verification failed: {e}")
    
    return 0

def display_progress(stats, final=False):
    """
    Display progress information for the import process.
    
    Args:
        stats: Dictionary containing progress statistics
        final: Whether this is the final progress display
    """
    # Calculate progress percentage
    if stats["total"] > 0:
        progress_pct = (stats["processed"] / stats["total"]) * 100
    else:
        progress_pct = 0
        
    # Calculate elapsed time
    elapsed = datetime.now() - stats["start_time"]
    
    # Calculate estimated time remaining
    if stats["processed"] > 0 and not final:
        items_per_second = stats["processed"] / elapsed.total_seconds()
        remaining_items = stats["total"] - stats["processed"]
        eta_seconds = remaining_items / items_per_second if items_per_second > 0 else 0
        eta = str(timedelta(seconds=int(eta_seconds)))
    else:
        eta = "completed" if final else "unknown"
        
    # Create progress bar
    bar_length = 30
    filled_length = int(bar_length * stats["processed"] // stats["total"]) if stats["total"] > 0 else 0
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    # Calculate success rate
    success_rate = (stats["success"] / stats["processed"]) * 100 if stats["processed"] > 0 else 0
    
    # Log progress
    if final:
        logger.info(f"FINAL RESULTS:")
        logger.info(
            f"Progress: [{bar}] {stats['processed']}/{stats['total']} "
            f"({progress_pct:.1f}%) | Success: {stats['success']} ({success_rate:.1f}%) | "
            f"Errors: {stats['error']} | Total time: {elapsed}"
        )
    else:
        logger.info(
            f"Progress: [{bar}] {stats['processed']}/{stats['total']} "
            f"({progress_pct:.1f}%) | Success: {stats['success']} | "
            f"Errors: {stats['error']} | Elapsed: {elapsed} | ETA: {eta}"
        )

def generate_summary_report(stats, mode):
    """
    Generate a summary report of the import process.
    
    Args:
        stats: Dictionary containing import statistics
        mode: Execution mode
        
    Returns:
        Path to the generated report file
    """
    # Create reports directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Generate report filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = f"reports/chembl_import_{timestamp}.json"
    
    # Calculate additional metrics
    elapsed_seconds = (datetime.now() - stats["start_time"]).total_seconds()
    compounds_per_second = stats["processed"] / elapsed_seconds if elapsed_seconds > 0 else 0
    average_batch_time = sum(stats["batch_times"]) / len(stats["batch_times"]) if stats["batch_times"] else 0
    success_rate = stats["success"] / stats["processed"] * 100 if stats["processed"] > 0 else 0
    
    # Add additional report data with enhanced strategy information
    report_data = {
        "timestamp": datetime.now().isoformat(),
        "mode": mode,
        "stats": {
            "total": stats["total"],
            "processed": stats["processed"],
            "success": stats["success"],
            "error": stats["error"],
            "success_rate_percent": success_rate,
            "errors_by_category": stats["errors"]
        },
        "performance": {
            "elapsed_seconds": elapsed_seconds,
            "compounds_per_second": compounds_per_second,
            "average_batch_time": average_batch_time
        },
        "discovery_strategy": {
            "search_terms": stats["search_terms"],
            "fetch_all": stats.get("fetch_all", False),
            "use_reference_seeds": stats.get("use_reference_seeds", False),
            "similarity_threshold": stats.get("similarity", 70),
            "expand_similar": stats.get("expand_similar", False)
        },
        "filtering_strategy": {
            "filter_properties": stats.get("filter_properties", False),
            "property_filters": stats.get("property_filters", {}),
            "use_smarts": stats.get("use_smarts", False),
            "smarts_classes": stats.get("smarts_classes", []),
            "rdkit_available": stats.get("rdkit_available", False)
        },
        "start_time": stats["start_time"].isoformat(),
        "end_time": datetime.now().isoformat()
    }
    
    # Save report
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2, default=str)
    
    logger.info(f"Summary report saved to {report_file}")
    
    # Also log summary to console with enhanced strategy information
    logger.info(f"Import Summary:")
    logger.info(f"  Mode: {mode}")
    logger.info(f"  Total compounds: {stats['total']}")
    logger.info(f"  Processed: {stats['processed']}")
    logger.info(f"  Success: {stats['success']} ({success_rate:.1f}%)")
    logger.info(f"  Errors: {stats['error']}")
    logger.info(f"  Total time: {timedelta(seconds=int(elapsed_seconds))}")
    logger.info(f"  Performance: {compounds_per_second:.2f} compounds/second")
    
    # Log discovery strategy
    logger.info(f"  Discovery Strategy:")
    if stats.get("use_reference_seeds", False):
        logger.info(f"    Used reference compounds as seeds for similarity search")
        logger.info(f"    Similarity threshold: {stats.get('similarity', 70)}%")
    elif stats.get("fetch_all", False):
        logger.info(f"    Fetched all ChEMBL compounds")
    else:
        logger.info(f"    Search terms: {stats['search_terms']}")
        if stats.get("expand_similar", False):
            logger.info(f"    Expanded with similar compounds (similarity >= {stats.get('similarity', 70)}%)")
    
    # Log filtering strategy
    logger.info(f"  Filtering Strategy:")
    if stats.get("filter_properties", False) and stats.get("rdkit_available", False):
        prop_filters = stats.get("property_filters", {})
        logger.info(f"    Molecular property filtering: Enabled")
        logger.info(f"    MW range: {prop_filters.get('min_mw', 30)}-{prop_filters.get('max_mw', 500)}")
        logger.info(f"    Min HBA: {prop_filters.get('min_hba', 2)}")
        logger.info(f"    Min HBD: {prop_filters.get('min_hbd', 1)}")
        logger.info(f"    Max ALogP: {prop_filters.get('max_alogp', 3.0)}")
    else:
        logger.info(f"    Molecular property filtering: Disabled")
    
    if stats.get("use_smarts", False) and stats.get("rdkit_available", False):
        logger.info(f"    SMARTS pattern matching: Enabled")
        logger.info(f"    Chemical classes: {stats.get('smarts_classes', ['all'])}")
    else:
        logger.info(f"    SMARTS pattern matching: Disabled")
    
    # Log error categories if any
    if stats["errors"]:
        logger.info("  Error categories:")
        for category, count in stats["errors"].items():
            logger.info(f"    {category}: {count}")
    
    return report_file

def filter_by_molecular_properties(chembl_client, compound_ids,
                                   min_mw=30, max_mw=500,
                                   min_hba=2, min_hbd=1,
                                   max_alogp=3.0):
    """
    Filter compounds based on molecular properties common to cryoprotectants.
    
    This function now uses the enhanced chembl_search_utils module for more
    reliable property-based filtering with improved error handling and caching.
    
    Args:
        chembl_client: ChEMBL client instance
        compound_ids: List of ChEMBL IDs to filter
        min_mw: Minimum molecular weight
        max_mw: Maximum molecular weight
        min_hba: Minimum number of hydrogen bond acceptors
        min_hbd: Minimum number of hydrogen bond donors
        max_alogp: Maximum ALogP value
        
    Returns:
        List of filtered ChEMBL IDs
    """
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available. Skipping molecular property filtering.")
        return compound_ids
    
    logger.info(f"Filtering {len(compound_ids)} compounds based on molecular properties")
    
    try:
        # Use the enhanced find_potential_cryoprotectants function from chembl_search_utils
        # This provides better caching, error handling, and resilience
        potential_cryoprotectants = find_potential_cryoprotectants(
            limit=len(compound_ids) * 2,  # Request more to ensure we get enough matches
            min_mw=min_mw,
            max_mw=max_mw,
            min_hba=min_hba,
            max_hba=20,  # Set a reasonable upper limit
            min_hbd=min_hbd,
            max_hbd=10,  # Set a reasonable upper limit
            min_logp=min_alogp * -1,  # Convert to range
            max_logp=max_alogp,
            use_cache=True,
            fallback_to_cache=True
        )
        
        # Extract ChEMBL IDs from the results
        potential_ids = set(comp.get("molecule_chembl_id") for comp in potential_cryoprotectants
                          if comp.get("molecule_chembl_id"))
        
        # Filter the original list to only include compounds that match our criteria
        filtered_ids = [cid for cid in compound_ids if cid in potential_ids]
        
        logger.info(f"Property filtering using chembl_search_utils: {len(filtered_ids)}/{len(compound_ids)} compounds passed")
        
        # If we didn't get any matches or very few, fall back to the original method
        if len(filtered_ids) < len(compound_ids) * 0.1:  # Less than 10% match rate
            logger.warning(f"Low match rate with enhanced filtering ({len(filtered_ids)}/{len(compound_ids)}). Falling back to direct filtering.")
            
            # Fallback to original implementation
            filtered_ids = []
            for chembl_id in compound_ids:
                try:
                    # Get compound data
                    compound_data = chembl_client.get_compound(chembl_id)
                    
                    # Extract SMILES
                    smiles = None
                    if 'molecule_structures' in compound_data and 'canonical_smiles' in compound_data['molecule_structures']:
                        smiles = compound_data['molecule_structures']['canonical_smiles']
                    
                    if not smiles:
                        logger.debug(f"No SMILES available for {chembl_id}, skipping property filtering")
                        continue
                        
                    # Create RDKit molecule
                    mol = Chem.MolFromSmiles(smiles)
                    if not mol:
                        logger.debug(f"Could not create RDKit molecule for {chembl_id}, skipping property filtering")
                        continue
                        
                    # Calculate properties
                    mw = Descriptors.MolWt(mol)
                    hba = Lipinski.NumHAcceptors(mol)
                    hbd = Lipinski.NumHDonors(mol)
                    alogp = Descriptors.MolLogP(mol)
                    
                    # Apply filters
                    if (min_mw <= mw <= max_mw and
                        hba >= min_hba and
                        hbd >= min_hbd and
                        alogp <= max_alogp):
                        filtered_ids.append(chembl_id)
                        logger.debug(f"Compound {chembl_id} passed property filters: MW={mw:.1f}, HBA={hba}, HBD={hbd}, ALogP={alogp:.1f}")
                    else:
                        logger.debug(f"Compound {chembl_id} failed property filters: MW={mw:.1f}, HBA={hba}, HBD={hbd}, ALogP={alogp:.1f}")
                        
                except Exception as e:
                    logger.warning(f"Error filtering compound {chembl_id}: {e}")
                    # Continue with next compound
                    continue
            
            logger.info(f"Property filtering with fallback method: {len(filtered_ids)}/{len(compound_ids)} compounds passed")
        
        return filtered_ids
        
    except Exception as e:
        logger.error(f"Error using enhanced property filtering: {e}")
        logger.warning("Falling back to original compound list due to filtering error")
        return compound_ids

def match_smarts_patterns(chembl_client, compound_ids):
    """
    Match compounds against SMARTS patterns for specific chemical classes
    relevant to cryoprotectants.
    
    This function now uses the enhanced chembl_search_utils module for more
    reliable chemical class filtering with improved error handling and caching.
    
    Args:
        chembl_client: ChEMBL client instance
        compound_ids: List of ChEMBL IDs to match
        
    Returns:
        Dictionary mapping chemical classes to lists of matching ChEMBL IDs
    """
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available. Skipping SMARTS pattern matching.")
        return {}
    
    logger.info(f"Matching {len(compound_ids)} compounds against SMARTS patterns")
    
    try:
        # Use the enhanced identify_compounds_by_chemical_class function from chembl_search_utils
        # This provides better caching, error handling, and resilience
        compounds_by_class = identify_compounds_by_chemical_class(
            use_cache=True,
            fallback_to_cache=True,
            limit_per_class=len(compound_ids)  # Set limit to ensure we get enough matches
        )
        
        # Filter the results to only include compounds from our original list
        filtered_matches = {}
        for class_name, compounds in compounds_by_class.items():
            # Extract ChEMBL IDs from the results
            class_ids = [comp.get("molecule_chembl_id") for comp in compounds if comp.get("molecule_chembl_id")]
            # Filter to only include compounds from our original list
            filtered_matches[class_name] = [cid for cid in class_ids if cid in compound_ids]
            logger.info(f"Chemical class {class_name}: {len(filtered_matches[class_name])}/{len(compound_ids)} compounds match")
        
        # If we didn't get any matches or very few, fall back to the original method
        total_matches = sum(len(matches) for matches in filtered_matches.values())
        if total_matches < len(compound_ids) * 0.1:  # Less than 10% match rate
            logger.warning(f"Low match rate with enhanced filtering ({total_matches} total matches). Falling back to direct matching.")
            
            # Define SMARTS patterns for relevant chemical classes
            patterns = {
                "polyol": "[OX2H][CX4][CX4][OX2H]", # Diol pattern
                "sugar": "[OX2H][CX4][CX4;$([CX4][OX2H])][CX4;$([CX4][OX2H])][CX4][OX2H]", # Simple sugar pattern
                "amide": "[NX3][CX3]=[OX1]", # Amide pattern
                "amine": "[NX3;H2,H1,H0;!$(NC=O)]", # Amine pattern
                "alcohol": "[OX2H]", # Alcohol pattern
                "ether": "[OX2]([CX4])[CX4]" # Ether pattern
            }
            
            # Compile patterns
            compiled_patterns = {}
            for name, smarts in patterns.items():
                try:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern:
                        compiled_patterns[name] = pattern
                    else:
                        logger.warning(f"Could not compile SMARTS pattern for {name}: {smarts}")
                except Exception as e:
                    logger.warning(f"Error compiling SMARTS pattern for {name}: {e}")
            
            # Initialize results
            matches = {name: [] for name in patterns.keys()}
            
            for chembl_id in compound_ids:
                try:
                    # Get compound data
                    compound_data = chembl_client.get_compound(chembl_id)
                    
                    # Extract SMILES
                    smiles = None
                    if 'molecule_structures' in compound_data and 'canonical_smiles' in compound_data['molecule_structures']:
                        smiles = compound_data['molecule_structures']['canonical_smiles']
                    
                    if not smiles:
                        logger.debug(f"No SMILES available for {chembl_id}, skipping SMARTS matching")
                        continue
                        
                    # Create RDKit molecule
                    mol = Chem.MolFromSmiles(smiles)
                    if not mol:
                        logger.debug(f"Could not create RDKit molecule for {chembl_id}, skipping SMARTS matching")
                        continue
                        
                    # Match against patterns
                    for name, pattern in compiled_patterns.items():
                        if mol.HasSubstructMatch(pattern):
                            matches[name].append(chembl_id)
                            logger.debug(f"Compound {chembl_id} matches {name} pattern")
                            
                except Exception as e:
                    logger.warning(f"Error matching SMARTS patterns for compound {chembl_id}: {e}")
                    # Continue with next compound
                    continue
            
            # Log summary
            for name, matched_ids in matches.items():
                logger.info(f"SMARTS matching with fallback method: {len(matched_ids)} compounds match {name} pattern")
            
            return matches
        
        return filtered_matches
        
    except Exception as e:
        logger.error(f"Error using enhanced chemical class filtering: {e}")
        logger.warning("Falling back to empty matches dictionary due to filtering error")
        return {
            "polyol": [],
            "sugar": [],
            "amide": [],
            "amine": [],
            "alcohol": [],
            "ether": []
        }

def expand_search_results(chembl_client, search_results, similarity=70, limit=500):
    """
    Expand search results by adding similar compounds.
    Uses the enhanced chembl_search_utils module for more reliable similarity search.
    
    Args:
        chembl_client: ChEMBL client instance
        search_results: Initial search results (ChEMBL IDs)
        similarity: Minimum similarity percentage (0-100)
        limit: Maximum number of compounds to return
        
    Returns:
        Expanded list of ChEMBL IDs
    """
    expanded_results = set(search_results)
    
    try:
        # Get SMILES for each compound in search_results
        reference_smiles = []
        for chembl_id in search_results:
            try:
                compound_data = chembl_client.get_compound(chembl_id)
                smiles = compound_data.get('molecule_structures', {}).get('canonical_smiles')
                if smiles:
                    reference_smiles.append(smiles)
            except Exception as e:
                logger.warning(f"Error getting SMILES for compound {chembl_id}: {e}")
                continue
        
        if reference_smiles:
            # Use the enhanced find_similar_compounds function from chembl_search_utils
            similar_compounds = find_similar_compounds(
                reference_smiles=reference_smiles,
                similarity_threshold=similarity,
                limit=limit,
                use_cache=True,
                fallback_to_cache=True
            )
            
            # Extract ChEMBL IDs from the results
            for compound in similar_compounds:
                chembl_id = compound.get("molecule_chembl_id")
                if chembl_id:
                    expanded_results.add(chembl_id)
            
            logger.info(f"Found {len(expanded_results) - len(search_results)} additional compounds using similarity search")
        else:
            logger.warning("No valid SMILES found for similarity search")
            
        # Fallback to original method if we didn't get many new compounds
        if len(expanded_results) < len(search_results) * 1.2:  # Less than 20% new compounds
            logger.warning(f"Low expansion rate with enhanced similarity search. Falling back to direct method.")
            
            # For each compound in the initial results
            for chembl_id in search_results:
                try:
                    # Find similar compounds using the client directly
                    similar_compounds = chembl_client.get_similar_compounds(chembl_id, similarity=similarity)
                    
                    # Add to expanded results
                    expanded_results.update(similar_compounds)
                    
                    # Stop if limit reached
                    if len(expanded_results) >= limit:
                        break
                except Exception as e:
                    logger.warning(f"Error expanding results for compound {chembl_id}: {e}")
                    # Continue with next compound
                    continue
            
            logger.info(f"Found {len(expanded_results) - len(search_results)} additional compounds using fallback similarity search")
    
    except Exception as e:
        logger.error(f"Error in similarity search expansion: {e}")
        logger.warning("Using original search results due to expansion error")
    
    return list(expanded_results)[:limit]

def get_compounds_from_reference_seeds(chembl_client, similarity=70, limit=500):
    """
    Get compounds similar to reference cryoprotectant compounds.
    
    Args:
        chembl_client: ChEMBL client instance
        similarity: Minimum similarity percentage (0-100)
        limit: Maximum number of compounds to return
        
    Returns:
        List of ChEMBL IDs
    """
    from chembl.reference_compounds import get_reference_compound_ids, get_extended_reference_compound_ids
    
    # Get reference compounds
    reference_ids = get_extended_reference_compound_ids()
    logger.info(f"Using {len(reference_ids)} reference compounds as seeds for similarity search")
    
    # Expand with similar compounds
    similar_compounds = expand_search_results(chembl_client, reference_ids, similarity=similarity, limit=limit)
    logger.info(f"Found {len(similar_compounds)} compounds similar to reference compounds")
    
    return similar_compounds

if __name__ == "__main__":
    sys.exit(main())
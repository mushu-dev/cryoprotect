#!/usr/bin/env python3
"""
ChEMBL Data Import Script

This script imports molecular data from ChEMBL, focusing on compounds similar to
the reference cryoprotectants. It ensures property completeness during import
and focuses on quality over quantity.

Enhanced with:
- Resilient ChEMBL API client with adaptive rate limiting and circuit breaking
- Comprehensive error handling and retry logic
- Checkpointing system for resumable operations
- Supabase MCP integration for database operations
- Progress tracking and reporting

Usage:
    python chembl_import.py

Author: Roo Agent
Date: 2025-05-01
"""

import json
import logging
import sys
import os
import time
import requests
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
import random
import traceback
import hashlib
import threading
from pathlib import Path
from queue import Queue, Empty

# Import ChEMBL modules
from chembl.client import ResilientChEMBLClient
from chembl.checkpoint import CheckpointManager
from chembl.error_handler import ErrorCategory, classify_error, get_recovery_strategy
from chembl.worker import ChEMBLWorker

# Configure logging
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/chembl_import.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Constants
BATCH_SIZE = 50
MAX_RETRIES = 5
CHECKPOINT_INTERVAL = 50  # molecules
MAX_WORKERS = 4  # Number of parallel workers

# Search terms for cryoprotectants
SEARCH_TERMS = [
    "cryoprotectant", 
    "cryoprotection", 
    "cryopreservation", 
    "antifreeze", 
    "freeze protection",
    "polyol", 
    "glycol",
    "glycerol",
    "dimethyl sulfoxide",
    "dmso",
    "trehalose",
    "sucrose",
    "glucose",
    "mannitol",
    "sorbitol",
    "ethylene glycol",
    "propylene glycol"
]

# Similarity search SMILES for reference compounds
REFERENCE_SMILES = [
    "C(C(CO)O)O",  # Glycerol
    "CS(=O)C",  # DMSO
    "C(CC(=O)O)N",  # beta-Alanine
    "CC(C)(C)O",  # tert-Butanol
    "C(=O)(N)N",  # Urea
    "C(CO)O",  # Ethylene glycol
    "CC(CO)O",  # Propylene glycol
    "C(C1C(C(C(C(O1)OC2C(OC(C2O)O)C(O)C(O)C(O)CO)O)O)O)O",  # Trehalose
    "C(C(=O)O)N"  # Glycine
]

def execute_sql(project_id: str, query: str, params: Optional[List[Any]] = None) -> List[Dict[str, Any]]:
    """
    Execute SQL query using Supabase MCP tool.
    
    Args:
        project_id: The Supabase project ID
        query: The SQL query to execute
        params: Optional parameters for the query
        
    Returns:
        A list of dictionaries representing the query results
    """
    try:
        from use_mcp_tool import use_mcp_tool
        
        arguments = {
            "project_id": project_id,
            "query": query
        }
        
        if params:
            arguments["params"] = params
            
        result = use_mcp_tool("supabase", "execute_sql", arguments)
        return result
    except Exception as e:
        error_msg = f"Error executing SQL query: {str(e)}"
        logger.error(error_msg)
        logger.debug(traceback.format_exc())
        return []

def execute_sql_with_retry(project_id: str, query: str, params: Optional[List[Any]] = None,
                          max_retries: int = 3, retry_delay: float = 1.0) -> List[Dict[str, Any]]:
    """
    Execute SQL query with retry logic.
    
    Args:
        project_id: The Supabase project ID
        query: The SQL query to execute
        params: Optional parameters for the query
        max_retries: Maximum number of retry attempts
        retry_delay: Initial delay between retries (will increase exponentially)
        
    Returns:
        A list of dictionaries representing the query results
    """
    for attempt in range(max_retries + 1):
        try:
            result = execute_sql(project_id, query, params)
            return result
        except Exception as e:
            if attempt < max_retries:
                # Calculate backoff with jitter
                backoff = retry_delay * (2 ** attempt) * (0.5 + random.random())
                logger.warning(f"SQL query failed (attempt {attempt + 1}/{max_retries + 1}): {str(e)}")
                logger.info(f"Retrying in {backoff:.2f} seconds...")
                time.sleep(backoff)
            else:
                logger.error(f"SQL query failed after {max_retries + 1} attempts: {str(e)}")
                raise

def get_chembl_client() -> ResilientChEMBLClient:
    """
    Get a configured ChEMBL client instance.
    
    Returns:
        A configured ResilientChEMBLClient instance
    """
    # Create cache directory
    os.makedirs("cache/chembl", exist_ok=True)
    
    # Configure client with appropriate rate limits
    client = ResilientChEMBLClient(
        cache_dir="cache/chembl",
        weekday_requests_per_second=3.0,
        weekend_requests_per_second=4.5,
        monday_requests_per_second=2.0,
        max_retries=MAX_RETRIES,
        failure_threshold=3,
        recovery_timeout=60,
        cache_ttl=86400 * 30,  # 30 days
        error_cache_ttl=3600,  # 1 hour
        memory_cache_size=1000
    )
    
    return client

def search_compounds_by_term(client: ResilientChEMBLClient, term: str, limit: int = 1000) -> List[Dict[str, Any]]:
    """
    Search for compounds in ChEMBL using a search term.
    
    Args:
        client: ResilientChEMBLClient instance
        term: The search term
        limit: Maximum number of results to return
        
    Returns:
        A list of compound data dictionaries
    """
    logger.info(f"Searching for compounds with term: {term}")
    
    try:
        # Use the resilient client to search for molecules
        results = client.search_molecules(
            query=term,
            limit=limit,
            use_cache=True,
            fallback_to_cache=True
        )
        
        molecules = results.get("Molecules", [])
        logger.info(f"Found {len(molecules)} compounds for term: {term}")
        return molecules
    except Exception as e:
        logger.error(f"Error searching for compounds with term '{term}': {str(e)}")
        logger.debug(traceback.format_exc())
        return []

def search_compounds_by_similarity(client: ResilientChEMBLClient, smiles: str, similarity: float = 70, limit: int = 500) -> List[Dict[str, Any]]:
    """
    Search for compounds in ChEMBL by structural similarity.
    
    Args:
        client: ResilientChEMBLClient instance
        smiles: The SMILES string to use as a reference
        similarity: The minimum similarity percentage (0-100)
        limit: Maximum number of results to return
        
    Returns:
        A list of compound data dictionaries
    """
    logger.info(f"Searching for compounds similar to: {smiles}")
    
    try:
        # Get similar compound IDs
        similar_ids = client.get_similar_compounds(
            chembl_id=smiles,  # This can be a SMILES string or ChEMBL ID
            similarity=int(similarity),
            limit=limit,
            use_cache=True,
            fallback_to_cache=True
        )
        
        # Fetch full compound data for each ID
        molecules = []
        for chembl_id in similar_ids:
            try:
                molecule_data = client.get_molecule_by_chembl_id(
                    chembl_id=chembl_id,
                    use_cache=True,
                    fallback_to_cache=True
                )
                
                if "Error" not in molecule_data:
                    molecules.append(molecule_data)
            except Exception as e:
                logger.warning(f"Error fetching data for similar compound {chembl_id}: {str(e)}")
                continue
        
        logger.info(f"Found {len(molecules)} compounds similar to: {smiles}")
        return molecules
    except Exception as e:
        logger.error(f"Error searching for compounds similar to '{smiles}': {str(e)}")
        logger.debug(traceback.format_exc())
        return []

def extract_molecule_data(molecule_data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Extract relevant data from a ChEMBL molecule record.
    
    Args:
        molecule_data: The ChEMBL molecule data
        
    Returns:
        A dictionary with the extracted data, or None if required fields are missing
    """
    try:
        # Skip if missing required fields
        chembl_id = molecule_data.get("ChEMBL ID") or molecule_data.get("molecule_chembl_id")
        if not chembl_id:
            return None
            
        # Get SMILES from different possible locations
        smiles = None
        if "molecule_structures" in molecule_data and molecule_data["molecule_structures"]:
            smiles = molecule_data["molecule_structures"].get("canonical_smiles")
        elif "SMILES" in molecule_data:
            smiles = molecule_data["SMILES"]
        elif "smiles" in molecule_data:
            smiles = molecule_data["smiles"]
        elif "canonical_smiles" in molecule_data:
            smiles = molecule_data["canonical_smiles"]
            
        if not smiles:
            return None
            
        # Extract basic properties
        molecule = {
            "chembl_id": chembl_id,
            "name": molecule_data.get("Name") or molecule_data.get("pref_name", chembl_id),
            "smiles": smiles,
            "inchi": molecule_data.get("InChI") or
                    (molecule_data.get("molecule_structures", {}) or {}).get("standard_inchi", ""),
            "inchikey": molecule_data.get("InChIKey") or
                       (molecule_data.get("molecule_structures", {}) or {}).get("standard_inchi_key", ""),
            "formula": molecule_data.get("Molecular Formula") or
                      (molecule_data.get("molecule_properties", {}) or {}).get("full_molformula", ""),
            "molecular_weight": molecule_data.get("Molecular Weight") or
                               (molecule_data.get("molecule_properties", {}) or {}).get("full_mwt"),
            "data_source": "ChEMBL"
        }
        
        # Extract properties
        properties = {}
        
        # Try to get properties from different possible locations
        mol_props = {}
        if "molecule_properties" in molecule_data and molecule_data["molecule_properties"]:
            mol_props = molecule_data["molecule_properties"]
        
        # Map ChEMBL properties to our property names
        property_mappings = {
            "alogp": "LogP",
            "LogP": "LogP",
            "hba": "Hydrogen Bond Acceptor Count",
            "H-Bond Acceptors": "Hydrogen Bond Acceptor Count",
            "hbd": "Hydrogen Bond Donor Count",
            "H-Bond Donors": "Hydrogen Bond Donor Count",
            "rtb": "Rotatable Bond Count",
            "Rotatable Bonds": "Rotatable Bond Count",
            "psa": "Polar Surface Area",
            "TPSA": "Polar Surface Area"
        }
        
        # Extract properties from molecule_properties or top-level
        for source_prop, target_prop in property_mappings.items():
            if source_prop in mol_props:
                properties[target_prop] = mol_props[source_prop]
            elif source_prop in molecule_data:
                properties[target_prop] = molecule_data[source_prop]
        
        # Only include molecules with all required properties
        required_properties = ["LogP", "Hydrogen Bond Acceptor Count", "Hydrogen Bond Donor Count", "Rotatable Bond Count"]
        if all(prop in properties for prop in required_properties):
            molecule["properties"] = properties
            return molecule
        else:
            missing_props = [prop for prop in required_properties if prop not in properties]
            logger.debug(f"Molecule {chembl_id} missing required properties: {missing_props}")
            return None
            
    except Exception as e:
        logger.warning(f"Error extracting molecule data: {str(e)}")
        return None

def insert_molecule(project_id: str, molecule: Dict[str, Any]) -> Optional[str]:
    """
    Insert a molecule into the database.
    
    Args:
        project_id: The Supabase project ID
        molecule: The molecule data
        
    Returns:
        The ID of the inserted molecule, or None if insertion failed
    """
    query = """
    INSERT INTO molecules 
        (chembl_id, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
    VALUES 
        ($1, $2, $3, $4, $5, $6, $7, $8)
    ON CONFLICT (chembl_id) DO UPDATE SET
        name = EXCLUDED.name,
        smiles = EXCLUDED.smiles,
        inchi = EXCLUDED.inchi,
        inchikey = EXCLUDED.inchikey,
        formula = EXCLUDED.formula,
        molecular_weight = EXCLUDED.molecular_weight,
        updated_at = NOW()
    RETURNING id;
    """
    
    params = [
        molecule["chembl_id"],
        molecule["name"],
        molecule["smiles"],
        molecule["inchi"],
        molecule["inchikey"],
        molecule["formula"],
        molecule["molecular_weight"],
        molecule["data_source"]
    ]
    
    try:
        result = execute_sql_with_retry(project_id, query, params)
        
        if result and len(result) > 0:
            return result[0]["id"]
        else:
            logger.error(f"Failed to insert molecule {molecule['chembl_id']}")
            return None
    except Exception as e:
        logger.error(f"Error inserting molecule {molecule['chembl_id']}: {str(e)}")
        logger.debug(traceback.format_exc())
        return None

def insert_property(project_id: str, molecule_id: str, property_name: str, property_value: Any) -> bool:
    """
    Insert a property for a molecule.
    
    Args:
        project_id: The Supabase project ID
        molecule_id: The ID of the molecule
        property_name: The name of the property
        property_value: The value of the property
        
    Returns:
        True if the property was inserted successfully, False otherwise
    """
    try:
        # Get the property type ID
        prop_type_query = """
        SELECT id FROM property_types
        WHERE name = $1
        """
        
        prop_type_result = execute_sql_with_retry(project_id, prop_type_query, [property_name])
        
        if not prop_type_result or len(prop_type_result) == 0:
            logger.error(f"Property type '{property_name}' not found")
            return False
        
        property_type_id = prop_type_result[0]["id"]
        
        # Determine the value type
        if isinstance(property_value, (int, float)):
            value_column = "numeric_value"
        elif isinstance(property_value, bool):
            value_column = "boolean_value"
        else:
            value_column = "text_value"
        
        query = f"""
        INSERT INTO molecular_properties
            (molecule_id, property_type_id, {value_column}, source)
        VALUES
            ($1, $2, $3, $4)
        ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
            {value_column} = EXCLUDED.{value_column},
            updated_at = NOW()
        RETURNING id;
        """
        
        params = [
            molecule_id,
            property_type_id,
            property_value,
            "ChEMBL"
        ]
        
        result = execute_sql_with_retry(project_id, query, params)
        
        if result and len(result) > 0:
            return True
        else:
            logger.error(f"Failed to insert property {property_name} for molecule {molecule_id}")
            return False
    except Exception as e:
        logger.error(f"Error inserting property {property_name} for molecule {molecule_id}: {str(e)}")
        logger.debug(traceback.format_exc())
        return False

def process_molecule_batch(project_id: str, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Process a batch of molecules and import them into the database.
    
    Args:
        project_id: The Supabase project ID
        molecules: List of molecule data from ChEMBL API
        
    Returns:
        Statistics about the batch processing
    """
    results = {
        "molecules_imported": 0,
        "properties_imported": 0,
        "errors": []
    }
    
    for molecule_data in molecules:
        try:
            # Extract molecule data
            molecule = extract_molecule_data(molecule_data)
            
            if not molecule:
                continue
                
            # Insert molecule
            molecule_id = insert_molecule(project_id, molecule)
            
            if not molecule_id:
                continue
                
            results["molecules_imported"] += 1
            
            # Insert properties
            for prop_name, value in molecule["properties"].items():
                if insert_property(project_id, molecule_id, prop_name, value):
                    results["properties_imported"] += 1
                else:
                    results["errors"].append({
                        "molecule_id": molecule_id,
                        "property": prop_name,
                        "error": "Failed to insert property"
                    })
                    
        except Exception as e:
            logger.error(f"Error processing molecule {molecule_data.get('molecule_chembl_id', 'unknown')}: {str(e)}")
            results["errors"].append({
                "molecule": molecule_data.get("molecule_chembl_id", "unknown"),
                "error": str(e)
            })
    
    return results

def create_checkpoint(project_id: str, data: Dict[str, Any], checkpoint_manager: Optional[CheckpointManager] = None) -> None:
    """
    Create a checkpoint to track import progress.
    
    Args:
        project_id: The Supabase project ID
        data: Checkpoint data to store
        checkpoint_manager: Optional CheckpointManager instance
    """
    # Save to database
    try:
        checkpoint_sql = """
        INSERT INTO import_checkpoints
            (import_type, checkpoint_data, created_at)
        VALUES
            ('chembl', $1, NOW())
        RETURNING id;
        """
        
        execute_sql_with_retry(project_id, checkpoint_sql, [json.dumps(data)])
        
        # Also save to a local file
        os.makedirs("checkpoints", exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        checkpoint_file = f"checkpoints/chembl_import_{timestamp}.json"
        
        with open(checkpoint_file, "w") as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Created checkpoint: {checkpoint_file}")
        
        # If checkpoint manager is provided, save using it as well
        if checkpoint_manager:
            checkpoint_manager.state.update(data)
            checkpoint_manager.save_checkpoint()
            
    except Exception as e:
        logger.error(f"Error creating checkpoint: {str(e)}")
        logger.debug(traceback.format_exc())
        
        # Try to save locally even if database save failed
        try:
            os.makedirs("checkpoints", exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            checkpoint_file = f"checkpoints/chembl_import_fallback_{timestamp}.json"
            
            with open(checkpoint_file, "w") as f:
                json.dump(data, f, indent=2)
                
            logger.info(f"Created fallback checkpoint: {checkpoint_file}")
        except Exception as e2:
            logger.error(f"Error creating fallback checkpoint: {str(e2)}")

def load_latest_checkpoint() -> Optional[Dict[str, Any]]:
    """
    Load the latest checkpoint file.
    
    Returns:
        Checkpoint data or None if no checkpoint exists
    """
    try:
        checkpoint_dir = Path("checkpoints")
        if not checkpoint_dir.exists():
            return None
            
        # Find all checkpoint files
        checkpoint_files = list(checkpoint_dir.glob("chembl_import_*.json"))
        if not checkpoint_files:
            return None
            
        # Sort by modification time (newest first)
        checkpoint_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        
        # Load the latest checkpoint
        with open(checkpoint_files[0], "r") as f:
            data = json.load(f)
            
        logger.info(f"Loaded checkpoint: {checkpoint_files[0]}")
        return data
    except Exception as e:
        logger.error(f"Error loading checkpoint: {str(e)}")
        return None

def import_chembl_data(project_id: str, target_count: int = 5000, resume: bool = True) -> Dict[str, Any]:
    """
    Import chemical data from ChEMBL focusing on cryoprotectant compounds.
    Enhanced with parallel processing, checkpointing, and resilient API client.
    
    Args:
        project_id: The Supabase project ID
        target_count: Target number of compounds to import
        resume: Whether to resume from the latest checkpoint
        
    Returns:
        Statistics about the import operation
    """
    logger.info(f"Starting ChEMBL data import with target of {target_count} molecules")
    
    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager("checkpoints/chembl")
    
    # Initialize results
    results = {
        "molecules_imported": 0,
        "properties_imported": 0,
        "errors": [],
        "batches_processed": 0,
        "search_terms_used": [],
        "similarity_searches": [],
        "start_time": datetime.now().isoformat(),
        "processed_compounds": set()
    }
    
    # Try to load checkpoint if resume is enabled
    if resume:
        loaded_checkpoint = False
        
        # First try to load from checkpoint manager
        if checkpoint_manager.load_checkpoint():
            logger.info("Resuming from checkpoint manager")
            
            # Update results from checkpoint
            results["molecules_imported"] = checkpoint_manager.state.get("total_processed", 0)
            results["properties_imported"] = checkpoint_manager.state.get("properties_imported", 0)
            results["batches_processed"] = checkpoint_manager.state.get("batches_processed", 0)
            results["search_terms_used"] = checkpoint_manager.state.get("search_terms_used", [])
            results["similarity_searches"] = checkpoint_manager.state.get("similarity_searches", [])
            results["processed_compounds"] = set(checkpoint_manager.state.get("processed_compounds", []))
            
            loaded_checkpoint = True
        else:
            # Try to load from file
            checkpoint_data = load_latest_checkpoint()
            if checkpoint_data:
                logger.info("Resuming from checkpoint file")
                
                # Update results from checkpoint
                results["molecules_imported"] = checkpoint_data.get("molecules_imported", 0)
                results["properties_imported"] = checkpoint_data.get("properties_imported", 0)
                results["batches_processed"] = checkpoint_data.get("batches_processed", 0)
                results["search_terms_used"] = checkpoint_data.get("search_terms_used", [])
                results["similarity_searches"] = checkpoint_data.get("similarity_searches", [])
                results["processed_compounds"] = set(checkpoint_data.get("processed_compounds", []))
                
                # Update checkpoint manager with loaded data
                checkpoint_manager.state.update(checkpoint_data)
                
                loaded_checkpoint = True
        
        if loaded_checkpoint:
            logger.info(f"Resumed from checkpoint: {results['molecules_imported']} molecules already imported")
    
    # Get ChEMBL client
    client = get_chembl_client()
    
    # Set up worker pool
    task_queue: Queue = Queue()
    result_queue: Queue = Queue()
    workers: List[ChEMBLWorker] = []
    
    # Start workers
    for i in range(MAX_WORKERS):
        worker = ChEMBLWorker(i, task_queue, result_queue)
        worker.start()
        workers.append(worker)
        logger.info(f"Started worker {i}")
    
    # Main processing block with proper try-finally structure
    try:
        # First, try similarity searches based on reference compounds
        for smiles in REFERENCE_SMILES:
            if results["molecules_imported"] >= target_count:
                break
                
            # Try different similarity thresholds
            for similarity in [80, 70, 60]:
                if results["molecules_imported"] >= target_count:
                    break
                    
                try:
                    # Search for compounds similar to the reference
                    similar_compounds = search_compounds_by_similarity(client, smiles, similarity)
                    
                    if not similar_compounds:
                        continue
                    
                    results["similarity_searches"].append({
                        "smiles": smiles,
                        "similarity": similarity,
                        "compounds_found": len(similar_compounds)
                    })
                    
                    # Process compounds in parallel
                    compounds_to_process = []
                    for compound in similar_compounds:
                        # Skip already processed compounds
                        compound_id = compound.get("ChEMBL ID") or compound.get("molecule_chembl_id")
                        if not compound_id or compound_id in results["processed_compounds"]:
                            continue
                            
                        compounds_to_process.append(compound)
                    
                    # Process in batches
                    for i in range(0, len(compounds_to_process), BATCH_SIZE):
                        if results["molecules_imported"] >= target_count:
                            break
                            
                        batch = compounds_to_process[i:i+BATCH_SIZE]
                        
                        # Start batch in checkpoint manager
                        batch_num = results["batches_processed"] + 1
                        checkpoint_manager.start_batch(batch_num, len(batch))
                        
                        # Submit tasks to workers
                        for compound in batch:
                            compound_id = compound.get("ChEMBL ID") or compound.get("molecule_chembl_id")
                            
                            # Skip if already processed
                            if compound_id in results["processed_compounds"]:
                                continue
                                
                            # Create task
                            task = {
                                "compound_data": compound,
                                "compound_id": compound_id,
                                "checkpoint_manager": checkpoint_manager,
                                "max_retries": MAX_RETRIES
                            }
                            
                            # Add to task queue
                            task_queue.put(task)
                        
                        # Process results
                        batch_success = 0
                        batch_errors = 0
                        
                        # Wait for all tasks to complete
                        for _ in range(len(batch)):
                            try:
                                result = result_queue.get(timeout=60)
                                
                                # Update results
                                if result["status"] == "success":
                                    results["molecules_imported"] += 1
                                    results["properties_imported"] += result.get("property_count", 0)
                                    results["processed_compounds"].add(result["compound_id"])
                                    batch_success += 1
                                else:
                                    results["errors"].append({
                                        "compound_id": result["compound_id"],
                                        "error": result.get("error", "Unknown error"),
                                        "error_type": result.get("error_type", "Unknown")
                                    })
                                    results["processed_compounds"].add(result["compound_id"])
                                    batch_errors += 1
                                    
                                # Mark task as done
                                result_queue.task_done()
                                
                            except Empty:
                                logger.warning("Timeout waiting for worker result")
                                batch_errors += 1
                        
                        # End batch in checkpoint manager
                        checkpoint_manager.end_batch(
                            processed=len(batch),
                            success=batch_success,
                            errors=batch_errors
                        )
                        
                        # Update batch count
                        results["batches_processed"] += 1
                        
                        # Create checkpoint
                        create_checkpoint(project_id, {
                            "smiles": smiles,
                            "similarity": similarity,
                            "batch": results["batches_processed"],
                            "molecules_imported": results["molecules_imported"],
                            "properties_imported": results["properties_imported"],
                            "processed_compounds": list(results["processed_compounds"]),
                            "batches_processed": results["batches_processed"],
                            "search_terms_used": results["search_terms_used"],
                            "similarity_searches": results["similarity_searches"]
                        }, checkpoint_manager)
                        
                        # Log progress
                        logger.info(f"Progress: {results['molecules_imported']}/{target_count} molecules imported")
                except Exception as e:
                    logger.error(f"Error in similarity search for {smiles}: {str(e)}")
                    logger.debug(traceback.format_exc())
                    results["errors"].append({
                        "smiles": smiles,
                        "similarity": similarity,
                        "error": str(e)
                    })
    
        # If we still need more compounds, try keyword searches
        if results["molecules_imported"] < target_count:
            # Shuffle the search terms to get a diverse set of compounds
            search_terms = SEARCH_TERMS.copy()
            random.shuffle(search_terms)
            
            for term in search_terms:
                if results["molecules_imported"] >= target_count:
                    break
                    
                try:
                    # Search for compounds matching the term using the ChEMBL API
                    logger.info(f"Searching ChEMBL for term: {term}")
                    search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties__full_mwt__lte=500&limit=1000&q={term}"
                    
                    # Use requests library with retry logic
                    for attempt in range(MAX_RETRIES):
                        try:
                            response = requests.get(search_url, timeout=30)
                            response.raise_for_status()  # Raise exception for 4XX/5XX responses
                            break
                        except (requests.RequestException, requests.Timeout) as req_err:
                            if attempt < MAX_RETRIES - 1:
                                backoff = (2 ** attempt) * (0.5 + random.random())
                                logger.warning(f"API request failed (attempt {attempt + 1}/{MAX_RETRIES}): {str(req_err)}")
                                logger.info(f"Retrying in {backoff:.2f} seconds...")
                                time.sleep(backoff)
                            else:
                                raise
                    
                    # Process the response
                    data = response.json()
                    compounds = data.get("molecules", [])
                    
                    if not compounds:
                        logger.info(f"No compounds found for term: {term}")
                        continue
                        
                    logger.info(f"Found {len(compounds)} compounds for term: {term}")
                    
                    results["search_terms_used"].append({
                        "term": term,
                        "compounds_found": len(compounds)
                    })
                    
                    # Process in batches using the process_chembl_batch function
                    for i in range(0, len(compounds), BATCH_SIZE):
                        if results["molecules_imported"] >= target_count:
                            break
                            
                        batch = compounds[i:i+BATCH_SIZE]
                        logger.info(f"Processing batch {i//BATCH_SIZE + 1} of {(len(compounds) + BATCH_SIZE - 1)//BATCH_SIZE} for term: {term}")
                        
                        # Process the batch using the function from task-005
                        batch_results = process_chembl_batch(project_id, batch)
                        
                        # Update results
                        results["molecules_imported"] += batch_results["molecules_imported"]
                        results["properties_imported"] += batch_results["properties_imported"]
                        results["errors"].extend(batch_results["errors"])
                        results["batches_processed"] += 1
                        
                        # Create checkpoint at intervals
                        if results["batches_processed"] % 5 == 0 or results["molecules_imported"] >= target_count:
                            create_checkpoint(project_id, {
                                "term": term,
                                "batch": i // BATCH_SIZE,
                                "molecules_imported": results["molecules_imported"],
                                "properties_imported": results["properties_imported"],
                                "processed_compounds": list(results["processed_compounds"]),
                                "batches_processed": results["batches_processed"],
                                "search_terms_used": results["search_terms_used"],
                                "similarity_searches": results["similarity_searches"]
                            }, checkpoint_manager)
                        
                        # Log progress
                        logger.info(f"Progress: {results['molecules_imported']}/{target_count} molecules imported")
                        
                except Exception as e:
                    logger.error(f"Error in term search for {term}: {str(e)}")
                    logger.debug(traceback.format_exc())
                    results["errors"].append({
                        "term": term,
                        "error": str(e)
                    })
    
    finally:
        # Clean up workers
        # Signal workers to stop
        for _ in workers:
            task_queue.put(None)
            
        # Wait for workers to finish
        for worker in workers:
            worker.join(timeout=5)
            
        logger.info("All workers stopped")
        
        # Create final checkpoint
        create_checkpoint(project_id, {
            "molecules_imported": results["molecules_imported"],
            "properties_imported": results["properties_imported"],
            "batches_processed": results["batches_processed"],
            "search_terms_used": results["search_terms_used"],
            "similarity_searches": results["similarity_searches"],
            "processed_compounds": list(results["processed_compounds"]),
            "status": "completed",
            "end_time": datetime.now().isoformat()
        }, checkpoint_manager)
    
    return results

def verify_imported_data(project_id: str) -> Dict[str, Any]:
    """
    Verify that the imported data meets the success criteria.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        Verification results
    """
    logger.info("Verifying imported data")
    
    results = {
        "total_molecules": 0,
        "molecules_with_properties": 0,
        "property_completeness_percentage": 0,
        "reference_compounds": 0,
        "reference_compounds_with_properties": 0
    }
    
    # Count total molecules
    total_query = "SELECT COUNT(*) FROM molecules;"
    total_result = execute_sql(project_id, total_query)
    
    if total_result and len(total_result) > 0:
        results["total_molecules"] = total_result[0]["count"]
    
    # Count molecules with properties
    with_props_query = "SELECT COUNT(DISTINCT molecule_id) FROM molecular_properties;"
    with_props_result = execute_sql(project_id, with_props_query)
    
    if with_props_result and len(with_props_result) > 0:
        results["molecules_with_properties"] = with_props_result[0]["count"]
    
    # Calculate property completeness percentage
    if results["total_molecules"] > 0:
        results["property_completeness_percentage"] = (results["molecules_with_properties"] / results["total_molecules"]) * 100
    
    # Count reference compounds
    ref_query = "SELECT COUNT(*) FROM molecules WHERE data_source = 'reference';"
    ref_result = execute_sql(project_id, ref_query)
    
    if ref_result and len(ref_result) > 0:
        results["reference_compounds"] = ref_result[0]["count"]
    
    # Count reference compounds with complete properties
    ref_props_query = """
    SELECT COUNT(DISTINCT m.id) 
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE m.data_source = 'reference'
    GROUP BY m.id
    HAVING COUNT(mp.id) >= 4;
    """
    ref_props_result = execute_sql(project_id, ref_props_query)
    
    if ref_props_result:
        results["reference_compounds_with_properties"] = len(ref_props_result)
    
    # Check if success criteria are met
    results["success_criteria"] = {
        "molecule_count": {
            "target": 5000,
            "actual": results["total_molecules"],
            "met": results["total_molecules"] >= 5000
        },
        "reference_compounds": {
            "target": 9,
            "actual": results["reference_compounds_with_properties"],
            "met": results["reference_compounds_with_properties"] >= 9
        },
        "property_completeness": {
            "target": 90,
            "actual": results["property_completeness_percentage"],
            "met": results["property_completeness_percentage"] >= 90
        }
    }
    
    return results

def update_project_state(results: Dict[str, Any], verification_results: Dict[str, Any]) -> None:
    """
    Update the project state file with the results of the operation.
    
    Args:
        results: The results of the operation
        verification_results: The verification results
    """
    try:
        # Read the current project state
        with open("project_state.json", "r") as f:
            project_state = json.load(f)
        
        # Update the high-level plan
        for phase in project_state["highLevelPlan"]:
            if phase["phase"] == "ChEMBL Data Import":
                phase["status"] = "Completed"
                phase["findings"] = f"Imported {results['molecules_imported']} molecules with {results['properties_imported']} properties. Property completeness: {verification_results['property_completeness_percentage']:.1f}%."
        
        # Update the current phase
        project_state["currentPhase"] = "Performance Optimization"
        
        # Update the database stats
        project_state["databaseStats"]["totalMolecules"] = verification_results["total_molecules"]
        project_state["databaseStats"]["moleculesWithProperties"] = verification_results["molecules_with_properties"]
        project_state["databaseStats"]["propertyCompleteness"] = verification_results["property_completeness_percentage"]
        
        # Update the task status
        task_id = "task-004"
        if task_id in project_state["tasks"]:
            project_state["tasks"][task_id]["status"] = "Completed"
            project_state["tasks"][task_id]["log"].append({
                "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z"),
                "message": f"Task completed - Imported {results['molecules_imported']} molecules with {results['properties_imported']} properties"
            })
        
        # Add a new task for Performance Optimization
        project_state["tasks"]["task-005"] = {
            "description": "Performance Optimization",
            "assignedTo": "master-orchestrator",
            "status": "Pending",
            "dependsOn": ["task-004"],
            "outputs": [],
            "log": []
        }
        
        # Add a journal entry
        project_state["journal"].append({
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z"),
            "entry": f"Phase 3 (ChEMBL Data Import) completed. Imported {results['molecules_imported']} molecules with {results['properties_imported']} properties. Property completeness: {verification_results['property_completeness_percentage']:.1f}%. Moving to Phase 4 (Performance Optimization)."
        })
        
        # Write the updated project state
        with open("project_state.json", "w") as f:
            json.dump(project_state, f, indent=2)
        
        logger.info("Updated project state")
    except Exception as e:
        logger.error(f"Error updating project state: {str(e)}")

def main():
    """Main function."""
    try:
        # Ensure logs directory exists
        os.makedirs("logs", exist_ok=True)
        
        logger.info("Starting ChEMBL data import")
        
        # Get the Supabase project ID
        project_id = "tsdlmynydfuypiugmkev"
        
        # Import ChEMBL data
        results = import_chembl_data(project_id)
        
        logger.info(f"Imported {results['molecules_imported']} molecules with {results['properties_imported']} properties")
        
        if results["errors"]:
            logger.warning(f"Encountered {len(results['errors'])} errors")
            for error in results["errors"][:10]:  # Log only the first 10 errors
                logger.warning(f"Error: {error}")
        
        # Verify imported data
        verification_results = verify_imported_data(project_id)
        
        logger.info(f"Total molecules: {verification_results['total_molecules']}")
        logger.info(f"Molecules with properties: {verification_results['molecules_with_properties']}")
        logger.info(f"Property completeness: {verification_results['property_completeness_percentage']:.1f}%")
        
        # Check if success criteria are met
        success_criteria = verification_results["success_criteria"]
        for criterion, details in success_criteria.items():
            status = "OK - Met" if details["met"] else "FAIL - Not met"
            logger.info(f"{criterion}: {details['actual']} / {details['target']} - {status}")
        
        # Update project state
        update_project_state(results, verification_results)
        
        logger.info("ChEMBL data import completed successfully")
        
        return 0
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1


def process_chembl_batch(project_id: str, batch: List[Dict[str, Any]], conn=None) -> Dict[str, Any]:
    """
    Process a batch of ChEMBL molecules and import into database.
    
    Args:
        project_id: The Supabase project ID
        batch: List of molecule data from ChEMBL API
        
    Returns:
        dict: Statistics about the batch processing
    """
    results = {
        "molecules_imported": 0,
        "properties_imported": 0,
        "errors": []
    }
    
    for molecule_data in batch:
        try:
            # Skip if missing required fields
            if not all(field in molecule_data for field in ["molecule_chembl_id", "molecule_structures"]):
                continue
                
            if not molecule_data.get("molecule_structures", {}).get("canonical_smiles"):
                continue
                
            # Extract basic properties
            molecule = {
                "chembl_id": molecule_data["molecule_chembl_id"],
                "name": molecule_data.get("pref_name", molecule_data["molecule_chembl_id"]),
                "smiles": molecule_data["molecule_structures"]["canonical_smiles"],
                "inchi": molecule_data["molecule_structures"].get("standard_inchi", ""),
                "inchikey": molecule_data["molecule_structures"].get("standard_inchi_key", ""),
                "formula": molecule_data.get("molecule_properties", {}).get("full_molformula", ""),
                "molecular_weight": molecule_data.get("molecule_properties", {}).get("full_mwt"),
                "data_source": "ChEMBL"
            }
            
            # Insert molecule
            molecule_sql = """
            INSERT INTO molecules
                (chembl_id, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
            VALUES
                ($1, $2, $3, $4, $5, $6, $7, $8)
            ON CONFLICT (chembl_id) DO UPDATE SET
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi,
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                molecular_weight = EXCLUDED.molecular_weight,
                updated_at = NOW()
            RETURNING id;
            """
            
            molecule_result = execute_sql_with_retry(
                project_id=project_id,
                query=molecule_sql,
                params=[
                    molecule["chembl_id"],
                    molecule["name"],
                    molecule["smiles"],
                    molecule["inchi"],
                    molecule["inchikey"],
                    molecule["formula"],
                    molecule["molecular_weight"],
                    molecule["data_source"]
                ]
            )
            
            if not molecule_result or len(molecule_result) == 0:
                continue
                
            molecule_id = molecule_result[0]["id"]
            results["molecules_imported"] += 1
            
            # Extract and insert properties
            properties = {}
            mol_props = molecule_data.get("molecule_properties", {})
            
            # Map ChEMBL properties to our property names
            if "alogp" in mol_props:
                properties["LogP"] = mol_props["alogp"]
            if "hba" in mol_props:
                properties["Hydrogen Bond Acceptor Count"] = mol_props["hba"]
            if "hbd" in mol_props:
                properties["Hydrogen Bond Donor Count"] = mol_props["hbd"]
            if "rtb" in mol_props:
                properties["Rotatable Bond Count"] = mol_props["rtb"]
            if "psa" in mol_props:
                properties["Polar Surface Area"] = mol_props["psa"]
                
            # Add properties
            for prop_name, value in properties.items():
                try:
                    # Get or create property type
                    prop_type_sql = """
                    INSERT INTO property_types (name, data_type)
                    VALUES ($1, 'numeric')
                    ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
                    RETURNING id;
                    """
                    
                    prop_type_result = execute_sql_with_retry(
                        project_id=project_id,
                        query=prop_type_sql,
                        params=[prop_name]
                    )
                    
                    if prop_type_result and len(prop_type_result) > 0:
                        property_type_id = prop_type_result[0]["id"]
                        
                        # Insert property value
                        property_sql = """
                        INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value, source)
                        VALUES ($1, $2, $3, $4)
                        ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                            numeric_value = EXCLUDED.numeric_value,
                            updated_at = NOW();
                        """
                        
                        execute_sql_with_retry(
                            project_id=project_id,
                            query=property_sql,
                            params=[molecule_id, property_type_id, value, "ChEMBL"]
                        )
                        
                        results["properties_imported"] += 1
                        
                except Exception as e:
                    results["errors"].append({
                        "molecule_id": molecule_id,
                        "property": prop_name,
                        "error": str(e)
                    })
                    
        except Exception as e:
            results["errors"].append({
                "molecule": molecule_data.get("molecule_chembl_id", "unknown"),
                "error": str(e)
            })
    
    return results

if __name__ == "__main__":
    sys.exit(main())
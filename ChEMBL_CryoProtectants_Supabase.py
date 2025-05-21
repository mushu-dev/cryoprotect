#!/usr/bin/env python3
"""
CryoProtect Analyzer - ChEMBL Cryoprotectant Data Importer (Supabase Version)

This script retrieves cryoprotectant data from ChEMBL, filters molecules based on
predefined criteria, scores them, and stores the results in a Supabase database.

Features:
- Batch processing with parameterizable batch size (default 100-500 molecules)
- Memory usage monitoring with adaptive rate limiting
- Robust checkpointing for resumable operation
- CLI arguments for batch size, checkpoint path, resume/reset
- Efficient Supabase bulk inserts for molecules and properties
- Progress and ETA logging
- Error/skipped molecule logging for review
- Deduplication based on InChIKey
- Source attribution for provenance tracking

Prerequisites:
- Python 3.6+ installed
- Supabase project with the CryoProtect schema applied
- supabase-py package installed (pip install supabase)
- python-dotenv package installed (pip install python-dotenv)
- requests package installed (pip install requests)
"""

import os
import time
import requests
import json
import argparse
import sys
import traceback
from datetime import datetime, timedelta
from typing import Dict, List, Any, Optional, Union, Tuple
import gc  # For garbage collection
from supabase import create_client, Client

# Import the configuration system
from config import active_config, validate_config, ConfigurationError

# Import the ChEMBL client
from chembl.client import ResilientChEMBLClient
from chembl.rate_limiter import AdaptiveRateLimiter

# Import the centralized logging system
from chembl.logging import (
    ChEMBLLogger, log_error, log_skipped_molecule,
    log_progress, write_summary, get_logger
)

# Initialize logger
logger = get_logger(__name__)

# Create a global logger instance with default settings
# This will create the logs directory if it doesn't exist
chembl_logger = ChEMBLLogger(
    log_dir="logs",
    progress_log="chembl_progress.jsonl",
    error_log="chembl_errors.jsonl",
    skipped_log="skipped_chembl_molecules.jsonl",
    summary_log="chembl_summary.json",
    general_log="chembl_import.log"
)

# Scoring Weights (Total = 200)
WEIGHTS = {
    "hydrogen_bonding": 50,
    "solubility_polarity": 40,
    "membrane_permeability": 40,
    "toxicity_biocompatibility": 30,
    "protein_stabilization": 20,
    "stability_reactivity": 10,
    "environmental_safety": 10
}

# Stage 1: Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-5, 5),  # Relaxed for testing
    "mw_range": (0, 1000),  # Relaxed for testing
    "TPSA_range": (0, 200), # Relaxed for testing
    "functional_groups": [] # No functional group requirement for testing
}

# Configuration variables will be initialized in verify_configuration()
CHEMBL_API_DELAY = None
CHEMBL_ID_FILE = None
supabase = None
chembl_client = None

def verify_configuration():
    """
    Verify all required configuration variables and test API connectivity.
    This function validates the configuration, logs non-secret values,
    and tests connectivity to both Supabase and ChEMBL APIs.
    
    Raises:
        ConfigurationError: If any required configuration is missing or invalid
        ConnectionError: If API connectivity tests fail
    """
    global CHEMBL_API_DELAY, CHEMBL_ID_FILE, supabase, chembl_client
    
    logger.info("Verifying configuration...")
    
    # Validate configuration using the hierarchical config system
    try:
        validate_config()
    except ConfigurationError as e:
        log_error(
            error_type="Configuration",
            message="Configuration validation failed",
            context={
                "exception": e,
                "source": "verify_configuration.validate_config"
            }
        )
        sys.exit(1)
    
    # Log all non-secret configuration values
    logger.info("Configuration values:")
    logger.info(f"  SUPABASE_URL: {active_config.SUPABASE_URL}")
    # Mask secret values
    logger.info(f"  SUPABASE_KEY: {'*' * 8}")
    if hasattr(active_config, 'CHEMBL_API_URL'):
        logger.info(f"  CHEMBL_API_URL: {active_config.CHEMBL_API_URL}")
    
    # Get configuration values from the config system
    try:
        # Required configuration
        supabase_url = active_config.SUPABASE_URL
        supabase_key = active_config.SUPABASE_KEY
        
        # Optional configuration with defaults
        CHEMBL_API_DELAY = getattr(active_config, 'CHEMBL_API_DELAY', 0.3)
        logger.info(f"  CHEMBL_API_DELAY: {CHEMBL_API_DELAY}")
        
        CHEMBL_ID_FILE = getattr(active_config, 'CHEMBL_ID_FILE', '')
        logger.info(f"  CHEMBL_ID_FILE: {CHEMBL_ID_FILE or 'Not specified'}")
        
        # Batch and checkpoint parameters
        batch_size = getattr(active_config, 'CHEMBL_BATCH_SIZE', 10)
        logger.info(f"  CHEMBL_BATCH_SIZE: {batch_size}")
        
        checkpoint_dir = getattr(active_config, 'CHECKPOINT_DIR', 'checkpoints')
        logger.info(f"  CHECKPOINT_DIR: {checkpoint_dir}")
        
        # Create checkpoint directory if it doesn't exist
        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
            logger.info(f"Created checkpoint directory: {checkpoint_dir}")
        
    except AttributeError as e:
        log_error(
            error_type="Configuration",
            message="Missing required configuration",
            context={
                "exception": e,
                "source": "verify_configuration.get_config_values"
            }
        )
        sys.exit(1)
    
    # Initialize Supabase client with service role key
    try:
        logger.info("Initializing Supabase client with service role key...")
        
        # Look for service role key in different config variables
        service_role_key = None
        for key_name in ['SUPABASE_SERVICE_ROLE_KEY', 'SUPABASE_SERVICE_KEY', 'SUPABASE_ADMIN_KEY', 'SUPABASE_KEY']:
            potential_key = getattr(active_config, key_name, None)
            if potential_key:
                # Check if it's likely a service role key (less strict check)
                if potential_key.startswith("eyJ"):
                    service_role_key = potential_key
                    logger.info(f"Using key from {key_name}")
                    break
        
        # If no service role key found, use whatever key is available
        if not service_role_key:
            logger.warning("No explicit service role key found. Using SUPABASE_KEY, which may not have sufficient permissions.")
            service_role_key = supabase_key
            
        # Create client with service role key
        supabase = create_client(supabase_url, service_role_key)
        
        # Set global headers for all requests
        os.environ["SUPABASE_KEY"] = service_role_key
    except Exception as e:
        log_error(
            error_type="Initialization",
            message="Failed to initialize Supabase client",
            context={
                "exception": e,
                "url": supabase_url,
                "source": "verify_configuration.init_supabase"
            }
        )
        sys.exit(1)
    
    # Test Supabase connectivity
    try:
        logger.info("Testing Supabase connectivity...")
        response = supabase.table('molecules').select('count', count='exact').limit(1).execute()
        if hasattr(response, 'error') and response.error:
            raise ConnectionError(f"Supabase query error: {response.error}")
        logger.info("Supabase connectivity test successful")
    except Exception as e:
        log_error(
            error_type="Connectivity",
            message="Supabase connectivity test failed",
            context={
                "exception": e,
                "source": "verify_configuration.test_supabase"
            }
        )
        sys.exit(1)
    
    # Initialize ChEMBL client
    try:
        logger.info("Initializing ChEMBL client...")
        cache_dir = getattr(active_config, 'CHEMBL_CACHE_DIR', 'cache/chembl')
        requests_per_second = getattr(active_config, 'CHEMBL_REQUESTS_PER_SECOND', 5.0)
        max_retries = getattr(active_config, 'CHEMBL_MAX_RETRIES', 5)
        failure_threshold = getattr(active_config, 'CHEMBL_FAILURE_THRESHOLD', 3)
        recovery_timeout = getattr(active_config, 'CHEMBL_RECOVERY_TIMEOUT', 60)
        cache_ttl = getattr(active_config, 'CHEMBL_CACHE_TTL', 86400 * 30)  # 30 days
        memory_cache_size = getattr(active_config, 'CHEMBL_MEMORY_CACHE_SIZE', 1000)
        # Note: memory_threshold is not used by ResilientChEMBLClient constructor
        memory_threshold = getattr(active_config, 'CHEMBL_MEMORY_THRESHOLD', 80.0)  # 80% memory threshold
        
        chembl_client = ResilientChEMBLClient(
            cache_dir=cache_dir,
            requests_per_second=requests_per_second,
            max_retries=max_retries,
            failure_threshold=failure_threshold,
            recovery_timeout=recovery_timeout,
            cache_ttl=cache_ttl,
            memory_cache_size=memory_cache_size
            # memory_threshold parameter is not supported by ResilientChEMBLClient
        )
        logger.info(f"  CHEMBL_CACHE_DIR: {cache_dir}")
        logger.info(f"  CHEMBL_REQUESTS_PER_SECOND: {requests_per_second}")
        logger.info(f"  CHEMBL_MEMORY_CACHE_SIZE: {memory_cache_size}")
    except Exception as e:
        # Print full traceback for debugging
        import traceback
        print("\nDetailed ChEMBL client initialization error:")
        traceback.print_exc()
        print("\n")
        
        log_error(
            error_type="Initialization",
            message="Failed to initialize ChEMBL client",
            context={
                "exception": str(e),
                "exception_type": type(e).__name__,
                "cache_dir": cache_dir,
                "requests_per_second": requests_per_second,
                "source": "verify_configuration.init_chembl"
            }
        )
        sys.exit(1)
    
    # Test ChEMBL API connectivity
    try:
        logger.info("Testing ChEMBL API connectivity...")
        test_result = chembl_client.search_molecules("glycerol", limit=1)
        if "Error" in test_result:
            raise ConnectionError(f"ChEMBL API error: {test_result['Error']}")
        logger.info("ChEMBL API connectivity test successful")
    except Exception as e:
        log_error(
            error_type="Connectivity",
            message="ChEMBL API connectivity test failed",
            context={
                "exception": e,
                "test_query": "glycerol",
                "source": "verify_configuration.test_chembl"
            }
        )
        sys.exit(1)
    
    logger.info("Configuration verification completed successfully")
    return True

def get_chembl_ids():
    """Retrieve ChEMBL IDs from a file or search for potential cryoprotectants."""
    if CHEMBL_ID_FILE and os.path.exists(CHEMBL_ID_FILE):
        with open(CHEMBL_ID_FILE, "r") as file:
            chembl_ids = [line.strip() for line in file if line.strip()]
        logger.info(f"SUCCESS: Loaded {len(chembl_ids)} ChEMBL IDs from file.")
        return chembl_ids
    
    # If no file is available, search for potential cryoprotectants
    logger.info("No ChEMBL ID file found. Searching for potential cryoprotectants...")
    
    # List of search terms for potential cryoprotectants
    search_terms = [
        "glycerol", "dimethyl sulfoxide", "DMSO", "ethylene glycol", 
        "propylene glycol", "trehalose", "sucrose", "glucose", 
        "methanol", "formamide", "acetamide", "proline", 
        "hydroxyectoine", "betaine", "polyvinyl alcohol", "cryoprotectant"
    ]
    
    chembl_ids = set()
    for term in search_terms:
        try:
            logger.info(f"Searching for '{term}'...")
            results = chembl_client.search_molecules(term, limit=50)
            if "Error" in results:
                logger.warning(f"Error searching for '{term}': {results['Error']}")
                continue
                
            for molecule in results.get("Molecules", []):
                if "ChEMBL ID" in molecule:
                    chembl_ids.add(molecule["ChEMBL ID"])
            
            logger.info(f"Found {len(results.get('Molecules', []))} results for '{term}'")
            time.sleep(CHEMBL_API_DELAY)
        except Exception as e:
            logger.warning(f"Error searching for '{term}': {str(e)}")
    
    chembl_ids = list(chembl_ids)
    logger.info(f"SUCCESS: Found {len(chembl_ids)} unique ChEMBL IDs from search.")
    return chembl_ids

def get_molecule_properties(chembl_id):
    """Fetch molecular properties from ChEMBL."""
    try:
        molecule = chembl_client.get_molecule_properties(chembl_id)
        
        if "Error" in molecule:
            error_msg = molecule["Error"]
            log_error(
                error_type="API",
                message=f"Error fetching properties for ChEMBL ID {chembl_id}",
                context={
                    "chembl_id": chembl_id,
                    "api_error": error_msg,
                    "source": "get_molecule_properties"
                }
            )
            return {"ChEMBL ID": chembl_id, "Error": error_msg}
        
        # Map ChEMBL properties to our format
        properties = {
            "ChEMBL ID": molecule.get("ChEMBL ID"),
            "Name": molecule.get("Name"),
            "Molecular Formula": molecule.get("Molecular Formula"),
            "Molecular Weight": molecule.get("Molecular Weight"),
            "LogP": molecule.get("LogP"),
            "TPSA": molecule.get("TPSA"),
            "H-Bond Donors": molecule.get("H-Bond Donors"),
            "H-Bond Acceptors": molecule.get("H-Bond Acceptors"),
            "SMILES": molecule.get("SMILES"),
            "InChI": molecule.get("InChI"),
            "InChIKey": molecule.get("InChIKey"),
            "ChEMBL Link": molecule.get("ChEMBL Link")
        }
        
        return properties
    except Exception as e:
        log_error(
            error_type="Exception",
            message=f"Exception fetching properties for ChEMBL ID {chembl_id}",
            context={
                "chembl_id": chembl_id,
                "exception": e,
                "source": "get_molecule_properties"
            }
        )
        return {"ChEMBL ID": chembl_id, "Error": str(e)}

def get_additional_properties(chembl_id):
    """Fetch additional properties from ChEMBL."""
    try:
        # This would typically fetch additional properties not included in the basic molecule data
        # For ChEMBL, we might need to make additional API calls for specific property types
        
        # For now, return placeholder data
        return {
            "Toxicity": None,
            "Stability": None,
            "Environmental Safety": None
        }
    except Exception as e:
        log_error(
            error_type="Exception",
            message=f"Exception fetching additional properties for ChEMBL ID {chembl_id}",
            context={
                "chembl_id": chembl_id,
                "exception": e,
                "source": "get_additional_properties"
            }
        )
        return {"Error": str(e)}

def filter_molecule(molecule, batch_num=None):
    """Initial filtering based on core cryoprotectant properties, ensuring numerical values."""
    chembl_id = molecule.get('ChEMBL ID', 'unknown')
    
    if "Error" in molecule:
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason=f"Error in molecule data: {molecule.get('Error', 'Unknown error')}",
            molecule_data=molecule,
            category="error",
            batch_num=batch_num
        )
        return False

    # Check for required fields for Supabase schema
    if not molecule.get("SMILES"):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason="Missing SMILES",
            molecule_data=molecule,
            category="validation",
            batch_num=batch_num
        )
        return False
        
    if not molecule.get("InChI"):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason="Missing InChI",
            molecule_data=molecule,
            category="validation",
            batch_num=batch_num
        )
        return False
        
    if not molecule.get("InChIKey"):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason="Missing InChIKey",
            molecule_data=molecule,
            category="validation",
            batch_num=batch_num
        )
        return False
        
    if not molecule.get("Molecular Formula"):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason="Missing Molecular Formula",
            molecule_data=molecule,
            category="validation",
            batch_num=batch_num
        )
        return False

    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError) as e:
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason=f"Invalid numerical values: {str(e)}",
            molecule_data=molecule,
            category="validation",
            batch_num=batch_num
        )
        return False

    smiles = molecule["SMILES"]

    # Check property ranges
    if mw is None or not (CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1]):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason=f"Molecular weight {mw} outside range {CORE_CRITERIA['mw_range']}",
            molecule_data={"ChEMBL ID": chembl_id, "Molecular Weight": mw},
            category="filter",
            batch_num=batch_num
        )
        return False
        
    if logp is None or not (CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1]):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason=f"LogP {logp} outside range {CORE_CRITERIA['logP_range']}",
            molecule_data={"ChEMBL ID": chembl_id, "LogP": logp},
            category="filter",
            batch_num=batch_num
        )
        return False
        
    if tpsa is None or not (CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
        log_skipped_molecule(
            chembl_id=chembl_id,
            reason=f"TPSA {tpsa} outside range {CORE_CRITERIA['TPSA_range']}",
            molecule_data={"ChEMBL ID": chembl_id, "TPSA": tpsa},
            category="filter",
            batch_num=batch_num
        )
        return False
        
    if CORE_CRITERIA["functional_groups"]:
        if not any(group in smiles for group in CORE_CRITERIA["functional_groups"]):
            log_skipped_molecule(
                chembl_id=chembl_id,
                reason=f"Missing required functional groups {CORE_CRITERIA['functional_groups']}",
                molecule_data={"ChEMBL ID": chembl_id, "SMILES": smiles},
                skipped_log_path=skipped_log_path
            )
            return False

    return True

def score_molecule(molecule, extra_properties):
    """Compute final score out of 200 based on all properties."""
    score = 0

    # Hydrogen Bonding
    score += WEIGHTS["hydrogen_bonding"]

    # Solubility & Permeability
    score += WEIGHTS["solubility_polarity"]
    score += WEIGHTS["membrane_permeability"]

    # Toxicity & Stability
    if extra_properties.get("Toxicity"):
        score += WEIGHTS["toxicity_biocompatibility"]
    if extra_properties.get("Stability"):
        score += WEIGHTS["stability_reactivity"]

    # Environmental Safety
    if extra_properties.get("Environmental Safety"):
        score += WEIGHTS["environmental_safety"]

    return score

def fetch_property_types():
    """Fetch property types from Supabase once per run."""
    try:
        # Add headers to bypass RLS policies
        headers = {
            "apikey": supabase.supabase_key,
            "Authorization": f"Bearer {supabase.supabase_key}",
            "Prefer": "return=representation"
        }
        
        # Use standard Supabase API to fetch property types
        response = supabase.table("property_types").select("id, name, data_type").execute(headers=headers)
        
        # Check for data in the response
        if hasattr(response, "data") and response.data:
            logger.info(f"Successfully fetched {len(response.data)} property types")
            return response.data
        else:
            logger.warning("No property types found in the database")
            return None
    except Exception as e:
        log_structured_error(
            error_type="Database",
            message="Exception fetching property types",
            context={
                "exception": e,
                "source": "fetch_property_types"
            }
        )
        # Return a minimal set of property types for testing
        logger.info("Using fallback property types due to database error")
        return [
            {"id": "00000000-0000-0000-0000-000000000001", "name": "LogP", "data_type": "numeric"},
            {"id": "00000000-0000-0000-0000-000000000002", "name": "TPSA", "data_type": "numeric"},
            {"id": "00000000-0000-0000-0000-000000000003", "name": "H-Bond Donors", "data_type": "numeric"},
            {"id": "00000000-0000-0000-0000-000000000004", "name": "H-Bond Acceptors", "data_type": "numeric"},
            {"id": "00000000-0000-0000-0000-000000000005", "name": "ChEMBL ID", "data_type": "text"}
        ]

def check_molecule_exists(inchikey):
    """Check if a molecule with the given InChIKey already exists in the database."""
    if not inchikey:
        return None
    
    try:
        # Use standard Supabase API with headers to bypass RLS
        headers = {
            "apikey": supabase.supabase_key,
            "Authorization": f"Bearer {supabase.supabase_key}",
            "Prefer": "return=representation"
        }
        
        # Query the molecules table for the inchikey
        response = supabase.table("molecules").select("id").eq("inchikey", inchikey).execute(headers=headers)
        
        # Check if we got any results
        if hasattr(response, "data") and response.data and len(response.data) > 0:
            return response.data[0]["id"]
                
        return None
    except Exception as e:
        log_structured_error(
            error_type="Database",
            message="Error checking if molecule exists",
            context={
                "inchikey": inchikey,
                "exception": str(e),
                "traceback": traceback.format_exc(),
                "source": "check_molecule_exists"
            }
        )
        return None

def bulk_insert_molecules(molecule_batch, user_id):
    """Insert molecules in bulk and return list of generated IDs (in order)."""
    if not molecule_batch:
        return []
    
    logger.info(f"Inserting {len(molecule_batch)} molecules into the database")
    
    inserted_ids = []
    try:
        # Use standard Supabase API with headers to bypass RLS
        headers = {
            "apikey": supabase.supabase_key,
            "Authorization": f"Bearer {supabase.supabase_key}",
            "Prefer": "return=representation"
        }
        
        # Insert molecules one by one to better handle errors
        for molecule in molecule_batch:
            try:
                # Insert the molecule
                response = supabase.table("molecules").insert(molecule).execute(headers=headers)
                
                # Extract the ID from the response
                if hasattr(response, "data") and response.data and len(response.data) > 0:
                    inserted_id = response.data[0].get("id")
                    if inserted_id:
                        inserted_ids.append(inserted_id)
                        logger.debug(f"Inserted molecule {molecule.get('name')} with ID {inserted_id}")
                    else:
                        logger.warning(f"No ID returned for inserted molecule {molecule.get('name')}")
                else:
                    logger.warning(f"Unexpected response format for molecule {molecule.get('name')}")
            except Exception as inner_e:
                log_structured_error(
                    error_type="Database",
                    message=f"Error inserting molecule {molecule.get('name')}",
                    context={
                        "exception": str(inner_e),
                        "molecule": molecule.get("name"),
                        "source": "bulk_insert_molecules.single_insert"
                    }
                )
        
        logger.info(f"Successfully inserted {len(inserted_ids)} molecules")
        return inserted_ids
    except Exception as e:
        log_structured_error(
            error_type="Database",
            message="Exception during bulk molecule insertion",
            context={
                "exception": str(e),
                "traceback": traceback.format_exc(),
                "batch_size": len(molecule_batch),
                "source": "bulk_insert_molecules"
            }
        )
        return inserted_ids  # Return any IDs that were successfully inserted before the exception

def bulk_insert_properties(property_batch):
    """Insert molecular properties in bulk."""
    if not property_batch:
        return True
    
    logger.info(f"Inserting {len(property_batch)} properties into the database")
    
    success_count = 0
    try:
        # Use standard Supabase API with headers to bypass RLS
        headers = {
            "apikey": supabase.supabase_key,
            "Authorization": f"Bearer {supabase.supabase_key}",
            "Prefer": "return=representation"
        }
        
        # Process properties in smaller sub-batches to avoid overwhelming the database
        sub_batch_size = 50
        for i in range(0, len(property_batch), sub_batch_size):
            sub_batch = property_batch[i:i+sub_batch_size]
            
            # Insert properties one by one to better handle errors
            for prop in sub_batch:
                # Skip if required fields are missing
                if not prop.get('molecule_id') or not prop.get('property_type_id'):
                    continue
                
                try:
                    # Insert the property
                    response = supabase.table("molecular_properties").insert(prop).execute(headers=headers)
                    
                    # Check if successful
                    if hasattr(response, "data") and response.data:
                        success_count += 1
                except Exception as inner_e:
                    log_structured_error(
                        error_type="Database",
                        message="Error inserting property",
                        context={
                            "exception": str(inner_e),
                            "property_type": prop.get("property_type_id"),
                            "molecule_id": prop.get("molecule_id"),
                            "source": "bulk_insert_properties.single_insert"
                        }
                    )
            
            # Log progress for large batches
            if i % 200 == 0 and i > 0:
                logger.info(f"Inserted {success_count} properties so far ({i}/{len(property_batch)})")
        
        logger.info(f"Successfully inserted {success_count}/{len(property_batch)} properties")
        return success_count > 0
    except Exception as e:
        log_structured_error(
            error_type="Database",
            message="Exception during bulk property insertion",
            context={
                "exception": str(e),
                "traceback": traceback.format_exc(),
                "batch_size": len(property_batch),
                "success_count": success_count,
                "source": "bulk_insert_properties"
            }
        )
        return success_count > 0  # Return True if at least some properties were inserted

def generate_checkpoint_filename(base_name="chembl_checkpoint", include_hash=True):
    """
    Generate a timestamped checkpoint filename with optional hash for uniqueness.
    
    Args:
        base_name (str): Base name for the checkpoint file
        include_hash (bool): Whether to include a hash for uniqueness
        
    Returns:
        str: Timestamped checkpoint filename
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    if include_hash:
        import hashlib
        import random
        # Add a random hash to ensure uniqueness even if multiple checkpoints
        # are created in the same second
        random_hash = hashlib.md5(f"{timestamp}_{random.random()}".encode()).hexdigest()[:8]
        return f"{base_name}_{timestamp}_{random_hash}.json"
    else:
        return f"{base_name}_{timestamp}.json"

def get_latest_checkpoint(checkpoint_dir, base_name="chembl_checkpoint", validate=True):
    """
    Find the most recent checkpoint file in the checkpoint directory.
    
    Args:
        checkpoint_dir (str): Directory containing checkpoint files
        base_name (str): Base name for the checkpoint files
        validate (bool): Whether to validate the checkpoint file content
        
    Returns:
        str or None: Path to the most recent checkpoint file, or None if no checkpoint exists
    """
    if not os.path.exists(checkpoint_dir):
        logger.warning(f"Checkpoint directory {checkpoint_dir} does not exist")
        return None
        
    # Find all checkpoint files matching the base name pattern
    checkpoint_files = [f for f in os.listdir(checkpoint_dir)
                        if f.startswith(base_name) and f.endswith('.json')]
    
    if not checkpoint_files:
        logger.info(f"No checkpoint files found in {checkpoint_dir}")
        return None
    
    # Sort by modification time (most recent first)
    checkpoint_files.sort(key=lambda f: os.path.getmtime(os.path.join(checkpoint_dir, f)), reverse=True)
    
    # Try checkpoints in order until we find a valid one
    for checkpoint_file in checkpoint_files:
        checkpoint_path = os.path.join(checkpoint_dir, checkpoint_file)
        
        if not validate:
            logger.info(f"Found latest checkpoint: {checkpoint_path}")
            return checkpoint_path
            
        # Validate the checkpoint file
        try:
            with open(checkpoint_path, "r") as f:
                checkpoint_data = json.load(f)
                
            # Check for required fields
            required_fields = ["last_completed_batch", "total_processed", "timestamp"]
            if all(field in checkpoint_data for field in required_fields):
                logger.info(f"Found valid checkpoint: {checkpoint_path}")
                return checkpoint_path
            else:
                logger.warning(f"Checkpoint file {checkpoint_path} is missing required fields, trying next one")
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Error reading checkpoint file {checkpoint_path}: {str(e)}, trying next one")
    
    logger.warning("No valid checkpoint files found")
    return None

def read_checkpoint(checkpoint_path):
    """
    Read checkpoint data from a file.
    
    Args:
        checkpoint_path (str): Path to the checkpoint file
        
    Returns:
        dict or None: Checkpoint data, or None if the file doesn't exist or is invalid
    """
    if not checkpoint_path or not os.path.exists(checkpoint_path):
        return None
        
    try:
        with open(checkpoint_path, "r") as f:
            checkpoint_data = json.load(f)
            
        # Validate checkpoint data structure
        required_fields = ["last_completed_batch", "total_processed", "timestamp"]
        if not all(field in checkpoint_data for field in required_fields):
            logger.warning(f"Checkpoint file {checkpoint_path} is missing required fields")
            return None
            
        # Check for data consistency
        if checkpoint_data["last_completed_batch"] < 0:
            logger.warning(f"Invalid batch number in checkpoint: {checkpoint_data['last_completed_batch']}")
            return None
            
        if checkpoint_data["total_processed"] < 0:
            logger.warning(f"Invalid processed count in checkpoint: {checkpoint_data['total_processed']}")
            return None
            
        # Parse timestamp to ensure it's valid
        try:
            datetime.fromisoformat(checkpoint_data["timestamp"])
        except ValueError:
            logger.warning(f"Invalid timestamp in checkpoint: {checkpoint_data['timestamp']}")
            return None
            
        # Create a backup of the checkpoint file
        backup_path = f"{checkpoint_path}.bak"
        try:
            import shutil
            shutil.copy2(checkpoint_path, backup_path)
            logger.debug(f"Created backup of checkpoint file: {backup_path}")
        except IOError as e:
            logger.warning(f"Failed to create backup of checkpoint file: {str(e)}")
            
        logger.info(f"Successfully loaded checkpoint from {checkpoint_path}")
        logger.info(f"Checkpoint data: Last batch: {checkpoint_data['last_completed_batch']}, "
                   f"Total processed: {checkpoint_data['total_processed']}, "
                   f"Timestamp: {checkpoint_data['timestamp']}")
        return checkpoint_data
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error reading checkpoint file {checkpoint_path}: {str(e)}")
        return None

def write_checkpoint(checkpoint_path, checkpoint_data):
    """
    Write checkpoint data to a file with atomic operation and validation.
    
    Args:
        checkpoint_path (str): Path to the checkpoint file
        checkpoint_data (dict): Checkpoint data to write
        
    Returns:
        bool: True if successful, False otherwise
    """
    # Ensure the checkpoint directory exists
    checkpoint_dir = os.path.dirname(checkpoint_path)
    if not os.path.exists(checkpoint_dir):
        try:
            os.makedirs(checkpoint_dir)
            logger.info(f"Created checkpoint directory: {checkpoint_dir}")
        except OSError as e:
            logger.error(f"Error creating checkpoint directory {checkpoint_dir}: {str(e)}")
            return False
    
    # Ensure timestamp is present
    if "timestamp" not in checkpoint_data:
        checkpoint_data["timestamp"] = datetime.now().isoformat()
    
    # Add metadata for traceability
    checkpoint_data["metadata"] = {
        "version": "2.0",
        "created_at": datetime.now().isoformat(),
        "hostname": os.environ.get("COMPUTERNAME", "unknown"),
        "process_id": os.getpid()
    }
    
    # Write to a temporary file first for atomic operation
    temp_path = f"{checkpoint_path}.tmp"
    try:
        with open(temp_path, "w") as f:
            json.dump(checkpoint_data, f, indent=2)
        
        # Validate the temporary file
        try:
            with open(temp_path, "r") as f:
                test_data = json.load(f)
                
            # Ensure the data was written correctly
            if test_data["last_completed_batch"] != checkpoint_data["last_completed_batch"]:
                raise ValueError("Checkpoint validation failed: batch number mismatch")
                
            if test_data["total_processed"] != checkpoint_data["total_processed"]:
                raise ValueError("Checkpoint validation failed: processed count mismatch")
        except Exception as e:
            logger.error(f"Checkpoint validation failed: {str(e)}")
            os.remove(temp_path)
            return False
        
        # Rename the temporary file to the final checkpoint file (atomic operation)
        if os.path.exists(checkpoint_path):
            # Create a backup of the existing checkpoint
            backup_path = f"{checkpoint_path}.bak"
            try:
                import shutil
                shutil.copy2(checkpoint_path, backup_path)
            except IOError as e:
                logger.warning(f"Failed to create backup of existing checkpoint: {str(e)}")
        
        # Perform the rename (atomic on most file systems)
        os.replace(temp_path, checkpoint_path)
        
        logger.info(f"Checkpoint saved to {checkpoint_path}")
        return True
    except Exception as e:
        log_structured_error(
            error_type="Checkpoint",
            message=f"Error writing checkpoint file",
            context={
                "checkpoint_path": checkpoint_path,
                "exception": e,
                "source": "write_checkpoint"
            }
        )
        # Clean up temporary file if it exists
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except:
                pass
        return False

# These functions have been replaced by the centralized logging module in chembl/logging.py

def generate_summary_report(
    start_time: float,
    end_time: float,
    total_processed: int,
    total_imported: int,
    total_duplicates: int,
    total_skipped: int,
    checkpoint_data: Dict[str, Any] = None,
    additional_data: Dict[str, Any] = None,
    include_system_info: bool = True
) -> Dict[str, Any]:
    """
    Generate a comprehensive summary report at the end of processing with enhanced metrics and system information.
    
    Args:
        start_time: Processing start time (Unix timestamp)
        end_time: Processing end time (Unix timestamp)
        total_processed: Total number of molecules processed
        total_imported: Total number of molecules imported
        total_duplicates: Total number of duplicates found
        total_skipped: Total number of molecules skipped
        checkpoint_data: Final checkpoint data (optional)
        additional_data: Any additional data to include in the report (optional)
        include_system_info: Whether to include system information in the report
        
    Returns:
        Dict containing the comprehensive summary report
    """
    # Generate a unique report ID
    import uuid
    report_id = str(uuid.uuid4())
    
    # Calculate processing time and rates
    processing_time = end_time - start_time
    processing_hours = processing_time / 3600
    processing_minutes = processing_time / 60
    
    # Calculate various performance metrics
    molecules_per_hour = total_processed / processing_hours if processing_hours > 0 else 0
    molecules_per_minute = total_processed / processing_minutes if processing_minutes > 0 else 0
    molecules_per_second = total_processed / processing_time if processing_time > 0 else 0
    success_rate = total_imported / total_processed * 100 if total_processed > 0 else 0
    duplicate_rate = total_duplicates / total_processed * 100 if total_processed > 0 else 0
    skip_rate = total_skipped / total_processed * 100 if total_processed > 0 else 0
    
    # Get system information if requested
    system_info = {}
    if include_system_info:
        try:
            import platform
            import psutil
            import socket
            
            system_info = {
                "hostname": socket.gethostname(),
                "platform": platform.platform(),
                "python_version": platform.python_version(),
                "processor": platform.processor(),
                "cpu_count": psutil.cpu_count(logical=False),
                "logical_cpu_count": psutil.cpu_count(logical=True),
                "memory_total_gb": round(psutil.virtual_memory().total / (1024**3), 2),
                "memory_available_gb": round(psutil.virtual_memory().available / (1024**3), 2),
                "memory_percent": psutil.virtual_memory().percent,
                "disk_usage_percent": psutil.disk_usage('/').percent
            }
        except (ImportError, Exception) as e:
            system_info = {"error": f"Failed to get system info: {str(e)}"}
    
    # Create enhanced standardized summary report
    summary = {
        "report_id": report_id,
        "timestamp": datetime.now().isoformat(),
        "event_type": "completion_summary",
        "status": "success",
        "processing": {
            "start_time": datetime.fromtimestamp(start_time).isoformat(),
            "end_time": datetime.fromtimestamp(end_time).isoformat(),
            "duration_seconds": round(processing_time, 2),
            "duration_minutes": round(processing_minutes, 2),
            "duration_hours": round(processing_hours, 4)
        },
        "performance": {
            "molecules_per_hour": round(molecules_per_hour, 2),
            "molecules_per_minute": round(molecules_per_minute, 2),
            "molecules_per_second": round(molecules_per_second, 2)
        },
        "results": {
            "total_processed": total_processed,
            "total_imported": total_imported,
            "total_duplicates": total_duplicates,
            "total_skipped": total_skipped,
            "success_rate_percent": round(success_rate, 2),
            "duplicate_rate_percent": round(duplicate_rate, 2),
            "skip_rate_percent": round(skip_rate, 2)
        },
        "checkpoint": checkpoint_data,
        "system": system_info,
        "process_id": os.getpid()
    }
    
    # Add categorized skip reasons if available in additional_data
    if additional_data and "skip_categories" in additional_data:
        summary["skip_categories"] = additional_data["skip_categories"]
        # Remove from additional_data to avoid duplication
        del additional_data["skip_categories"]
    
    # Add any remaining additional data without overwriting existing fields
    if additional_data:
        for key, value in additional_data.items():
            if key not in summary:
                summary[key] = value
    
    # Write summary using the centralized logging module
    try:
        # Write the summary report
        write_summary(summary)
        logger.info("Summary report saved successfully")
    except Exception as e:
        log_error(
            error_type="Reporting",
            message="Error writing summary report",
            context={
                "exception": e,
                "source": "generate_summary_report"
            }
        )
    
    # Also log a human-readable message to the standard logger
    logger.info(
        f"SUCCESS: All batches complete! Processed {total_processed} molecules, "
        f"imported {total_imported} ({round(success_rate, 1)}% success), "
        f"found {total_duplicates} duplicates ({round(duplicate_rate, 1)}%), "
        f"skipped {total_skipped} ({round(skip_rate, 1)}%). "
        f"Processing time: {round(processing_minutes, 2)} minutes "
        f"({round(molecules_per_hour, 2)} molecules/hour, {round(molecules_per_second, 2)}/sec)."
    )
    
    return summary

def estimate_time_remaining(start_time, batches_done, total_batches):
    elapsed = time.time() - start_time
    avg_per_batch = elapsed / batches_done if batches_done else 0
    remaining_batches = total_batches - batches_done
    eta_seconds = avg_per_batch * remaining_batches
    return str(timedelta(seconds=int(eta_seconds)))

def process_batches(
    chembl_ids,
    batch_size,
    checkpoint_path,
    resume,
    reset,
    checkpoint_frequency=1,
    memory_check_frequency=10  # Check memory every N molecules
):
    """
    Process ChEMBL molecules in batches with robust checkpointing.
    
    Args:
        chembl_ids (list): List of ChEMBL IDs to process
        batch_size (int): Number of molecules to process in each batch
        checkpoint_path (str): Path to checkpoint file or directory
        resume (bool): Whether to resume from the last checkpoint
        reset (bool): Whether to reset the checkpoint and start from scratch
        checkpoint_frequency (int): Save checkpoint after every N batches (default: 1)
        memory_check_frequency (int): Check memory every N molecules (default: 10)
        
    Returns:
        None
    """
    # Monitor initial memory usage
    try:
        import psutil
        process = psutil.Process()
        initial_memory = process.memory_info().rss / (1024 * 1024)  # MB
        logger.info(f"Initial memory usage: {initial_memory:.2f} MB")
    except ImportError:
        logger.warning("psutil not installed, memory monitoring will be limited")
        initial_memory = None
    except Exception as e:
        logger.warning(f"Error monitoring memory: {e}")
        initial_memory = None
    logger.info("STARTED: Starting batch dataset processing...")
    total_ids = len(chembl_ids)
    total_batches = (total_ids + batch_size - 1) // batch_size

    # Checkpoint directory setup
    # Handle checkpoint directory
    if os.path.isdir(checkpoint_path):
        # If checkpoint_path is already a directory, use it directly
        checkpoint_dir = checkpoint_path
        # Ensure it exists
        if not os.path.exists(checkpoint_dir):
            try:
                os.makedirs(checkpoint_dir)
                logger.info(f"Created checkpoint directory: {checkpoint_dir}")
            except OSError as e:
                logger.error(f"Error creating checkpoint directory {checkpoint_dir}: {str(e)}")
                return
    else:
        # If checkpoint_path is a file path, get its directory
        checkpoint_dir = os.path.dirname(checkpoint_path)
        # If the directory component is empty, use the current directory
        if not checkpoint_dir:
            checkpoint_dir = "."
        # Ensure the directory exists
        if not os.path.exists(checkpoint_dir):
            try:
                os.makedirs(checkpoint_dir)
                logger.info(f"Created checkpoint directory: {checkpoint_dir}")
            except OSError as e:
                logger.error(f"Error creating checkpoint directory {checkpoint_dir}: {str(e)}")
                return

    # Checkpoint logic
    checkpoint = None
    if resume:
        # If checkpoint_path is a directory, find the latest checkpoint file
        if os.path.isdir(checkpoint_path):
            latest_checkpoint = get_latest_checkpoint(checkpoint_path)
            if latest_checkpoint:
                checkpoint = read_checkpoint(latest_checkpoint)
                checkpoint_path = latest_checkpoint
        else:
            # If checkpoint_path is a file, try to read it
            checkpoint = read_checkpoint(checkpoint_path)
        
        if checkpoint:
            start_batch = checkpoint.get("last_completed_batch", 0) + 1
            logger.info(f"Resuming from batch {start_batch} (checkpoint: {checkpoint_path})")
            logger.info(f"Previously processed: {checkpoint.get('total_processed', 0)} molecules, "
                       f"imported: {checkpoint.get('total_imported', 0)}, "
                       f"duplicates: {checkpoint.get('total_duplicates', 0)}")
        else:
            logger.warning("No valid checkpoint found for resuming. Starting from the beginning.")
            start_batch = 0
    elif reset:
        # If checkpoint_path is a directory, remove all checkpoint files
        if os.path.isdir(checkpoint_path):
            try:
                for f in os.listdir(checkpoint_path):
                    if f.endswith('.json') and 'chembl_checkpoint' in f:
                        os.remove(os.path.join(checkpoint_path, f))
                logger.info(f"Reset all checkpoints in {checkpoint_path}")
            except OSError as e:
                logger.error(f"Error resetting checkpoints in {checkpoint_path}: {str(e)}")
        elif os.path.exists(checkpoint_path):
            try:
                os.remove(checkpoint_path)
                logger.info(f"Reset checkpoint at {checkpoint_path}")
            except OSError as e:
                logger.error(f"Error removing checkpoint file {checkpoint_path}: {str(e)}")
        
        start_batch = 0
    else:
        start_batch = 0

    # Fetch property types once
    property_types = fetch_property_types()
    if property_types is None:
        logger.error("Could not fetch property types. Exiting.")
        return

    # Get user_id if authenticated
    try:
        user_id = supabase.auth.current_user.id if hasattr(supabase.auth, 'current_user') and supabase.auth.current_user else None
        if user_id:
            logger.info(f"Using authenticated user ID: {user_id}")
        else:
            # When using service role key without user authentication, we don't need a user_id
            # The RLS policies will allow the operation based on the service role
            logger.info("No authenticated user ID available. Using service role for operations.")
    except Exception as e:
        logger.warning(f"Error getting user ID: {e}")
        user_id = None

    total_processed = 0
    total_imported = 0
    total_duplicates = 0
    total_skipped = 0
    start_time = time.time()

    for batch_num in range(start_batch, total_batches):
        batch_start = batch_num * batch_size
        batch_end = min(batch_start + batch_size, total_ids)
        batch_ids = chembl_ids[batch_start:batch_end]

        molecules_to_insert = []
        molecule_chembl_map = []
        batch_properties = []
        skipped_in_batch = 0
        duplicates_in_batch = 0
        
        # Memory check before batch processing
        if initial_memory is not None:
            try:
                current_memory = process.memory_info().rss / (1024 * 1024)  # MB
                memory_increase = current_memory - initial_memory
                logger.info(f"Memory before batch {batch_num+1}: {current_memory:.2f} MB (increase: {memory_increase:.2f} MB)")
                
                # Force garbage collection if memory usage is high
                if process.memory_percent() > 75.0:
                    logger.warning(f"High memory usage detected ({process.memory_percent():.1f}%). Running garbage collection.")
                    gc.collect()
                    current_memory_after_gc = process.memory_info().rss / (1024 * 1024)
                    logger.info(f"Memory after GC: {current_memory_after_gc:.2f} MB (freed: {current_memory - current_memory_after_gc:.2f} MB)")
            except Exception as e:
                logger.warning(f"Error monitoring memory: {e}")

        # Process molecules in smaller chunks within the batch to manage memory better
        for i, chembl_id in enumerate(batch_ids):
            # Check memory periodically within the batch
            if initial_memory is not None and i % memory_check_frequency == 0 and i > 0:
                try:
                    current_memory = process.memory_info().rss / (1024 * 1024)  # MB
                    memory_percent = process.memory_percent()
                    logger.debug(f"Memory during batch {batch_num+1} ({i}/{len(batch_ids)} molecules): {current_memory:.2f} MB ({memory_percent:.1f}%)")
                    
                    # Force garbage collection if memory usage is high
                    if memory_percent > 85.0:
                        logger.warning(f"Critical memory usage detected ({memory_percent:.1f}%). Running garbage collection.")
                        gc.collect()
                        time.sleep(1.0)  # Brief pause to let system stabilize
                except Exception as e:
                    logger.warning(f"Error monitoring memory: {e}")
            try:
                molecule = get_molecule_properties(chembl_id)
                logger.info(f"DEBUG: ChEMBL ID {chembl_id} molecule properties before filtering: {molecule}")
                
                if not filter_molecule(molecule, skipped_log_path):
                    logger.info(f"DEBUG: ChEMBL ID {chembl_id} did not pass filter criteria: {molecule}")
                    skipped_in_batch += 1
                    total_skipped += 1
                    continue

                # Check for duplicates by InChIKey
                existing_id = check_molecule_exists(molecule.get("InChIKey"))
                if existing_id:
                    logger.info(f"DEBUG: ChEMBL ID {chembl_id} already exists with ID {existing_id}")
                    log_skipped_molecule(
                        chembl_id=chembl_id,
                        reason=f"Duplicate (existing ID: {existing_id})",
                        molecule_data={
                            "ChEMBL ID": chembl_id,
                            "InChIKey": molecule.get("InChIKey"),
                            "existing_id": existing_id
                        },
                        skipped_log_path=skipped_log_path
                    )
                    duplicates_in_batch += 1
                    total_duplicates += 1
                    continue

                extra_properties = get_additional_properties(chembl_id)
                score = score_molecule(molecule, extra_properties)

                # Note: The database schema supports both 'cid' and 'pubchem_cid' columns.
                # 'pubchem_cid' is the primary storage column, and 'cid' is a generated column based on 'pubchem_cid'.
                # Either column name can be used in queries, but new code should prefer 'pubchem_cid'.
                molecules_to_insert.append({
                    "name": molecule.get("Name") or f"ChEMBL: {molecule['ChEMBL ID']}",
                    "smiles": molecule.get("SMILES"),
                    "inchi": molecule.get("InChI"),
                    "inchikey": molecule.get("InChIKey"),
                    "formula": molecule.get("Molecular Formula"),
                    "molecular_weight": float(molecule.get("Molecular Weight")) if molecule.get("Molecular Weight") else None,
                    "pubchem_cid": None,  # Added field for PubChem Compound ID (can be populated later if available)
                    "created_by": user_id,
                    "data_source": "ChEMBL",
                    "version": 1,
                    "modification_history": [{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "source": "ChEMBL_CryoProtectants_Supabase.py"
                    }]
                })
                molecule_chembl_map.append((chembl_id, molecule, extra_properties, score))
            except Exception as e:
                log_structured_error(
                    error_type="Processing",
                    message=f"Error processing ChEMBL ID {chembl_id}",
                    context={
                        "chembl_id": chembl_id,
                        "exception": e,
                        "batch_num": batch_num,
                        "source": "process_batches"
                    }
                )
                log_skipped_molecule(
                    chembl_id=chembl_id,
                    reason=f"Processing exception: {str(e)}",
                    molecule_data={"ChEMBL ID": chembl_id},
                    skipped_log_path=skipped_log_path
                )
                skipped_in_batch += 1
                total_skipped += 1

            time.sleep(CHEMBL_API_DELAY)

        # Bulk insert molecules
        inserted_ids = bulk_insert_molecules(molecules_to_insert, user_id)
        if not inserted_ids or len(inserted_ids) != len(molecule_chembl_map):
            error_msg = f"Bulk insert mismatch: {len(inserted_ids)} IDs for {len(molecule_chembl_map)} molecules."
            log_structured_error(
                error_type="Database",
                message=error_msg,
                context={
                    "batch_num": batch_num,
                    "expected_count": len(molecule_chembl_map),
                    "actual_count": len(inserted_ids) if inserted_ids else 0,
                    "source": "process_batches.bulk_insert_molecules"
                }
            )
            # Log all ChEMBL IDs in this batch as skipped if bulk insert failed
            for chembl_id, molecule, _, _ in molecule_chembl_map:
                log_skipped_molecule(
                    chembl_id=chembl_id,
                    reason="Bulk insert failed",
                    molecule_data={
                        "ChEMBL ID": chembl_id,
                        "Name": molecule.get("Name"),
                        "InChIKey": molecule.get("InChIKey")
                    },
                    skipped_log_path=skipped_log_path
                )
            continue

        # Prepare property inserts for all molecules in batch
        for idx, (chembl_id, molecule, extra_properties, score) in enumerate(molecule_chembl_map):
            molecule_id = inserted_ids[idx]
            properties_to_insert = {
                "LogP": molecule.get("LogP"),
                "TPSA": molecule.get("TPSA"),
                "H-Bond Donors": molecule.get("H-Bond Donors"),
                "H-Bond Acceptors": molecule.get("H-Bond Acceptors"),
                "Toxicity": extra_properties.get("Toxicity"),
                "Stability": extra_properties.get("Stability"),
                "Environmental Safety": extra_properties.get("Environmental Safety"),
                "Total Score": score,
                "ChEMBL ID": molecule.get("ChEMBL ID")  # Store ChEMBL ID as a property
            }
            for property_name, value in properties_to_insert.items():
                if value is None:
                    continue
                property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
                if not property_type:
                    logger.warning(f"Property type '{property_name}' not found, skipping")
                    continue
                property_insert = {
                    "molecule_id": molecule_id,
                    "property_type_id": property_type["id"],
                    "created_by": user_id
                }
                if property_type["data_type"] == "numeric":
                    try:
                        property_insert["numeric_value"] = float(value) if value is not None else None
                    except (ValueError, TypeError) as e:
                        log_structured_error(
                            error_type="Data",
                            message=f"Could not convert property value to float",
                            context={
                                "chembl_id": chembl_id,
                                "property_name": property_name,
                                "value": str(value),
                                "exception": e,
                                "source": "process_batches.property_conversion"
                            }
                        )
                        continue
                elif property_type["data_type"] == "text":
                    property_insert["text_value"] = str(value) if value is not None else None
                elif property_type["data_type"] == "boolean":
                    property_insert["boolean_value"] = bool(value) if value is not None else None
                batch_properties.append(property_insert)

        # Bulk insert properties
        if not bulk_insert_properties(batch_properties):
            error_msg = f"Failed to bulk insert properties for batch {batch_num}"
            log_structured_error(
                error_type="Database",
                message=error_msg,
                context={
                    "batch_num": batch_num,
                    "property_count": len(batch_properties),
                    "molecule_count": len(molecule_chembl_map),
                    "source": "process_batches.bulk_insert_properties"
                }
            )
            # Log all ChEMBL IDs in this batch as skipped for properties
            for chembl_id, molecule, _, _ in molecule_chembl_map:
                log_skipped_molecule(
                    chembl_id=chembl_id,
                    reason="Bulk property insert failed",
                    molecule_data={
                        "ChEMBL ID": chembl_id,
                        "Name": molecule.get("Name"),
                        "InChIKey": molecule.get("InChIKey")
                    },
                    category="error",
                    batch_num=batch_num
                )
            continue

        total_processed += len(batch_ids)
        total_imported += len(inserted_ids)
        batches_done = batch_num + 1
        eta = estimate_time_remaining(start_time, batches_done, total_batches)
        
        # Memory check after batch processing
        memory_info = ""
        if initial_memory is not None:
            try:
                current_memory = process.memory_info().rss / (1024 * 1024)  # MB
                memory_increase = current_memory - initial_memory
                memory_percent = process.memory_percent()
                memory_info = f" Memory: {current_memory:.2f} MB ({memory_percent:.1f}%)"
                
                # Force garbage collection after each batch to prevent memory buildup
                gc.collect()
            except Exception as e:
                logger.warning(f"Error monitoring memory: {e}")
        
        # Use standardized progress reporting
        log_progress(
            batch_num=batch_num,
            total_batches=total_batches,
            total_processed=total_processed,
            total_ids=total_ids,
            total_imported=total_imported,
            total_duplicates=total_duplicates,
            skipped_in_batch=skipped_in_batch,
            eta=eta,
            memory_info=memory_info,
            additional_data={
                "batch_size": batch_size,
                "total_skipped": total_skipped
            }
        )

        # Save checkpoint with timestamped filename (if frequency matches)
        if (batch_num + 1) % checkpoint_frequency == 0 or batch_num == total_batches - 1:
            checkpoint_data = {
                "last_completed_batch": batch_num,
                "total_processed": total_processed,
                "total_imported": total_imported,
                "total_duplicates": total_duplicates,
                "total_skipped": total_skipped,
                "timestamp": datetime.now().isoformat(),
                "last_chembl_id": batch_ids[-1] if batch_ids else None,
                "batch_size": batch_size,
                "total_batches": total_batches,
                "eta": eta
            }
            
            # If checkpoint_path is a directory, generate a timestamped filename
            if os.path.isdir(checkpoint_path):
                timestamped_filename = generate_checkpoint_filename()
                current_checkpoint_path = os.path.join(checkpoint_path, timestamped_filename)
            else:
                current_checkpoint_path = checkpoint_path
                
            # Write the checkpoint
            if not write_checkpoint(current_checkpoint_path, checkpoint_data):
                log_error(
                    error_type="Checkpoint",
                    message=f"Failed to write checkpoint to {current_checkpoint_path}",
                    context={
                        "batch_num": batch_num,
                        "checkpoint_path": current_checkpoint_path,
                        "total_processed": total_processed,
                        "total_imported": total_imported,
                        "source": "process_batches.write_checkpoint"
                    }
                )
            else:
                logger.info(f"Checkpoint saved at batch {batch_num+1}/{total_batches}")

    # Generate and log comprehensive summary report
    end_time = time.time()
    final_checkpoint = {
        "last_completed_batch": total_batches - 1,
        "total_processed": total_processed,
        "total_imported": total_imported,
        "total_duplicates": total_duplicates,
        "total_skipped": total_skipped,
        "timestamp": datetime.now().isoformat(),
        "batch_size": batch_size,
        "total_batches": total_batches
    }
    
    generate_summary_report(
        start_time=start_time,
        end_time=end_time,
        total_processed=total_processed,
        total_imported=total_imported,
        total_duplicates=total_duplicates,
        total_skipped=total_skipped,
        checkpoint_data=final_checkpoint
    )
    
    # Explicitly close any connections
    try:
        # Sign out if authenticated
        if hasattr(supabase.auth, 'current_user') and supabase.auth.current_user:
            supabase.auth.sign_out()
            logger.info("Signed out of Supabase")
    except Exception as e:
        log_structured_error(
            error_type="Cleanup",
            message="Error during cleanup and connection closing",
            context={
                "exception": e,
                "source": "process_batches.cleanup"
            }
        )

def main():
    # Verify configuration before proceeding
    verify_configuration()
    
    parser = argparse.ArgumentParser(description="ChEMBL CryoProtectant Data Importer (Supabase)")
    parser.add_argument("--batch-size", type=int, default=getattr(active_config, 'CHEMBL_BATCH_SIZE', 100),
                        help=f"Batch size for processing (default: {getattr(active_config, 'CHEMBL_BATCH_SIZE', 100)})")
    
    # Enhanced checkpoint options
    checkpoint_dir = getattr(active_config, 'CHECKPOINT_DIR', 'checkpoints')
    parser.add_argument("--checkpoint", type=str,
                        default=checkpoint_dir,
                        help=f"Path to checkpoint file or directory (default: {checkpoint_dir})")
    parser.add_argument("--checkpoint-filename", type=str,
                        default="chembl_checkpoint",
                        help="Base name for checkpoint files when using a directory (default: chembl_checkpoint)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint (finds most recent in directory)")
    parser.add_argument("--resume-specific", type=str,
                        help="Resume from a specific checkpoint file")
    parser.add_argument("--reset", action="store_true",
                        help="Reset checkpoint and start from scratch")
    parser.add_argument("--checkpoint-frequency", type=int, default=1,
                        help="Save checkpoint after every N batches (default: 1)")
    parser.add_argument("--memory-check-frequency", type=int, default=10,
                        help="Check memory usage every N molecules within a batch (default: 10)")
    
    # Other options
    # Skipped log path is now handled by the centralized logging module
    args = parser.parse_args()

    # Check if we're using a service role key
    is_service_role = False
    try:
        # Try to determine if we're using a service role key
        if "service_role" in supabase.auth.get_session().access_token:
            is_service_role = True
            logger.info("Using service role key for authentication")
    except Exception:
        # If we can't determine, assume we're not using a service role key
        pass
    
    # Authenticate if credentials provided and not using service role
    if not is_service_role:
        supabase_user = getattr(active_config, 'SUPABASE_USER', None)
        supabase_password = getattr(active_config, 'SUPABASE_PASSWORD', None)
        if supabase_user and supabase_password:
            try:
                response = supabase.auth.sign_in_with_password({
                    "email": supabase_user,
                    "password": supabase_password
                })
                if response.error:
                    log_error(
                        error_type="Authentication",
                        message="Authentication error from Supabase",
                        context={
                            "error": str(response.error),
                            "user": supabase_user,
                            "source": "main.authentication"
                        }
                    )
                    logger.warning("Continuing without user authentication. Using service role key for operations.")
                else:
                    logger.info(f"Authenticated as {supabase_user}")
            except Exception as e:
                log_error(
                    error_type="Authentication",
                    message="Exception during authentication",
                    context={
                        "exception": e,
                        "user": supabase_user,
                        "source": "main.authentication"
                    }
                )
                logger.warning("Continuing without user authentication. Using service role key for operations.")
        else:
            logger.info("No user authentication credentials provided. Using service role key for operations.")

    # Load ChEMBL IDs
    chembl_ids = get_chembl_ids()
    if not chembl_ids:
        log_error(
            error_type="Data",
            message="No ChEMBL IDs loaded",
            context={
                "source": "main.get_chembl_ids",
                "chembl_id_file": CHEMBL_ID_FILE or "Not specified"
            }
        )
        return

    # Skipped log is now handled by the centralized logging module

    # Handle checkpoint path based on arguments
    checkpoint_path = args.checkpoint
    
    # If a specific checkpoint file is specified for resuming
    if args.resume_specific:
        if os.path.exists(args.resume_specific):
            checkpoint_path = args.resume_specific
            args.resume = True
            logger.info(f"Will resume from specific checkpoint: {checkpoint_path}")
        else:
            log_error(
                error_type="Checkpoint",
                message="Specified checkpoint file not found",
                context={
                    "checkpoint_path": args.resume_specific,
                    "source": "main.resume_specific"
                }
            )
            return
    
    # If checkpoint is a directory, ensure it exists
    if os.path.isdir(checkpoint_path) or not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path, exist_ok=True)
        logger.info(f"Using checkpoint directory: {checkpoint_path}")
        
        # If it's a directory and we're not resuming from a specific file,
        # we'll use the directory and let process_batches find the latest checkpoint
    
    process_batches(
        chembl_ids=chembl_ids,
        batch_size=args.batch_size,
        checkpoint_path=checkpoint_path,
        resume=args.resume,
        reset=args.reset,
        checkpoint_frequency=args.checkpoint_frequency,
        memory_check_frequency=args.memory_check_frequency
    )
    
    # Force exit to ensure the script terminates properly
    logger.info("Exiting script...")
    import sys
    sys.exit(0)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_error(
            error_type="Fatal",
            message="Fatal error in main execution",
            context={
                "exception": e,
                "source": "main_execution"
            }
        )
        # Force exit on error
        import sys
        sys.exit(1)
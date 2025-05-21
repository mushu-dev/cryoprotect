#!/usr/bin/env python3
"""
ChEMBL Simplified Import Script

This script imports cryoprotectant data from ChEMBL using the official
chembl_webresource_client package with the new simplified database module.
Features:

- Streamlined database operations
- Official ChEMBL client integration
- Structured logging
- Checkpointing for resumable operations
- Batch processing with transaction support
- Comprehensive error handling and reporting
"""

import os
import sys
import json
import time
import argparse
import traceback
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Dict, Any, Optional

# Import the official ChEMBL client if available, otherwise mock it for testing
# Import the official ChEMBL client if available, otherwise mock it for testing
CHEMBL_CLIENT_AVAILABLE = False
TEST_MODE = True  # Set to True to use mock client for testing

if TEST_MODE:
    print("Running in test mode with mock ChEMBL client")
    # Create a simple mock client for testing
    class MockChemblClient:
        def filter(self, *args, **kwargs):
            return self
        def only(self, *args):
            return self
        def __getitem__(self, key):
            # Create a list of mock molecule objects
            mock_data = [
                {
                    "molecule_chembl_id": "CHEMBL25",
                    "pref_name": "ASPIRIN",
                    "molecular_formula": "C9H8O4",
                    "molecule_properties": {
                        "full_mwt": 180.16,
                        "alogp": 1.23,
                        "psa": 63.6,
                        "hba": 4,
                        "hbd": 1,
                        "num_ro5_violations": 0,
                        "med_chem_friendly": "Yes"
                    },
                    "molecule_structures": {
                        "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                        "standard_inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                        "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
                    },
                    "get": lambda x, default=None: {
                        "molecule_chembl_id": "CHEMBL25",
                        "pref_name": "ASPIRIN",
                        "molecular_formula": "C9H8O4"
                    }.get(x, default)
                },
                {
                    "molecule_chembl_id": "CHEMBL122",
                    "pref_name": "GLUCOSE",
                    "molecular_formula": "C6H12O6",
                    "molecule_properties": {
                        "full_mwt": 180.16,
                        "alogp": -2.93,
                        "psa": 110.38,
                        "hba": 6,
                        "hbd": 5,
                        "num_ro5_violations": 1,
                        "med_chem_friendly": "Yes"
                    },
                    "molecule_structures": {
                        "canonical_smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
                        "standard_inchi": "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
                        "standard_inchi_key": "WQZGKKKJIJFFOK-GASJEMHNSA-N"
                    },
                    "get": lambda x, default=None: {
                        "molecule_chembl_id": "CHEMBL122",
                        "pref_name": "GLUCOSE",
                        "molecular_formula": "C6H12O6"
                    }.get(x, default)
                }
            ]
            
            # If key is a slice, return appropriate sliced list
            if isinstance(key, slice):
                return mock_data[key]
            else:
                # If key is an integer, return the single item
                return mock_data[key]

    def new_client(resource_name):
        return MockChemblClient()
else:
    try:
        from chembl_webresource_client.new_client import new_client
        CHEMBL_CLIENT_AVAILABLE = True
    except ImportError:
        print("Warning: chembl_webresource_client not found. Using mock client for testing.")
        # Create a simple mock client for testing
        class MockChemblClient:
            def filter(self, *args, **kwargs):
                return self
            def only(self, *args):
                return []

        def new_client(resource_name):
            return MockChemblClient()

        CHEMBL_CLIENT_AVAILABLE = False

# Import our new simplified database module
from database import db, utils
from dotenv import load_dotenv

# Constants
DEFAULT_BATCH_SIZE = 25
DEFAULT_TARGET_LIMIT = 500
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 5
DEFAULT_THROTTLE_DELAY = 0.5

# Properties we're interested in
# Property mapping is split into molecule data and property data
MOLECULE_MAPPING = {
    "molecule_chembl_id": "chembl_id",
    "pref_name": "name",
    "molecular_formula": "formula",
    "molecule_structures.canonical_smiles": "smiles",
    "molecule_structures.standard_inchi": "inchi",
    "molecule_structures.standard_inchi_key": "inchikey"
}

PROPERTY_MAPPING = {
    "molecule_properties.full_mwt": "molecular_weight",
    "molecule_properties.alogp": "logp",
    "molecule_properties.psa": "tpsa",
    "molecule_properties.hba": "h_bond_acceptors",
    "molecule_properties.hbd": "h_bond_donors",
    "molecule_properties.num_ro5_violations": "ro5_violations",
    "molecule_properties.med_chem_friendly": "med_chem_friendly"
}

# Combined mapping for fetching from ChEMBL
ALL_PROPERTIES = {**MOLECULE_MAPPING, **PROPERTY_MAPPING}

# Filtering criteria for cryoprotectants
FILTER_CRITERIA = {
    "molecule_properties.full_mwt__range": (0, 1000),
    "molecule_properties.alogp__range": (-5, 5),
    "molecule_properties.psa__range": (0, 200),
    "molecule_properties.hbd__lte": 10,
    "molecule_properties.hba__lte": 10
}

# Set up logging directory
LOGS_DIR = Path("logs")
LOGS_DIR.mkdir(exist_ok=True)

# Set up checkpoints directory
CHECKPOINTS_DIR = Path("checkpoints")
CHECKPOINTS_DIR.mkdir(exist_ok=True)

# Set up logging
import logging
LOG_FILE = LOGS_DIR / "chembl_import_simplified.log"
SKIPPED_LOG = LOGS_DIR / "chembl_skipped_simplified.log"
ERROR_LOG = LOGS_DIR / "chembl_errors_simplified.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("chembl_import")

def log_skipped_molecule(chembl_id, reason):
    """Log a skipped molecule to a separate file."""
    with open(SKIPPED_LOG, "a") as f:
        f.write(f"{datetime.now().isoformat()}: {chembl_id} - {reason}\n")
    logger.warning(f"Skipped {chembl_id}: {reason}")

def log_error(chembl_id, error):
    """Log an error to a separate file."""
    with open(ERROR_LOG, "a") as f:
        f.write(f"{datetime.now().isoformat()}: {chembl_id} - {error}\n")
    logger.error(f"Error processing {chembl_id}: {error}")

def init_database():
    """Initialize database connection using configuration from environment variables."""
    # Load environment variables
    load_dotenv()
    
    # Create config from Supabase settings
    config = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
        'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    logger.info(f"Initializing database connection to {config['host']}:{config['port']}/{config['database']}")
    
    # Initialize database pool
    return db.init_connection_pool(config=config)

def get_checkpoint_path(name="chembl_import"):
    """Get the path to the checkpoint file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return CHECKPOINTS_DIR / f"{name}_{timestamp}.json"

def save_checkpoint(data, checkpoint_path=None):
    """Save checkpoint data to a file."""
    if checkpoint_path is None:
        checkpoint_path = get_checkpoint_path()
        
    # Add timestamp to data
    data["timestamp"] = datetime.now().isoformat()
    
    try:
        with open(checkpoint_path, "w") as f:
            json.dump(data, f, indent=2)
            
        logger.info(f"Checkpoint saved to {checkpoint_path}")
        return True
    except Exception as e:
        logger.error(f"Error saving checkpoint: {str(e)}")
        return False

def load_checkpoint(checkpoint_path):
    """Load checkpoint data from a file."""
    try:
        if not os.path.exists(checkpoint_path):
            logger.info(f"No checkpoint file found at {checkpoint_path}")
            return None
            
        with open(checkpoint_path, "r") as f:
            data = json.load(f)
            
        logger.info(f"Checkpoint loaded from {checkpoint_path}")
        return data
    except Exception as e:
        logger.error(f"Error loading checkpoint: {str(e)}")
        return None

def fetch_chembl_data(target_limit=DEFAULT_TARGET_LIMIT, batch_size=DEFAULT_BATCH_SIZE,
                      throttle_delay=DEFAULT_THROTTLE_DELAY, max_retries=DEFAULT_MAX_RETRIES,
                      retry_delay=DEFAULT_RETRY_DELAY, offset=0):
    """
    Fetch cryoprotectant data from ChEMBL.

    Args:
        target_limit: Maximum number of compounds to fetch
        batch_size: Number of compounds to fetch in each batch
        throttle_delay: Delay between API calls to avoid rate limiting
        max_retries: Maximum number of retry attempts
        retry_delay: Delay between retry attempts
        offset: Starting offset for pagination

    Returns:
        List of fetched compounds
    """
    if TEST_MODE:
        logger.info("Running in test mode with mock ChEMBL data")
        # Return mock data from the test client
        compound_client = new_client("molecule")
        mock_compounds = compound_client[0:10]
        return mock_compounds
    
    logger.info(f"Fetching data from ChEMBL API (limit={target_limit}, batch={batch_size}, offset={offset})")
    
    molecules = []
    current_offset = offset
    compound_client = new_client("molecule")
    
    while len(molecules) < target_limit:
        retry_count = 0
        batch_size_current = min(batch_size, target_limit - len(molecules))
        
        while retry_count <= max_retries:
            try:
                logger.info(f"Fetching batch of {batch_size_current} compounds from offset {current_offset}")
                
                # Apply filters to find potential cryoprotectants
                batch = compound_client.filter(
                    **FILTER_CRITERIA
                ).only(
                    list(ALL_PROPERTIES.keys())
                )[current_offset:current_offset + batch_size_current]
                
                # If no more results, break out of the loop
                if not batch:
                    logger.info("No more compounds available from ChEMBL")
                    break
                    
                # Add batch to molecules list
                molecules.extend(batch)
                current_offset += len(batch)
                
                # Log progress
                logger.info(f"Fetched {len(batch)} compounds ({len(molecules)}/{target_limit})")
                
                # Add artificial delay to avoid rate limiting
                time.sleep(throttle_delay)
                
                # Break out of retry loop on success
                break
                
            except Exception as e:
                retry_count += 1
                logger.warning(f"Error fetching from ChEMBL API: {str(e)}. Retry {retry_count}/{max_retries}")
                
                if retry_count > max_retries:
                    logger.error(f"Failed to fetch data after {max_retries} retries. Last error: {str(e)}")
                    break
                    
                # Wait before retrying with exponential backoff
                backoff_time = retry_delay * (2 ** (retry_count - 1))
                time.sleep(backoff_time)
        
        # If we've reached the target limit or no more results, break out of the main loop
        if len(molecules) >= target_limit or not batch:
            break
    
    logger.info(f"Fetched a total of {len(molecules)} compounds from ChEMBL")
    return molecules

def transform_chembl_data(molecule):
    """
    Transform ChEMBL molecule data into our database format.
    
    Args:
        molecule: ChEMBL molecule data
        
    Returns:
        Tuple of (molecule_data, property_data)
    """
    # Create molecule data dictionary (only molecule fields)
    molecule_data = {
        "chembl_id": None,
        "name": None,
        "formula": None,
        "smiles": None,
        "inchi": None,
        "inchikey": None,
        "data_source": "ChEMBL"
    }
    
    # Create property data dictionary
    property_data = {}
    
    # Extract basic molecule properties at the root level
    for chembl_key, our_key in MOLECULE_MAPPING.items():
        if "." not in chembl_key:
            molecule_data[our_key] = molecule.get(chembl_key)
            
    # Extract nested molecule properties (structures)
    if "molecule_structures" in molecule and molecule["molecule_structures"]:
        for chembl_key, our_key in MOLECULE_MAPPING.items():
            if chembl_key.startswith("molecule_structures."):
                struct_name = chembl_key.split(".")[1]
                molecule_data[our_key] = molecule["molecule_structures"].get(struct_name)
    
    # Extract property values (these will go to molecular_properties table)
    if "molecule_properties" in molecule and molecule["molecule_properties"]:
        for chembl_key, our_key in PROPERTY_MAPPING.items():
            if chembl_key.startswith("molecule_properties."):
                prop_name = chembl_key.split(".")[1]
                property_value = molecule["molecule_properties"].get(prop_name)
                if property_value is not None:
                    property_data[our_key] = property_value
    
    return molecule_data, property_data

def validate_chembl_molecule(molecule_data, property_data):
    """
    Validate that a ChEMBL molecule has the required properties.
    
    Args:
        molecule_data: Transformed ChEMBL molecule data
        property_data: Transformed ChEMBL property data
        
    Returns:
        (is_valid, reason): Tuple of boolean and reason string
    """
    # Check for required fields
    if not molecule_data.get("chembl_id"):
        return False, "Missing ChEMBL ID"
        
    if not molecule_data.get("smiles"):
        return False, "Missing SMILES"
        
    if not molecule_data.get("inchi"):
        return False, "Missing InChI"
        
    if not molecule_data.get("inchikey"):
        return False, "Missing InChI Key"
        
    # Validate numerical properties if available
    try:
        # Check molecular weight
        mw = float(property_data.get("molecular_weight", 0))
        if mw > 0 and not (0 <= mw <= 1000):
            return False, f"Molecular weight {mw} outside range (0-1000)"
            
        # Check logP if available
        if "logp" in property_data:
            logp = float(property_data.get("logp", 0))
            if not (-5 <= logp <= 5):
                return False, f"LogP {logp} outside range (-5-5)"
            
        # Check TPSA if available
        if "tpsa" in property_data:
            tpsa = float(property_data.get("tpsa", 0))
            if not (0 <= tpsa <= 200):
                return False, f"TPSA {tpsa} outside range (0-200)"
            
    except (ValueError, TypeError):
        return False, "Invalid numerical property values"
        
    return True, ""

def insert_molecule(molecule_data):
    """
    Insert a ChEMBL molecule into the database.
    
    Args:
        molecule_data: Transformed and validated ChEMBL molecule data
        
    Returns:
        (success, molecule_id): Tuple of success boolean and molecule ID
    """
    try:
        chembl_id = molecule_data.get("chembl_id")
        
        # First check if molecule already exists by ChEMBL ID
        existing = db.execute_query(
            "SELECT id FROM molecules WHERE chembl_id = %s",
            (chembl_id,)
        )
        
        if existing and len(existing) > 0:
            # Molecule exists, update it
            logger.info(f"Updating existing molecule {chembl_id}")
            
            update_fields = []
            update_values = []
            
            for key, value in molecule_data.items():
                if key != "chembl_id" and value is not None:
                    update_fields.append(f"{key} = %s")
                    update_values.append(value)
                    
            if not update_fields:
                logger.warning(f"No fields to update for {chembl_id}")
                return True, existing[0]["id"]
                
            # Add timestamp and ChEMBL ID to parameters
            update_values.append(datetime.now().isoformat())
            update_values.append(chembl_id)
            
            # Execute update
            query = f"""
                UPDATE molecules
                SET {", ".join(update_fields)}, updated_at = %s
                WHERE chembl_id = %s
                RETURNING id
            """
            
            result = db.execute_query(query, tuple(update_values))
            
            if result and len(result) > 0:
                return True, result[0]["id"]
            else:
                log_error(chembl_id, "Failed to update molecule")
                return False, None
        else:
            # Molecule doesn't exist, insert it
            logger.info(f"Inserting new molecule {chembl_id}")
            
            # Prepare fields and values
            fields = []
            values = []
            placeholders = []
            
            for key, value in molecule_data.items():
                if value is not None:
                    fields.append(key)
                    values.append(value)
                    placeholders.append("%s")
                    
            # Add timestamps
            fields.extend(["created_at", "updated_at"])
            now = datetime.now().isoformat()
            values.extend([now, now])
            placeholders.extend(["%s", "%s"])
            
            # Execute insert
            query = f"""
                INSERT INTO molecules
                ({", ".join(fields)})
                VALUES ({", ".join(placeholders)})
                RETURNING id
            """
            
            result = db.execute_query(query, tuple(values))
            
            if result and len(result) > 0:
                return True, result[0]["id"]
            else:
                log_error(chembl_id, "Failed to insert molecule")
                return False, None
                
    except Exception as e:
        chembl_id = molecule_data.get("chembl_id", "unknown")
        log_error(chembl_id, f"Error inserting molecule: {str(e)}")
        logger.error(traceback.format_exc())
        return False, None

def insert_molecular_properties(molecule_id, property_data):
    """
    Insert molecular properties for a molecule.
    
    Args:
        molecule_id: Molecule ID
        property_data: Dictionary of property name to value
        
    Returns:
        Number of properties inserted
    """
    if not property_data:
        return 0
    
    # Define property types and their data types
    property_types = {
        "molecular_weight": "numeric",
        "logp": "numeric",
        "tpsa": "numeric",
        "h_bond_donors": "numeric",
        "h_bond_acceptors": "numeric",
        "ro5_violations": "numeric",
        "med_chem_friendly": "text"
    }
    
    with db.transaction() as cursor:
        inserted_count = 0
        
        for property_name, property_value in property_data.items():
            if property_value is None:
                continue
                
            property_type = property_types.get(property_name, "text")
            
            # Find or create property type
            cursor.execute(
                "SELECT id FROM property_types WHERE name = %s",
                (property_name,)
            )
            result = cursor.fetchone()
            
            if result:
                property_type_id = result["id"]
            else:
                # Create new property type
                cursor.execute(
                    """
                    INSERT INTO property_types
                    (name, data_type, created_at, updated_at)
                    VALUES (%s, %s, NOW(), NOW())
                    RETURNING id
                    """,
                    (property_name, property_type)
                )
                result = cursor.fetchone()
                property_type_id = result["id"] if result else None
                
            if not property_type_id:
                logger.warning(f"Failed to create property type {property_name}")
                continue
               
            # Check if the property already exists 
            cursor.execute(
                """
                SELECT id FROM molecular_properties 
                WHERE molecule_id = %s AND property_type_id = %s
                """,
                (molecule_id, property_type_id)
            )
            existing_prop = cursor.fetchone()
            
            if property_type == "numeric":
                if existing_prop:
                    # Update existing property
                    cursor.execute(
                        """
                        UPDATE molecular_properties
                        SET numeric_value = %s, updated_at = NOW()
                        WHERE molecule_id = %s AND property_type_id = %s
                        """,
                        (property_value, molecule_id, property_type_id)
                    )
                else:
                    # Insert new property
                    cursor.execute(
                        """
                        INSERT INTO molecular_properties
                        (molecule_id, property_type_id, numeric_value, source, created_at, updated_at)
                        VALUES (%s, %s, %s, %s, NOW(), NOW())
                        """,
                        (molecule_id, property_type_id, property_value, "ChEMBL")
                    )
            else:
                if existing_prop:
                    # Update existing property
                    cursor.execute(
                        """
                        UPDATE molecular_properties
                        SET text_value = %s, updated_at = NOW()
                        WHERE molecule_id = %s AND property_type_id = %s
                        """,
                        (str(property_value), molecule_id, property_type_id)
                    )
                else:
                    # Insert new property
                    cursor.execute(
                        """
                        INSERT INTO molecular_properties
                        (molecule_id, property_type_id, text_value, source, created_at, updated_at)
                        VALUES (%s, %s, %s, %s, NOW(), NOW())
                        """,
                        (molecule_id, property_type_id, str(property_value), "ChEMBL")
                    )
                
            inserted_count += 1
            
    return inserted_count

def process_chembl_molecules(molecules, batch_size=DEFAULT_BATCH_SIZE):
    """
    Process and import ChEMBL molecules in batches.
    
    Args:
        molecules: List of ChEMBL molecules
        batch_size: Number of molecules to process in each batch
        
    Returns:
        Dict with processing statistics
    """
    stats = {
        "total": len(molecules),
        "processed": 0,
        "inserted": 0,
        "updated": 0,
        "skipped": 0,
        "errors": 0,
        "property_count": 0
    }
    
    # Process in batches to avoid memory issues
    for i in range(0, len(molecules), batch_size):
        batch = molecules[i:i+batch_size]
        
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules) + batch_size - 1)//batch_size} ({len(batch)} molecules)")
        
        for molecule in batch:
            stats["processed"] += 1
            chembl_id = molecule.get("molecule_chembl_id", "unknown")
            
            try:
                # Transform ChEMBL data to our format (separating molecule and property data)
                molecule_data, property_data = transform_chembl_data(molecule)
                
                # Validate molecule data
                is_valid, reason = validate_chembl_molecule(molecule_data, property_data)
                
                if not is_valid:
                    log_skipped_molecule(chembl_id, reason)
                    stats["skipped"] += 1
                    continue
                    
                # Insert or update molecule
                success, molecule_id = insert_molecule(molecule_data)
                
                if not success:
                    stats["errors"] += 1
                    continue
                    
                # Determine if this was an insert or update
                existing = db.execute_query(
                    "SELECT created_at, updated_at FROM molecules WHERE id = %s",
                    (molecule_id,)
                )
                
                if existing and len(existing) > 0:
                    created_at = existing[0].get("created_at")
                    updated_at = existing[0].get("updated_at")
                    
                    if created_at and updated_at and (updated_at - created_at).total_seconds() < 5:
                        # If created_at and updated_at are within 5 seconds, treat as new insert
                        stats["inserted"] += 1
                    else:
                        # Otherwise treat as update
                        stats["updated"] += 1
                        
                # Insert molecular properties
                property_count = insert_molecular_properties(molecule_id, property_data)
                stats["property_count"] += property_count
                
                # Log progress
                if stats["processed"] % 10 == 0 or stats["processed"] == stats["total"]:
                    progress = stats["processed"] / stats["total"] * 100
                    logger.info(f"Progress: {progress:.1f}% ({stats['processed']}/{stats['total']})")
                    
            except Exception as e:
                log_error(chembl_id, f"Unhandled error: {str(e)}")
                logger.error(traceback.format_exc())
                stats["errors"] += 1
                
        # Save checkpoint after each batch
        checkpoint_data = {
            "stats": stats,
            "batch_idx": i // batch_size + 1,
            "total_batches": (len(molecules) + batch_size - 1) // batch_size
        }
        save_checkpoint(checkpoint_data)
        
    return stats

def main():
    """Main function to run the ChEMBL import."""
    parser = argparse.ArgumentParser(description="Import cryoprotectant data from ChEMBL")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help=f"Batch size for processing (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--limit", type=int, default=DEFAULT_TARGET_LIMIT, help=f"Maximum number of compounds to import (default: {DEFAULT_TARGET_LIMIT})")
    parser.add_argument("--throttle", type=float, default=DEFAULT_THROTTLE_DELAY, help=f"Delay between API calls in seconds (default: {DEFAULT_THROTTLE_DELAY})")
    parser.add_argument("--retry", type=int, default=DEFAULT_MAX_RETRIES, help=f"Maximum number of API retry attempts (default: {DEFAULT_MAX_RETRIES})")
    parser.add_argument("--retry-delay", type=float, default=DEFAULT_RETRY_DELAY, help=f"Initial delay between retries in seconds (default: {DEFAULT_RETRY_DELAY})")
    parser.add_argument("--checkpoint", type=str, help="Path to checkpoint file to resume from")
    args = parser.parse_args()
    
    # Initialize database connection
    if not init_database():
        logger.error("Failed to initialize database connection. Exiting.")
        return 1
        
    logger.info("Database connection initialized successfully.")
    
    try:
        start_time = time.time()
        
        # Load checkpoint if specified
        offset = 0
        if args.checkpoint:
            checkpoint = load_checkpoint(args.checkpoint)
            if checkpoint and "stats" in checkpoint:
                offset = checkpoint["stats"]["processed"]
                logger.info(f"Resuming from checkpoint with offset {offset}")
                
        # Fetch data from ChEMBL
        logger.info(f"Fetching up to {args.limit} compounds from ChEMBL")
        molecules = fetch_chembl_data(
            target_limit=args.limit,
            batch_size=args.batch_size,
            throttle_delay=args.throttle,
            max_retries=args.retry,
            retry_delay=args.retry_delay,
            offset=offset
        )
        
        if not molecules:
            logger.error("No molecules fetched from ChEMBL. Exiting.")
            return 1
            
        # Process molecules
        logger.info(f"Processing {len(molecules)} ChEMBL molecules")
        stats = process_chembl_molecules(molecules, batch_size=args.batch_size)
        
        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        elapsed_str = str(timedelta(seconds=int(elapsed_time)))
        
        # Generate final report
        report = {
            "timestamp": datetime.now().isoformat(),
            "elapsed_time": elapsed_str,
            "elapsed_seconds": elapsed_time,
            "stats": stats,
            "molecules_per_second": stats["processed"] / elapsed_time if elapsed_time > 0 else 0
        }
        
        # Save final report
        report_path = Path("reports") / f"chembl_import_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        Path("reports").mkdir(exist_ok=True)
        
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
            
        # Log summary
        logger.info(f"Import completed in {elapsed_str}")
        logger.info(f"Processed: {stats['processed']}/{stats['total']} molecules")
        logger.info(f"Results: {stats['inserted']} inserted, {stats['updated']} updated, {stats['skipped']} skipped, {stats['errors']} errors")
        logger.info(f"Properties: {stats['property_count']} properties added")
        logger.info(f"Report saved to {report_path}")
        
        return 0
        
    except KeyboardInterrupt:
        logger.info("Import interrupted by user")
        return 130
        
    except Exception as e:
        logger.error(f"Unhandled error: {str(e)}")
        logger.error(traceback.format_exc())
        return 1
        
    finally:
        # Close database connections
        db.close_all_connections()
        logger.info("Database connections closed")

if __name__ == "__main__":
    sys.exit(main())
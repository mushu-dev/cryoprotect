#!/usr/bin/env python3
"""
ChEMBL Integrated Import Script

This script implements an improved approach for importing cryoprotectant data from ChEMBL
using the official chembl_webresource_client package. It features:

- Official ChEMBL client integration
- Structured, centralized logging
- Checkpointing for resumable operations
- Direct Supabase connection for database operations
- Batch processing with transaction support
- Comprehensive error handling and reporting
"""

import os
import sys
import json
import time
import uuid
import signal
import atexit
import argparse
import traceback
import shutil
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Dict, Any, Optional, Union, Tuple

# Import the official ChEMBL client
from chembl_webresource_client.new_client import new_client

# Import configuration system
from config import active_config, validate_config, ConfigurationError

# Import the centralized logging system
from chembl.logging import (
    ChEMBLLogger, log_error, log_skipped_molecule,
    log_progress, write_summary, get_logger
)

# Import service role helper for user authentication
from service_role_helper import get_user_id, ensure_user_profile

# Import RLS verification and remediation utilities
from rls_utils import ensure_rls_restored

# Import direct Supabase connection
from supabase_direct import SupabaseDirectConnection

# Initialize logger
logger = get_logger(__name__)

# Create a global logger instance with default settings
# This will create the logs directory if it doesn't exist
chembl_logger = ChEMBLLogger(
    log_dir="logs",
    progress_log="chembl_integrated_import_progress.jsonl",
    error_log="chembl_integrated_import_errors.jsonl",
    skipped_log="chembl_integrated_import_skipped.jsonl",
    summary_log="chembl_integrated_import_summary.json",
    general_log="chembl_integrated_import.log"
)


# Ensure checkpoints directory exists
Path(active_config.CHECKPOINT_DIR).mkdir(parents=True, exist_ok=True)

# Default checkpoint file path
DEFAULT_CHECKPOINT_FILE = Path(active_config.CHECKPOINT_DIR) / "chembl_import_checkpoint.json"

# Initialize Supabase client
supabase = None

# Standard reference compounds to ensure they are always imported
standard_reference_ids = [
    "CHEMBL25",    # Aspirin
    "CHEMBL1118",  # Caffeine
    "CHEMBL1234",  # Glycerol (common cryoprotectant)
    "CHEMBL444",   # Glucose
    "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
    "CHEMBL9335",  # Dimethyl sulfoxide (DMSO, common cryoprotectant)
    "CHEMBL15151"  # Trehalose (common cryoprotectant)
]

# Initialize database connection
db = None


class ChEMBLProgressTracker:
    """Class to track and manage progress statistics for the ChEMBL import process."""
    
    def __init__(self, total_compounds, batch_size, checkpoint_file=DEFAULT_CHECKPOINT_FILE):
        self.total_compounds = total_compounds
        self.batch_size = batch_size
        self.total_batches = (total_compounds + batch_size - 1) // batch_size
        self.checkpoint_file = checkpoint_file
        
        # Progress statistics
        self.start_time = time.time()
        self.current_batch = 0
        self.total_processed = 0
        self.total_imported = 0
        self.total_skipped = 0
        self.total_errors = 0
        self.compounds_in_current_batch = 0
        self.batch_start_time = time.time()
        self.batch_times = []
        self.status = "Running"
        
        # Recent log messages (circular buffer)
        self.max_log_entries = 50
        self.recent_logs = []
        
        # Load from checkpoint if exists
        self._load_from_checkpoint()
        
        # Register signal handlers for graceful shutdown
        signal.signal(signal.SIGINT, self._handle_exit)
        signal.signal(signal.SIGTERM, self._handle_exit)
        
        # Register exit handler
        atexit.register(self._handle_exit)
    
    def _handle_exit(self, *args):
        """Handle script exit (save checkpoint and set status)."""
        if self.status == "Running":
            self.status = "Paused"
            self.save_checkpoint()
        # Don't exit here as it might interfere with other exit handlers
    
    def _load_from_checkpoint(self):
        """Load progress data from checkpoint file if it exists."""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, "r") as f:
                    checkpoint = json.load(f)
                
                # Basic checkpoint data
                self.current_batch = checkpoint.get("last_completed_batch", 0) + 1
                self.total_processed = checkpoint.get("total_processed", 0)
                self.total_imported = checkpoint.get("total_imported", 0)
                
                # Enhanced progress data
                self.total_skipped = checkpoint.get("total_skipped", 0)
                self.total_errors = checkpoint.get("total_errors", 0)
                self.batch_times = checkpoint.get("batch_times", [])
                
                # Adjust start time to maintain accurate elapsed time
                if "elapsed_seconds" in checkpoint:
                    self.start_time = time.time() - checkpoint["elapsed_seconds"]
                
                logger.info(f"Loaded progress from checkpoint: {self.total_processed}/{self.total_compounds} compounds processed")
            except Exception as e:
                logger.error(f"Error loading checkpoint: {str(e)}")
    
    def save_checkpoint(self):
        """Save progress data to checkpoint file."""
        checkpoint_data = {
            "last_completed_batch": self.current_batch - 1 if self.current_batch > 0 else 0,
            "total_processed": self.total_processed,
            "total_imported": self.total_imported,
            "total_skipped": self.total_skipped,
            "total_errors": self.total_errors,
            "batch_times": self.batch_times[-20:],  # Keep only the last 20 batch times
            "elapsed_seconds": time.time() - self.start_time,
            "timestamp": datetime.now().isoformat(),
            "status": self.status
        }
        
        with open(self.checkpoint_file, "w") as f:
            json.dump(checkpoint_data, f, indent=2)
    
    def start_batch(self, batch_num, batch_size):
        """Record the start of a new batch."""
        self.current_batch = batch_num + 1  # 1-indexed for display
        self.compounds_in_current_batch = batch_size
        self.batch_start_time = time.time()
        self.add_log_message(f"Starting batch {self.current_batch}/{self.total_batches} with {batch_size} compounds")
    
    def end_batch(self, processed, imported, skipped, errors):
        """Record the end of a batch."""
        batch_time = time.time() - self.batch_start_time
        self.batch_times.append(batch_time)
        
        self.total_processed += processed
        self.total_imported += imported
        self.total_skipped += skipped
        self.total_errors += errors
        
        # Calculate progress metrics
        progress_percent = round((self.total_processed / self.total_compounds) * 100, 1)
        eta = self.estimate_time_remaining()
        
        self.add_log_message(
            f"Completed batch {self.current_batch}/{self.total_batches}: "
            f"{processed} processed, {imported} imported, {skipped} skipped, {errors} errors. "
            f"Progress: {progress_percent}%, ETA: {eta}"
        )
        
        # Save checkpoint after each batch
        self.save_checkpoint()
    
    def add_error(self, error_message):
        """Record an error."""
        self.total_errors += 1
        self.add_log_message(f"ERROR: {error_message}")
    
    def add_skipped(self, chembl_id, reason):
        """Record a skipped compound."""
        self.total_skipped += 1
        self.add_log_message(f"SKIPPED: ChEMBL ID {chembl_id} - {reason}")
    
    def add_log_message(self, message):
        """Add a log message to the recent logs buffer."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_entry = f"{timestamp} - {message}"
        
        # Add to circular buffer
        self.recent_logs.append(log_entry)
        if len(self.recent_logs) > self.max_log_entries:
            self.recent_logs.pop(0)
    
    def estimate_time_remaining(self):
        """Estimate time remaining based on average batch processing time."""
        if not self.batch_times:
            return "Unknown"
        
        # Use the last 10 batch times for a more accurate recent average
        recent_batch_times = self.batch_times[-10:] if len(self.batch_times) >= 10 else self.batch_times
        avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
        
        remaining_batches = self.total_batches - self.current_batch
        eta_seconds = avg_batch_time * remaining_batches
        
        return str(timedelta(seconds=int(eta_seconds)))
    
    def get_avg_time_per_batch(self):
        """Get the average time per batch."""
        if not self.batch_times:
            return "00:00:00"
        
        avg_seconds = sum(self.batch_times) / len(self.batch_times)
        return str(timedelta(seconds=int(avg_seconds)))
    
    def get_elapsed_time(self):
        """Get the total elapsed time."""
        elapsed_seconds = time.time() - self.start_time
        return str(timedelta(seconds=int(elapsed_seconds)))
    
    def get_progress_percentage(self):
        """Get the progress percentage."""
        if self.total_compounds == 0:
            return 0
        return round((self.total_processed / self.total_compounds) * 100, 1)
    
    def set_status(self, status):
        """Set the current status of the import process."""
        self.status = status
        self.add_log_message(f"Status changed to: {status}")
        self.save_checkpoint()
    
    def get_progress_data(self):
        """Get all progress data as a dictionary for the dashboard."""
        return {
            "total_compounds": self.total_compounds,
            "total_processed": self.total_processed,
            "total_imported": self.total_imported,
            "total_skipped": self.total_skipped,
            "total_errors": self.total_errors,
            "current_batch": self.current_batch,
            "total_batches": self.total_batches,
            "compounds_in_current_batch": self.compounds_in_current_batch,
            "elapsed_time": self.get_elapsed_time(),
            "estimated_time_remaining": self.estimate_time_remaining(),
            "avg_time_per_batch": self.get_avg_time_per_batch(),
            "progress_percentage": self.get_progress_percentage(),
            "status": self.status,
            "recent_logs": self.recent_logs
        }

def get_db_connection() -> SupabaseDirectConnection:
    """
    Get and cache the Supabase direct connection.
    
    Returns:
        SupabaseDirectConnection: Database connection instance
    """
    global db
    
    if db is None:
        db = SupabaseDirectConnection.get_instance()
        logger.info("Using direct Supabase connection for database operations")
    
    return db

def fetch_compound_by_id(chembl_id):
    """
    Fetch a compound by ChEMBL ID.
    
    Args:
        chembl_id: ChEMBL ID to fetch
        
    Returns:
        dict: Compound data or None if not found
    """
    try:
        logger.info(f"Fetching standard reference compound: {chembl_id}")
        molecule = new_client.molecule
        compound = molecule.get(chembl_id)
        
        # Add properties to compound
        if 'molecule_properties' in compound:
            mol_props = compound.get('molecule_properties', {})
            properties = []
            for prop_key, prop_value in mol_props.items():
                if prop_value is not None:
                    properties.append({
                        'property_name': prop_key,
                        'value': prop_value
                    })
            compound['properties'] = properties
        
        return compound
    except Exception as e:
        log_error(
            error_type="API",
            message=f"Error fetching reference compound {chembl_id}",
            context={
                "exception": e,
                "compound_id": chembl_id,
                "source": "fetch_compound_by_id"
            }
        )
        return None

def insert_property_type(name: str, data_type: str = "text", description: str = None, units: str = None) -> Optional[str]:
    """
    Insert a new property type into the database using direct Supabase connection.
    
    Args:
        name: The name of the property type
        data_type: The data type of the property (numeric, text, boolean)
        description: Optional description of the property
        units: Optional units for the property
        
    Returns:
        str: Property type ID if inserted successfully, None otherwise
    """
    try:
        # Prepare property type data
        property_type_data = {
            "name": name,
            "data_type": data_type,
            "description": description or f"Auto-added by ChEMBL import on {datetime.now().isoformat()}",
            "units": units
        }
        
        # Get database connection
        db = get_db_connection()
        
        # Construct SQL for insertion with parameterized query
        columns = ", ".join(property_type_data.keys())
        placeholders = ", ".join([f"%({k})s" for k in property_type_data.keys()])
        
        sql = f"""
        INSERT INTO property_types ({columns})
        VALUES ({placeholders})
        RETURNING id;
        """
        
        # Execute SQL via direct connection
        result = db.execute_query(sql, property_type_data)
        
        # Extract data from result
        if result and len(result) > 0 and 'id' in result[0]:
            logger.info(f"Successfully inserted new property type: {name}")
            return result[0]['id']
        else:
            logger.error(f"Failed to insert property type: {name}")
            return None
            
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error inserting property type: {str(e)}",
            context={
                "exception": e,
                "property_name": name,
                "source": "insert_property_type"
            }
        )
        return None

def get_property_types() -> Dict[str, str]:
    """
    Get property types from the database using direct Supabase connection.
    
    Returns:
        Dict[str, str]: Dictionary mapping property names to property type IDs
    """
    try:
        # Get database connection
        db = get_db_connection()
        
        # Construct SQL query
        sql = "SELECT id, name FROM property_types;"
        
        # Execute SQL via direct connection
        result = db.execute_query(sql)
        
        # Extract data from result
        if result:
            return {pt['name'].lower(): pt['id'] for pt in result}
        else:
            logger.error("Failed to get property types: empty result")
            return {}
            
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error getting property types: {str(e)}",
            context={
                "exception": e,
                "source": "get_property_types"
            }
        )
        raise

def check_molecule_exists(inchikey: str) -> Optional[str]:
    """
    Check if a molecule with the given InChIKey already exists in the database using direct Supabase connection.
    
    Args:
        inchikey: InChIKey to check
        
    Returns:
        str: Molecule ID if exists, None otherwise
    """
    try:
        # Get database connection
        db = get_db_connection()
        
        # Construct SQL query with parameterized query
        sql = "SELECT id FROM molecules WHERE inchikey = %(inchikey)s;"
        
        # Execute SQL via direct connection
        result = db.execute_query(sql, {"inchikey": inchikey})
        
        # Extract data from result
        if result and len(result) > 0:
            return result[0]['id']
        else:
            return None
            
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error checking if molecule exists: {str(e)}",
            context={
                "exception": e,
                "inchikey": inchikey,
                "source": "check_molecule_exists"
            }
        )
        return None

def insert_molecule(molecule_data: Dict[str, Any]) -> Optional[str]:
    """
    Insert a molecule into the database using direct Supabase connection.
    
    Args:
        molecule_data: Molecule data to insert
        
    Returns:
        str: Molecule ID if inserted successfully, None otherwise
    """
    try:
        # Get database connection
        db = get_db_connection()
        
        # Filter out None values
        filtered_data = {k: v for k, v in molecule_data.items() if v is not None}
        
        # Construct SQL for insertion with parameterized query
        columns = list(filtered_data.keys())
        placeholders = [f"%({k})s" for k in columns]
        
        sql = f"""
        INSERT INTO molecules ({', '.join(columns)})
        VALUES ({', '.join(placeholders)})
        RETURNING id;
        """
        
        # Execute SQL via direct connection
        result = db.execute_query(sql, filtered_data)
        
        # Extract data from result
        if result and len(result) > 0 and 'id' in result[0]:
            return result[0]['id']
        else:
            logger.error(f"Failed to insert molecule: {molecule_data.get('name')}")
            return None
            
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error inserting molecule: {str(e)}",
            context={
                "exception": e,
                "molecule_name": molecule_data.get('name'),
                "source": "insert_molecule"
            }
        )
        return None

def insert_property(property_data: Dict[str, Any]) -> bool:
    """
    Insert a molecular property into the database using direct Supabase connection.
    
    Args:
        property_data: Property data to insert
        
    Returns:
        bool: True if inserted successfully, False otherwise
    """
    try:
        # Get database connection
        db = get_db_connection()
        
        # Filter out None values
        filtered_data = {k: v for k, v in property_data.items() if v is not None}
        
        # Construct SQL for insertion with parameterized query
        columns = list(filtered_data.keys())
        placeholders = [f"%({k})s" for k in columns]
        
        sql = f"""
        INSERT INTO molecular_properties ({', '.join(columns)})
        VALUES ({', '.join(placeholders)})
        RETURNING id;
        """
        
        # Execute SQL via direct connection
        result = db.execute_query(sql, filtered_data)
        
        # Check if successful
        return result and len(result) > 0 and 'id' in result[0]
            
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error inserting property: {str(e)}",
            context={
                "exception": e,
                "property_type": property_data.get('property_type_id'),
                "molecule_id": property_data.get('molecule_id'),
                "source": "insert_property"
            }
        )
        return False


class ChEMBLCheckpointManager:
    """
    Manages checkpoints for resumable ChEMBL data import.
    
    This class handles saving and loading checkpoint data, allowing the import process
    to be resumed after interruption. It provides robust error handling, checkpoint
    verification, and detailed status information.
    """
    
    def __init__(self, checkpoint_file=DEFAULT_CHECKPOINT_FILE):
        """
        Initialize the checkpoint manager.
        
        Args:
            checkpoint_file: Path to the checkpoint file
        """
        self.checkpoint_file = checkpoint_file
        self.backup_file = f"{checkpoint_file}.bak"
        self.last_save_time = 0
        self.save_count = 0
    
    def save(self, compounds, next_index, status="Running"):
        """
        Save the current import state to a checkpoint file.
        
        Args:
            compounds: List of compounds that have been processed
            next_index: Index to start from on next run
            status: Current status of the import
            
        Returns:
            True if checkpoint was saved successfully, False otherwise
        """
        # Create checkpoint data
        checkpoint_data = {
            "compounds": compounds,
            "next_index": next_index,
            "timestamp": datetime.now().isoformat(),
            "status": status,
            "version": "1.0"  # For future compatibility
        }
        
        try:
            # Backup existing checkpoint if it exists
            if os.path.exists(self.checkpoint_file):
                try:
                    shutil.copy2(self.checkpoint_file, self.backup_file)
                except Exception as e:
                    logger.warning(f"Failed to create backup of checkpoint file: {str(e)}")
            
            # Write new checkpoint
            with open(self.checkpoint_file, "w") as f:
                json.dump(checkpoint_data, f, indent=2)
            
            self.last_save_time = time.time()
            self.save_count += 1
            
            logger.info(f"Checkpoint saved: {len(compounds)} compounds processed, "
                      f"next index: {next_index}, status: {status}")
            return True
            
        except Exception as e:
            logger.error(f"Error saving checkpoint: {str(e)}")
            return False
    
    def load(self):
        """
        Load the last checkpoint if it exists.
        
        Returns:
            Dictionary with checkpoint data if found, None otherwise
        """
        # Check if checkpoint file exists
        if not os.path.exists(self.checkpoint_file):
            # Try backup file if main file doesn't exist
            if os.path.exists(self.backup_file):
                logger.warning("Main checkpoint file not found, trying backup...")
                return self._load_file(self.backup_file)
            
            logger.info("No checkpoint found. Starting from the beginning.")
            return None
        
        # Try to load the main checkpoint file
        checkpoint = self._load_file(self.checkpoint_file)
        
        # If main file is corrupted, try backup
        if checkpoint is None and os.path.exists(self.backup_file):
            logger.warning("Main checkpoint file corrupted, trying backup...")
            checkpoint = self._load_file(self.backup_file)
            
            # If backup loaded successfully, restore it as the main checkpoint
            if checkpoint is not None:
                try:
                    shutil.copy2(self.backup_file, self.checkpoint_file)
                    logger.info("Restored checkpoint from backup file")
                except Exception as e:
                    logger.warning(f"Failed to restore checkpoint from backup: {str(e)}")
        
        return checkpoint
    
    def _load_file(self, file_path):
        """
        Load checkpoint data from a file.
        
        Args:
            file_path: Path to the checkpoint file
            
        Returns:
            Dictionary with checkpoint data if loaded successfully, None otherwise
        """
        try:
            with open(file_path, "r") as f:
                checkpoint = json.load(f)
            
            # Validate checkpoint data
            required_keys = ["compounds", "next_index", "timestamp"]
            
            if not all(key in checkpoint for key in required_keys):
                logger.error(f"Checkpoint file {file_path} is missing required keys")
                return None
            
            # Log checkpoint details
            logger.info(f"Checkpoint loaded from {file_path}: "
                      f"{len(checkpoint.get('compounds', []))} compounds processed, "
                      f"next index: {checkpoint.get('next_index', 0)}, "
                      f"timestamp: {checkpoint.get('timestamp', 'Unknown')}")
            
            return checkpoint
            
        except json.JSONDecodeError as e:
            logger.error(f"Error parsing checkpoint file {file_path}: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error loading checkpoint from {file_path}: {str(e)}")
            return None


def fetch_cryoprotectant_compounds(limit: int = 1000, checkpoint_interval: int = 100) -> List[Dict[str, Any]]:
    """
    Fetch cryoprotectant compounds from ChEMBL.
    
    Args:
        limit: Maximum number of compounds to fetch
        checkpoint_interval: Save checkpoint after this many compounds
        
    Returns:
        List of compounds with properties
    """
    logger.info(f"Fetching up to {limit} cryoprotectant compounds from ChEMBL")
    
    # Initialize ChEMBL client
    molecule = new_client.molecule
    
    # Construct a query for potential cryoprotectants
    # Cryoprotectants often have specific properties:
    # - Multiple hydrogen bond donors/acceptors
    # - Moderate LogP values (water solubility)
    # - Specific molecular weight range
    
    # Start with a broad search for common cryoprotectants
    query_terms = [
        "cryoprotect",
        "glycerol",
        "dmso",
        "dimethyl sulfoxide",
        "ethylene glycol",
        "propylene glycol",
        "trehalose",
        "sucrose",
        "glucose",
        "formamide",
        "acetamide",
        "methanol",
        "polyvinyl alcohol"
    ]
    
    all_compounds = []
    checkpoint_file = Path(active_config.CHECKPOINT_DIR) / "chembl_integrated_checkpoint.json"
    
    # Initialize checkpoint manager
    checkpoint_manager = ChEMBLCheckpointManager(checkpoint_file)
    
    # Load checkpoint if it exists
    start_index = 0
    checkpoint = checkpoint_manager.load()
    if checkpoint:
        all_compounds = checkpoint.get("compounds", [])
        start_index = checkpoint.get("next_index", 0)
        logger.info(f"Loaded checkpoint with {len(all_compounds)} compounds, starting from index {start_index}")
    
    # First, fetch standard reference compounds
    reference_compounds = []
    for chembl_id in standard_reference_ids:
        compound = fetch_compound_by_id(chembl_id)
        if compound:
            reference_compounds.append(compound)
            # Add to all_compounds
            all_compounds.append(compound)
            logger.info(f"Added reference compound: {chembl_id}")
            # Slight delay to be gentle on the API
            time.sleep(0.2)
    
    logger.info(f"Fetched {len(reference_compounds)} reference compounds")
    
    # Initialize progress tracker for the query terms
    progress_tracker = ChEMBLProgressTracker(
        total_compounds=len(query_terms) - start_index,
        batch_size=1,  # Process one query term at a time
        checkpoint_file=str(checkpoint_file) + ".progress"
    )
    
    # Process each query term
    for i, term in enumerate(query_terms[start_index:], start=start_index):
        try:
            # Update progress tracker
            progress_tracker.start_batch(i - start_index, 1)
            logger.info(f"Searching for '{term}' ({i+1}/{len(query_terms)})")
            
            # Search for compounds matching the term
            results = molecule.filter(
                pref_name__icontains=term
            ).only(
                'molecule_chembl_id',
                'pref_name',
                'molecule_structures'
            )
            
            # Track statistics for this term
            term_processed = 0
            term_imported = 0
            term_skipped = 0
            term_errors = 0
            
            # Process results
            for compound in results:
                # Skip if we already have enough compounds
                if len(all_compounds) >= limit:
                    break
                    
                # Skip if no structures available
                if not compound.get('molecule_structures'):
                    log_skipped_molecule(
                        chembl_id=compound.get('molecule_chembl_id', 'unknown'),
                        reason="Missing molecular structures",
                        molecule_data={"name": compound.get('pref_name')},
                        category="validation"
                    )
                    term_skipped += 1
                    progress_tracker.add_skipped(
                        compound.get('molecule_chembl_id', 'unknown'),
                        "Missing molecular structures"
                    )
                    continue
                    
                # Get compound details
                compound_id = compound.get('molecule_chembl_id')
                try:
                    # Get full compound details
                    full_details = molecule.get(compound_id)
                    
                    # Properties are already included in the molecule_properties field
                    # We'll create a properties list for compatibility with the rest of the code
                    properties = []
                    if 'molecule_properties' in full_details:
                        mol_props = full_details.get('molecule_properties', {})
                        for prop_key, prop_value in mol_props.items():
                            if prop_value is not None:
                                properties.append({
                                    'property_name': prop_key,
                                    'value': prop_value
                                })
                    
                    # Add properties to compound
                    full_details['properties'] = properties
                    
                    # Add to results
                    all_compounds.append(full_details)
                    term_imported += 1
                    term_processed += 1
                    
                    # Log progress
                    if len(all_compounds) % 10 == 0:
                        logger.info(f"Fetched {len(all_compounds)}/{limit} compounds")
                    
                    # Save checkpoint if needed
                    if len(all_compounds) % checkpoint_interval == 0:
                        checkpoint_manager.save(all_compounds, i, "Running")
                        
                except Exception as e:
                    term_errors += 1
                    progress_tracker.add_error(f"Error fetching details for {compound_id}: {str(e)}")
                    log_error(
                        error_type="API",
                        message=f"Error fetching details for {compound_id}",
                        context={
                            "exception": e,
                            "compound_id": compound_id,
                            "source": "fetch_cryoprotectant_compounds.get_details"
                        }
                    )
                    continue
                    
                # Slight delay to be gentle on the API
                time.sleep(0.2)
            
            # Update progress tracker for this term
            progress_tracker.end_batch(term_processed, term_imported, term_skipped, term_errors)
            
            # Break if we have enough compounds
            if len(all_compounds) >= limit:
                logger.info(f"Reached limit of {limit} compounds.")
                break
                
        except Exception as e:
            progress_tracker.add_error(f"Error processing term '{term}': {str(e)}")
            log_error(
                error_type="API",
                message=f"Error processing term '{term}'",
                context={
                    "exception": e,
                    "term": term,
                    "source": "fetch_cryoprotectant_compounds.process_term"
                }
            )
            
            # Save checkpoint before continuing to next term
            checkpoint_manager.save(all_compounds, i, "Error")
            
            # Continue with next term
            continue
    
    # Save final checkpoint
    checkpoint_manager.save(all_compounds, len(query_terms), "Completed")
    
    # Update progress tracker status
    progress_tracker.set_status("Completed")
    progress_tracker.add_log_message(f"Import completed with {len(all_compounds)} compounds")
    
    # Log final statistics
    logger.info(f"ChEMBL import completed:")
    logger.info(f"  - Total compounds: {len(all_compounds)}")
    logger.info(f"  - Total processed: {progress_tracker.total_processed}")
    logger.info(f"  - Total imported: {progress_tracker.total_imported}")
    logger.info(f"  - Total skipped: {progress_tracker.total_skipped}")
    logger.info(f"  - Total errors: {progress_tracker.total_errors}")
    logger.info(f"  - Elapsed time: {progress_tracker.get_elapsed_time()}")
    
    return all_compounds


def transform_chembl_to_molecule(compound: Dict[str, Any], user_profile_id: str) -> Dict[str, Any]:
    """
    Transform ChEMBL compound data to match our molecule table schema.
    
    Args:
        compound: ChEMBL compound data
        user_profile_id: User profile ID for created_by field
        
    Returns:
        Dictionary matching our molecule table schema
    """
    # Extract structures
    structures = compound.get('molecule_structures', {}) or {}
    
    # Extract molecule properties
    mol_props = compound.get('molecule_properties', {}) or {}
    
    # Get ChEMBL version if available
    chembl_version = getattr(compound, 'chembl_version', None) or "Unknown"
    
    # Get the preferred name or a fallback
    name = compound.get('pref_name')
    if not name:
        # Try alternative name fields
        name = (compound.get('molecule_synonyms', [{}])[0].get('synonym') if compound.get('molecule_synonyms')
                else compound.get('molecule_chembl_id', 'Unknown Compound'))
    
    # Transform to our schema
    return {
        "name": name,
        "smiles": structures.get('canonical_smiles'),
        "inchi": structures.get('standard_inchi'),
        "inchikey": structures.get('standard_inchi_key'),
        "formula": mol_props.get('full_molformula'),
        "molecular_weight": mol_props.get('full_mwt'),
        "pubchem_cid": None,  # Added field for PubChem Compound ID (can be populated later if available)
        "chembl_id": compound.get('molecule_chembl_id'),  # Store ChEMBL ID in dedicated column
        "created_by": user_profile_id,
        "data_source": f"ChEMBL v{chembl_version} ID: {compound.get('molecule_chembl_id')}",
        "version": 1,
        "modification_history": json.dumps([{
            "timestamp": datetime.now().isoformat(),
            "action": "created",
            "user_id": user_profile_id
        }])
    }


def transform_chembl_to_properties(compound: Dict[str, Any], molecule_id: str, user_profile_id: str, property_type_map: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Transform ChEMBL compound properties to match our molecular_properties table schema.
    
    Args:
        compound: ChEMBL compound data
        molecule_id: Molecule ID from our database
        user_profile_id: User profile ID for created_by field
        property_type_map: Map of property names to property type IDs
        
    Returns:
        List of dictionaries matching our molecular_properties table schema
    """
    properties = []
    
    # Process standard molecule properties
    mol_props = compound.get('molecule_properties', {}) or {}
    
    # Enhanced property mappings with units and descriptions
    property_mappings = {
        'alogp': {
            'name': 'LogP',
            'description': 'Calculated octanol/water partition coefficient',
            'data_type': 'numeric',
            'units': None
        },
        'full_mwt': {
            'name': 'Molecular Weight',
            'description': 'Molecular weight of the compound',
            'data_type': 'numeric',
            'units': 'g/mol'
        },
        'hba': {
            'name': 'Hydrogen Bond Acceptor Count',
            'description': 'Number of hydrogen bond acceptors',
            'data_type': 'numeric',
            'units': None
        },
        'hbd': {
            'name': 'Hydrogen Bond Donor Count',
            'description': 'Number of hydrogen bond donors',
            'data_type': 'numeric',
            'units': None
        },
        'psa': {
            'name': 'Topological Polar Surface Area',
            'description': 'Topological polar surface area',
            'data_type': 'numeric',
            'units': 'Å²'
        },
        'rtb': {
            'name': 'Rotatable Bond Count',
            'description': 'Number of rotatable bonds',
            'data_type': 'numeric',
            'units': None
        },
        'cx_logp': {
            'name': 'LogP (ChemAxon)',
            'description': 'ChemAxon calculated octanol/water partition coefficient',
            'data_type': 'numeric',
            'units': None
        },
        'cx_logd': {
            'name': 'LogD',
            'description': 'Distribution coefficient at pH 7.4',
            'data_type': 'numeric',
            'units': None
        },
        'aromatic_rings': {
            'name': 'Aromatic Ring Count',
            'description': 'Number of aromatic rings',
            'data_type': 'numeric',
            'units': None
        },
        'heavy_atoms': {
            'name': 'Heavy Atom Count',
            'description': 'Number of non-hydrogen atoms',
            'data_type': 'numeric',
            'units': None
        },
        'num_ro5_violations': {
            'name': 'Rule of Five Violations',
            'description': 'Number of Lipinski Rule of Five violations',
            'data_type': 'numeric',
            'units': None
        },
        'qed_weighted': {
            'name': 'QED Weighted',
            'description': 'Weighted quantitative estimate of drug-likeness',
            'data_type': 'numeric',
            'units': None
        },
        'mw_freebase': {
            'name': 'Molecular Weight (Freebase)',
            'description': 'Molecular weight of the compound as freebase',
            'data_type': 'numeric',
            'units': 'g/mol'
        },
        'mw_monoisotopic': {
            'name': 'Molecular Weight (Monoisotopic)',
            'description': 'Monoisotopic mass of the compound',
            'data_type': 'numeric',
            'units': 'g/mol'
        },
        'full_molformula': {
            'name': 'Molecular Formula',
            'description': 'Molecular formula of the compound',
            'data_type': 'text',
            'units': None
        }
    }
    
    # Get the Unknown Property Type ID for fallback
    unknown_property_type_id = property_type_map.get('unknown property type')
    
    # Add properties from molecule_properties
    for prop_key, prop_info in property_mappings.items():
        if prop_key in mol_props and mol_props[prop_key] is not None:
            prop_name = prop_info['name']
            
            # Try to get property type ID
            property_type_id = property_type_map.get(prop_name.lower())
            
            # If property type doesn't exist, try to insert it
            if not property_type_id:
                logger.info(f"Property type '{prop_name}' not found in database. Attempting to insert it.")
                
                # Try to insert the new property type
                property_type_id = insert_property_type(
                    name=prop_name,
                    data_type=prop_info['data_type'],
                    description=prop_info['description'],
                    units=prop_info['units']
                )
                
                # Log the insertion attempt
                if property_type_id:
                    logger.info(f"Successfully inserted new property type: {prop_name}")
                    # Update the property_type_map with the new property type
                    property_type_map[prop_name.lower()] = property_type_id
                else:
                    logger.warning(f"Failed to insert property type: {prop_name}. Using fallback.")
                    # Use Unknown Property Type as fallback
                    if unknown_property_type_id:
                        property_type_id = unknown_property_type_id
                        logger.info(f"Using 'Unknown Property Type' as fallback for {prop_name}")
                    else:
                        logger.error(f"Unknown Property Type not found in database. Skipping property {prop_name}.")
                        log_error(
                            error_type="Property",
                            message=f"Unknown Property Type not found in database",
                            context={
                                "property_name": prop_name,
                                "molecule_id": molecule_id,
                                "source": "transform_chembl_to_properties"
                            }
                        )
                        continue
            
            # If we have a property type ID (either original, newly inserted, or fallback), add the property
            if property_type_id:
                # Determine value type and convert accordingly
                value = mol_props[prop_key]
                numeric_value = None
                text_value = None
                boolean_value = None
                
                if prop_info['data_type'] == 'numeric':
                    try:
                        numeric_value = float(value)
                    except (ValueError, TypeError):
                        # If conversion fails, store as text
                        text_value = str(value)
                elif prop_info['data_type'] == 'boolean':
                    if isinstance(value, bool):
                        boolean_value = value
                    else:
                        # Try to convert to boolean if possible
                        try:
                            boolean_value = bool(value)
                        except (ValueError, TypeError):
                            text_value = str(value)
                else:
                    # Default to text
                    text_value = str(value)
                
                properties.append({
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,
                    "numeric_value": numeric_value,
                    "text_value": text_value,
                    "boolean_value": boolean_value,
                    "created_by": user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_key}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": user_profile_id
                    }])
                })
    
    # Process additional properties from ChEMBL API
    # These might include bioactivity data, drug classifications, etc.
    additional_properties = []
    
    # Extract from compound.molecule_hierarchy if available
    if 'molecule_hierarchy' in compound and compound['molecule_hierarchy']:
        hierarchy = compound['molecule_hierarchy']
        if 'parent_chembl_id' in hierarchy and hierarchy['parent_chembl_id']:
            additional_properties.append({
                'property_name': 'Parent ChEMBL ID',
                'value': hierarchy['parent_chembl_id'],
                'data_type': 'text'
            })
    
    # Extract from compound.molecule_synonyms if available
    if 'molecule_synonyms' in compound and compound['molecule_synonyms']:
        for i, synonym in enumerate(compound['molecule_synonyms']):
            if i < 5:  # Limit to 5 synonyms to avoid excessive properties
                syn_type = synonym.get('syn_type', 'Synonym')
                additional_properties.append({
                    'property_name': f"{syn_type} Synonym",
                    'value': synonym.get('synonym', ''),
                    'data_type': 'text'
                })
    
    # Extract from compound.biotherapeutic if available
    if 'biotherapeutic' in compound and compound['biotherapeutic']:
        bio = compound['biotherapeutic']
        if 'description' in bio and bio['description']:
            additional_properties.append({
                'property_name': 'Biotherapeutic Description',
                'value': bio['description'],
                'data_type': 'text'
            })
    
    # Process additional properties
    for prop in additional_properties:
        prop_name = prop['property_name']
        if prop_name and prop['value'] is not None:
            # Try to get property type ID
            property_type_id = property_type_map.get(prop_name.lower())
            
            # If property type doesn't exist, try to insert it
            if not property_type_id:
                logger.info(f"Property type '{prop_name}' not found in database. Attempting to insert it.")
                
                # Try to insert the new property type
                property_type_id = insert_property_type(
                    name=prop_name,
                    data_type=prop.get('data_type', 'text'),
                    description=f"Auto-added by ChEMBL import - {prop_name}"
                )
                
                # Log the insertion attempt
                if property_type_id:
                    logger.info(f"Successfully inserted new property type: {prop_name}")
                    # Update the property_type_map with the new property type
                    property_type_map[prop_name.lower()] = property_type_id
                else:
                    logger.warning(f"Failed to insert property type: {prop_name}. Using fallback.")
                    # Use Unknown Property Type as fallback
                    if unknown_property_type_id:
                        property_type_id = unknown_property_type_id
                        logger.info(f"Using 'Unknown Property Type' as fallback for {prop_name}")
                    else:
                        logger.error(f"Unknown Property Type not found in database. Skipping property {prop_name}.")
                        continue
            
            # If we have a property type ID (either original, newly inserted, or fallback), add the property
            if property_type_id:
                # Determine value type
                value = prop['value']
                numeric_value = None
                text_value = None
                boolean_value = None
                
                if prop.get('data_type') == 'numeric':
                    try:
                        numeric_value = float(value)
                    except (ValueError, TypeError):
                        text_value = str(value)
                elif prop.get('data_type') == 'boolean':
                    if isinstance(value, bool):
                        boolean_value = value
                    else:
                        try:
                            boolean_value = bool(value)
                        except (ValueError, TypeError):
                            text_value = str(value)
                else:
                    text_value = str(value)
                
                properties.append({
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,
                    "numeric_value": numeric_value,
                    "text_value": text_value,
                    "boolean_value": boolean_value,
                    "created_by": user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_name}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": user_profile_id
                    }])
                })
    
    # Process any properties from the compound.properties list (if present)
    for prop in compound.get('properties', []):
        prop_name = prop.get('property_name')
        if prop_name and prop.get('value') is not None:
            # Try to get property type ID
            property_type_id = property_type_map.get(prop_name.lower())
            
            # If property type doesn't exist, try to insert it
            if not property_type_id:
                logger.info(f"Property type '{prop_name}' not found in database. Attempting to insert it.")
                
                # Determine value type
                value = prop.get('value')
                data_type = "text"
                try:
                    float(value)
                    data_type = "numeric"
                except (ValueError, TypeError):
                    if isinstance(value, bool):
                        data_type = "boolean"
                
                # Try to insert the new property type
                property_type_id = insert_property_type(
                    name=prop_name,
                    data_type=data_type,
                    description=f"Auto-added by ChEMBL import"
                )
                
                # Log the insertion attempt
                if property_type_id:
                    logger.info(f"Successfully inserted new property type: {prop_name}")
                    # Update the property_type_map with the new property type
                    property_type_map[prop_name.lower()] = property_type_id
                else:
                    logger.warning(f"Failed to insert property type: {prop_name}. Using fallback.")
                    # Use Unknown Property Type as fallback
                    if unknown_property_type_id:
                        property_type_id = unknown_property_type_id
                        logger.info(f"Using 'Unknown Property Type' as fallback for {prop_name}")
                    else:
                        logger.error(f"Unknown Property Type not found in database. Skipping property {prop_name}.")
                        continue
            
            # If we have a property type ID (either original, newly inserted, or fallback), add the property
            if property_type_id:
                # Determine value type
                value = prop.get('value')
                numeric_value = None
                text_value = None
                boolean_value = None
                
                try:
                    numeric_value = float(value)
                except (ValueError, TypeError):
                    if isinstance(value, bool):
                        boolean_value = value
                    else:
                        text_value = str(value)
                
                properties.append({
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,
                    "numeric_value": numeric_value,
                    "text_value": text_value,
                    "boolean_value": boolean_value,
                    "created_by": user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_name}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": user_profile_id
                    }])
                })
    
    return properties


@ensure_rls_restored
def import_compounds_to_database(compounds: List[Dict[str, Any]], batch_size: int = 50, dry_run: bool = False,
                                checkpoint_file: str = None) -> Dict[str, int]:
    """
    Import compounds to database with proper error handling and batch processing.
    
    Args:
        compounds: List of ChEMBL compounds
        batch_size: Number of compounds to process in each batch
        dry_run: If True, don't actually insert data, just simulate
        
    Returns:
        Dictionary with import statistics
    """
    logger.info(f"Importing {len(compounds)} compounds to database")
    
    # Get authenticated user
    user_id = get_user_id()
    user_profile_id = ensure_user_profile()
    
    if not user_profile_id:
        log_error(
            error_type="Authentication",
            message="Failed to get user profile ID",
            context={
                "user_id": user_id,
                "source": "import_compounds_to_database.get_user"
            }
        )
        return {"error": "Failed to get user profile ID"}
    
    # Get property type map
    property_type_map = {}
    
    # Skip database operations in dry run mode
    if not dry_run:
        try:
            property_type_map = get_property_types()
        except Exception as e:
            log_error(
                error_type="Database",
                message="Error getting property types",
                context={
                    "exception": e,
                    "source": "import_compounds_to_database.get_property_types"
                }
            )
            return {"error": f"Failed to get property types: {str(e)}"}
    
    # Prepare statistics
    stats = {
        "total_compounds": len(compounds),
        "processed": 0,
        "molecules_inserted": 0,
        "molecules_skipped": 0,
        "properties_inserted": 0,
        "errors": 0,
        "batches_processed": 0,
        "batches_failed": 0
    }
    
    # Initialize progress tracker
    if not checkpoint_file:
        checkpoint_file = Path(active_config.CHECKPOINT_DIR) / "chembl_import_progress.json"
    
    progress_tracker = ChEMBLProgressTracker(
        total_compounds=len(compounds),
        batch_size=batch_size,
        checkpoint_file=checkpoint_file
    )
    
    if dry_run:
        logger.info("DRY RUN MODE: No data will be inserted into the database")
        
        # Create a mock property type map for dry run
        mock_property_types = {
            "logp": "mock-property-type-logp",
            "molecular weight": "mock-property-type-mw",
            "hydrogen bond acceptor count": "mock-property-type-hba",
            "hydrogen bond donor count": "mock-property-type-hbd",
            "topological polar surface area": "mock-property-type-psa",
            "rotatable bond count": "mock-property-type-rtb",
            "logd": "mock-property-type-logd",
            "aromatic ring count": "mock-property-type-arc",
            "heavy atom count": "mock-property-type-hac",
            "rule of five violations": "mock-property-type-ro5",
            "unknown property type": "mock-property-type-unknown"  # Add Unknown Property Type for fallback
        }
        
        # Track dynamic property type insertions and fallback usages for dry run reporting
        dynamic_property_insertions = []
        fallback_usages = []
        
        # Process compounds to validate transformation
        for i, compound in enumerate(compounds[:min(10, len(compounds))]):
            try:
                # Transform compound to molecule
                molecule_data = transform_chembl_to_molecule(compound, user_profile_id)
                logger.info(f"DRY RUN: Would insert molecule: {molecule_data['name']} ({molecule_data['inchikey']})")
                
                # Simulate property type fallback logic for dry run
                # Add some unknown property types to test fallback logic
                test_unknown_props = {
                    'test_unknown_prop_1': 'numeric',
                    'test_unknown_prop_2': 'text'
                }
                
                # Add unknown properties to the compound for testing fallback
                if 'properties' not in compound:
                    compound['properties'] = []
                
                for prop_name, data_type in test_unknown_props.items():
                    # Add a test property that's not in mock_property_types
                    test_value = 123.45 if data_type == 'numeric' else 'test value'
                    compound['properties'].append({
                        'property_name': prop_name,
                        'value': test_value
                    })
                    
                    # Simulate dynamic insertion (50% success rate for testing)
                    if len(dynamic_property_insertions) % 2 == 0:
                        # Simulate successful insertion
                        mock_property_types[prop_name.lower()] = f"mock-property-type-{prop_name.lower()}"
                        dynamic_property_insertions.append({
                            'name': prop_name,
                            'data_type': data_type,
                            'status': 'success'
                        })
                    else:
                        # Simulate failed insertion (will use fallback)
                        fallback_usages.append({
                            'name': prop_name,
                            'fallback_type': 'unknown property type'
                        })
                
                # Transform properties
                properties = transform_chembl_to_properties(
                    compound, "dry-run-molecule-id", user_profile_id, mock_property_types
                )
                
                # Log detailed property information
                logger.info(f"DRY RUN: Would insert {len(properties)} properties for molecule {molecule_data['name']}")
                for prop in properties[:5]:  # Show first 5 properties only to avoid log clutter
                    prop_type = next((k for k, v in mock_property_types.items() if v == prop['property_type_id']), "unknown")
                    value = prop['numeric_value'] or prop['text_value'] or prop['boolean_value']
                    logger.info(f"  - Property: {prop_type}, Value: {value}")
                
                if len(properties) > 5:
                    logger.info(f"  - ... and {len(properties) - 5} more properties")
                
                stats["processed"] += 1
                stats["molecules_inserted"] += 1
                stats["properties_inserted"] += len(properties)
            except Exception as e:
                log_error(
                    error_type="Transformation",
                    message=f"Error transforming compound in dry run",
                    context={
                        "exception": e,
                        "compound_id": compound.get('molecule_chembl_id', 'unknown'),
                        "source": "import_compounds_to_database.dry_run"
                    }
                )
                stats["errors"] += 1
        
        # Log summary of what would happen in a real run
        logger.info("\nDRY RUN SUMMARY:")
        logger.info(f"Total compounds processed: {stats['processed']}")
        logger.info(f"Molecules that would be inserted: {stats['molecules_inserted']}")
        logger.info(f"Properties that would be inserted: {stats['properties_inserted']}")
        logger.info(f"Errors encountered: {stats['errors']}")
        
        # Log property type fallback information
        logger.info("\nProperty Type Fallback Logic:")
        logger.info(f"Dynamic property type insertions: {len(dynamic_property_insertions)}")
        for insertion in dynamic_property_insertions:
            logger.info(f"  - {insertion['name']} ({insertion['data_type']}): {insertion['status']}")
        
        logger.info(f"Fallback to 'Unknown Property Type': {len(fallback_usages)}")
        for fallback in fallback_usages:
            logger.info(f"  - {fallback['name']} → {fallback['fallback_type']}")
        
        logger.info("\nDatabase operations that would be performed:")
        logger.info("1. For each batch of compounds:")
        logger.info("   a. Check which molecules already exist by InChIKey")
        logger.info("   b. Insert new molecules in a batch transaction")
        logger.info("   c. Insert properties for all molecules in a batch transaction")
        logger.info("   d. For unknown property types:")
        logger.info("      i. Attempt to insert new property type")
        logger.info("      ii. If insertion fails, use 'Unknown Property Type' as fallback")
        
        return stats
    
    try:
        # Get database connection
        db = get_db_connection()
        
        # Track previous values for progress tracking
        previous_molecules_inserted = 0
        
        # Process compounds in batches
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            batch_num = i//batch_size + 1
            total_batches = (len(compounds)-1)//batch_size + 1
            
            # Update progress tracker
            progress_tracker.start_batch(batch_num-1, len(batch))
            logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} compounds)")
            
            # Lists to store batch data
            batch_molecules = []
            batch_molecule_data = []
            batch_properties = []
            batch_skipped = []
            
            # Track batch statistics
            batch_processed = 0
            batch_imported = 0
            batch_skipped_count = 0
            batch_errors = 0
            batch_errors = []
            
            # First pass: transform compounds and check which molecules already exist
            existing_molecules = {}
            new_molecules = []
            
            # Collect all InChIKeys for batch lookup
            inchikeys_to_check = []
            inchikey_to_molecule_data = {}
            inchikey_to_compound = {}
            
            for compound in batch:
                try:
                    # Transform compound to molecule
                    molecule_data = transform_chembl_to_molecule(compound, user_profile_id)
                    
                    # Check if molecule has InChIKey
                    inchikey = molecule_data.get('inchikey')
                    if not inchikey:
                        log_skipped_molecule(
                            chembl_id=compound.get('molecule_chembl_id', 'unknown'),
                            reason="Missing InChIKey",
                            molecule_data={"name": molecule_data.get('name')},
                            category="validation"
                        )
                        batch_skipped.append({
                            "chembl_id": compound.get('molecule_chembl_id', 'unknown'),
                            "reason": "Missing InChIKey"
                        })
                        stats["molecules_skipped"] += 1
                        continue
                    
                    # Store for batch lookup
                    inchikeys_to_check.append(inchikey)
                    inchikey_to_molecule_data[inchikey] = molecule_data
                    inchikey_to_compound[inchikey] = compound
                    
                except Exception as e:
                    log_error(
                        error_type="Transformation",
                        message=f"Error transforming compound",
                        context={
                            "exception": e,
                            "compound_id": compound.get('molecule_chembl_id', 'unknown'),
                            "source": "import_compounds_to_database.transform"
                        }
                    )
                    batch_errors.append({
                        "chembl_id": compound.get('molecule_chembl_id', 'unknown'),
                        "error": str(e),
                        "stage": "transform"
                    })
                    stats["errors"] += 1
            
            # Skip empty batches
            if not inchikeys_to_check:
                logger.info(f"Batch {batch_num} has no valid compounds to process, skipping")
                continue
            
            # Batch check which molecules already exist
            try:
                # Construct SQL for batch lookup
                placeholders = ', '.join([f"'{inchikey}'" for inchikey in inchikeys_to_check])
                sql = f"SELECT id, inchikey FROM molecules WHERE inchikey IN ({placeholders});"
                
                # Execute SQL via direct connection
                result = db.execute_query(sql)
                
                # Process results
                if result:
                    for row in result:
                        existing_molecules[row['inchikey']] = row['id']
                
                logger.debug(f"Found {len(existing_molecules)} existing molecules in batch {batch_num}")
                
                # Determine which molecules need to be inserted
                for inchikey in inchikeys_to_check:
                    if inchikey in existing_molecules:
                        # Molecule already exists
                        stats["molecules_skipped"] += 1
                    else:
                        # New molecule to insert
                        new_molecules.append(inchikey)
                
                logger.debug(f"Will insert {len(new_molecules)} new molecules in batch {batch_num}")
                
            except Exception as e:
                log_error(
                    error_type="Database",
                    message=f"Error checking existing molecules for batch {batch_num}",
                    context={
                        "exception": e,
                        "source": "import_compounds_to_database.check_existing"
                    }
                )
                batch_errors.append({
                    "error": str(e),
                    "stage": "check_existing",
                    "batch": batch_num
                })
                stats["errors"] += 1
                stats["batches_failed"] += 1
                continue
            
            # Prepare batch insert for new molecules
            if new_molecules:
                try:
                    # Construct SQL for batch insert
                    molecule_insert_queries = []
                    molecule_ids = {}
                    
                    for inchikey in new_molecules:
                        molecule_data = inchikey_to_molecule_data[inchikey]
                        
                        # Filter out None values
                        filtered_data = {k: v for k, v in molecule_data.items() if v is not None}
                        
                        # Generate a UUID for the molecule
                        molecule_id = str(uuid.uuid4())
                        filtered_data['id'] = molecule_id
                        molecule_ids[inchikey] = molecule_id
                        
                        # Construct SQL for insertion
                        columns = list(filtered_data.keys())
                        values = []
                        
                        for k in columns:
                            v = filtered_data[k]
                            if isinstance(v, str):
                                # Escape single quotes in strings
                                escaped_value = v.replace("'", "''")
                                values.append(f"'{escaped_value}'")
                            elif v is None:
                                values.append("NULL")
                            else:
                                values.append(str(v))
                        
                        sql = f"""
                        INSERT INTO molecules ({', '.join(columns)})
                        VALUES ({', '.join(values)});
                        """
                        
                        molecule_insert_queries.append(sql)
                        batch_molecule_data.append(filtered_data)
                    
                    # Execute batch insert within a transaction
                    if molecule_insert_queries:
                        success = db.execute_batch(molecule_insert_queries, transaction=True)
                        
                        if success:
                            logger.info(f"Successfully inserted {len(molecule_insert_queries)} molecules in batch {batch_num}")
                            stats["molecules_inserted"] += len(molecule_insert_queries)
                            
                            # Update existing_molecules with newly inserted molecules
                            for inchikey, molecule_id in molecule_ids.items():
                                existing_molecules[inchikey] = molecule_id
                        else:
                            log_error(
                                error_type="Database",
                                message=f"Failed to insert molecules for batch {batch_num}",
                                context={
                                    "source": "import_compounds_to_database.insert_molecules"
                                }
                            )
                            batch_errors.append({
                                "error": "Failed to insert molecules",
                                "stage": "insert_molecules",
                                "batch": batch_num
                            })
                            stats["errors"] += len(molecule_insert_queries)
                            stats["batches_failed"] += 1
                            continue
                
                except Exception as e:
                    log_error(
                        error_type="Database",
                        message=f"Error inserting molecules for batch {batch_num}",
                        context={
                            "exception": e,
                            "source": "import_compounds_to_database.insert_molecules"
                        }
                    )
                    batch_errors.append({
                        "error": str(e),
                        "stage": "insert_molecules",
                        "batch": batch_num
                    })
                    stats["errors"] += 1
                    stats["batches_failed"] += 1
                    continue
            
            # Now process properties for all molecules in the batch
            try:
                property_insert_queries = []
                
                for inchikey in inchikeys_to_check:
                    if inchikey in existing_molecules:
                        molecule_id = existing_molecules[inchikey]
                        compound = inchikey_to_compound[inchikey]
                        
                        # Transform properties
                        properties = transform_chembl_to_properties(
                            compound, molecule_id, user_profile_id, property_type_map
                        )
                        
                        # Prepare property insert queries
                        for prop in properties:
                            # Filter out None values
                            filtered_data = {k: v for k, v in prop.items() if v is not None}
                            
                            # Construct SQL for insertion
                            columns = list(filtered_data.keys())
                            values = []
                            
                            for k in columns:
                                v = filtered_data[k]
                                if isinstance(v, str):
                                    # Escape single quotes in strings
                                    escaped_value = v.replace("'", "''")
                                    values.append(f"'{escaped_value}'")
                                elif v is None:
                                    values.append("NULL")
                                else:
                                    values.append(str(v))
                            
                            sql = f"""
                            INSERT INTO molecular_properties ({', '.join(columns)})
                            VALUES ({', '.join(values)});
                            """
                            
                            property_insert_queries.append(sql)
                            batch_properties.append(filtered_data)
                
                # Execute batch insert for properties within a transaction
                if property_insert_queries:
                    success = db.execute_batch(property_insert_queries, transaction=True)
                    
                    if success:
                        logger.info(f"Successfully inserted {len(property_insert_queries)} properties in batch {batch_num}")
                        stats["properties_inserted"] += len(property_insert_queries)
                    else:
                        log_error(
                            error_type="Database",
                            message=f"Failed to insert properties for batch {batch_num}",
                            context={
                                "source": "import_compounds_to_database.insert_properties"
                            }
                        )
                        batch_errors.append({
                            "error": "Failed to insert properties",
                            "stage": "insert_properties",
                            "batch": batch_num
                        })
                        stats["errors"] += len(property_insert_queries)
                        stats["batches_failed"] += 1
                        continue
            
            except Exception as e:
                log_error(
                    error_type="Database",
                    message=f"Error inserting properties for batch {batch_num}",
                    context={
                        "exception": e,
                        "source": "import_compounds_to_database.insert_properties"
                    }
                )
                batch_errors.append({
                    "error": str(e),
                    "stage": "insert_properties",
                    "batch": batch_num
                })
                stats["errors"] += 1
                stats["batches_failed"] += 1
                continue
            
            # Update statistics
            stats["processed"] += len(batch)
            stats["batches_processed"] += 1
            batch_processed = len(batch)
            batch_imported = stats["molecules_inserted"] - previous_molecules_inserted
            batch_skipped_count = len(batch_skipped)
            batch_errors_count = len(batch_errors)
            
            # Update progress tracker
            progress_tracker.end_batch(
                processed=batch_processed,
                imported=batch_imported,
                skipped=batch_skipped_count,
                errors=batch_errors_count
            )
            
            logger.info(f"Successfully processed batch {batch_num}/{total_batches}")
            
            # Log any skipped items or errors
            if batch_skipped:
                logger.debug(f"Skipped {len(batch_skipped)} compounds in batch {batch_num}")
                for skipped in batch_skipped[:5]:  # Log first 5 only
                    logger.debug(f"  - Skipped {skipped['chembl_id']}: {skipped['reason']}")
                    progress_tracker.add_skipped(skipped['chembl_id'], skipped['reason'])
            
            if batch_errors:
                logger.warning(f"Encountered {len(batch_errors)} errors in batch {batch_num}")
                for error in batch_errors[:5]:  # Log first 5 only
                    progress_tracker.add_error(f"Error in batch {batch_num}: {error['error']} ({error['stage']})")
                    logger.warning(f"  - Error in {error.get('stage', 'unknown')}: {error.get('error', 'unknown')}")
            
            # Sleep to avoid overloading the database
            time.sleep(1)
            
            # Update previous values for next batch
            previous_molecules_inserted = stats["molecules_inserted"]
    
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error processing compounds",
            context={
                "exception": e,
                "source": "import_compounds_to_database"
            }
        )
        stats["errors"] += 1
        progress_tracker.add_error(f"Fatal error: {str(e)}")
        progress_tracker.set_status("Error")
    else:
        # Set progress tracker status to completed
        progress_tracker.set_status("Completed")
        
        # Log final statistics
        logger.info(f"ChEMBL import completed:")
        logger.info(f"  - Total compounds: {stats['total_compounds']}")
        logger.info(f"  - Total processed: {stats['processed']}")
        logger.info(f"  - Total molecules inserted: {stats['molecules_inserted']}")
        logger.info(f"  - Total properties inserted: {stats['properties_inserted']}")
        logger.info(f"  - Total skipped: {stats['molecules_skipped']}")
        logger.info(f"  - Total errors: {stats['errors']}")
        logger.info(f"  - Elapsed time: {progress_tracker.get_elapsed_time()}")
    
    # Add progress data to stats
    stats["progress_data"] = progress_tracker.get_progress_data()
    
    return stats


def main(args=None):
    """Main entry point."""
    if args is None:
        parser = argparse.ArgumentParser(description="Import data from ChEMBL to CryoProtect")
        parser.add_argument("--limit", type=int, default=2500, help="Maximum number of compounds to import")
        parser.add_argument("--batch-size", type=int, default=50, help="Batch size for database operations")
        parser.add_argument("--checkpoint-interval", type=int, default=100, help="Interval for saving checkpoints")
        parser.add_argument("--dry-run", action="store_true", help="Don't actually insert data, just simulate")
        parser.add_argument("--apply-fixes", action="store_true", help="Apply fixes for properties and cross-references after import")
        parser.add_argument("--fix-batch-size", type=int, default=50, help="Batch size for fix operations")
        parser.add_argument("--fix-limit", type=int, default=None, help="Maximum number of molecules to fix")
        args = parser.parse_args()

    start_time = time.time()

    try:
        # Validate configuration
        logger.info("Validating configuration...")
        validate_config()

        # Fetch compounds from ChEMBL
        logger.info(f"Fetching compounds from ChEMBL (limit: {args.limit})...")
        compounds = fetch_cryoprotectant_compounds(limit=args.limit, checkpoint_interval=args.checkpoint_interval)

        # Import to database if not a dry run
        if not args.dry_run:
            logger.info(f"Importing {len(compounds)} compounds to database (batch size: {args.batch_size})...")
            checkpoint_file = Path(active_config.CHECKPOINT_DIR) / "chembl_db_import_progress.json"
            stats = import_compounds_to_database(
                compounds,
                batch_size=args.batch_size,
                dry_run=False,
                checkpoint_file=str(checkpoint_file)
            )
            logger.info(f"Import statistics: {json.dumps(stats, indent=2)}")
        else:
            logger.info(f"Dry run completed. {len(compounds)} compounds would be imported.")
            stats = import_compounds_to_database(compounds, batch_size=args.batch_size, dry_run=True)

        # Apply fixes if requested
        if args.apply_fixes and not args.dry_run:
            try:
                # Import the enhanced property calculator and ChEMBL to PubChem resolver
                from integrated_chembl_import_fix import process_import_results

                logger.info("Applying fixes for properties and cross-references...")

                # Process the import results
                fix_stats = process_import_results(
                    stats,
                    fix_immediately=True,
                    batch_size=args.fix_batch_size,
                    limit=args.fix_limit
                )

                # Add fix stats to the main stats
                stats["enhancements"] = fix_stats["enhancements"]

                logger.info("Fixes applied successfully.")
            except ImportError as e:
                logger.warning(f"Could not import fixes module: {str(e)}")
                logger.warning("Skipping fixes. Run integrated_chembl_import_fix.py manually if needed.")
            except Exception as e:
                logger.error(f"Error applying fixes: {str(e)}")
                logger.warning("Import was successful, but fixes failed. Run integrated_chembl_import_fix.py manually if needed.")

        # Generate summary report
        end_time = time.time()
        summary = {
            "start_time": datetime.fromtimestamp(start_time).isoformat(),
            "end_time": datetime.fromtimestamp(end_time).isoformat(),
            "duration_seconds": round(end_time - start_time, 2),
            "compounds": {
                "total_fetched": len(compounds),
                "total_processed": stats.get("processed", 0),
                "inserted": stats.get("molecules_inserted", 0),
                "skipped": stats.get("molecules_skipped", 0)
            },
            "properties": {
                "total_inserted": stats.get("properties_inserted", 0)
            },
            "errors": stats.get("errors", 0),
            "args": vars(args)
        }

        # Add enhancement statistics if available
        if "enhancements" in stats:
            summary["enhancements"] = stats["enhancements"]

        # Write summary to file
        write_summary(summary)

        logger.info("Done.")

        # Return stats for integration with other scripts
        return stats

    except ConfigurationError as e:
        logger.error(f"Configuration error: {str(e)}")
        return 1
    except Exception as e:
        log_error(
            error_type="Application",
            message=f"Unhandled error in main: {str(e)}",
            context={
                "exception": e,
                "traceback": traceback.format_exc(),
                "source": "main"
            }
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
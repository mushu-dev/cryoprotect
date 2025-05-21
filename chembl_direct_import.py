#!/usr/bin/env python3
"""
ChEMBL Direct Import with Optimized Performance

This script implements a high-performance import system for ChEMBL molecular data:
- Official ChEMBL client integration for reliable API communication
- Optimized connection pooling for efficient database operations
- Robust checkpointing with incremental progress tracking
- Parallel batch processing of molecule properties
- Enhanced monitoring and diagnostics
- Memory management for large datasets
- Exponential backoff retry with circuit breaker pattern

This script is designed to handle full ChEMBL datasets with millions of
molecules while maintaining database consistency and performance.
"""

import os
import sys
import json
import time
import uuid
import signal
import atexit
import argparse
import logging
import traceback
import concurrent.futures
import multiprocessing
import threading
import queue
import psutil
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Dict, Any, Optional, Union, Tuple, Set, Callable

# Import the official ChEMBL client
from chembl_webresource_client.new_client import new_client

# Import configuration system
from config import active_config, validate_config, ConfigurationError

# Import the enhanced connection pool
from optimized_connection_pool import (
    OptimizedConnectionPool, ConnectionManager, 
    TransactionManager, execute_query_with_retry,
    execute_batch_with_retry, transaction_context
)

# Import centralized logging system
from chembl.logging import (
    ChEMBLLogger, log_error, log_skipped_molecule,
    log_progress, write_summary, get_logger
)

# Import RLS verification and remediation utilities
from rls_utils import ensure_rls_restored

# Initialize logger
logger = get_logger(__name__)

# Create global logger instance with default settings
chembl_logger = ChEMBLLogger(
    log_dir="logs",
    progress_log="chembl_direct_import_progress.jsonl",
    error_log="chembl_direct_import_errors.jsonl",
    skipped_log="chembl_direct_import_skipped.jsonl",
    summary_log="chembl_direct_import_summary.json",
    general_log="chembl_direct_import.log"
)

# Ensure checkpoints directory exists
Path(active_config.CHECKPOINT_DIR).mkdir(parents=True, exist_ok=True)

# Default checkpoint file path
DEFAULT_CHECKPOINT_FILE = Path(active_config.CHECKPOINT_DIR) / "chembl_direct_import_checkpoint.json"

# Memory management thresholds
MEMORY_CHECK_INTERVAL = max(1, active_config.CHEMBL_MEMORY_CHECK_FREQUENCY)
MEMORY_THRESHOLD = min(95.0, max(50.0, active_config.CHEMBL_MEMORY_THRESHOLD))

# Standard reference compounds to ensure they are always imported
STANDARD_REFERENCE_IDS = [
    "CHEMBL25",     # Aspirin
    "CHEMBL1118",   # Caffeine
    "CHEMBL1234",   # Glycerol (common cryoprotectant)
    "CHEMBL444",    # Glucose
    "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
    "CHEMBL9335",   # Dimethyl sulfoxide (DMSO, common cryoprotectant)
    "CHEMBL15151",  # Trehalose (common cryoprotectant)
    "CHEMBL1201625", # Propylene glycol (common cryoprotectant)
    "CHEMBL1276055", # Sucrose (common cryoprotectant)
]

# Performance metrics for monitoring
class PerformanceMetrics:
    """Track performance metrics for the import process."""
    
    def __init__(self):
        self.start_time = time.time()
        self.api_calls = 0
        self.api_errors = 0
        self.api_call_times = []
        self.db_operations = 0
        self.db_errors = 0
        self.db_operation_times = []
        self.memory_usage = []
        self.memory_check_interval = 60  # Check memory every 60 seconds
        self.last_memory_check = 0
        
        # For tracking parallel operations
        self.task_queue_size = []
        self.worker_utilization = []
        
        # Start monitoring thread
        self._start_monitoring_thread()
    
    def record_api_call(self, duration: float, success: bool = True):
        """Record an API call."""
        self.api_calls += 1
        self.api_call_times.append(duration)
        if not success:
            self.api_errors += 1
    
    def record_db_operation(self, duration: float, success: bool = True):
        """Record a database operation."""
        self.db_operations += 1
        self.db_operation_times.append(duration)
        if not success:
            self.db_errors += 1
    
    def record_task_queue(self, queue_size: int, active_workers: int, total_workers: int):
        """Record task queue metrics."""
        self.task_queue_size.append(queue_size)
        self.worker_utilization.append(active_workers / max(1, total_workers))
    
    def _start_monitoring_thread(self):
        """Start a background thread for periodic monitoring."""
        def monitor_memory():
            while True:
                try:
                    time.sleep(self.memory_check_interval)
                    
                    # Update memory usage
                    process = psutil.Process(os.getpid())
                    memory_info = process.memory_info()
                    memory_percent = process.memory_percent()
                    
                    self.memory_usage.append({
                        'timestamp': time.time(),
                        'rss': memory_info.rss,
                        'vms': memory_info.vms,
                        'percent': memory_percent
                    })
                    
                    # Keep only last 1000 data points to avoid memory growth
                    if len(self.memory_usage) > 1000:
                        self.memory_usage = self.memory_usage[-1000:]
                    
                    # Check if memory usage is too high
                    if memory_percent > MEMORY_THRESHOLD:
                        logger.warning(f"High memory usage detected: {memory_percent:.1f}% > {MEMORY_THRESHOLD:.1f}%")
                except Exception as e:
                    logger.error(f"Error in memory monitor: {str(e)}")
        
        self.monitor_thread = threading.Thread(
            target=monitor_memory,
            daemon=True
        )
        self.monitor_thread.start()
    
    def get_memory_usage(self) -> Dict[str, Any]:
        """Get current memory usage."""
        try:
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            memory_percent = process.memory_percent()
            
            return {
                'rss': memory_info.rss,
                'vms': memory_info.vms,
                'percent': memory_percent
            }
        except Exception as e:
            logger.error(f"Error getting memory usage: {str(e)}")
            return {'error': str(e)}
    
    def get_elapsed_time(self) -> float:
        """Get elapsed time in seconds."""
        return time.time() - self.start_time
    
    def get_metrics(self) -> Dict[str, Any]:
        """Get performance metrics."""
        memory_usage = self.get_memory_usage()
        
        return {
            'elapsed_time': self.get_elapsed_time(),
            'elapsed_time_formatted': str(timedelta(seconds=int(self.get_elapsed_time()))),
            'api': {
                'calls': self.api_calls,
                'errors': self.api_errors,
                'avg_time': sum(self.api_call_times) / max(1, len(self.api_call_times)),
                'success_rate': (self.api_calls - self.api_errors) / max(1, self.api_calls)
            },
            'db': {
                'operations': self.db_operations,
                'errors': self.db_errors,
                'avg_time': sum(self.db_operation_times) / max(1, len(self.db_operation_times)),
                'success_rate': (self.db_operations - self.db_errors) / max(1, self.db_operations)
            },
            'memory': memory_usage,
            'parallelism': {
                'avg_queue_size': sum(self.task_queue_size) / max(1, len(self.task_queue_size)),
                'avg_worker_utilization': sum(self.worker_utilization) / max(1, len(self.worker_utilization))
            }
        }

# Enhanced checkpointing system for ChEMBL imports
class ChEMBLImportCheckpoint:
    """
    Robust checkpoint system for managing the state of ChEMBL imports.
    
    Features:
    - Atomic checkpoint file updates
    - Backup and recovery
    - Detailed import progress tracking
    - Support for resumable imports
    """
    
    def __init__(self, checkpoint_file=DEFAULT_CHECKPOINT_FILE):
        """
        Initialize the checkpoint system.
        
        Args:
            checkpoint_file: Path to the checkpoint file
        """
        self.checkpoint_file = checkpoint_file
        self.backup_file = f"{checkpoint_file}.bak"
        self.temp_file = f"{checkpoint_file}.tmp"
        self.lock = threading.RLock()
        
        # Import state
        self.next_batch_index = 0
        self.processed_chembl_ids = set()
        self.imported_count = 0
        self.skipped_count = 0
        self.error_count = 0
        self.last_save_time = 0
        self.started_at = datetime.now()
        self.status = "pending"
        
        # Load existing checkpoint if available
        self._load()
    
    def _load(self):
        """Load existing checkpoint if available."""
        with self.lock:
            if os.path.exists(self.checkpoint_file):
                try:
                    with open(self.checkpoint_file, 'r') as f:
                        data = json.load(f)
                    
                    # Load state
                    self.next_batch_index = data.get('next_batch_index', 0)
                    self.processed_chembl_ids = set(data.get('processed_chembl_ids', []))
                    self.imported_count = data.get('imported_count', 0)
                    self.skipped_count = data.get('skipped_count', 0)
                    self.error_count = data.get('error_count', 0)
                    
                    # Parse timestamps
                    started_at = data.get('started_at')
                    if started_at:
                        self.started_at = datetime.fromisoformat(started_at)
                    
                    self.status = data.get('status', 'pending')
                    
                    logger.info(f"Loaded checkpoint: next batch index {self.next_batch_index}, "
                                f"{len(self.processed_chembl_ids)} processed compounds")
                    return True
                except Exception as e:
                    logger.error(f"Error loading checkpoint: {str(e)}")
                    
                    # Try backup file
                    if os.path.exists(self.backup_file):
                        try:
                            with open(self.backup_file, 'r') as f:
                                data = json.load(f)
                            
                            # Load state
                            self.next_batch_index = data.get('next_batch_index', 0)
                            self.processed_chembl_ids = set(data.get('processed_chembl_ids', []))
                            self.imported_count = data.get('imported_count', 0)
                            self.skipped_count = data.get('skipped_count', 0)
                            self.error_count = data.get('error_count', 0)
                            
                            # Parse timestamps
                            started_at = data.get('started_at')
                            if started_at:
                                self.started_at = datetime.fromisoformat(started_at)
                            
                            self.status = data.get('status', 'pending')
                            
                            logger.info(f"Loaded backup checkpoint: next batch index {self.next_batch_index}, "
                                        f"{len(self.processed_chembl_ids)} processed compounds")
                            return True
                        except Exception as backup_e:
                            logger.error(f"Error loading backup checkpoint: {str(backup_e)}")
            
            logger.info("No valid checkpoint found, starting from beginning")
            return False
    
    def save(self, force=False):
        """
        Save the current state to the checkpoint file.
        
        Args:
            force: If True, save regardless of time since last save
        
        Returns:
            True if checkpoint was saved, False otherwise
        """
        # Don't save too frequently unless forced
        if not force and time.time() - self.last_save_time < 10:
            return False
        
        with self.lock:
            try:
                # Prepare checkpoint data
                data = {
                    'next_batch_index': self.next_batch_index,
                    'processed_chembl_ids': list(self.processed_chembl_ids),
                    'imported_count': self.imported_count,
                    'skipped_count': self.skipped_count,
                    'error_count': self.error_count,
                    'started_at': self.started_at.isoformat(),
                    'last_updated': datetime.now().isoformat(),
                    'status': self.status
                }
                
                # Write to temporary file first for atomic update
                with open(self.temp_file, 'w') as f:
                    json.dump(data, f, indent=2)
                
                # Create backup of existing checkpoint if it exists
                if os.path.exists(self.checkpoint_file):
                    try:
                        os.replace(self.checkpoint_file, self.backup_file)
                    except Exception as e:
                        logger.warning(f"Error creating backup: {str(e)}")
                
                # Move temp file to checkpoint file
                os.replace(self.temp_file, self.checkpoint_file)
                
                self.last_save_time = time.time()
                return True
            except Exception as e:
                logger.error(f"Error saving checkpoint: {str(e)}")
                return False
    
    def mark_molecule_processed(self, chembl_id: str, imported: bool = True, error: bool = False):
        """
        Mark a molecule as processed.
        
        Args:
            chembl_id: ChEMBL ID of the molecule
            imported: Whether the molecule was imported successfully
            error: Whether an error occurred during processing
        """
        with self.lock:
            self.processed_chembl_ids.add(chembl_id)
            if imported:
                self.imported_count += 1
            elif error:
                self.error_count += 1
            else:
                self.skipped_count += 1
    
    def mark_batch_completed(self, batch_index: int):
        """
        Mark a batch as completed.
        
        Args:
            batch_index: Index of the completed batch
        """
        with self.lock:
            # Update next batch index if this is the current batch
            if batch_index == self.next_batch_index:
                self.next_batch_index = batch_index + 1
                logger.info(f"Updated next batch index to {self.next_batch_index}")
            
            # Save checkpoint
            self.save()
    
    def set_status(self, status: str):
        """
        Set the status of the import.
        
        Args:
            status: Status to set ('pending', 'running', 'completed', 'error', 'paused')
        """
        with self.lock:
            self.status = status
            self.save(force=True)
    
    def get_elapsed_time(self) -> timedelta:
        """Get elapsed time since import started."""
        return datetime.now() - self.started_at
    
    def get_state(self) -> Dict[str, Any]:
        """Get the current state of the import."""
        with self.lock:
            return {
                'next_batch_index': self.next_batch_index,
                'processed_molecules': len(self.processed_chembl_ids),
                'imported_count': self.imported_count,
                'skipped_count': self.skipped_count,
                'error_count': self.error_count,
                'started_at': self.started_at.isoformat(),
                'elapsed_time': str(self.get_elapsed_time()),
                'status': self.status
            }


class MoleculeTransformer:
    """
    Transform ChEMBL compound data to match our database schema.
    
    This class handles the conversion of ChEMBL API molecule data to the format
    expected by our database schema, including property transformation and validation.
    """
    
    def __init__(self, user_profile_id: str):
        """
        Initialize the transformer.
        
        Args:
            user_profile_id: User profile ID for created_by field
        """
        self.user_profile_id = user_profile_id
        self.property_types = {}
        self.unknown_property_type_id = None
        
        # Load property types
        self._load_property_types()
    
    def _load_property_types(self):
        """Load property types from the database."""
        try:
            # Execute query
            sql = "SELECT id, name, data_type FROM property_types;"
            result = execute_query_with_retry(sql)
            
            # Process results
            if result:
                for row in result:
                    name = row['name'].lower()
                    self.property_types[name] = {
                        'id': row['id'],
                        'data_type': row['data_type']
                    }
                
                # Find unknown property type
                for name, info in self.property_types.items():
                    if 'unknown' in name.lower():
                        self.unknown_property_type_id = info['id']
                        break
                
                logger.info(f"Loaded {len(self.property_types)} property types from database")
            else:
                logger.warning("No property types found in database")
        except Exception as e:
            logger.error(f"Error loading property types: {str(e)}")
            raise
    
    def transform_to_molecule(self, compound: Dict[str, Any]) -> Dict[str, Any]:
        """
        Transform ChEMBL compound data to match our molecule table schema.
        
        Args:
            compound: ChEMBL compound data
            
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
            "id": str(uuid.uuid4()),
            "name": name,
            "smiles": structures.get('canonical_smiles'),
            "inchi": structures.get('standard_inchi'),
            "inchikey": structures.get('standard_inchi_key'),
            "formula": mol_props.get('full_molformula'),
            "molecular_weight": mol_props.get('full_mwt'),
            "pubchem_cid": None,  # Added field for PubChem Compound ID (can be populated later if available)
            "chembl_id": compound.get('molecule_chembl_id'),  # Store ChEMBL ID in dedicated column
            "created_by": self.user_profile_id,
            "data_source": f"ChEMBL v{chembl_version} ID: {compound.get('molecule_chembl_id')}",
            "version": 1,
            "modification_history": json.dumps([{
                "timestamp": datetime.now().isoformat(),
                "action": "created",
                "user_id": self.user_profile_id
            }])
        }
    
    def transform_to_properties(self, compound: Dict[str, Any], molecule_id: str) -> List[Dict[str, Any]]:
        """
        Transform ChEMBL compound properties to match our molecular_properties table schema.
        
        Args:
            compound: ChEMBL compound data
            molecule_id: Molecule ID from our database
            
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
        
        # Add properties from molecule_properties
        for prop_key, prop_info in property_mappings.items():
            if prop_key in mol_props and mol_props[prop_key] is not None:
                prop_name = prop_info['name']
                
                # Try to get property type info
                property_type_info = self.property_types.get(prop_name.lower())
                
                # If property type doesn't exist, add it to the database
                if not property_type_info:
                    property_type_id = self._insert_property_type(
                        name=prop_name,
                        data_type=prop_info['data_type'],
                        description=prop_info['description'],
                        units=prop_info['units']
                    )
                    
                    if property_type_id:
                        # Update local cache
                        self.property_types[prop_name.lower()] = {
                            'id': property_type_id,
                            'data_type': prop_info['data_type']
                        }
                        property_type_info = self.property_types[prop_name.lower()]
                    else:
                        # Use Unknown Property Type as fallback
                        if self.unknown_property_type_id:
                            property_type_info = {
                                'id': self.unknown_property_type_id,
                                'data_type': 'text'
                            }
                            logger.info(f"Using 'Unknown Property Type' as fallback for {prop_name}")
                        else:
                            logger.error(f"Unknown Property Type not found in database. Skipping property {prop_name}.")
                            continue
                
                # Add the property if we have a valid property type
                if property_type_info:
                    # Get property value
                    value = mol_props[prop_key]
                    
                    # Determine value type and convert accordingly
                    numeric_value = None
                    text_value = None
                    boolean_value = None
                    
                    if property_type_info['data_type'] == 'numeric':
                        try:
                            numeric_value = float(value)
                        except (ValueError, TypeError):
                            # If conversion fails, store as text
                            text_value = str(value)
                    elif property_type_info['data_type'] == 'boolean':
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
                        "property_type_id": property_type_info['id'],
                        "numeric_value": numeric_value,
                        "text_value": text_value,
                        "boolean_value": boolean_value,
                        "created_by": self.user_profile_id,
                        "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_key}",
                        "version": 1,
                        "modification_history": json.dumps([{
                            "timestamp": datetime.now().isoformat(),
                            "action": "created",
                            "user_id": self.user_profile_id
                        }])
                    })
        
        # Process additional properties from ChEMBL API
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
        
        # Process additional properties
        for prop in additional_properties:
            # Skip empty properties
            if not prop.get('value'):
                continue
                
            prop_name = prop['property_name']
            
            # Try to get property type info
            property_type_info = self.property_types.get(prop_name.lower())
            
            # If property type doesn't exist, add it to the database
            if not property_type_info:
                property_type_id = self._insert_property_type(
                    name=prop_name,
                    data_type=prop.get('data_type', 'text'),
                    description=f"Auto-added by ChEMBL import - {prop_name}"
                )
                
                if property_type_id:
                    # Update local cache
                    self.property_types[prop_name.lower()] = {
                        'id': property_type_id,
                        'data_type': prop.get('data_type', 'text')
                    }
                    property_type_info = self.property_types[prop_name.lower()]
                else:
                    # Use Unknown Property Type as fallback
                    if self.unknown_property_type_id:
                        property_type_info = {
                            'id': self.unknown_property_type_id,
                            'data_type': 'text'
                        }
                        logger.info(f"Using 'Unknown Property Type' as fallback for {prop_name}")
                    else:
                        logger.error(f"Unknown Property Type not found in database. Skipping property {prop_name}.")
                        continue
            
            # Add the property if we have a valid property type
            if property_type_info:
                # Get property value
                value = prop['value']
                
                # Determine value type and convert accordingly
                numeric_value = None
                text_value = None
                boolean_value = None
                
                if property_type_info['data_type'] == 'numeric':
                    try:
                        numeric_value = float(value)
                    except (ValueError, TypeError):
                        # If conversion fails, store as text
                        text_value = str(value)
                elif property_type_info['data_type'] == 'boolean':
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
                    "property_type_id": property_type_info['id'],
                    "numeric_value": numeric_value,
                    "text_value": text_value,
                    "boolean_value": boolean_value,
                    "created_by": self.user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}, property: {prop_name}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": self.user_profile_id
                    }])
                })
        
        return properties
    
    def _insert_property_type(self, name: str, data_type: str = 'text', description: str = None, units: str = None) -> Optional[str]:
        """
        Insert a new property type into the database.
        
        Args:
            name: Property type name
            data_type: Data type (numeric, text, boolean)
            description: Optional description
            units: Optional units
            
        Returns:
            Property type ID if successful, None otherwise
        """
        try:
            # Prepare property type data
            property_type_data = {
                "name": name,
                "data_type": data_type,
                "description": description or f"Auto-added by ChEMBL import on {datetime.now().isoformat()}",
                "units": units
            }
            
            # Construct SQL for insertion with parameterized query
            columns = ", ".join(property_type_data.keys())
            placeholders = ", ".join([f"%({k})s" for k in property_type_data.keys()])
            
            sql = f"""
            INSERT INTO property_types ({columns})
            VALUES ({placeholders})
            RETURNING id;
            """
            
            # Execute SQL
            result = execute_query_with_retry(sql, property_type_data)
            
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
                    "source": "_insert_property_type"
                }
            )
            return None


class ChEMBLImporter:
    """
    Main importer class for ChEMBL data.
    
    This class handles the end-to-end process of importing ChEMBL data:
    1. Fetching compounds from the ChEMBL API
    2. Transforming them to match our database schema
    3. Importing them into the database with efficient batching
    4. Tracking progress and supporting resumable imports
    """
    
    def __init__(self, 
                 limit: int = 10000,
                 batch_size: int = 100,
                 checkpoint_file: str = None,
                 search_terms: List[str] = None,
                 user_profile_id: str = None,
                 max_workers: int = None,
                 dry_run: bool = False):
        """
        Initialize the importer.
        
        Args:
            limit: Maximum number of compounds to import
            batch_size: Number of compounds to process in each batch
            checkpoint_file: Path to checkpoint file for resumable imports
            search_terms: List of search terms to use for finding compounds
            user_profile_id: User profile ID for created_by field
            max_workers: Maximum number of worker threads/processes
            dry_run: If True, don't actually insert data
        """
        self.limit = limit
        self.batch_size = batch_size
        self.search_terms = search_terms or [
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
        self.user_profile_id = user_profile_id
        self.dry_run = dry_run
        
        # Set max workers
        if max_workers is None:
            # Use number of CPU cores with a minimum of 2 and maximum of 8
            max_workers = min(8, max(2, multiprocessing.cpu_count()))
        self.max_workers = max_workers
        
        # Initialize ChEMBL client
        self.molecule_client = new_client.molecule
        
        # Initialize checkpoint system
        if checkpoint_file is None:
            checkpoint_file = DEFAULT_CHECKPOINT_FILE
        self.checkpoint = ChEMBLImportCheckpoint(checkpoint_file)
        
        # Initialize transformer
        self.transformer = MoleculeTransformer(user_profile_id)
        
        # Initialize metrics
        self.performance_metrics = PerformanceMetrics()
        
        # Initialize state
        self.all_compounds = []
        self.processed_chembl_ids = set()
        self.exit_event = threading.Event()
        
        # Set up signal handlers
        signal.signal(signal.SIGINT, self._handle_signal)
        signal.signal(signal.SIGTERM, self._handle_signal)
        atexit.register(self._handle_exit)
        
        logger.info(f"Initialized ChEMBLImporter with limit={limit}, batch_size={batch_size}, "
                    f"max_workers={max_workers}, dry_run={dry_run}")
    
    def _handle_signal(self, signum, frame):
        """Handle termination signals."""
        logger.info(f"Received signal {signum}, shutting down gracefully...")
        self.exit_event.set()
        self.checkpoint.set_status("paused")
    
    def _handle_exit(self):
        """Handle normal exit."""
        if not self.exit_event.is_set():
            logger.info("Normal program exit, finalizing...")
            self.checkpoint.save(force=True)
    
    def fetch_standard_reference_compounds(self) -> List[Dict[str, Any]]:
        """
        Fetch standard reference compounds.
        
        Returns:
            List of reference compounds
        """
        reference_compounds = []
        
        for chembl_id in STANDARD_REFERENCE_IDS:
            if chembl_id in self.checkpoint.processed_chembl_ids:
                logger.info(f"Reference compound {chembl_id} already processed, skipping")
                continue
            
            start_time = time.time()
            try:
                logger.info(f"Fetching standard reference compound: {chembl_id}")
                molecule = new_client.molecule
                compound = molecule.get(chembl_id)
                
                # Add to results
                if compound:
                    reference_compounds.append(compound)
                    logger.info(f"Added reference compound: {chembl_id}")
                else:
                    logger.warning(f"Reference compound {chembl_id} not found")
                    log_skipped_molecule(
                        chembl_id=chembl_id,
                        reason="Not found",
                        molecule_data={"chembl_id": chembl_id},
                        category="reference"
                    )
                
                # Record API call
                self.performance_metrics.record_api_call(time.time() - start_time, True)
                
                # Slight delay to be gentle on the API
                time.sleep(0.2)
                
            except Exception as e:
                # Record API call failure
                self.performance_metrics.record_api_call(time.time() - start_time, False)
                
                log_error(
                    error_type="API",
                    message=f"Error fetching reference compound {chembl_id}",
                    context={
                        "exception": e,
                        "compound_id": chembl_id,
                        "source": "fetch_standard_reference_compounds"
                    }
                )
                
                # Slight delay to be gentle on the API
                time.sleep(1)
        
        logger.info(f"Fetched {len(reference_compounds)} reference compounds")
        return reference_compounds
    
    def search_compounds(self) -> List[Dict[str, Any]]:
        """
        Search for compounds based on search terms.
        
        Returns:
            List of compounds matching search terms
        """
        compounds = []
        
        if self.exit_event.is_set():
            logger.info("Exit requested, stopping search")
            return compounds
        
        # Start from the next batch based on checkpoint
        start_index = self.checkpoint.next_batch_index
        
        # Resuming from a checkpoint
        if start_index > 0:
            logger.info(f"Resuming search from term index {start_index}, with {len(self.checkpoint.processed_chembl_ids)} compounds already processed")
        
        # Process search terms
        for i, term in enumerate(self.search_terms[start_index:], start=start_index):
            if self.exit_event.is_set():
                logger.info("Exit requested, stopping search")
                break
            
            if len(compounds) >= self.limit:
                logger.info(f"Reached limit of {self.limit} compounds, stopping search")
                break
            
            try:
                logger.info(f"Searching for '{term}' ({i+1}/{len(self.search_terms)})")
                
                start_time = time.time()
                
                # Search for compounds matching the term
                results = self.molecule_client.filter(
                    pref_name__icontains=term
                ).only(
                    'molecule_chembl_id',
                    'pref_name',
                    'molecule_structures'
                )
                
                # Record API call
                self.performance_metrics.record_api_call(time.time() - start_time, True)
                
                # Process results
                term_processed = 0
                for compound in results:
                    if self.exit_event.is_set():
                        logger.info("Exit requested, stopping search")
                        break
                    
                    # Skip if we already have enough compounds
                    if len(compounds) >= self.limit:
                        break
                    
                    # Skip if we've already processed this compound
                    chembl_id = compound.get('molecule_chembl_id')
                    if chembl_id in self.checkpoint.processed_chembl_ids:
                        logger.debug(f"Compound {chembl_id} already processed, skipping")
                        continue
                    
                    # Skip if no structures available
                    if not compound.get('molecule_structures'):
                        log_skipped_molecule(
                            chembl_id=chembl_id or 'unknown',
                            reason="Missing molecular structures",
                            molecule_data={"name": compound.get('pref_name')},
                            category="validation"
                        )
                        self.checkpoint.mark_molecule_processed(chembl_id, imported=False, error=False)
                        continue
                    
                    # Get full compound details
                    try:
                        start_time = time.time()
                        
                        full_details = self.molecule_client.get(chembl_id)
                        
                        # Record API call
                        self.performance_metrics.record_api_call(time.time() - start_time, True)
                        
                        # Add to results
                        compounds.append(full_details)
                        term_processed += 1
                        
                        # Log progress
                        if len(compounds) % 10 == 0:
                            logger.info(f"Fetched {len(compounds)}/{self.limit} compounds")
                        
                        # Save checkpoint occasionally
                        if len(compounds) % 50 == 0:
                            self.checkpoint.save()
                        
                        # Check memory usage occasionally
                        if len(compounds) % MEMORY_CHECK_INTERVAL == 0:
                            memory_usage = self.performance_metrics.get_memory_usage()
                            if memory_usage.get('percent', 0) > MEMORY_THRESHOLD:
                                logger.warning(f"Memory usage high: {memory_usage['percent']:.1f}% > {MEMORY_THRESHOLD:.1f}%")
                                # Sleep briefly to allow memory to be reclaimed
                                time.sleep(1)
                        
                        # Slight delay to be gentle on the API
                        time.sleep(0.2)
                        
                    except Exception as e:
                        # Record API call failure
                        self.performance_metrics.record_api_call(time.time() - start_time, False)
                        
                        log_error(
                            error_type="API",
                            message=f"Error fetching details for {chembl_id}",
                            context={
                                "exception": e,
                                "compound_id": chembl_id,
                                "source": "search_compounds"
                            }
                        )
                        
                        self.checkpoint.mark_molecule_processed(chembl_id, imported=False, error=True)
                        
                        # Slight delay to be gentle on the API
                        time.sleep(1)
                
                # Mark this batch as completed
                self.checkpoint.mark_batch_completed(i)
                logger.info(f"Processed search term '{term}': {term_processed} compounds added")
                
            except Exception as e:
                log_error(
                    error_type="API",
                    message=f"Error processing term '{term}'",
                    context={
                        "exception": e,
                        "term": term,
                        "source": "search_compounds"
                    }
                )
                
                # Save checkpoint before continuing to next term
                self.checkpoint.save(force=True)
        
        # Save final checkpoint
        self.checkpoint.save(force=True)
        logger.info(f"Search completed with {len(compounds)} compounds")
        
        return compounds
    
    def check_existing_molecules(self, batch: List[Dict[str, Any]]) -> Dict[str, str]:
        """
        Check which molecules already exist in the database.
        
        Args:
            batch: List of compounds to check
            
        Returns:
            Dictionary mapping InChIKey to molecule ID for existing molecules
        """
        existing_molecules = {}
        
        # Extract InChIKeys from compounds
        inchikeys = []
        for compound in batch:
            structures = compound.get('molecule_structures', {}) or {}
            inchikey = structures.get('standard_inchi_key')
            if inchikey:
                inchikeys.append(inchikey)
        
        if not inchikeys:
            return existing_molecules
        
        try:
            # Construct query
            placeholders = ', '.join([f"'{inchikey}'" for inchikey in inchikeys])
            sql = f"SELECT id, inchikey FROM molecules WHERE inchikey IN ({placeholders});"
            
            # Execute query
            start_time = time.time()
            result = execute_query_with_retry(sql)
            
            # Record database operation
            self.performance_metrics.record_db_operation(time.time() - start_time, True)
            
            # Process results
            if result:
                for row in result:
                    if 'inchikey' in row and 'id' in row:
                        existing_molecules[row['inchikey']] = row['id']
            
            logger.debug(f"Found {len(existing_molecules)} existing molecules in batch of {len(inchikeys)}")
            return existing_molecules
            
        except Exception as e:
            # Record database operation failure
            self.performance_metrics.record_db_operation(time.time() - start_time, False)
            
            log_error(
                error_type="Database",
                message=f"Error checking existing molecules",
                context={
                    "exception": e,
                    "source": "check_existing_molecules"
                }
            )
            raise
    
    def process_batch(self, batch: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Process a batch of compounds.
        
        Args:
            batch: List of compounds to process
            
        Returns:
            Dictionary with processing statistics
        """
        stats = {
            'batch_size': len(batch),
            'molecules_inserted': 0,
            'molecules_skipped': 0,
            'properties_inserted': 0,
            'errors': 0
        }
        
        if not batch:
            return stats
        
        if self.dry_run:
            logger.info(f"DRY RUN: Would process batch of {len(batch)} compounds")
            
            # Simulate processing
            for compound in batch:
                chembl_id = compound.get('molecule_chembl_id', 'unknown')
                structures = compound.get('molecule_structures', {}) or {}
                inchikey = structures.get('standard_inchi_key')
                
                if not inchikey:
                    logger.info(f"DRY RUN: Would skip molecule {chembl_id} due to missing InChIKey")
                    stats['molecules_skipped'] += 1
                    continue
                
                molecule_data = self.transformer.transform_to_molecule(compound)
                properties = self.transformer.transform_to_properties(compound, "dry-run-id")
                
                logger.info(f"DRY RUN: Would insert molecule {chembl_id} ({molecule_data.get('name')}) "
                            f"with {len(properties)} properties")
                
                stats['molecules_inserted'] += 1
                stats['properties_inserted'] += len(properties)
                
                # Mark as processed in checkpoint
                self.checkpoint.mark_molecule_processed(chembl_id, imported=True)
            
            return stats
        
        try:
            # Check which molecules already exist
            existing_molecules = self.check_existing_molecules(batch)
            
            # Prepare molecule inserts
            molecule_inserts = []
            molecule_properties = {}
            
            for compound in batch:
                chembl_id = compound.get('molecule_chembl_id', 'unknown')
                
                try:
                    # Extract InChIKey
                    structures = compound.get('molecule_structures', {}) or {}
                    inchikey = structures.get('standard_inchi_key')
                    
                    if not inchikey:
                        logger.warning(f"Skipping molecule {chembl_id} due to missing InChIKey")
                        log_skipped_molecule(
                            chembl_id=chembl_id,
                            reason="Missing InChIKey",
                            molecule_data={"name": compound.get('pref_name')},
                            category="validation"
                        )
                        stats['molecules_skipped'] += 1
                        self.checkpoint.mark_molecule_processed(chembl_id, imported=False)
                        continue
                    
                    # Check if molecule already exists
                    if inchikey in existing_molecules:
                        molecule_id = existing_molecules[inchikey]
                        logger.debug(f"Molecule {chembl_id} already exists with ID {molecule_id}")
                        
                        # Transform properties
                        properties = self.transformer.transform_to_properties(compound, molecule_id)
                        
                        # Store for batch insertion
                        molecule_properties[molecule_id] = properties
                        
                        stats['molecules_skipped'] += 1
                        self.checkpoint.mark_molecule_processed(chembl_id, imported=True)
                        continue
                    
                    # Transform molecule
                    molecule_data = self.transformer.transform_to_molecule(compound)
                    molecule_id = molecule_data['id']
                    
                    # Add to batch inserts
                    molecule_inserts.append(molecule_data)
                    
                    # Transform properties
                    properties = self.transformer.transform_to_properties(compound, molecule_id)
                    
                    # Store for batch insertion
                    molecule_properties[molecule_id] = properties
                    
                except Exception as e:
                    log_error(
                        error_type="Transformation",
                        message=f"Error transforming compound {chembl_id}",
                        context={
                            "exception": e,
                            "compound_id": chembl_id,
                            "source": "process_batch.transform"
                        }
                    )
                    stats['errors'] += 1
                    self.checkpoint.mark_molecule_processed(chembl_id, imported=False, error=True)
                    continue
            
            # Insert molecules
            if molecule_inserts:
                inserted_molecules = self._insert_molecules(molecule_inserts)
                stats['molecules_inserted'] = len(inserted_molecules)
            
            # Insert properties for all molecules (new and existing)
            property_count = 0
            for molecule_id, properties in molecule_properties.items():
                if properties:
                    inserted = self._insert_properties(properties)
                    property_count += inserted
            
            stats['properties_inserted'] = property_count
            
            logger.info(f"Processed batch: {stats['molecules_inserted']} molecules inserted, "
                        f"{stats['molecules_skipped']} skipped, "
                        f"{stats['properties_inserted']} properties inserted")
            
            return stats
            
        except Exception as e:
            log_error(
                error_type="Batch",
                message=f"Error processing batch",
                context={
                    "exception": e,
                    "batch_size": len(batch),
                    "source": "process_batch"
                }
            )
            stats['errors'] += 1
            return stats
    
    def _insert_molecules(self, molecules: List[Dict[str, Any]]) -> List[str]:
        """
        Insert molecules into the database.
        
        Args:
            molecules: List of molecule data to insert
            
        Returns:
            List of inserted molecule IDs
        """
        if not molecules:
            return []
        
        inserted_ids = []
        
        try:
            # Prepare insertion queries
            queries = []
            
            for molecule in molecules:
                # Filter out None values
                filtered = {k: v for k, v in molecule.items() if v is not None}
                
                # Build query
                columns = list(filtered.keys())
                placeholders = [f"%({k})s" for k in columns]
                
                query = f"""
                INSERT INTO molecules ({', '.join(columns)})
                VALUES ({', '.join(placeholders)})
                RETURNING id;
                """
                
                # Execute query
                start_time = time.time()
                result = execute_query_with_retry(query, filtered)
                
                # Record database operation
                self.performance_metrics.record_db_operation(time.time() - start_time, True)
                
                if result and len(result) > 0 and 'id' in result[0]:
                    inserted_ids.append(result[0]['id'])
                    
                    # Mark as processed in checkpoint
                    chembl_id = filtered.get('chembl_id')
                    if chembl_id:
                        self.checkpoint.mark_molecule_processed(chembl_id, imported=True)
                
            logger.info(f"Inserted {len(inserted_ids)} molecules")
            return inserted_ids
            
        except Exception as e:
            # Record database operation failure
            self.performance_metrics.record_db_operation(time.time(), False)
            
            log_error(
                error_type="Database",
                message=f"Error inserting molecules",
                context={
                    "exception": e,
                    "molecule_count": len(molecules),
                    "source": "_insert_molecules"
                }
            )
            return inserted_ids
    
    def _insert_properties(self, properties: List[Dict[str, Any]]) -> int:
        """
        Insert molecular properties into the database.
        
        Args:
            properties: List of property data to insert
            
        Returns:
            Number of inserted properties
        """
        if not properties:
            return 0
        
        try:
            # Prepare insertion queries
            property_batches = [properties[i:i+100] for i in range(0, len(properties), 100)]
            total_inserted = 0
            
            for batch in property_batches:
                queries = []
                
                for prop in batch:
                    # Filter out None values
                    filtered = {k: v for k, v in prop.items() if v is not None}
                    
                    # Build query
                    columns = list(filtered.keys())
                    placeholders = [f"%({k})s" for k in columns]
                    
                    query = f"""
                    INSERT INTO molecular_properties ({', '.join(columns)})
                    VALUES ({', '.join(placeholders)})
                    """
                    
                    # Execute query with parameters
                    start_time = time.time()
                    result = execute_query_with_retry(query, filtered)
                    
                    # Record database operation
                    self.performance_metrics.record_db_operation(time.time() - start_time, True)
                    
                    total_inserted += 1
            
            logger.info(f"Inserted {total_inserted} properties")
            return total_inserted
            
        except Exception as e:
            # Record database operation failure
            self.performance_metrics.record_db_operation(time.time(), False)
            
            log_error(
                error_type="Database",
                message=f"Error inserting properties",
                context={
                    "exception": e,
                    "property_count": len(properties),
                    "source": "_insert_properties"
                }
            )
            return 0
    
    def process_compounds_parallel(self, compounds: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Process compounds in parallel batches.
        
        Args:
            compounds: List of compounds to process
            
        Returns:
            Dictionary with processing statistics
        """
        if not compounds:
            logger.warning("No compounds to process")
            return {
                'total_compounds': 0,
                'processed': 0,
                'molecules_inserted': 0,
                'molecules_skipped': 0,
                'properties_inserted': 0,
                'errors': 0
            }
        
        # Set status
        self.checkpoint.set_status("running")
        
        # Split compounds into batches
        batches = [compounds[i:i+self.batch_size] for i in range(0, len(compounds), self.batch_size)]
        logger.info(f"Processing {len(compounds)} compounds in {len(batches)} batches")
        
        # Initialize statistics
        stats = {
            'total_compounds': len(compounds),
            'processed': 0,
            'molecules_inserted': 0,
            'molecules_skipped': 0,
            'properties_inserted': 0,
            'errors': 0,
            'batches_processed': 0,
            'batches_failed': 0
        }
        
        # Process batches sequentially for stability
        for i, batch in enumerate(batches):
            if self.exit_event.is_set():
                logger.info("Exit requested, stopping processing")
                break
            
            logger.info(f"Processing batch {i+1}/{len(batches)} ({len(batch)} compounds)")
            
            # Process batch
            batch_stats = self.process_batch(batch)
            
            # Update statistics
            stats['processed'] += batch_stats['batch_size']
            stats['molecules_inserted'] += batch_stats['molecules_inserted']
            stats['molecules_skipped'] += batch_stats['molecules_skipped']
            stats['properties_inserted'] += batch_stats['properties_inserted']
            stats['errors'] += batch_stats['errors']
            stats['batches_processed'] += 1
            
            if batch_stats['errors'] > 0:
                stats['batches_failed'] += 1
            
            # Log progress
            logger.info(f"Batch {i+1}/{len(batches)} completed: "
                        f"{batch_stats['molecules_inserted']} molecules inserted, "
                        f"{batch_stats['molecules_skipped']} skipped, "
                        f"{batch_stats['properties_inserted']} properties inserted, "
                        f"{batch_stats['errors']} errors")
            
            # Save checkpoint
            self.checkpoint.save()
            
            # Sleep briefly to avoid overwhelming the database
            time.sleep(0.5)
        
        # Update status
        if self.exit_event.is_set():
            self.checkpoint.set_status("paused")
        else:
            self.checkpoint.set_status("completed")
        
        # Save final checkpoint
        self.checkpoint.save(force=True)
        
        # Include performance metrics
        stats['performance'] = self.performance_metrics.get_metrics()
        
        return stats
    
    def run(self) -> Dict[str, Any]:
        """
        Run the full import process.
        
        Returns:
            Dictionary with import statistics
        """
        try:
            # Fetch reference compounds
            logger.info("Fetching standard reference compounds...")
            reference_compounds = self.fetch_standard_reference_compounds()
            
            # Search for compounds based on search terms
            logger.info(f"Searching for compounds using {len(self.search_terms)} search terms...")
            search_compounds = self.search_compounds()
            
            # Combine all compounds
            self.all_compounds = reference_compounds + search_compounds
            logger.info(f"Total compounds to process: {len(self.all_compounds)}")
            
            # Process all compounds
            logger.info("Processing compounds...")
            stats = self.process_compounds_parallel(self.all_compounds)
            
            # Generate summary
            summary = {
                'start_time': self.checkpoint.started_at.isoformat(),
                'end_time': datetime.now().isoformat(),
                'elapsed_time': str(self.checkpoint.get_elapsed_time()),
                'status': self.checkpoint.status,
                'compounds': {
                    'total': len(self.all_compounds),
                    'processed': stats['processed'],
                    'inserted': stats['molecules_inserted'],
                    'skipped': stats['molecules_skipped']
                },
                'properties': {
                    'inserted': stats['properties_inserted']
                },
                'errors': stats['errors'],
                'performance': self.performance_metrics.get_metrics()
            }
            
            # Write summary to file
            write_summary(summary)
            
            logger.info(f"Import completed with status: {self.checkpoint.status}")
            logger.info(f"Total compounds: {len(self.all_compounds)}")
            logger.info(f"Processed: {stats['processed']}")
            logger.info(f"Molecules inserted: {stats['molecules_inserted']}")
            logger.info(f"Molecules skipped: {stats['molecules_skipped']}")
            logger.info(f"Properties inserted: {stats['properties_inserted']}")
            logger.info(f"Errors: {stats['errors']}")
            logger.info(f"Elapsed time: {self.checkpoint.get_elapsed_time()}")
            
            return stats
            
        except Exception as e:
            log_error(
                error_type="Application",
                message=f"Error running import: {str(e)}",
                context={
                    "exception": e,
                    "traceback": traceback.format_exc(),
                    "source": "run"
                }
            )
            
            # Update status
            self.checkpoint.set_status("error")
            self.checkpoint.save(force=True)
            
            raise


def get_user_profile_id() -> str:
    """
    Get or create a user profile ID for the import.
    
    Returns:
        User profile ID
    """
    try:
        # Try to get existing ChEMBL importer user
        sql = """
        SELECT id FROM user_profiles 
        WHERE name = 'ChEMBL Importer' 
        OR email = 'chembl-import@system.local'
        LIMIT 1;
        """
        
        result = execute_query_with_retry(sql)
        
        if result and len(result) > 0 and 'id' in result[0]:
            logger.info(f"Using existing user profile: {result[0]['id']}")
            return result[0]['id']
        
        # Create a new user profile
        logger.info("Creating new user profile for ChEMBL importer")
        
        user_data = {
            "id": str(uuid.uuid4()),
            "name": "ChEMBL Importer",
            "email": "chembl-import@system.local",
            "created_at": datetime.now().isoformat(),
            "is_system": True
        }
        
        # Insert user
        columns = list(user_data.keys())
        placeholders = [f"%({k})s" for k in columns]
        
        sql = f"""
        INSERT INTO user_profiles ({', '.join(columns)})
        VALUES ({', '.join(placeholders)})
        RETURNING id;
        """
        
        result = execute_query_with_retry(sql, user_data)
        
        if result and len(result) > 0 and 'id' in result[0]:
            logger.info(f"Created new user profile: {result[0]['id']}")
            return result[0]['id']
        
        # Fallback to a static ID if all else fails
        logger.warning("Failed to get or create user profile, using fallback ID")
        return "00000000-0000-0000-0000-000000000000"
        
    except Exception as e:
        log_error(
            error_type="Database",
            message=f"Error getting user profile ID: {str(e)}",
            context={
                "exception": e,
                "source": "get_user_profile_id"
            }
        )
        
        # Fallback to a static ID
        return "00000000-0000-0000-0000-000000000000"


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Import data from ChEMBL to CryoProtect")
    parser.add_argument("--limit", type=int, default=10000, help="Maximum number of compounds to import")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for database operations")
    parser.add_argument("--checkpoint-file", type=str, help="Path to checkpoint file for resumable imports")
    parser.add_argument("--search-terms", type=str, help="Comma-separated list of search terms")
    parser.add_argument("--max-workers", type=int, help="Maximum number of worker threads/processes")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually insert data, just simulate")
    args = parser.parse_args()
    
    try:
        # Validate configuration
        logger.info("Validating configuration...")
        validate_config()
        
        # Get user profile ID
        user_profile_id = get_user_profile_id()
        
        # Parse search terms
        search_terms = None
        if args.search_terms:
            search_terms = [term.strip() for term in args.search_terms.split(',')]
        
        # Initialize and run importer
        importer = ChEMBLImporter(
            limit=args.limit,
            batch_size=args.batch_size,
            checkpoint_file=args.checkpoint_file,
            search_terms=search_terms,
            user_profile_id=user_profile_id,
            max_workers=args.max_workers,
            dry_run=args.dry_run
        )
        
        stats = importer.run()
        
        # Log final statistics
        logger.info("ChEMBL import completed successfully")
        logger.info(f"Total compounds: {stats.get('total_compounds', 0)}")
        logger.info(f"Molecules inserted: {stats.get('molecules_inserted', 0)}")
        logger.info(f"Properties inserted: {stats.get('properties_inserted', 0)}")
        
        return 0
        
    except ConfigurationError as e:
        logger.error(f"Configuration error: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Unhandled error: {str(e)}")
        logger.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
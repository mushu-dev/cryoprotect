#!/usr/bin/env python3
"""
Unified ChEMBL Import with Integrated Fixes

This script combines the ChEMBL import process with property calculation and
PubChem cross-reference resolution in a single, efficient workflow.

It provides:
- Official ChEMBL client integration
- Comprehensive molecular property calculation
- PubChem cross-reference resolution
- Direct Supabase connection
- Batch processing with transaction support
- Detailed reporting and checkpointing
"""

import os
import sys
import json
import time
import uuid
import signal
import atexit
import logging
import argparse
import traceback
import requests
import hashlib
import pickle
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple, Set
from urllib.parse import quote

# Import the official ChEMBL client
try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_CLIENT_AVAILABLE = True
except ImportError:
    CHEMBL_CLIENT_AVAILABLE = False
    print("WARNING: chembl_webresource_client not available. Install with: pip install chembl_webresource_client")

# Try importing RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, MolSurf, AllChem, DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    # Try to use mock_rdkit
    try:
        import mock_rdkit
        mock_rdkit.create_mock_rdkit()
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, MolSurf, AllChem, DataStructs
        RDKIT_AVAILABLE = True
        print("Using mock RDKit for basic property calculations")
    except ImportError:
        RDKIT_AVAILABLE = False
        print("WARNING: RDKit not available. Property calculation will be disabled.")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/unified_chembl_import.log')
    ]
)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)
os.makedirs('reports', exist_ok=True)
os.makedirs('checkpoints', exist_ok=True)
os.makedirs('cache', exist_ok=True)

logger = logging.getLogger(__name__)

# Database connection parameters (from environment variables)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

# Default checkpoint file
DEFAULT_CHECKPOINT_FILE = "checkpoints/unified_chembl_import_checkpoint.json"

# Cache settings
CACHE_DIR = "cache"
CACHE_EXPIRY = 7 * 24 * 60 * 60  # 7 days in seconds

# API Configuration
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/compound/"
REQUEST_DELAY = 0.5  # Time between API requests (seconds)
MAX_RETRIES = 3      # Maximum number of retries for API requests
CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"

# Property definitions
PROPERTY_DEFINITIONS = {
    'LogP': {
        'description': 'Octanol-water partition coefficient',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.MolLogP if RDKIT_AVAILABLE else None
    },
    'TPSA': {
        'description': 'Topological polar surface area',
        'data_type': 'numeric',
        'unit': 'Å²',
        'rdkit_func': MolSurf.TPSA if RDKIT_AVAILABLE else None
    },
    'Molecular Weight': {
        'description': 'The molecular weight of the compound',
        'data_type': 'numeric',
        'unit': 'g/mol',
        'rdkit_func': Descriptors.MolWt if RDKIT_AVAILABLE else None
    },
    'Heavy Atom Count': {
        'description': 'Number of non-hydrogen atoms',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.HeavyAtomCount if RDKIT_AVAILABLE else None
    },
    'Hydrogen Bond Donor Count': {
        'description': 'Number of hydrogen bond donors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHDonors if RDKIT_AVAILABLE else None
    },
    'Hydrogen Bond Acceptor Count': {
        'description': 'Number of hydrogen bond acceptors',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Lipinski.NumHAcceptors if RDKIT_AVAILABLE else None
    },
    'Rotatable Bond Count': {
        'description': 'Number of rotatable bonds',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.NumRotatableBonds if RDKIT_AVAILABLE else None
    },
    'Ring Count': {
        'description': 'Number of rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': Descriptors.RingCount if RDKIT_AVAILABLE else None
    },
    'Aromatic Ring Count': {
        'description': 'Number of aromatic rings',
        'data_type': 'numeric',
        'unit': '',
        'rdkit_func': lambda mol: rdMolDescriptors.CalcNumAromaticRings(mol) if RDKIT_AVAILABLE else None
    }
}

# Standard reference compounds to ensure they are always imported
STANDARD_REFERENCE_IDS = [
    "CHEMBL25",    # Aspirin
    "CHEMBL1118",  # Caffeine
    "CHEMBL1234",  # Glycerol (common cryoprotectant)
    "CHEMBL444",   # Glucose
    "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
    "CHEMBL9335",  # Dimethyl sulfoxide (DMSO, common cryoprotectant)
    "CHEMBL15151"  # Trehalose (common cryoprotectant)
]

# Cache management functions
def get_cache_key(url, params=None):
    """
    Generate a cache key from a URL and parameters.

    Args:
        url: The URL
        params: Optional parameters

    Returns:
        str: Cache key
    """
    key = url
    if params:
        if isinstance(params, dict):
            # Sort dictionary items for consistent keys
            key += "?" + "&".join(f"{k}={v}" for k, v in sorted(params.items()))
        else:
            key += str(params)

    # Create a hash of the key
    hashed_key = hashlib.md5(key.encode('utf-8')).hexdigest()
    return hashed_key

def save_to_cache(key, data):
    """
    Save data to cache.

    Args:
        key: Cache key
        data: Data to cache
    """
    cache_file = os.path.join(CACHE_DIR, f"{key}.pkl")

    try:
        with open(cache_file, 'wb') as f:
            pickle.dump({
                'data': data,
                'timestamp': time.time()
            }, f)
    except Exception as e:
        logger.warning(f"Error saving to cache: {str(e)}")

def load_from_cache(key):
    """
    Load data from cache if available and not expired.

    Args:
        key: Cache key

    Returns:
        The cached data or None if not available or expired
    """
    cache_file = os.path.join(CACHE_DIR, f"{key}.pkl")

    if not os.path.exists(cache_file):
        return None

    try:
        with open(cache_file, 'rb') as f:
            cache_entry = pickle.load(f)

            # Check if cache entry is expired
            if time.time() - cache_entry['timestamp'] > CACHE_EXPIRY:
                logger.debug(f"Cache expired for key {key}")
                return None

            return cache_entry['data']
    except Exception as e:
        logger.warning(f"Error loading from cache: {str(e)}")
        return None

def clear_expired_cache():
    """Clear expired cache entries."""
    try:
        now = time.time()
        for filename in os.listdir(CACHE_DIR):
            if filename.endswith('.pkl'):
                cache_file = os.path.join(CACHE_DIR, filename)

                try:
                    with open(cache_file, 'rb') as f:
                        cache_entry = pickle.load(f)

                        # Check if cache entry is expired
                        if now - cache_entry['timestamp'] > CACHE_EXPIRY:
                            os.remove(cache_file)
                            logger.debug(f"Removed expired cache file: {filename}")
                except Exception as e:
                    logger.warning(f"Error processing cache file {filename}: {str(e)}")
    except Exception as e:
        logger.warning(f"Error clearing expired cache: {str(e)}")

# Connect to database
def get_db_connection():
    """Get database connection using psycopg2."""
    import psycopg2
    from psycopg2.extras import RealDictCursor

    try:
        conn = psycopg2.connect(**DB_PARAMS)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {str(e)}")
        raise

# Common search terms for cryoprotectants (for backward compatibility)
DEFAULT_SEARCH_TERMS = [
    "cryoprotect",
    "cryopreservation",
    "antifreeze",
    "freeze protection",
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

# Molecular property filters for identifying potential cryoprotectants
CRYOPROTECTANT_PROPERTY_FILTERS = [
    # Small molecules with hydrogen bonding capability
    {
        'full_mwt__lte': 200.0,  # Molecular weight less than or equal to 200
        'hba__gte': 2,           # At least 2 hydrogen bond acceptors
        'hbd__gte': 1            # At least 1 hydrogen bond donor
    },
    # Medium-sized molecules with high hydrogen bonding capability
    {
        'full_mwt__lte': 350.0,  # Molecular weight less than or equal to 350
        'full_mwt__gte': 150.0,  # Molecular weight greater than or equal to 150
        'hba__gte': 3,           # At least 3 hydrogen bond acceptors
        'hbd__gte': 2            # At least 2 hydrogen bond donors
    },
    # Polyols and similar compounds
    {
        'full_mwt__lte': 400.0,  # Molecular weight less than or equal to 400
        'hba_lipinski__gte': 4,  # At least 4 hydrogen bond acceptors (Lipinski)
        'hbd_lipinski__gte': 3   # At least 3 hydrogen bond donors (Lipinski)
    },
    # Compounds with good water solubility (based on LogP)
    {
        'alogp__lte': 1.0,       # LogP less than or equal to 1.0
        'full_mwt__lte': 300.0,  # Molecular weight less than or equal to 300
        'hba__gte': 2            # At least 2 hydrogen bond acceptors
    }
]

class ChEMBLImportProgressTracker:
    """
    Tracks and manages progress for the ChEMBL import process.
    Provides checkpoint functionality for resumable imports.
    """
    
    def __init__(self, total_compounds, batch_size, checkpoint_file=DEFAULT_CHECKPOINT_FILE):
        """
        Initialize the progress tracker.
        
        Args:
            total_compounds: Total number of compounds to process
            batch_size: Batch size for processing
            checkpoint_file: File to save checkpoints to
        """
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

class PubChemResolver:
    """
    Resolves ChEMBL IDs to PubChem CIDs using multiple methods.
    """

    def __init__(self, conn):
        """
        Initialize the resolver.

        Args:
            conn: Database connection
        """
        self.conn = conn
        self.session = requests.Session()
        self.last_request_time = 0

        # Clean up expired cache entries
        clear_expired_cache()

    def _api_request(self, url: str, params: Dict[str, str] = None) -> Optional[Dict[str, Any]]:
        """
        Make a request to the PubChem API with rate limiting and retries.

        Args:
            url: URL to request
            params: Dictionary of URL parameters

        Returns:
            JSON response as dictionary or None if request failed
        """
        # Check cache first
        cache_key = get_cache_key(url, params)
        cached_data = load_from_cache(cache_key)

        if cached_data:
            logger.debug(f"Using cached data for URL: {url}")
            return cached_data

        # Rate limiting
        now = time.time()
        time_since_last_request = now - self.last_request_time
        if time_since_last_request < REQUEST_DELAY:
            time.sleep(REQUEST_DELAY - time_since_last_request)

        self.last_request_time = time.time()

        # Make request with retries
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()

                # Parse the response
                data = response.json()

                # Save to cache
                save_to_cache(cache_key, data)

                return data
            except (requests.RequestException, json.JSONDecodeError) as e:
                logger.warning(f"API request failed (attempt {attempt}/{MAX_RETRIES}): {str(e)}")

                if attempt < MAX_RETRIES:
                    # Exponential backoff
                    sleep_time = REQUEST_DELAY * (2 ** (attempt - 1))
                    time.sleep(sleep_time)
                else:
                    logger.error(f"API request failed after {MAX_RETRIES} attempts: {str(e)}")
                    return None
    
    def find_pubchem_by_chembl_id(self, chembl_id: str) -> Optional[int]:
        """
        Find PubChem CID by ChEMBL ID using direct PubChem lookup.

        Args:
            chembl_id: ChEMBL ID

        Returns:
            PubChem CID as integer or None if not found
        """
        # Format changed to match PubChem API expectations - ChEMBL requires "ChEMBL:" prefix
        formatted_chembl_id = chembl_id.replace("CHEMBL", "ChEMBL:")
        url = f"{PUBCHEM_API_BASE}/compound/xref/RegistryID/{quote(formatted_chembl_id)}/cids/JSON"

        try:
            # Log the URL for debugging purposes
            logger.debug(f"Trying PubChem lookup with URL: {url}")

            response = self._api_request(url)

            if not response:
                # Try alternative endpoint format as fallback
                fallback_url = f"{PUBCHEM_API_BASE}/compound/xref/SourceName/ChEMBL/SourceID/{quote(chembl_id.replace('CHEMBL', ''))}/cids/JSON"
                logger.debug(f"Trying fallback URL: {fallback_url}")
                response = self._api_request(fallback_url)

            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID

            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for ChEMBL ID {chembl_id}: {str(e)}")
            return None
    
    def find_pubchem_by_inchikey(self, inchikey: str) -> Optional[int]:
        """
        Find PubChem CID by InChIKey.
        
        Args:
            inchikey: InChIKey
            
        Returns:
            PubChem CID as integer or None if not found
        """
        url = f"{PUBCHEM_API_BASE}/compound/inchikey/{quote(inchikey)}/cids/JSON"
        
        try:
            response = self._api_request(url)
            
            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID
            
            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for InChIKey {inchikey}: {str(e)}")
            return None
    
    def find_pubchem_by_smiles(self, smiles: str) -> Optional[int]:
        """
        Find PubChem CID by SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            PubChem CID as integer or None if not found
        """
        url = f"{PUBCHEM_API_BASE}/compound/smiles/cids/JSON"
        
        try:
            response = self._api_request(url, params={'smiles': smiles})
            
            if response and 'IdentifierList' in response and 'CID' in response['IdentifierList']:
                cids = response['IdentifierList']['CID']
                if cids:
                    return cids[0]  # Return the first CID
            
            return None
        except Exception as e:
            logger.error(f"Error finding PubChem CID for SMILES {smiles}: {str(e)}")
            return None
    
    def update_pubchem_cid(self, molecule_id: str, pubchem_cid: int, method: str) -> bool:
        """
        Update the PubChem CID for a molecule.
        
        Args:
            molecule_id: Molecule ID
            pubchem_cid: PubChem CID
            method: Method used to find the CID
            
        Returns:
            True if update was successful, False otherwise
        """
        try:
            # Construct query - don't try to update any generated columns
            query = """
            UPDATE molecules
            SET pubchem_cid = %s,
                updated_at = NOW(),
                modification_history = COALESCE(modification_history, '[]'::jsonb) || %s::jsonb
            WHERE id = %s;
            """

            # Create modification history entry
            history_entry = json.dumps([{
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "operation": "update",
                "field": "pubchem_cid",
                "value": pubchem_cid,
                "method": method
            }])

            # Execute query
            with self.conn.cursor() as cursor:
                cursor.execute(
                    query,
                    (
                        pubchem_cid,
                        history_entry,
                        molecule_id
                    )
                )
            
            # Commit changes
            self.conn.commit()
            
            logger.info(f"Updated PubChem CID for molecule {molecule_id} to {pubchem_cid} (method: {method})")
            return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error updating PubChem CID for molecule {molecule_id}: {str(e)}")
            return False
    
    def resolve_pubchem_cid(self, molecule_id: str, chembl_id: str, inchikey: str, smiles: str) -> Optional[Tuple[int, str]]:
        """
        Resolve PubChem CID for a molecule using all available methods.

        Args:
            molecule_id: Molecule ID
            chembl_id: ChEMBL ID
            inchikey: InChIKey
            smiles: SMILES string

        Returns:
            Tuple of (PubChem CID, method) or None if not found
        """
        # Changed order to prioritize more reliable methods first

        # Method 1: InChIKey lookup - most reliable
        if inchikey:
            logger.info(f"Trying to find PubChem CID for {molecule_id} using InChIKey: {inchikey}")
            pubchem_cid = self.find_pubchem_by_inchikey(inchikey)
            if pubchem_cid:
                return (pubchem_cid, "inchikey")

        # Method 2: SMILES lookup
        if smiles:
            logger.info(f"Trying to find PubChem CID for {molecule_id} using SMILES")
            pubchem_cid = self.find_pubchem_by_smiles(smiles)
            if pubchem_cid:
                return (pubchem_cid, "smiles")

        # Method 3: Direct ChEMBL ID lookup - least reliable based on our tests
        if chembl_id:
            logger.info(f"Trying to find PubChem CID for {molecule_id} using ChEMBL ID: {chembl_id}")
            pubchem_cid = self.find_pubchem_by_chembl_id(chembl_id)
            if pubchem_cid:
                return (pubchem_cid, "chembl_id")

        logger.warning(f"Could not find PubChem CID for molecule {molecule_id} using any method")
        return None

class PropertyCalculator:
    """
    Calculates molecular properties using RDKit.
    """
    
    def __init__(self, conn):
        """
        Initialize the property calculator.
        
        Args:
            conn: Database connection
        """
        self.conn = conn
        self.property_types = {}
        
        # Load property types
        self._load_property_types()
    
    def _load_property_types(self):
        """Load property types from the database."""
        try:
            query = "SELECT id, name, data_type FROM property_types;"
            with self.conn.cursor() as cursor:
                cursor.execute(query)
                results = cursor.fetchall()
                
                for row in results:
                    self.property_types[row[1]] = {
                        'id': row[0],
                        'data_type': row[2]
                    }
                
                logger.info(f"Loaded {len(self.property_types)} property types from database")
                
                # Check if all properties exist, create missing ones
                self._ensure_property_types_exist()
        except Exception as e:
            logger.error(f"Error loading property types: {str(e)}")
            raise
    
    def _ensure_property_types_exist(self):
        """Ensure all required property types exist in the database."""
        for prop_name, definition in PROPERTY_DEFINITIONS.items():
            if prop_name not in self.property_types:
                logger.info(f"Creating missing property type: {prop_name}")
                self._create_property_type(
                    prop_name, 
                    definition['description'],
                    definition['data_type'],
                    definition['unit']
                )
    
    def _create_property_type(self, name: str, description: str, data_type: str, unit: str):
        """
        Create a new property type in the database.
        
        Args:
            name: Property type name
            description: Property type description
            data_type: Data type (numeric, text, boolean)
            unit: Unit of measurement
        """
        try:
            query = """
            INSERT INTO property_types 
            (id, name, description, data_type, units, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            RETURNING id, name, data_type;
            """
            
            with self.conn.cursor() as cursor:
                prop_id = str(uuid.uuid4())
                cursor.execute(query, (prop_id, name, description, data_type, unit))
                result = cursor.fetchone()
                self.conn.commit()
                
                # Add to property types dict
                self.property_types[name] = {
                    'id': result[0],
                    'data_type': result[2]
                }
                
                logger.info(f"Created property type: {name} (ID: {prop_id})")
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error creating property type {name}: {str(e)}")
            raise
    
    def calculate_rdkit_properties(self, smiles: str, inchi: str = None) -> Optional[Dict[str, Any]]:
        """
        Calculate molecular properties using RDKit.
        
        Args:
            smiles: SMILES string
            inchi: Optional InChI string
            
        Returns:
            Dictionary of calculated properties or None if calculation failed
        """
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. Cannot calculate properties.")
            return None
        
        try:
            # Try to get a valid RDKit molecule
            mol = None
            
            # First try from SMILES
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
            
            # If SMILES failed, try from InChI
            if mol is None and inchi:
                mol = Chem.MolFromInchi(inchi)
            
            # If both failed, return None
            if mol is None:
                logger.warning(f"Could not convert to RDKit molecule: SMILES='{smiles}', InChI='{inchi}'")
                return None
            
            # Calculate properties
            properties = {}
            
            # Add standard molecular properties
            for prop_name, definition in PROPERTY_DEFINITIONS.items():
                try:
                    if definition['rdkit_func'] is not None:
                        value = definition['rdkit_func'](mol)
                        properties[prop_name] = value
                except Exception as e:
                    logger.warning(f"Error calculating {prop_name}: {str(e)}")
            
            # Add calculated molecular formula
            try:
                properties['Molecular Formula'] = AllChem.CalcMolFormula(mol)
            except Exception as e:
                logger.warning(f"Error calculating molecular formula: {str(e)}")
            
            return properties
        except Exception as e:
            logger.error(f"Error calculating properties: {str(e)}")
            return None
    
    def get_existing_properties(self, molecule_id: str) -> Dict[str, Any]:
        """
        Get existing properties for a molecule.
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary mapping property types to values
        """
        query = """
        SELECT property_type, numeric_value, text_value, boolean_value
        FROM molecular_properties
        WHERE molecule_id = %s;
        """
        
        try:
            with self.conn.cursor() as cursor:
                cursor.execute(query, (molecule_id,))
                results = cursor.fetchall()
                
                properties = {}
                for row in results:
                    prop_type = row[0]
                    
                    # Get the appropriate value based on data type
                    if row[1] is not None:
                        properties[prop_type] = row[1]
                    elif row[2] is not None:
                        properties[prop_type] = row[2]
                    elif row[3] is not None:
                        properties[prop_type] = row[3]
                
                return properties
        except Exception as e:
            logger.error(f"Error fetching properties for molecule {molecule_id}: {str(e)}")
            return {}
    
    def update_molecule_properties(self, molecule_id: str, properties: Dict[str, Any], existing_properties: Dict[str, Any] = None, user_id: str = None) -> bool:
        """
        Update molecular properties for a molecule.
        
        Args:
            molecule_id: Molecule ID
            properties: Dictionary of property name to value
            existing_properties: Optional dictionary of existing properties
            user_id: Optional user ID for created_by field
            
        Returns:
            True if update was successful, False otherwise
        """
        if not properties:
            logger.warning(f"No properties to update for molecule {molecule_id}")
            return False
        
        # Get existing properties if not provided
        if existing_properties is None:
            existing_properties = self.get_existing_properties(molecule_id)
        
        try:
            # Start transaction
            with self.conn:
                with self.conn.cursor() as cursor:
                    # Update molecular_properties table
                    for prop_name, value in properties.items():
                        # Skip if property already exists
                        if prop_name in existing_properties:
                            continue
                        
                        # Skip if property type not defined
                        if prop_name not in self.property_types:
                            logger.warning(f"Property type {prop_name} not found in database")
                            continue
                        
                        # Get property type details
                        prop_type_id = self.property_types[prop_name]['id']
                        data_type = self.property_types[prop_name]['data_type']
                        
                        # Determine which value field to use
                        numeric_value = None
                        text_value = None
                        boolean_value = None
                        
                        if data_type == 'numeric':
                            numeric_value = value
                        elif data_type == 'text':
                            text_value = str(value)
                        elif data_type == 'boolean':
                            boolean_value = bool(value)
                        
                        # Insert property
                        query = """
                        INSERT INTO molecular_properties 
                        (id, molecule_id, property_type_id, property_type, property_name,
                         numeric_value, text_value, boolean_value, 
                         source, created_at, updated_at, created_by)
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW(), %s);
                        """
                        
                        cursor.execute(
                            query, 
                            (
                                str(uuid.uuid4()),  # id
                                molecule_id,        # molecule_id
                                prop_type_id,       # property_type_id
                                prop_name,          # property_type
                                prop_name,          # property_name
                                numeric_value,      # numeric_value
                                text_value,         # text_value
                                boolean_value,      # boolean_value
                                'Calculated',       # source
                                user_id             # created_by
                            )
                        )
                    
                    # Update the JSONB properties field
                    jsonb_properties = {}
                    
                    # Combine existing and new properties
                    for prop_name, value in {**existing_properties, **properties}.items():
                        jsonb_properties[prop_name] = value
                    
                    # Convert to JSONB
                    jsonb_str = json.dumps(jsonb_properties)
                    
                    # Update molecule
                    query = """
                    UPDATE molecules
                    SET properties = %s, updated_at = NOW()
                    WHERE id = %s;
                    """
                    
                    cursor.execute(query, (jsonb_str, molecule_id))
                    
                    # Commit transaction
                    self.conn.commit()
                    
                    logger.info(f"Updated properties for molecule {molecule_id}")
                    return True
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Error updating properties for molecule {molecule_id}: {str(e)}")
            return False

def fetch_compound_by_id(chembl_id):
    """
    Fetch a compound by ChEMBL ID.

    Args:
        chembl_id: ChEMBL ID to fetch

    Returns:
        dict: Compound data or None if not found
    """
    if not CHEMBL_CLIENT_AVAILABLE:
        logger.error("ChEMBL client not available")
        return None

    # Check cache first
    cache_key = get_cache_key(f"{CHEMBL_API_BASE}/molecule/{chembl_id}")
    cached_compound = load_from_cache(cache_key)

    if cached_compound:
        logger.info(f"Using cached data for compound {chembl_id}")
        return cached_compound

    try:
        logger.info(f"Fetching compound with ChEMBL ID: {chembl_id}")
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

        # Save to cache
        save_to_cache(cache_key, compound)

        return compound
    except Exception as e:
        logger.error(f"Error fetching compound {chembl_id}: {str(e)}")
        return None

def fetch_compounds_by_property(property_filters=None, limit=1000, use_existing_implementation=True):
    """
    Fetch compounds from ChEMBL based on molecular property filters.

    Args:
        property_filters: List of property filter dictionaries (used only if use_existing_implementation=False)
        limit: Maximum number of compounds to fetch
        use_existing_implementation: Whether to use the existing implementation from chembl_search_utils.py

    Returns:
        list: List of compounds
    """
    if not CHEMBL_CLIENT_AVAILABLE:
        logger.error("ChEMBL client not available")
        return []

    logger.info(f"Fetching compounds from ChEMBL by molecular properties (limit: {limit})")

    # First check if we can use the existing implementation in chembl_search_utils.py
    if use_existing_implementation:
        try:
            # Import the function from chembl_search_utils.py
            from chembl_search_utils import find_potential_cryoprotectants

            logger.info("Using existing property-based search implementation from chembl_search_utils.py")

            # Fetch compounds using the existing implementation
            compounds = find_potential_cryoprotectants(limit=limit)
            logger.info(f"Found {len(compounds)} compounds using property-based search")

            # Transform the compounds to the format expected by the rest of the code
            transformed_compounds = []
            compounds_by_id = {}  # Used to deduplicate

            # First, fetch standard reference compounds
            for chembl_id in STANDARD_REFERENCE_IDS:
                compound = fetch_compound_by_id(chembl_id)
                if compound and chembl_id not in compounds_by_id:
                    compounds_by_id[chembl_id] = compound
                    transformed_compounds.append(compound)
                    logger.info(f"Added reference compound: {chembl_id}")

            # Add compounds from property-based search
            for compound in compounds:
                chembl_id = compound.get('molecule_chembl_id')

                # Skip if already in results or reached limit
                if chembl_id in compounds_by_id or len(compounds_by_id) >= limit:
                    continue

                # Transform to the format expected by the rest of the code
                full_details = fetch_compound_by_id(chembl_id)

                if full_details:
                    compounds_by_id[chembl_id] = full_details
                    transformed_compounds.append(full_details)

                    # Log progress
                    if len(compounds_by_id) % 10 == 0:
                        logger.info(f"Transformed {len(compounds_by_id)}/{limit} compounds")

            logger.info(f"Total compounds fetched by property criteria: {len(transformed_compounds)}")
            return transformed_compounds[:limit]

        except ImportError:
            logger.warning("Could not import find_potential_cryoprotectants from chembl_search_utils.py")
            logger.info("Falling back to built-in property-based search implementation")

    # If the import failed or use_existing_implementation=False, use the built-in implementation
    # Use default property filters if none provided
    if property_filters is None:
        property_filters = CRYOPROTECTANT_PROPERTY_FILTERS

    # Initialize ChEMBL client
    molecule = new_client.molecule

    all_compounds = []
    compounds_by_id = {}  # Used to deduplicate

    # First, fetch standard reference compounds
    for chembl_id in STANDARD_REFERENCE_IDS:
        compound = fetch_compound_by_id(chembl_id)
        if compound:
            compounds_by_id[chembl_id] = compound
            all_compounds.append(compound)
            logger.info(f"Added reference compound: {chembl_id}")

    logger.info(f"Fetched {len(all_compounds)} reference compounds")

    # Process each property filter
    for filter_index, property_filter in enumerate(property_filters):
        if len(compounds_by_id) >= limit:
            break

        try:
            logger.info(f"Searching with property filter #{filter_index+1}: {property_filter}")

            # Search for compounds matching the property filter
            results = molecule.filter(**property_filter).only(
                'molecule_chembl_id',
                'pref_name',
                'molecule_structures'
            )

            found_count = 0

            # Process results
            for compound in results:
                chembl_id = compound.get('molecule_chembl_id')

                # Skip if already in results or reached limit
                if chembl_id in compounds_by_id or len(compounds_by_id) >= limit:
                    continue

                # Skip if no structures available
                if not compound.get('molecule_structures'):
                    logger.debug(f"Skipping {chembl_id}: Missing molecular structures")
                    continue

                # Get full compound details
                try:
                    full_details = molecule.get(chembl_id)

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
                    compounds_by_id[chembl_id] = full_details
                    all_compounds.append(full_details)
                    found_count += 1

                    # Log progress
                    if len(compounds_by_id) % 10 == 0:
                        logger.info(f"Fetched {len(compounds_by_id)}/{limit} compounds")

                except Exception as e:
                    logger.warning(f"Error fetching details for {chembl_id}: {str(e)}")
                    continue

            logger.info(f"Found {found_count} compounds for property filter #{filter_index+1}")

        except Exception as e:
            logger.error(f"Error processing property filter #{filter_index+1}: {str(e)}")
            continue

    logger.info(f"Total compounds fetched by property criteria: {len(all_compounds)}")
    return all_compounds[:limit]

def fetch_cryoprotectant_compounds(search_terms=None, limit=1000, use_properties=True):
    """
    Fetch cryoprotectant compounds from ChEMBL.

    Args:
        search_terms: List of search terms for keyword search (if use_properties is False)
        limit: Maximum number of compounds to fetch
        use_properties: Whether to use property-based search (True) or keyword search (False)

    Returns:
        list: List of compounds
    """
    if not CHEMBL_CLIENT_AVAILABLE:
        logger.error("ChEMBL client not available")
        return []

    # Use property-based search by default (new method)
    if use_properties:
        return fetch_compounds_by_property(property_filters=CRYOPROTECTANT_PROPERTY_FILTERS, limit=limit)

    # Legacy keyword search method
    logger.info(f"Fetching cryoprotectant compounds from ChEMBL by keywords (limit: {limit})")

    # Use default search terms if none provided
    if search_terms is None:
        search_terms = DEFAULT_SEARCH_TERMS

    # Initialize ChEMBL client
    molecule = new_client.molecule

    all_compounds = []
    compounds_by_id = {}  # Used to deduplicate

    # First, fetch standard reference compounds
    for chembl_id in STANDARD_REFERENCE_IDS:
        compound = fetch_compound_by_id(chembl_id)
        if compound:
            compounds_by_id[chembl_id] = compound
            all_compounds.append(compound)
            logger.info(f"Added reference compound: {chembl_id}")

    logger.info(f"Fetched {len(all_compounds)} reference compounds")

    # Process each search term
    for term in search_terms:
        if len(compounds_by_id) >= limit:
            break

        try:
            logger.info(f"Searching for keyword '{term}'")

            # Search for compounds matching the term
            results = molecule.filter(
                pref_name__icontains=term
            ).only(
                'molecule_chembl_id',
                'pref_name',
                'molecule_structures'
            )

            found_count = 0

            # Process results
            for compound in results:
                chembl_id = compound.get('molecule_chembl_id')

                # Skip if already in results or reached limit
                if chembl_id in compounds_by_id or len(compounds_by_id) >= limit:
                    continue

                # Skip if no structures available
                if not compound.get('molecule_structures'):
                    logger.debug(f"Skipping {chembl_id}: Missing molecular structures")
                    continue

                # Get full compound details
                try:
                    full_details = molecule.get(chembl_id)

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
                    compounds_by_id[chembl_id] = full_details
                    all_compounds.append(full_details)
                    found_count += 1

                    # Log progress
                    if len(compounds_by_id) % 10 == 0:
                        logger.info(f"Fetched {len(compounds_by_id)}/{limit} compounds")

                except Exception as e:
                    logger.warning(f"Error fetching details for {chembl_id}: {str(e)}")
                    continue

            logger.info(f"Found {found_count} compounds for keyword '{term}'")

        except Exception as e:
            logger.error(f"Error processing search term '{term}': {str(e)}")
            continue

    logger.info(f"Total compounds fetched by keywords: {len(all_compounds)}")
    return all_compounds[:limit]

def transform_chembl_to_molecule(compound, user_id):
    """
    Transform ChEMBL compound data to match our molecule table schema.
    
    Args:
        compound: ChEMBL compound data
        user_id: User ID for created_by field
        
    Returns:
        dict: Molecule data for insertion
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
    
    # Calculate initial properties for JSONB field
    properties = {}
    for prop_key, prop_value in mol_props.items():
        if prop_value is not None:
            properties[prop_key] = prop_value

    # Ensure common properties are explicitly set if available in mol_props
    if mol_props.get('full_mwt'):
        properties['Molecular Weight'] = mol_props.get('full_mwt')
    if mol_props.get('alogp'):
        properties['LogP'] = mol_props.get('alogp')
    if mol_props.get('hba'):
        properties['Hydrogen Bond Acceptor Count'] = mol_props.get('hba')
    if mol_props.get('hbd'):
        properties['Hydrogen Bond Donor Count'] = mol_props.get('hbd')
    if mol_props.get('rtb'):
        properties['Rotatable Bond Count'] = mol_props.get('rtb')
    
    # Transform to our schema
    return {
        "name": name,
        "smiles": structures.get('canonical_smiles'),
        "inchi": structures.get('standard_inchi'),
        "inchikey": structures.get('standard_inchi_key'),
        "formula": mol_props.get('full_molformula'),
        "molecular_weight": mol_props.get('full_mwt'),
        "pubchem_cid": None,  # Will be populated later
        "chembl_id": compound.get('molecule_chembl_id'),
        "created_by": user_id,
        "data_source": f"ChEMBL v{chembl_version} ID: {compound.get('molecule_chembl_id')}",
        "version": 1,
        "properties": json.dumps(properties),
        "modification_history": json.dumps([{
            "timestamp": datetime.now().isoformat(),
            "action": "created",
            "user_id": user_id
        }])
    }

def upsert_reference_user(conn):
    """
    Get a reference user ID for creating molecules.
    Does not depend on user_profiles table existing.

    Args:
        conn: Database connection

    Returns:
        str: User ID
    """
    # Use a fixed UUID for the system user to ensure consistency
    # This approach doesn't require the user_profiles table to exist
    system_user_id = "00000000-0000-0000-0000-000000000001"
    logger.info(f"Using fixed system user ID: {system_user_id}")
    return system_user_id

def import_compounds_to_database(compounds, batch_size=50, dry_run=False):
    """
    Import compounds to database with property calculation and cross-reference resolution.
    
    Args:
        compounds: List of ChEMBL compounds
        batch_size: Number of compounds to process in each batch
        dry_run: If True, don't actually insert data, just simulate
        
    Returns:
        dict: Import statistics
    """
    logger.info(f"Importing {len(compounds)} compounds to database")
    
    # Initialize connections
    conn = None
    
    # Statistics
    stats = {
        "total_compounds": len(compounds),
        "processed": 0,
        "molecules_inserted": 0,
        "molecules_skipped": 0,
        "properties_inserted": 0,
        "errors": 0,
        "pubchem_cross_refs": {
            "resolved_count": 0,
            "resolved_by_method": {
                "chembl_id": 0,
                "inchikey": 0,
                "smiles": 0
            }
        }
    }
    
    # Exit early if there are no compounds to process
    if len(compounds) == 0:
        logger.info("No compounds to import")
        return stats
    
    # Process compounds
    try:
        # Get database connection
        conn = get_db_connection()
        
        # Get reference user
        user_id = upsert_reference_user(conn)
        
        # Initialize property calculator and PubChem resolver
        property_calculator = PropertyCalculator(conn)
        pubchem_resolver = PubChemResolver(conn)
        
        # Initialize progress tracker
        progress_tracker = ChEMBLImportProgressTracker(
            total_compounds=len(compounds),
            batch_size=batch_size
        )
        
        # Process compounds in batches
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            batch_num = i//batch_size + 1
            total_batches = (len(compounds)-1)//batch_size + 1
            
            # Update progress tracker
            progress_tracker.start_batch(batch_num-1, len(batch))
            logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} compounds)")
            
            # Track batch statistics
            batch_processed = 0
            batch_imported = 0
            batch_skipped = 0
            batch_errors = 0
            
            # Skip batch if in dry run mode
            if dry_run:
                # Simulate batch processing
                for j, compound in enumerate(batch):
                    logger.info(f"DRY RUN: Would process compound {j+1}/{len(batch)}: {compound.get('molecule_chembl_id')}")
                    batch_processed += 1
                    batch_imported += 1
                
                # Update progress tracker
                progress_tracker.end_batch(
                    processed=batch_processed,
                    imported=batch_imported,
                    skipped=batch_skipped,
                    errors=batch_errors
                )
                
                # Update statistics
                stats["processed"] += batch_processed
                stats["molecules_inserted"] += batch_imported
                stats["molecules_skipped"] += batch_skipped
                stats["errors"] += batch_errors
                
                continue
            
            # Process each compound in the batch
            for compound in batch:
                chembl_id = compound.get('molecule_chembl_id')
                
                try:
                    # Transform compound to molecule
                    molecule_data = transform_chembl_to_molecule(compound, user_id)
                    
                    # Check if molecule with InChIKey already exists
                    inchikey = molecule_data.get('inchikey')
                    if inchikey:
                        query = "SELECT id FROM molecules WHERE inchikey = %s;"
                        with conn.cursor() as cursor:
                            cursor.execute(query, (inchikey,))
                            result = cursor.fetchone()
                            
                            if result:
                                # Molecule already exists
                                logger.info(f"Skipping molecule with InChIKey {inchikey}: Already exists")
                                batch_skipped += 1
                                stats["molecules_skipped"] += 1
                                progress_tracker.add_skipped(chembl_id, "Already exists")
                                continue
                    
                    # Insert molecule
                    query = """
                    INSERT INTO molecules (
                        name, smiles, inchi, inchikey, formula, molecular_weight,
                        pubchem_cid, chembl_id, created_by, data_source, version,
                        properties, modification_history, created_at, updated_at
                    )
                    VALUES (
                        %s, %s, %s, %s, %s, %s, 
                        %s, %s, %s, %s, %s,
                        %s, %s, NOW(), NOW()
                    )
                    RETURNING id;
                    """
                    
                    with conn.cursor() as cursor:
                        cursor.execute(
                            query,
                            (
                                molecule_data.get('name'),
                                molecule_data.get('smiles'),
                                molecule_data.get('inchi'),
                                molecule_data.get('inchikey'),
                                molecule_data.get('formula'),
                                molecule_data.get('molecular_weight'),
                                molecule_data.get('pubchem_cid'),
                                molecule_data.get('chembl_id'),
                                molecule_data.get('created_by'),
                                molecule_data.get('data_source'),
                                molecule_data.get('version'),
                                molecule_data.get('properties'),
                                molecule_data.get('modification_history')
                            )
                        )
                        result = cursor.fetchone()
                        
                        if result:
                            molecule_id = result[0]
                            batch_imported += 1
                            stats["molecules_inserted"] += 1
                            
                            # Process properties
                            smiles = molecule_data.get('smiles')
                            inchi = molecule_data.get('inchi')
                            
                            # Calculate properties using RDKit
                            calculated_properties = property_calculator.calculate_rdkit_properties(smiles, inchi)
                            
                            if calculated_properties:
                                # Get existing properties
                                existing_properties = property_calculator.get_existing_properties(molecule_id)
                                
                                # Update properties
                                success = property_calculator.update_molecule_properties(
                                    molecule_id, 
                                    calculated_properties,
                                    existing_properties,
                                    user_id
                                )
                                
                                if success:
                                    stats["properties_inserted"] += len(calculated_properties)
                                    logger.info(f"Successfully updated {len(calculated_properties)} properties for molecule {molecule_data.get('name')}")
                            
                            # Resolve PubChem CID
                            pubchem_result = pubchem_resolver.resolve_pubchem_cid(
                                molecule_id,
                                molecule_data.get('chembl_id'),
                                molecule_data.get('inchikey'),
                                molecule_data.get('smiles')
                            )
                            
                            if pubchem_result:
                                pubchem_cid, method = pubchem_result
                                success = pubchem_resolver.update_pubchem_cid(
                                    molecule_id,
                                    pubchem_cid,
                                    method
                                )
                                
                                if success:
                                    stats["pubchem_cross_refs"]["resolved_count"] += 1
                                    stats["pubchem_cross_refs"]["resolved_by_method"][method] += 1
                                    logger.info(f"Successfully updated PubChem CID for molecule {molecule_data.get('name')} to {pubchem_cid} (method: {method})")
                        else:
                            logger.error(f"Failed to insert molecule: {molecule_data.get('name')}")
                            batch_errors += 1
                            stats["errors"] += 1
                except Exception as e:
                    logger.error(f"Error processing compound {chembl_id}: {str(e)}")
                    batch_errors += 1
                    stats["errors"] += 1
                    progress_tracker.add_error(f"Error processing compound {chembl_id}: {str(e)}")
                
                batch_processed += 1
            
            # Update progress tracker
            progress_tracker.end_batch(
                processed=batch_processed,
                imported=batch_imported,
                skipped=batch_skipped,
                errors=batch_errors
            )
            
            # Update statistics
            stats["processed"] += batch_processed
            conn.commit()
            
            # Sleep to avoid overloading the database
            time.sleep(0.5)
    
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error importing compounds: {str(e)}")
        stats["errors"] += 1
    finally:
        if conn:
            conn.close()
    
    return stats

def verify_imported_data(limit=None):
    """
    Verify the quality of imported data.
    
    Args:
        limit: Optional limit on number of molecules to check
        
    Returns:
        dict: Verification results
    """
    logger.info("Verifying imported data")
    
    # Initialize connections
    conn = None
    
    try:
        # Get database connection
        conn = get_db_connection()
        
        # Get molecules with ChEMBL IDs
        limit_clause = f"LIMIT {limit}" if limit else ""
        
        query = f"""
        SELECT id, name, chembl_id, pubchem_cid, smiles, inchi, inchikey,
               properties, updated_at
        FROM molecules
        WHERE chembl_id IS NOT NULL
        ORDER BY updated_at DESC
        {limit_clause};
        """
        
        with conn.cursor() as cursor:
            cursor.execute(query)
            molecules = cursor.fetchall()
        
        if not molecules:
            logger.error("No ChEMBL molecules found in the database")
            return {
                "error": "No ChEMBL molecules found",
                "timestamp": datetime.now().isoformat()
            }
        
        logger.info(f"Found {len(molecules)} ChEMBL molecules in the database")
        
        # Expected properties
        expected_properties = {
            'LogP', 'TPSA', 'Molecular Weight', 'Heavy Atom Count',
            'Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count',
            'Rotatable Bond Count', 'Ring Count', 'Aromatic Ring Count'
        }
        
        # Initialize results
        property_results = {
            "total_molecules": len(molecules),
            "molecules_with_properties": 0,
            "molecules_missing_properties": 0,
            "molecules_with_jsonb_properties": 0,
            "molecules_missing_jsonb_properties": 0,
            "property_coverage": {},
            "molecules_by_property_count": {}
        }
        
        cross_ref_results = {
            "total_molecules": len(molecules),
            "molecules_with_pubchem_cid": 0,
            "molecules_missing_pubchem_cid": 0,
            "molecules_with_inchikey": 0,
            "molecules_without_inchikey": 0
        }
        
        # Initialize property coverage
        for prop in expected_properties:
            property_results["property_coverage"][prop] = 0
        
        # Initialize molecules by property count
        for i in range(len(expected_properties) + 1):
            property_results["molecules_by_property_count"][i] = 0
        
        # Process each molecule
        for molecule in molecules:
            # Get properties from the JSONB field
            jsonb_properties = molecule[7] or {}
            
            if isinstance(jsonb_properties, str):
                try:
                    jsonb_properties = json.loads(jsonb_properties)
                except json.JSONDecodeError:
                    jsonb_properties = {}
            
            # Check properties in database tables
            query = """
            SELECT property_type
            FROM molecular_properties
            WHERE molecule_id = %s;
            """
            
            with conn.cursor() as cursor:
                cursor.execute(query, (molecule[0],))
                db_properties = {row[0] for row in cursor.fetchall()}
                
                # Count properties found in database tables
                db_property_count = sum(1 for prop in expected_properties if prop in db_properties)
                property_results["molecules_by_property_count"][db_property_count] += 1
                
                if db_property_count == len(expected_properties):
                    property_results["molecules_with_properties"] += 1
                else:
                    property_results["molecules_missing_properties"] += 1
                
                # Update property coverage
                for prop in expected_properties:
                    if prop in db_properties:
                        property_results["property_coverage"][prop] += 1
                
                # Check properties in JSONB field
                jsonb_property_count = sum(1 for prop in expected_properties if prop in jsonb_properties)
                
                if jsonb_property_count == len(expected_properties):
                    property_results["molecules_with_jsonb_properties"] += 1
                else:
                    property_results["molecules_missing_jsonb_properties"] += 1
            
            # Check cross-references
            if molecule[3]:  # pubchem_cid
                cross_ref_results["molecules_with_pubchem_cid"] += 1
            else:
                cross_ref_results["molecules_missing_pubchem_cid"] += 1
            
            if molecule[6]:  # inchikey
                cross_ref_results["molecules_with_inchikey"] += 1
            else:
                cross_ref_results["molecules_without_inchikey"] += 1
        
        # Calculate percentages
        for prop in expected_properties:
            count = property_results["property_coverage"][prop]
            percentage = (count / property_results["total_molecules"]) * 100
            property_results["property_coverage"][prop] = {
                "count": count,
                "percentage": percentage
            }
        
        property_results["property_completeness_percentage"] = (property_results["molecules_with_properties"] / property_results["total_molecules"]) * 100
        property_results["jsonb_property_completeness_percentage"] = (property_results["molecules_with_jsonb_properties"] / property_results["total_molecules"]) * 100
        
        cross_ref_results["pubchem_cross_reference_percentage"] = (cross_ref_results["molecules_with_pubchem_cid"] / cross_ref_results["total_molecules"]) * 100
        cross_ref_results["inchikey_coverage_percentage"] = (cross_ref_results["molecules_with_inchikey"] / cross_ref_results["total_molecules"]) * 100
        
        # Combine results
        results = {
            "timestamp": datetime.now().isoformat(),
            "molecule_count": len(molecules),
            "property_completeness": property_results,
            "pubchem_cross_references": cross_ref_results
        }
        
        return results
    
    except Exception as e:
        logger.error(f"Error verifying imported data: {str(e)}")
        return {
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }
    finally:
        if conn:
            conn.close()

def main():
    """Main entry point."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Unified ChEMBL Import with Integrated Fixes")
    parser.add_argument("--limit", type=int, default=2500, help="Maximum number of compounds to import")
    parser.add_argument("--batch-size", type=int, default=50, help="Batch size for database operations")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually insert data, just simulate")
    parser.add_argument("--verify-only", action="store_true", help="Only verify existing data, don't import")
    parser.add_argument("--verify-limit", type=int, default=None, help="Maximum number of molecules to verify")
    search_group = parser.add_argument_group('search options')
    search_group.add_argument("--use-keywords", action="store_true",
                      help="Use keyword search instead of property-based search")
    search_group.add_argument("--use-builtin-implementation", action="store_true",
                      help="When using property-based search, use the built-in implementation instead of the one from chembl_search_utils.py")
    args = parser.parse_args()
    
    # Create timestamps for reports
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Handle verify-only mode
    if args.verify_only:
        logger.info("Running in verification-only mode")
        
        # Verify imported data
        results = verify_imported_data(limit=args.verify_limit)
        
        # Save report
        report_file = f"reports/chembl_verification_{timestamp}.json"
        os.makedirs(os.path.dirname(report_file), exist_ok=True)
        
        with open(report_file, "w") as f:
            json.dump(results, f, indent=2)
        
        logger.info(f"Saved verification report to {report_file}")
        
        # Print summary
        print("\n=== ChEMBL Import Verification Summary ===")
        print(f"Total molecules checked: {results['molecule_count']}")
        
        prop_results = results["property_completeness"]
        print(f"Property completeness: {prop_results['property_completeness_percentage']:.2f}%")
        print(f"JSONB property completeness: {prop_results['jsonb_property_completeness_percentage']:.2f}%")
        
        xref_results = results["pubchem_cross_references"]
        print(f"PubChem cross-reference coverage: {xref_results['pubchem_cross_reference_percentage']:.2f}%")
        print(f"InChIKey coverage: {xref_results['inchikey_coverage_percentage']:.2f}%")
        
        print(f"\nDetailed report saved to: {report_file}")
        print("============================================\n")
        
        # Return code based on verification results
        prop_threshold = 95.0  # 95% property completeness
        xref_threshold = 90.0  # 90% cross-reference coverage
        
        if (prop_results['property_completeness_percentage'] >= prop_threshold and
            prop_results['jsonb_property_completeness_percentage'] >= prop_threshold and
            xref_results['pubchem_cross_reference_percentage'] >= xref_threshold):
            print("Verification PASSED: The import meets or exceeds quality thresholds.")
            return 0
        else:
            print("Verification WARNING: The import does not meet all quality thresholds.")
            print(f"Expected property completeness: >= {prop_threshold}%")
            print(f"Expected JSONB property completeness: >= {prop_threshold}%")
            print(f"Expected PubChem cross-reference coverage: >= {xref_threshold}%")
            return 1
    
    # Regular import mode
    start_time = time.time()
    
    try:
        # Verify required dependencies
        if not CHEMBL_CLIENT_AVAILABLE:
            logger.error("ChEMBL client not available. Please install chembl_webresource_client")
            return 1
        
        if not RDKIT_AVAILABLE:
            logger.warning("RDKit not available. Property calculation will be limited.")
        
        # Fetch compounds from ChEMBL
        logger.info(f"Fetching compounds from ChEMBL (limit: {args.limit})...")
        if args.use_keywords:
            logger.info("Using keyword-based search for cryoprotectants")
            compounds = fetch_cryoprotectant_compounds(limit=args.limit, use_properties=False)
        else:
            logger.info("Using property-based search for cryoprotectants")

            # When using property-based search, pass along the implementation preference
            if not args.use_builtin_implementation:
                logger.info("Using implementation from chembl_search_utils.py")
                compounds = fetch_compounds_by_property(
                    limit=args.limit,
                    use_existing_implementation=True
                )
            else:
                logger.info("Using built-in property-based search implementation")
                compounds = fetch_compounds_by_property(
                    limit=args.limit,
                    use_existing_implementation=False
                )
        
        # Import to database
        if args.dry_run:
            logger.info(f"Dry run completed. {len(compounds)} compounds would be imported.")
            stats = import_compounds_to_database(compounds, batch_size=args.batch_size, dry_run=True)
        else:
            logger.info(f"Importing {len(compounds)} compounds to database (batch size: {args.batch_size})...")
            stats = import_compounds_to_database(compounds, batch_size=args.batch_size, dry_run=False)
            logger.info(f"Import statistics: {json.dumps(stats, indent=2)}")
        
        # Verify imported data
        if not args.dry_run:
            logger.info("Verifying imported data...")
            verification_results = verify_imported_data()
            
            # Add verification results to stats
            stats["verification"] = verification_results
            
            # Log verification summary
            prop_results = verification_results["property_completeness"]
            xref_results = verification_results["pubchem_cross_references"]
            
            logger.info(f"Property completeness: {prop_results['property_completeness_percentage']:.2f}%")
            logger.info(f"JSONB property completeness: {prop_results['jsonb_property_completeness_percentage']:.2f}%")
            logger.info(f"PubChem cross-reference coverage: {xref_results['pubchem_cross_reference_percentage']:.2f}%")
        
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
            "pubchem_cross_references": stats.get("pubchem_cross_refs", {}),
            "errors": stats.get("errors", 0),
            "args": vars(args),
        }
        
        # Add verification results if available
        if "verification" in stats:
            summary["verification"] = stats["verification"]
        
        # Write summary to file
        report_file = f"reports/chembl_import_{timestamp}.json"
        os.makedirs(os.path.dirname(report_file), exist_ok=True)
        
        with open(report_file, "w") as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Saved import report to {report_file}")
        logger.info("Import completed successfully.")
        
        return 0
        
    except Exception as e:
        logger.error(f"Unhandled error in main: {str(e)}")
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())
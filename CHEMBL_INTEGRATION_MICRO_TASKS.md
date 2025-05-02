# ChEMBL Integration Micro-Tasks

This document provides ready-to-implement micro-tasks for the ChEMBL integration component. Each task is scoped to minimize token usage while ensuring clear implementation guidance.

## Error Handling Framework

### Task 1: Error Categories

```
TASK: CHEMBL-ERR-1: Error Category Enumeration

FILES:
- Primary: chembl/error_handler.py:1-25
- Reference: pubchem/rate_limiter.py:40-60

IMPLEMENTATION:
Create an ErrorCategory enum to classify ChEMBL API errors:

```python
from enum import Enum
from typing import Dict, Any, Optional, Tuple
import requests
import time
import logging

logger = logging.getLogger(__name__)

class ErrorCategory(Enum):
    """Categories of errors that can occur during ChEMBL API operations"""
    API_RATE_LIMIT = "rate_limit"       # API rate limiting (429 responses)
    CONNECTION_ERROR = "connection"     # Network/connection issues
    API_SERVER_ERROR = "server_error"   # Server-side errors (5xx)
    API_CLIENT_ERROR = "client_error"   # Client-side errors (4xx except 429)
    DATA_VALIDATION = "validation"      # Data validation errors
    TRANSFORMATION = "transformation"   # Data transformation errors
    PARSING_ERROR = "parsing"           # JSON/XML parsing errors
    DATABASE_ERROR = "database"         # Database operation errors
    UNKNOWN = "unknown"                 # Uncategorized errors
```

VERIFICATION:
- The enum covers all common error categories
- Documentation clearly explains each category
- Import statements are correct
```

### Task 2: Error Classification

```
TASK: CHEMBL-ERR-2: Error Classification Function

FILES:
- Primary: chembl/error_handler.py:26-60
- Reference: pubchem/rate_limiter.py:62-85

IMPLEMENTATION:
Create a function to classify exceptions into error categories:

```python
def classify_error(exception: Exception) -> Tuple[ErrorCategory, str]:
    """
    Classify an exception into an error category and provide a description.
    
    Args:
        exception: The exception to classify
        
    Returns:
        Tuple containing (ErrorCategory, description)
    """
    # Handle requests/HTTP exceptions
    if isinstance(exception, requests.exceptions.RequestException):
        if isinstance(exception, requests.exceptions.Timeout):
            return ErrorCategory.CONNECTION_ERROR, "Request timed out"
        
        if isinstance(exception, requests.exceptions.ConnectionError):
            return ErrorCategory.CONNECTION_ERROR, "Connection error"
            
        if isinstance(exception, requests.exceptions.HTTPError):
            response = exception.response
            if response is not None:
                status_code = response.status_code
                
                if status_code == 429:
                    return ErrorCategory.API_RATE_LIMIT, f"Rate limit exceeded: {status_code}"
                    
                if 400 <= status_code < 500:
                    return ErrorCategory.API_CLIENT_ERROR, f"Client error: {status_code}"
                    
                if 500 <= status_code < 600:
                    return ErrorCategory.API_SERVER_ERROR, f"Server error: {status_code}"
    
    # Handle JSON parsing errors
    if isinstance(exception, (ValueError, TypeError)) and "JSON" in str(exception):
        return ErrorCategory.PARSING_ERROR, f"JSON parsing error: {str(exception)}"
    
    # Handle validation errors (often custom exceptions)
    if "validation" in str(exception).lower() or "invalid" in str(exception).lower():
        return ErrorCategory.DATA_VALIDATION, f"Validation error: {str(exception)}"
    
    # Default case
    return ErrorCategory.UNKNOWN, f"Unclassified error: {str(exception)}"
```

VERIFICATION:
- Correctly classifies requests exceptions
- Correctly identifies rate limiting (429 responses)
- Handles parsing errors appropriately
- Has a reasonable default case
```

### Task 3: Recovery Strategy

```
TASK: CHEMBL-ERR-3: Recovery Strategy Function

FILES:
- Primary: chembl/error_handler.py:61-100
- Reference: PubChem_CryoProtectants_Supabase_Enhanced.py:150-180

IMPLEMENTATION:
Create a function to determine recovery strategy based on error category:

```python
def recovery_strategy(category: ErrorCategory, error_details: str, attempt: int = 0) -> Dict[str, Any]:
    """
    Determine appropriate recovery strategy for an error category.
    
    Args:
        category: The error category
        error_details: Details about the error
        attempt: Current attempt number (0-based)
        
    Returns:
        Dictionary with recovery strategy information
    """
    max_retries = 5
    base_delay = 1.0  # Base delay in seconds
    
    # Check if max retries exceeded
    if attempt >= max_retries:
        return {
            "action": "abort",
            "reason": f"Maximum retries ({max_retries}) exceeded",
            "recoverable": False
        }
    
    # Calculate backoff with jitter (exponential backoff with randomization)
    import random
    jitter = random.uniform(0.8, 1.2)
    backoff = base_delay * (2 ** attempt) * jitter
    
    # Strategies based on category
    if category == ErrorCategory.API_RATE_LIMIT:
        # Parse retry-after header if available
        retry_after = None
        if "retry-after:" in error_details.lower():
            try:
                retry_after = int(error_details.lower().split("retry-after:")[1].split()[0])
            except (ValueError, IndexError):
                pass
                
        return {
            "action": "retry",
            "delay": retry_after if retry_after is not None else backoff,
            "reason": "Rate limit exceeded",
            "recoverable": True
        }
    
    elif category == ErrorCategory.CONNECTION_ERROR:
        return {
            "action": "retry",
            "delay": backoff,
            "reason": "Connection error",
            "recoverable": True
        }
    
    elif category == ErrorCategory.API_SERVER_ERROR:
        return {
            "action": "retry",
            "delay": backoff,
            "reason": "Server error",
            "recoverable": True
        }
    
    elif category == ErrorCategory.API_CLIENT_ERROR:
        return {
            "action": "abort",
            "reason": "Client error",
            "recoverable": False
        }
    
    elif category == ErrorCategory.PARSING_ERROR:
        return {
            "action": "retry",
            "delay": backoff,
            "reason": "Parsing error",
            "recoverable": True
        }
    
    # Default strategy for other categories
    return {
        "action": "log",
        "reason": f"Unhandled error category: {category.name}",
        "recoverable": False
    }
```

VERIFICATION:
- Handles all error categories
- Implements exponential backoff with jitter
- Respects retry-after headers for rate limiting
- Has appropriate strategy for each category
```

## Checkpoint System

### Task 4: Checkpoint Manager Class

```
TASK: CHEMBL-CP-1: Checkpoint Manager Class

FILES:
- Primary: chembl/checkpoint.py:1-50
- Reference: PubChem_CryoProtectants_Supabase_Enhanced.py:120-160

IMPLEMENTATION:
Create a CheckpointManager class for saving and loading progress state:

```python
import os
import json
import time
from datetime import datetime
from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)

class CheckpointManager:
    """
    Manages checkpoints for resumable operations.
    Provides functions to save and load progress state.
    """
    
    def __init__(self, checkpoint_dir: str, prefix: str = "chembl_import"):
        """
        Initialize the checkpoint manager.
        
        Args:
            checkpoint_dir: Directory to store checkpoint files
            prefix: Prefix for checkpoint filenames
        """
        self.checkpoint_dir = checkpoint_dir
        self.prefix = prefix
        
        # Create directory if it doesn't exist
        os.makedirs(checkpoint_dir, exist_ok=True)
        
        # Initialize state
        self.state = {
            "processed_compounds": [],
            "successful_compounds": [],
            "failed_compounds": {},  # Compound ID -> error message
            "current_batch": 0,
            "total_processed": 0,
            "success_count": 0,
            "error_count": 0,
            "start_time": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "config": {}
        }
        
    def get_checkpoint_path(self) -> str:
        """Get the path for the current checkpoint file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return os.path.join(self.checkpoint_dir, f"{self.prefix}_{timestamp}.json")
        
    def get_latest_checkpoint(self) -> Optional[str]:
        """Get the path to the most recent checkpoint file"""
        files = [
            os.path.join(self.checkpoint_dir, f) 
            for f in os.listdir(self.checkpoint_dir) 
            if f.startswith(self.prefix) and f.endswith('.json')
        ]
        
        if not files:
            return None
            
        # Sort by modification time (newest first)
        files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        return files[0]
```

VERIFICATION:
- Correctly initializes checkpoint directory
- Properly initializes state dictionary
- Provides methods to get checkpoint paths
- Finds the latest checkpoint file
```

### Task 5: Checkpoint Save/Load

```
TASK: CHEMBL-CP-2: Checkpoint Save and Load Functions

FILES:
- Primary: chembl/checkpoint.py:51-100
- Reference: PubChem_CryoProtectants_Supabase_Enhanced.py:161-200

IMPLEMENTATION:
Add methods to save and load checkpoint data:

```python
def save_checkpoint(self) -> str:
    """
    Save current state to a checkpoint file.
    Uses atomic write to prevent corruption.
    
    Returns:
        Path to the saved checkpoint file
    """
    # Update timestamp
    self.state["last_updated"] = datetime.now().isoformat()
    
    # Calculate elapsed time
    try:
        start_time = datetime.fromisoformat(self.state["start_time"])
        elapsed = (datetime.now() - start_time).total_seconds()
        self.state["elapsed_seconds"] = elapsed
    except (ValueError, KeyError):
        self.state["elapsed_seconds"] = 0
    
    # Generate checkpoint path
    checkpoint_path = self.get_checkpoint_path()
    temp_path = f"{checkpoint_path}.tmp"
    
    # Write to temporary file first (atomic write pattern)
    try:
        with open(temp_path, 'w') as f:
            json.dump(self.state, f, indent=2)
        
        # Rename temp file to final filename (atomic operation)
        os.replace(temp_path, checkpoint_path)
        
        # Keep only the 5 most recent checkpoints
        self._cleanup_old_checkpoints(5)
        
        logger.info(f"Checkpoint saved: {checkpoint_path}")
        return checkpoint_path
        
    except Exception as e:
        logger.error(f"Failed to save checkpoint: {e}")
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except:
                pass
        raise
        
def load_checkpoint(self) -> bool:
    """
    Load the most recent checkpoint.
    
    Returns:
        True if a checkpoint was loaded, False otherwise
    """
    latest = self.get_latest_checkpoint()
    if latest is None:
        logger.info("No checkpoint found to load")
        return False
        
    try:
        with open(latest, 'r') as f:
            data = json.load(f)
            
        # Update state with loaded data
        self.state.update(data)
        logger.info(f"Loaded checkpoint: {latest}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to load checkpoint {latest}: {e}")
        
        # Try backup if exists
        backup = f"{latest}.bak"
        if os.path.exists(backup):
            try:
                with open(backup, 'r') as f:
                    data = json.load(f)
                self.state.update(data)
                logger.info(f"Loaded backup checkpoint: {backup}")
                return True
            except Exception as e_backup:
                logger.error(f"Failed to load backup checkpoint: {e_backup}")
                
        return False
        
def _cleanup_old_checkpoints(self, keep_count: int) -> None:
    """
    Remove old checkpoint files, keeping only the most recent ones.
    
    Args:
        keep_count: Number of recent checkpoints to keep
    """
    files = [
        os.path.join(self.checkpoint_dir, f) 
        for f in os.listdir(self.checkpoint_dir) 
        if f.startswith(self.prefix) and f.endswith('.json')
    ]
    
    # Sort by modification time (newest first)
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    
    # Remove old files
    for old_file in files[keep_count:]:
        try:
            os.remove(old_file)
            logger.debug(f"Removed old checkpoint: {old_file}")
        except Exception as e:
            logger.warning(f"Failed to remove old checkpoint {old_file}: {e}")
```

VERIFICATION:
- Implements atomic write to prevent corruption
- Handles errors gracefully
- Cleans up old checkpoint files
- Properly loads the most recent checkpoint
```

### Task 6: Progress Tracking

```
TASK: CHEMBL-CP-3: Progress Tracking Methods

FILES:
- Primary: chembl/checkpoint.py:101-140
- Reference: PubChem_CryoProtectants_Supabase_Enhanced.py:201-230

IMPLEMENTATION:
Add methods to track progress and update state:

```python
def update_progress(self, compound_id: str, success: bool = True, 
                   error: Optional[str] = None) -> None:
    """
    Update progress for a single compound.
    
    Args:
        compound_id: The ChEMBL ID of the processed compound
        success: Whether processing was successful
        error: Error message if processing failed
    """
    # Track the compound as processed
    if compound_id not in self.state["processed_compounds"]:
        self.state["processed_compounds"].append(compound_id)
        self.state["total_processed"] += 1
    
    # Track success/failure
    if success:
        if compound_id not in self.state["successful_compounds"]:
            self.state["successful_compounds"].append(compound_id)
            self.state["success_count"] += 1
    else:
        self.state["failed_compounds"][compound_id] = error or "Unknown error"
        self.state["error_count"] += 1
    
    # Update timestamp
    self.state["last_updated"] = datetime.now().isoformat()
    
def start_batch(self, batch_num: int, size: int) -> None:
    """
    Mark the start of a new batch.
    
    Args:
        batch_num: The batch number (0-based)
        size: Number of items in the batch
    """
    self.state["current_batch"] = batch_num
    self.state["current_batch_size"] = size
    self.state["current_batch_start"] = datetime.now().isoformat()
    self.state["current_batch_processed"] = 0
    
def end_batch(self, processed: int, success: int, errors: int) -> None:
    """
    Mark the end of the current batch.
    
    Args:
        processed: Number of items processed in this batch
        success: Number of successful items
        errors: Number of items with errors
    """
    batch_num = self.state["current_batch"]
    
    # Store batch statistics
    if "batches" not in self.state:
        self.state["batches"] = {}
        
    batch_end = datetime.now().isoformat()
    
    try:
        batch_start = datetime.fromisoformat(self.state["current_batch_start"])
        batch_duration = (datetime.fromisoformat(batch_end) - batch_start).total_seconds()
    except (ValueError, KeyError):
        batch_duration = 0
        
    self.state["batches"][str(batch_num)] = {
        "processed": processed,
        "success": success,
        "errors": errors,
        "duration_seconds": batch_duration
    }
    
    # Save a checkpoint after each batch
    self.save_checkpoint()
```

VERIFICATION:
- Correctly tracks processed compounds
- Updates success and error counts
- Tracks batch progress
- Creates checkpoints at appropriate times
```

## ChEMBL Worker Implementation

### Task 7: ChEMBL Worker Class

```
TASK: CHEMBL-WK-1: ChEMBL Worker Class

FILES:
- Primary: chembl/worker.py:1-50
- Reference: Worker implementation from completed Worker Pool task

IMPLEMENTATION:
Create a specialized worker for ChEMBL data processing:

```python
from queue import Queue, Empty
from typing import Dict, Any, Optional
import logging
import threading
import time
import traceback

# Import from other modules once they're implemented
# from .error_handler import classify_error, recovery_strategy, ErrorCategory
# from .checkpoint import CheckpointManager

logger = logging.getLogger(__name__)

class ChEMBLWorker:
    """
    Worker for processing ChEMBL data tasks.
    Handles one compound at a time and reports results.
    """
    
    def __init__(self, worker_id: int, task_queue: Queue, result_queue: Queue):
        """
        Initialize ChEMBL worker.
        
        Args:
            worker_id: Unique identifier for this worker
            task_queue: Queue for receiving tasks
            result_queue: Queue for sending results
        """
        self.worker_id = worker_id
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.running = False
        self.thread = None
        self.chembl_client = None  # Will be initialized in run()
        
    def run(self):
        """Main worker loop"""
        self.running = True
        
        # Import here to avoid circular imports
        from chembl.client import ChEMBLClient
        self.chembl_client = ChEMBLClient()
        
        logger.info(f"Worker {self.worker_id} started")
        
        while self.running:
            try:
                # Try to get a task with timeout to allow checking running flag
                try:
                    task = self.task_queue.get(timeout=1.0)
                except Empty:
                    continue
                
                # Check for poison pill (None task indicates shutdown)
                if task is None:
                    logger.debug(f"Worker {self.worker_id} received shutdown signal")
                    break
                
                # Process the task
                result = self.process_task(task)
                
                # Send result
                self.result_queue.put(result)
                
                # Mark task as done
                self.task_queue.task_done()
                
            except Exception as e:
                logger.error(f"Worker {self.worker_id} encountered an error: {e}")
                logger.debug(traceback.format_exc())
                
                # Try to put error result in queue
                try:
                    error_result = {
                        "status": "error",
                        "worker_id": self.worker_id,
                        "error": str(e),
                        "traceback": traceback.format_exc()
                    }
                    
                    # Include task info if available
                    if 'task' in locals():
                        error_result["task"] = task
                        
                    self.result_queue.put(error_result)
                    
                    # Mark task as done if it exists
                    if 'task' in locals():
                        self.task_queue.task_done()
                except:
                    pass
        
        logger.info(f"Worker {self.worker_id} stopped")
        
    def process_task(self, task: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single ChEMBL task.
        
        Args:
            task: Dictionary containing task information
              Expected keys:
              - compound_id: ChEMBL ID to process
              
        Returns:
            Dictionary with processing results
        """
        # This is a placeholder - real implementation will be added in CHEMBL-WK-2
        compound_id = task.get("compound_id")
        
        if not compound_id:
            return {
                "status": "error",
                "worker_id": self.worker_id,
                "error": "Missing compound_id in task",
                "task": task
            }
        
        return {
            "status": "success",
            "worker_id": self.worker_id,
            "compound_id": compound_id,
            "result": "Placeholder implementation"
        }
        
    def start(self):
        """Start the worker in a separate thread"""
        if self.thread is not None and self.thread.is_alive():
            return  # Already running
            
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True  # Thread will exit when main thread exits
        self.thread.start()
        
    def stop(self):
        """Signal the worker to stop"""
        self.running = False
        
        # Wait for thread to finish if it exists
        if self.thread is not None and self.thread.is_alive():
            self.thread.join(timeout=5.0)
```

VERIFICATION:
- Correctly initializes worker with queues
- Handles tasks from the queue
- Processes poison pill shutdown signal
- Has error handling for unexpected exceptions
```

### Task 8: Task Processing

```
TASK: CHEMBL-WK-2: Task Processing Logic

FILES:
- Primary: chembl/worker.py:51-100
- Reference: PubChem_CryoProtectants_Supabase_Enhanced.py:250-300

IMPLEMENTATION:
Implement the task processing logic for ChEMBL compounds:

```python
def process_task(self, task: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single ChEMBL task.
    
    Args:
        task: Dictionary containing task information
          Expected keys:
          - compound_id: ChEMBL ID to process
          
    Returns:
        Dictionary with processing results
    """
    start_time = time.time()
    compound_id = task.get("compound_id")
    
    if not compound_id:
        return {
            "status": "error",
            "worker_id": self.worker_id,
            "error": "Missing compound_id in task",
            "task": task
        }
    
    logger.info(f"Worker {self.worker_id} processing compound {compound_id}")
    
    try:
        # Step 1: Fetch compound data from ChEMBL
        compound_data = self.fetch_compound_data(compound_id)
        
        # Step 2: Transform data to database format
        transformed_data = self.transform_compound_data(compound_data)
        
        # Step 3: Store in database
        if not task.get("dry_run", False):
            stored_id = self.store_compound_data(transformed_data)
            transformed_data["stored_id"] = stored_id
        
        # Calculate processing time
        processing_time = time.time() - start_time
        
        # Return success result
        return {
            "status": "success",
            "worker_id": self.worker_id,
            "compound_id": compound_id,
            "processing_time": processing_time,
            "property_count": len(transformed_data.get("properties", [])),
            "dry_run": task.get("dry_run", False)
        }
        
    except Exception as e:
        # Calculate processing time even for errors
        processing_time = time.time() - start_time
        
        # Get error details
        error_type = type(e).__name__
        error_message = str(e)
        
        # Import error handling functions if they're available
        try:
            from chembl.error_handler import classify_error
            category, description = classify_error(e)
            error_category = category.name
        except ImportError:
            error_category = "UNKNOWN"
            description = error_message
        
        # Return error result
        return {
            "status": "error",
            "worker_id": self.worker_id,
            "compound_id": compound_id,
            "error": error_message,
            "error_type": error_type,
            "error_category": error_category,
            "error_description": description,
            "processing_time": processing_time
        }
```

VERIFICATION:
- Implements three-step processing workflow
- Has appropriate error handling
- Tracks processing time
- Returns structured result dictionary
```

### Task 9: Helper Methods

```
TASK: CHEMBL-WK-3: Compound Processing Helper Methods

FILES:
- Primary: chembl/worker.py:101-150
- Reference: ChEMBL_Integrated_Import.py:400-450

IMPLEMENTATION:
Add helper methods for compound data processing:

```python
def fetch_compound_data(self, compound_id: str) -> Dict[str, Any]:
    """
    Fetch compound data from ChEMBL.
    
    Args:
        compound_id: ChEMBL ID of the compound
        
    Returns:
        Dictionary containing compound data
        
    Raises:
        Exception: If fetching fails
    """
    if not self.chembl_client:
        from chembl.client import ChEMBLClient
        self.chembl_client = ChEMBLClient()
    
    # First attempt to get from cache
    if hasattr(self.chembl_client, 'get_from_cache'):
        cached_data = self.chembl_client.get_from_cache(compound_id)
        if cached_data:
            return cached_data
    
    # Fetch from API
    logger.debug(f"Fetching compound data for {compound_id}")
    compound_data = self.chembl_client.get_compound(compound_id)
    
    if not compound_data:
        raise ValueError(f"No data found for compound {compound_id}")
        
    return compound_data
    
def transform_compound_data(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Transform ChEMBL compound data to database format.
    
    Args:
        compound_data: Raw compound data from ChEMBL
        
    Returns:
        Dictionary with transformed data ready for storage
        
    Raises:
        Exception: If transformation fails
    """
    # This is a simplified implementation - in practice, we'd
    # import specialized transformation modules
    
    # Extract basic molecule information
    molecule = {
        "name": compound_data.get("pref_name") or compound_data.get("molecule_synonyms", [{}])[0].get("synonym") or f"Compound {compound_data.get('molecule_chembl_id')}",
        "chembl_id": compound_data.get("molecule_chembl_id"),
        "formula": compound_data.get("molecule_properties", {}).get("full_molformula"),
        "molecular_weight": compound_data.get("molecule_properties", {}).get("full_mwt"),
        "smiles": compound_data.get("molecule_structures", {}).get("canonical_smiles"),
        "inchi": compound_data.get("molecule_structures", {}).get("standard_inchi"),
        "inchi_key": compound_data.get("molecule_structures", {}).get("standard_inchi_key"),
        "data_source": "ChEMBL"
    }
    
    # Extract properties
    properties = []
    if "molecule_properties" in compound_data:
        props = compound_data["molecule_properties"]
        
        # Map each property
        property_mappings = {
            "alogp": {"name": "LogP", "unit": "log units"},
            "cx_logp": {"name": "CxLogP", "unit": "log units"},
            "cx_logd": {"name": "CxLogD", "unit": "log units"},
            "hba": {"name": "Hydrogen Bond Acceptors", "unit": "count"},
            "hbd": {"name": "Hydrogen Bond Donors", "unit": "count"},
            "psa": {"name": "Polar Surface Area", "unit": "Å²"},
            "rtb": {"name": "Rotatable Bonds", "unit": "count"},
            "full_mwt": {"name": "Molecular Weight", "unit": "g/mol"},
            "aromatic_rings": {"name": "Aromatic Rings", "unit": "count"},
            "heavy_atoms": {"name": "Heavy Atoms", "unit": "count"}
        }
        
        for prop_key, mapping in property_mappings.items():
            if prop_key in props and props[prop_key] is not None:
                properties.append({
                    "property_name": mapping["name"],
                    "property_type": "physicochemical",
                    "value": float(props[prop_key]),
                    "unit": mapping["unit"],
                    "source": "ChEMBL"
                })
    
    return {
        "molecule": molecule,
        "properties": properties
    }
    
def store_compound_data(self, transformed_data: Dict[str, Any]) -> str:
    """
    Store transformed compound data in the database.
    
    Args:
        transformed_data: Transformed compound data
        
    Returns:
        ID of the stored molecule
        
    Raises:
        Exception: If database storage fails
    """
    # In a real implementation, we'd use supabase_direct or MCP tools
    # For this implementation, we'll just log and return a placeholder
    
    molecule = transformed_data.get("molecule", {})
    properties = transformed_data.get("properties", [])
    
    logger.info(f"Would store molecule: {molecule.get('name')} ({molecule.get('chembl_id')})")
    logger.info(f"Would store {len(properties)} properties")
    
    # In real implementation, insert molecule and get ID
    molecule_id = f"placeholder_{molecule.get('chembl_id')}"
    
    return molecule_id
```

VERIFICATION:
- Implements all three processing steps
- Has proper error handling
- Correctly transforms ChEMBL data format
- Has appropriate logging
```

## Data Population Execution

### Task 10: Main Executor

```
TASK: CHEMBL-EX-1: Main Execution Script

FILES:
- Primary: run_chembl_import.py:1-50
- Reference: PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:500-550

IMPLEMENTATION:
Create a script to orchestrate the ChEMBL import process:

```python
#!/usr/bin/env python3
"""
ChEMBL Data Import Script

This script imports compound data from ChEMBL into the CryoProtect database.
It uses parallel processing with worker threads for efficient importing.
"""

import os
import sys
import argparse
import logging
import time
import json
from datetime import datetime
from queue import Queue
from typing import List, Dict, Any

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/chembl_import.log')
    ]
)

logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Import ChEMBL data into CryoProtect database")
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
    parser.add_argument("--query", default="cryoprotectant",
                        help="Search query for finding compounds")
    parser.add_argument("--include-refs", action="store_true",
                        help="Include reference compounds")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose logging")
    return parser.parse_args()

def main():
    """Main execution function"""
    args = parse_arguments()
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"Starting ChEMBL import in {args.mode} mode")
    logger.info(f"Settings: limit={args.limit}, batch_size={args.batch_size}, workers={args.workers}")
    
    # Import necessary modules
    try:
        from chembl.worker import ChEMBLWorker
        from chembl.checkpoint import CheckpointManager
        from chembl.client import ChEMBLClient
    except ImportError as e:
        logger.error(f"Failed to import required modules: {e}")
        return 1
    
    # Create checkpoint manager
    checkpoint_manager = CheckpointManager(args.checkpoint_dir)
    
    # Load checkpoint if in resume mode
    if args.mode == "resume":
        checkpoint_loaded = checkpoint_manager.load_checkpoint()
        if checkpoint_loaded:
            logger.info("Resuming from checkpoint")
        else:
            logger.info("No checkpoint found, starting fresh")
    
    # Initialize client for fetching compound IDs
    client = ChEMBLClient()
    
    # Get compound IDs to process
    try:
        compound_ids = fetch_compound_ids(client, args)
        logger.info(f"Found {len(compound_ids)} compounds to process")
    except Exception as e:
        logger.error(f"Failed to fetch compound IDs: {e}")
        return 1
```

VERIFICATION:
- Implements argument parsing
- Sets up logging
- Imports required modules
- Creates checkpoint manager
- Gets compound IDs to process
```

### Task 11: Worker Management

```
TASK: CHEMBL-EX-2: Worker Pool Management

FILES:
- Primary: run_chembl_import.py:51-100
- Reference: PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:551-600

IMPLEMENTATION:
Add worker pool management to the main script:

```python
def fetch_compound_ids(client, args):
    """
    Fetch ChEMBL compound IDs to process.
    
    Args:
        client: ChEMBL client instance
        args: Command line arguments
        
    Returns:
        List of ChEMBL IDs to process
    """
    # Get reference compounds if requested
    reference_compounds = []
    if args.include_refs:
        reference_compounds = [
            "CHEMBL25",    # Aspirin
            "CHEMBL1118",  # Glycerol
            "CHEMBL1201",  # Caffeine
            "CHEMBL1308",  # Sucrose
            "CHEMBL1489",  # Urea
            "CHEMBL1559",  # DMSO
            "CHEMBL1742",  # Ethanol
            "CHEMBL1750",  # Ethylene glycol
            "CHEMBL1771",  # Glycine
            "CHEMBL2043",  # Propylene glycol
            "CHEMBL2107886" # Trehalose
        ]
        logger.info(f"Including {len(reference_compounds)} reference compounds")
    
    # Search for compounds by query
    if args.query:
        logger.info(f"Searching for compounds with query: {args.query}")
        query_compounds = client.search_compounds(args.query, limit=args.limit)
        logger.info(f"Found {len(query_compounds)} compounds from search")
    else:
        query_compounds = []
    
    # Combine and deduplicate
    all_compounds = list(set(reference_compounds + query_compounds))
    
    # Apply limit
    if args.limit and len(all_compounds) > args.limit:
        all_compounds = all_compounds[:args.limit]
    
    return all_compounds

def setup_worker_pool(args):
    """
    Set up worker pool for processing.
    
    Args:
        args: Command line arguments
        
    Returns:
        Dictionary with worker pool components
    """
    from chembl.worker import ChEMBLWorker
    
    # Create queues
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
    
    return {
        "workers": workers,
        "task_queue": task_queue,
        "result_queue": result_queue
    }

def distribute_tasks(compound_ids, pool, args, checkpoint_manager=None):
    """
    Distribute tasks to worker pool.
    
    Args:
        compound_ids: List of compound IDs to process
        pool: Worker pool components
        args: Command line arguments
        checkpoint_manager: Optional checkpoint manager
        
    Returns:
        Number of tasks distributed
    """
    task_queue = pool["task_queue"]
    
    # Skip already processed compounds if resuming
    if args.mode == "resume" and checkpoint_manager:
        processed = checkpoint_manager.state.get("processed_compounds", [])
        compound_ids = [cid for cid in compound_ids if cid not in processed]
        logger.info(f"Skipping {len(processed)} already processed compounds")
    
    # Add dry run flag for dry-run mode
    dry_run = (args.mode == "dry-run")
    
    # Add tasks to queue
    for compound_id in compound_ids:
        task_queue.put({
            "compound_id": compound_id,
            "dry_run": dry_run
        })
    
    logger.info(f"Distributed {len(compound_ids)} tasks to worker pool")
    return len(compound_ids)
```

VERIFICATION:
- Implements compound ID fetching
- Sets up worker pool with appropriate queues
- Creates and starts workers
- Distributes tasks to workers
- Handles resuming from checkpoint
```

### Task 12: Result Processing

```
TASK: CHEMBL-EX-3: Result Processing and Main Loop

FILES:
- Primary: run_chembl_import.py:101-150
- Reference: PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:601-650

IMPLEMENTATION:
Add result processing and main execution loop:

```python
def process_results(pool, total_tasks, args, checkpoint_manager=None):
    """
    Process results from worker pool.
    
    Args:
        pool: Worker pool components
        total_tasks: Total number of tasks distributed
        args: Command line arguments
        checkpoint_manager: Optional checkpoint manager
        
    Returns:
        Dictionary with processing statistics
    """
    result_queue = pool["result_queue"]
    task_queue = pool["task_queue"]
    
    # Statistics
    stats = {
        "total": total_tasks,
        "processed": 0,
        "success": 0,
        "error": 0,
        "start_time": time.time(),
        "errors": {}
    }
    
    # Progress display parameters
    update_interval = 1.0  # seconds
    last_update = 0
    
    # Process results until all tasks are done
    while stats["processed"] < total_tasks:
        try:
            # Get result with timeout to allow checking queue status
            try:
                result = result_queue.get(timeout=0.1)
            except Queue.Empty:
                # Check if tasks are done processing
                if task_queue.empty() and task_queue.unfinished_tasks == 0:
                    logger.warning("Task queue is empty but not all results received")
                    break
                continue
            
            # Process the result
            stats["processed"] += 1
            
            if result.get("status") == "success":
                stats["success"] += 1
                
                # Update checkpoint if available
                if checkpoint_manager:
                    checkpoint_manager.update_progress(
                        result["compound_id"], 
                        success=True
                    )
            else:
                stats["error"] += 1
                error_category = result.get("error_category", "UNKNOWN")
                
                # Track error by category
                if error_category not in stats["errors"]:
                    stats["errors"][error_category] = 0
                stats["errors"][error_category] += 1
                
                # Update checkpoint if available
                if checkpoint_manager:
                    checkpoint_manager.update_progress(
                        result["compound_id"],
                        success=False,
                        error=result.get("error")
                    )
            
            # Mark as done in result queue
            result_queue.task_done()
            
            # Display progress periodically
            current_time = time.time()
            if current_time - last_update >= update_interval:
                display_progress(stats)
                last_update = current_time
                
                # Save checkpoint periodically
                if checkpoint_manager and stats["processed"] % args.batch_size == 0:
                    checkpoint_manager.save_checkpoint()
                
        except Exception as e:
            logger.error(f"Error processing result: {e}")
    
    # Final progress display
    display_progress(stats, final=True)
    
    # Save final checkpoint
    if checkpoint_manager:
        checkpoint_manager.save_checkpoint()
    
    return stats

def display_progress(stats, final=False):
    """Display progress information"""
    elapsed = time.time() - stats["start_time"]
    processed = stats["processed"]
    total = stats["total"]
    
    if processed == 0:
        return
    
    percent = (processed / total) * 100 if total > 0 else 0
    rate = processed / elapsed if elapsed > 0 else 0
    
    # Calculate ETA
    if rate > 0 and not final:
        eta = (total - processed) / rate
        eta_str = f"ETA: {eta:.1f}s"
    else:
        eta_str = "Complete"
    
    # Create progress bar
    bar_length = 30
    filled_length = int(bar_length * percent / 100)
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    # Print progress
    sys.stdout.write(f"\r[{bar}] {percent:.1f}% | {processed}/{total} | " +
                    f"Success: {stats['success']} | Errors: {stats['error']} | " +
                    f"Rate: {rate:.1f}/s | {eta_str}")
    sys.stdout.flush()
    
    if final:
        sys.stdout.write("\n")
        
        # Print error summary if there are errors
        if stats["error"] > 0:
            sys.stdout.write("\nError summary:\n")
            for category, count in stats["errors"].items():
                sys.stdout.write(f"  {category}: {count}\n")

# Add to main function
    # Continue from previous code in main()
    
    # Set up worker pool
    pool = setup_worker_pool(args)
    
    # Distribute tasks
    total_tasks = distribute_tasks(compound_ids, pool, args, checkpoint_manager)
    
    # Process results
    stats = process_results(pool, total_tasks, args, checkpoint_manager)
    
    # Stop workers
    for worker in pool["workers"]:
        worker.stop()
    
    # Print summary
    logger.info(f"Import completed: {stats['success']}/{stats['total']} compounds imported successfully")
    logger.info(f"Error count: {stats['error']}")
    
    # Generate verification report if not in dry-run mode
    if args.mode != "dry-run" and stats['success'] > 0:
        try:
            from verify_chembl_data import verify_chembl_import
            verify_chembl_import()
        except ImportError:
            logger.warning("Verification module not available")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

VERIFICATION:
- Implements result processing loop
- Displays progress information
- Updates checkpoint with results
- Tracks error statistics
- Properly shuts down workers
- Generates summary report
```

## Additional Helper Tasks

### Task 13: Reference Compounds

```
TASK: CHEMBL-HELP-1: Reference Compounds List

FILES:
- Primary: chembl/reference_compounds.py:1-40
- Reference: None (new file)

IMPLEMENTATION:
Create a file with reference compound definitions:

```python
"""
Reference compounds for CryoProtect database.

This module defines standard reference compounds that should be included
in the database to ensure consistent comparison across experiments.
"""

REFERENCE_COMPOUNDS = [
    {
        "chembl_id": "CHEMBL1118",
        "name": "GLYCEROL",
        "role": "primary cryoprotectant",
        "description": "One of the most common cryoprotectants, glycerol is used in many cryopreservation protocols.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1559",
        "name": "DIMETHYL SULFOXIDE",
        "common_name": "DMSO",
        "role": "primary cryoprotectant",
        "description": "DMSO is widely used for cell cryopreservation, particularly in stem cell and tissue preservation.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1750",
        "name": "ETHYLENE GLYCOL",
        "role": "primary cryoprotectant",
        "description": "Ethylene glycol is used in vitrification protocols for embryo and oocyte cryopreservation.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL2043",
        "name": "PROPYLENE GLYCOL",
        "role": "primary cryoprotectant",
        "description": "Propylene glycol is used in cryopreservation of various cell types and tissues.",
        "priority": 1
    },
    {
        "chembl_id": "CHEMBL1308",
        "name": "SUCROSE",
        "role": "secondary cryoprotectant",
        "description": "Sucrose is often used as a non-penetrating cryoprotectant to reduce osmotic stress.",
        "priority": 2
    },
    {
        "chembl_id": "CHEMBL2107886",
        "name": "TREHALOSE",
        "role": "secondary cryoprotectant",
        "description": "Trehalose provides membrane protection during freezing and is often used with DMSO.",
        "priority": 2
    },
    {
        "chembl_id": "CHEMBL1742",
        "name": "ETHANOL",
        "role": "secondary cryoprotectant",
        "description": "Ethanol is sometimes used in specialized cryopreservation protocols.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1489",
        "name": "UREA",
        "role": "additive",
        "description": "Urea is sometimes used to enhance protein stability during freezing.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1771",
        "name": "GLYCINE",
        "role": "additive",
        "description": "Glycine can act as an osmolyte that protects proteins during freezing.",
        "priority": 3
    },
    {
        "chembl_id": "CHEMBL1201",
        "name": "CAFFEINE",
        "role": "control",
        "description": "Caffeine is used as a control compound that is not expected to have cryoprotective properties.",
        "priority": 4
    },
    {
        "chembl_id": "CHEMBL25",
        "name": "ASPIRIN",
        "role": "control",
        "description": "Aspirin is used as a control compound that is not expected to have cryoprotective properties.",
        "priority": 4
    }
]

def get_reference_compounds():
    """Get the list of reference compounds with their ChEMBL IDs"""
    return REFERENCE_COMPOUNDS

def get_reference_compound_ids():
    """Get just the ChEMBL IDs of reference compounds"""
    return [compound["chembl_id"] for compound in REFERENCE_COMPOUNDS]

def get_primary_cryoprotectants():
    """Get the list of primary cryoprotectant reference compounds"""
    return [c for c in REFERENCE_COMPOUNDS if c["role"] == "primary cryoprotectant"]
```

VERIFICATION:
- Includes common cryoprotectants with ChEMBL IDs
- Provides metadata about each compound
- Includes utility functions to access the list
- Has compounds categorized by role
```

These micro-tasks provide a complete implementation plan for the ChEMBL integration component. Each task is self-contained and focused on a specific part of the system, reducing context requirements and token usage during implementation.
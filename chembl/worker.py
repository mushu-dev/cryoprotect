from queue import Queue, Empty
from typing import Dict, Any, Optional, List
import logging
import threading
import time
import traceback
import uuid

# Import from other modules
from .error_handler import classify_error, get_recovery_strategy, ErrorCategory
from .checkpoint import CheckpointManager

# Import direct PostgreSQL connection utilities
from postgres_direct import PostgresDirectConnection
from sql_executor import (
    get_db,
    execute_query,
    bulk_insert,
    execute_batch,
    with_retry,
    with_transaction,
    process_in_batches
)

# Import enhanced property manager
from property_utils import PropertyManager

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
        Process a single ChEMBL task with enhanced retry logic and error handling.
        
        This method implements a retry loop with exponential backoff for recoverable errors.
        It uses the error classification and recovery strategy functions to determine
        appropriate actions for different error types. Progress is tracked via the
        checkpoint manager. Enhanced with direct PostgreSQL connection and batch processing.
        
        Args:
            task: Dictionary containing task information
              Expected keys:
              - compound_id: ChEMBL ID to process
              - checkpoint_manager: Optional CheckpointManager instance
              - max_retries: Optional maximum number of retries (default: 3)
              - dry_run: Optional flag to skip database operations (default: False)
              - batch_size: Optional batch size for property operations (default: 100)
              
        Returns:
            Dictionary with processing results
        """
        start_time = time.time()
        compound_id = task.get("compound_id")
        max_retries = task.get("max_retries", 3)
        checkpoint_manager = task.get("checkpoint_manager")
        batch_size = task.get("batch_size", 100)
        
        # Validate required parameters
        if not compound_id:
            return {
                "status": "error",
                "worker_id": self.worker_id,
                "error": "Missing compound_id in task",
                "task": task
            }
        
        logger.info(f"Worker {self.worker_id} processing compound {compound_id}")
        
        # Initialize attempt counter and result
        attempt = 0
        result = None
        
        # Create a unique operation ID for tracking
        operation_id = str(uuid.uuid4())
        
        # Retry loop with enhanced error handling
        while attempt <= max_retries:
            try:
                # Step 1: Fetch compound data from ChEMBL with optimized retry
                logger.debug(f"Worker {self.worker_id} fetching data for {compound_id} (attempt {attempt+1}, operation {operation_id})")
                compound_data = self._fetch_compound_data_with_retry(compound_id)
                
                # Step 2: Transform data to database format
                logger.debug(f"Worker {self.worker_id} transforming data for {compound_id}")
                transformed_data = self.transform_compound_data(compound_data)
                
                # Step 3: Store in database with transaction and batch processing
                if not task.get("dry_run", False):
                    logger.debug(f"Worker {self.worker_id} storing data for {compound_id}")
                    
                    # Use transaction to ensure atomicity
                    with get_db().transaction() as conn:
                        stored_id = self.store_compound_data(transformed_data, conn=conn)
                        transformed_data["stored_id"] = stored_id
                
                # Calculate processing time
                processing_time = time.time() - start_time
                
                # Update checkpoint with detailed information
                if checkpoint_manager:
                    checkpoint_data = {
                        "operation_id": operation_id,
                        "processing_time": processing_time,
                        "property_count": len(transformed_data.get("properties", [])),
                        "timestamp": time.time()
                    }
                    checkpoint_manager.update_progress(compound_id, success=True, data=checkpoint_data)
                
                # Return enhanced success result with more details
                result = {
                    "status": "success",
                    "worker_id": self.worker_id,
                    "compound_id": compound_id,
                    "processing_time": processing_time,
                    "property_count": len(transformed_data.get("properties", [])),
                    "dry_run": task.get("dry_run", False),
                    "attempts": attempt + 1,
                    "operation_id": operation_id,
                    "stored_id": transformed_data.get("stored_id")
                }
                
                # Success - exit retry loop
                break
                
            except Exception as e:
                # Calculate processing time even for errors
                processing_time = time.time() - start_time
                
                # Get error details with enhanced tracking
                error_type = type(e).__name__
                error_message = str(e)
                error_traceback = traceback.format_exc()
                
                # Classify the error with better error handling
                try:
                    category, description = classify_error(e)
                    error_category = category.name
                except Exception as classify_error:
                    logger.warning(f"Error during error classification: {classify_error}")
                    category = ErrorCategory.UNKNOWN
                    error_category = "UNKNOWN"
                    description = error_message
                
                # Get recovery strategy with fallback options
                try:
                    recovery_action = get_recovery_strategy(category, description, attempt)
                except Exception as recovery_error:
                    logger.warning(f"Error determining recovery strategy: {recovery_error}")
                    recovery_action = "ABORT"  # Default to abort on error
                
                logger.info(f"Worker {self.worker_id} encountered {error_category} error on attempt {attempt+1}: {error_message}")
                logger.info(f"Recovery action: {recovery_action}")
                
                # Enhanced error handling with more detailed logging
                if recovery_action == "RETRY":
                    # Increment attempt counter
                    attempt += 1
                    
                    if attempt <= max_retries:
                        # Calculate backoff with jitter and progressive delay
                        import random
                        base_delay = 2.0  # Increased base delay for better resilience
                        jitter = random.uniform(0.8, 1.2)
                        backoff = base_delay * (2 ** attempt) * jitter
                        
                        # Cap maximum backoff at 60 seconds
                        backoff = min(backoff, 60.0)
                        
                        logger.info(f"Retrying after {backoff:.2f} seconds (attempt {attempt+1}/{max_retries+1})")
                        time.sleep(backoff)
                        
                        # Log retry attempt with operation ID for tracking
                        logger.debug(f"Retry attempt {attempt+1} for operation {operation_id}, compound {compound_id}")
                        continue
                    else:
                        # Max retries exceeded with enhanced error reporting
                        logger.warning(f"Maximum retries ({max_retries}) exceeded for {compound_id}")
                        
                        # Update checkpoint with detailed error information
                        if checkpoint_manager:
                            error_data = {
                                "operation_id": operation_id,
                                "error_type": error_type,
                                "error_category": error_category,
                                "attempts": attempt,
                                "processing_time": processing_time,
                                "timestamp": time.time()
                            }
                            checkpoint_manager.update_progress(
                                compound_id,
                                success=False,
                                error=f"Max retries exceeded: {error_message}",
                                data=error_data
                            )
                        
                        # Create enhanced error result with more details
                        result = {
                            "status": "error",
                            "worker_id": self.worker_id,
                            "compound_id": compound_id,
                            "error": error_message,
                            "error_type": error_type,
                            "error_category": error_category,
                            "error_description": description,
                            "error_traceback": error_traceback,
                            "processing_time": processing_time,
                            "attempts": attempt,
                            "recovery_action": "MAX_RETRIES_EXCEEDED",
                            "operation_id": operation_id
                        }
                        break
                        
                elif recovery_action == "SKIP":
                    # Skip this compound with enhanced reporting
                    logger.info(f"Skipping compound {compound_id} due to {error_category}")
                    
                    # Update checkpoint with skip reason
                    if checkpoint_manager:
                        skip_data = {
                            "operation_id": operation_id,
                            "error_type": error_type,
                            "error_category": error_category,
                            "processing_time": processing_time,
                            "timestamp": time.time()
                        }
                        checkpoint_manager.update_progress(
                            compound_id,
                            success=False,
                            error=f"Skipped: {error_message}",
                            data=skip_data
                        )
                    
                    # Create enhanced skip result
                    result = {
                        "status": "skipped",
                        "worker_id": self.worker_id,
                        "compound_id": compound_id,
                        "error": error_message,
                        "error_type": error_type,
                        "error_category": error_category,
                        "error_description": description,
                        "processing_time": processing_time,
                        "attempts": attempt + 1,
                        "recovery_action": "SKIP",
                        "operation_id": operation_id
                    }
                    break
                    
                else:  # ABORT or any other action with enhanced reporting
                    # Abort processing for this compound
                    logger.warning(f"Aborting processing for {compound_id} due to {error_category}")
                    
                    # Update checkpoint with abort reason and details
                    if checkpoint_manager:
                        abort_data = {
                            "operation_id": operation_id,
                            "error_type": error_type,
                            "error_category": error_category,
                            "error_traceback": error_traceback,
                            "processing_time": processing_time,
                            "timestamp": time.time()
                        }
                        checkpoint_manager.update_progress(
                            compound_id,
                            success=False,
                            error=f"Aborted: {error_message}",
                            data=abort_data
                        )
                    
                    # Create enhanced abort result
                    result = {
                        "status": "error",
                        "worker_id": self.worker_id,
                        "compound_id": compound_id,
                        "error": error_message,
                        "error_type": error_type,
                        "error_category": error_category,
                        "error_description": description,
                        "error_traceback": error_traceback,
                        "processing_time": processing_time,
                        "attempts": attempt + 1,
                        "recovery_action": "ABORT",
                        "operation_id": operation_id
                    }
                    break
        
        return result
        
    def _fetch_compound_data_with_retry(self, compound_id: str) -> Dict[str, Any]:
        """
        Fetch compound data from ChEMBL with enhanced retry logic.
        
        This method uses the enhanced retry decorator from sql_executor for improved
        resilience with exponential backoff, jitter, and comprehensive error handling.
        
        Args:
            compound_id: ChEMBL ID of the compound
            
        Returns:
            Dictionary containing compound data
            
        Raises:
            Exception: If fetching fails after retries
        """
        # Define a retry-wrapped function for fetching compound data
        @with_retry(max_retries=5, retry_delay=2, retry_backoff=1.5)
        def fetch_with_enhanced_retry(cid):
            logger.debug(f"Fetching compound data for {cid} with enhanced retry")
            return self.fetch_compound_data(cid)
        
        # Call the retry-wrapped function
        try:
            return fetch_with_enhanced_retry(compound_id)
        except Exception as e:
            logger.error(f"Failed to fetch compound data for {compound_id} after multiple retries: {str(e)}")
            raise
        
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
    
    def join(self, timeout=None):
        """Join the worker thread"""
        if self.thread is not None and self.thread.is_alive():
            self.thread.join(timeout=timeout)
            
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
                logger.debug(f"Retrieved compound {compound_id} from cache")
                return cached_data
        
        # Fetch from API
        logger.debug(f"Fetching compound data for {compound_id} from API")
        compound_data = self.chembl_client.get_compound(compound_id)
        
        if not compound_data:
            raise ValueError(f"No data found for compound {compound_id}")
            
        # Store in cache if possible
        if hasattr(self.chembl_client, 'store_in_cache'):
            self.chembl_client.store_in_cache(compound_id, compound_data)
            
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
        
        # Extract cross-references if available
        cross_references = {}
        if "cross_references" in compound_data:
            for ref_type, ref_id in compound_data.get("cross_references", {}).items():
                cross_references[ref_type] = ref_id
        
        # Add cross-references to molecule data
        molecule["cross_references"] = cross_references
        
        # Extract properties using the enhanced extraction method
        properties = self._extract_properties(compound_data)
        
        return {
            "molecule": molecule,
            "properties": properties
        }
        
    def _extract_properties(self, compound_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Extract comprehensive property data from ChEMBL compound data.
        
        Args:
            compound_data: ChEMBL compound data
            
        Returns:
            List of extracted properties
        """
        properties = []
        
        # Basic molecular properties
        molecule_props = compound_data.get('molecule_properties', {})
        if molecule_props:
            # Map each property with more comprehensive coverage
            property_mappings = {
                "alogp": {"name": "LogP", "unit": "log units", "type": "physicochemical"},
                "cx_logp": {"name": "CxLogP", "unit": "log units", "type": "physicochemical"},
                "cx_logd": {"name": "CxLogD", "unit": "log units", "type": "physicochemical"},
                "hba": {"name": "Hydrogen Bond Acceptors", "unit": "count", "type": "physicochemical"},
                "hbd": {"name": "Hydrogen Bond Donors", "unit": "count", "type": "physicochemical"},
                "psa": {"name": "Polar Surface Area", "unit": "Å²", "type": "physicochemical"},
                "rtb": {"name": "Rotatable Bonds", "unit": "count", "type": "physicochemical"},
                "full_mwt": {"name": "Molecular Weight", "unit": "g/mol", "type": "physicochemical"},
                "mw_freebase": {"name": "Freebase Molecular Weight", "unit": "g/mol", "type": "physicochemical"},
                "aromatic_rings": {"name": "Aromatic Rings", "unit": "count", "type": "structural"},
                "heavy_atoms": {"name": "Heavy Atoms", "unit": "count", "type": "structural"},
                "qed_weighted": {"name": "QED Weighted", "unit": "score", "type": "druglikeness"},
                "ro3_pass": {"name": "Rule of Three Pass", "unit": "boolean", "type": "druglikeness"},
                "num_ro5_violations": {"name": "Lipinski Violations", "unit": "count", "type": "druglikeness"},
                "med_chem_friendly": {"name": "Med Chem Friendly", "unit": "boolean", "type": "druglikeness"},
                "full_molformula": {"name": "Molecular Formula", "unit": "formula", "type": "structural"}
            }
            
            for prop_key, mapping in property_mappings.items():
                if prop_key in molecule_props and molecule_props[prop_key] is not None:
                    # Handle boolean properties
                    if mapping["unit"] == "boolean":
                        value = molecule_props[prop_key] == 'Y'
                    else:
                        try:
                            value = float(molecule_props[prop_key])
                        except (ValueError, TypeError):
                            value = molecule_props[prop_key]
                    
                    properties.append({
                        "property_name": mapping["name"],
                        "property_type": mapping["type"],
                        "value": value,
                        "unit": mapping["unit"],
                        "source": "ChEMBL"
                    })
        
        # Physical properties (from bioactivities if available)
        bioactivities = compound_data.get('activities', [])
        if bioactivities:
            # Extract properties from bioactivities
            for activity in bioactivities:
                if activity.get('type') == 'LogP':
                    properties.append({
                        "property_name": "Experimental LogP",
                        "property_type": "physicochemical",
                        "value": float(activity.get('value')),
                        "unit": "log units",
                        "source": "ChEMBL-Bioactivity"
                    })
                elif activity.get('type') == 'Solubility':
                    properties.append({
                        "property_name": "Solubility",
                        "property_type": "physicochemical",
                        "value": float(activity.get('value')),
                        "unit": activity.get('units', 'mg/mL'),
                        "source": "ChEMBL-Bioactivity"
                    })
                elif activity.get('type') == 'Melting Point':
                    properties.append({
                        "property_name": "Melting Point",
                        "property_type": "physicochemical",
                        "value": float(activity.get('value')),
                        "unit": activity.get('units', '°C'),
                        "source": "ChEMBL-Bioactivity"
                    })
        
        # Extract additional properties from compound_data if available
        if 'molecule_type' in compound_data:
            properties.append({
                "property_name": "Molecule Type",
                "property_type": "classification",
                "value": compound_data['molecule_type'],
                "unit": "category",
                "source": "ChEMBL"
            })
            
        if 'max_phase' in compound_data:
            properties.append({
                "property_name": "Max Development Phase",
                "property_type": "clinical",
                "value": int(compound_data['max_phase']),
                "unit": "phase",
                "source": "ChEMBL"
            })
            
        # Add cryoprotectant-specific properties if available
        # This would be expanded based on domain knowledge
        
        return properties
        
    @with_retry(max_retries=3, retry_backoff=2.0)
    @with_transaction
    def store_compound_data(self, transformed_data: Dict[str, Any], conn=None) -> str:
        """
        Store transformed compound data in the database using direct PostgreSQL connection.
        Uses enhanced PropertyManager for efficient property handling.
        
        Args:
            transformed_data: Transformed compound data
            conn: Database connection object (injected by @with_transaction)
            
        Returns:
            ID of the stored molecule
            
        Raises:
            Exception: If database storage fails
        """
        molecule = transformed_data.get("molecule", {})
        properties = transformed_data.get("properties", [])
        
        logger.info(f"Storing molecule: {molecule.get('name')} ({molecule.get('chembl_id')})")
        logger.info(f"Storing {len(properties)} properties")
        
        try:
            # Insert molecule using direct SQL execution
            results = execute_query("""
                INSERT INTO molecules (
                    name, formula, molecular_weight, smiles, inchi, inchi_key,
                    chembl_id, pubchem_cid, data_source, created_at, updated_at
                ) VALUES (
                    %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
                ) RETURNING id
            """, (
                molecule.get('name'),
                molecule.get('formula'),
                molecule.get('molecular_weight'),
                molecule.get('smiles'),
                molecule.get('inchi'),
                molecule.get('inchi_key'),
                molecule.get('chembl_id'),
                molecule.get('pubchem_cid'),
                molecule.get('data_source', 'ChEMBL')
            ))
            
            # Get the first result if available
            result = results[0] if results and len(results) > 0 else None
            
            if not result:
                raise ValueError(f"Failed to insert molecule: {molecule.get('name')}")
                
            molecule_id = result['id']
            
            # Use enhanced PropertyManager for batch property operations
            property_manager = PropertyManager()
            
            # Convert properties to format expected by PropertyManager
            property_dict = {}
            for prop in properties:
                property_dict[prop['property_name']] = prop['value']
                
                # Also store unit as a separate property if available
                if prop.get('unit'):
                    property_dict[f"{prop['property_name']}_unit"] = prop['unit']
                
                # Store property type as metadata
                if prop.get('property_type'):
                    property_dict[f"{prop['property_name']}_type"] = prop['property_type']
                
                # Store source as metadata
                if prop.get('source'):
                    property_dict[f"{prop['property_name']}_source"] = prop['source']
            
            # Set all properties in a single batch operation
            success_count, total_count = property_manager.set_properties(
                molecule_id=molecule_id,
                properties=property_dict
            )
            
            logger.info(f"Added {success_count}/{total_count} properties for molecule {molecule_id}")
            
            return molecule_id
            
        except Exception as e:
            logger.error(f"Error storing compound data: {str(e)}")
            raise
            
    def _get_data_type(self, value: Any) -> str:
        """
        Determine the data type for a property value.
        
        Args:
            value: Property value
            
        Returns:
            Data type string ('numeric', 'text', or 'boolean')
        """
        if isinstance(value, bool):
            return 'boolean'
        elif isinstance(value, (int, float)):
            return 'numeric'
        else:
            return 'text'
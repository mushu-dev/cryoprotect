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
                
    def update_progress(self, compound_id: str, success: bool = True,
                       error: Optional[str] = None, data: Optional[Dict[str, Any]] = None) -> None:
        """
        Update progress for a single compound.
        
        Args:
            compound_id: The ChEMBL ID of the processed compound
            success: Whether processing was successful
            error: Error message if processing failed
            data: Optional additional data to store with the progress update
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
        
        # Store additional data if provided
        if data and isinstance(data, dict):
            if "compound_data" not in self.state:
                self.state["compound_data"] = {}
            self.state["compound_data"][compound_id] = data
            
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
            
        # Record batch statistics
        self.state["batches"][str(batch_num)] = {
            "processed": processed,
            "success": success,
            "errors": errors,
            "start_time": self.state.get("current_batch_start"),
            "end_time": batch_end,
            "duration_seconds": batch_duration,
            "items_per_second": processed / batch_duration if batch_duration > 0 else 0
        }
        
        # Update overall statistics
        self.state["current_batch_processed"] = processed
        
        # Calculate and store performance metrics
        if "performance" not in self.state:
            self.state["performance"] = {
                "avg_items_per_second": 0,
                "total_duration_seconds": 0,
                "batch_count": 0
            }
            
        perf = self.state["performance"]
        perf["batch_count"] += 1
        perf["total_duration_seconds"] += batch_duration
        
        # Recalculate average performance
        total_processed = sum(
            batch_data.get("processed", 0)
            for batch_data in self.state["batches"].values()
        )
        
        if perf["total_duration_seconds"] > 0:
            perf["avg_items_per_second"] = total_processed / perf["total_duration_seconds"]
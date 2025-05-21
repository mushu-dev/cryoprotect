"""
Checkpoint management for resumable imports.

This module provides robust checkpointing to allow import operations
to be resumed from the point of failure.
"""

import os
import json
import time
import shutil
import logging
from typing import Dict, Any, Optional, List, Set, Union


class CheckpointManager:
    """
    Manages checkpoints for resumable molecular data imports.
    
    This class provides functionality to:
    - Save progress at regular intervals
    - Track which items have been processed
    - Resume operations from the last successful point
    - Back up checkpoint files to prevent data loss
    """
    
    def __init__(
        self,
        checkpoint_file: str,
        backup_interval: int = 10,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the checkpoint manager.
        
        Args:
            checkpoint_file: Path to the checkpoint file
            backup_interval: How many save operations before creating a backup
            logger: Logger instance for recording checkpoint operations
        """
        self.checkpoint_file = checkpoint_file
        self.backup_interval = backup_interval
        self.logger = logger or logging.getLogger(__name__)
        
        # Initialize checkpoint data structure
        self.data: Dict[str, Any] = {
            'metadata': {
                'created_at': time.time(),
                'last_updated': time.time(),
                'save_count': 0,
                'version': '1.0'
            },
            'progress': {
                'total_items': 0,
                'processed_items': 0,
                'successful_items': 0,
                'failed_items': 0,
                'skipped_items': 0
            },
            'processed_ids': set(),
            'failed_ids': set(),
            'skipped_ids': set(),
            'last_batch': [],
            'current_position': None,
            'custom_state': {}
        }
        
        # Load existing checkpoint if available
        self._load()
    
    def _load(self) -> None:
        """Load checkpoint data from file if it exists."""
        if not os.path.exists(self.checkpoint_file):
            self.logger.info(f"No checkpoint file found at {self.checkpoint_file}")
            return
            
        try:
            with open(self.checkpoint_file, 'r') as f:
                data = json.load(f)
                
            # Convert lists back to sets for id tracking
            if 'processed_ids' in data:
                data['processed_ids'] = set(data['processed_ids'])
            if 'failed_ids' in data:
                data['failed_ids'] = set(data['failed_ids'])
            if 'skipped_ids' in data:
                data['skipped_ids'] = set(data['skipped_ids'])
                
            self.data.update(data)
            self.logger.info(
                f"Loaded checkpoint from {self.checkpoint_file}",
                extra={'structured_data': {'checkpoint_state': self.get_summary()}}
            )
        except (json.JSONDecodeError, IOError) as e:
            self.logger.error(f"Failed to load checkpoint: {str(e)}")
            
            # Try to load from backup if main file is corrupted
            backup_file = f"{self.checkpoint_file}.bak"
            if os.path.exists(backup_file):
                try:
                    with open(backup_file, 'r') as f:
                        data = json.load(f)
                    
                    # Convert lists back to sets for id tracking
                    if 'processed_ids' in data:
                        data['processed_ids'] = set(data['processed_ids'])
                    if 'failed_ids' in data:
                        data['failed_ids'] = set(data['failed_ids'])
                    if 'skipped_ids' in data:
                        data['skipped_ids'] = set(data['skipped_ids'])
                        
                    self.data.update(data)
                    self.logger.warning(
                        f"Loaded checkpoint from backup file {backup_file}"
                    )
                except (json.JSONDecodeError, IOError) as e2:
                    self.logger.error(f"Failed to load backup checkpoint: {str(e2)}")
    
    def save(self) -> None:
        """Save checkpoint data to file."""
        # Update metadata
        self.data['metadata']['last_updated'] = time.time()
        self.data['metadata']['save_count'] += 1
        
        # Create temp serializable copy with sets converted to lists
        temp_data = self.data.copy()
        temp_data['processed_ids'] = list(self.data['processed_ids'])
        temp_data['failed_ids'] = list(self.data['failed_ids'])
        temp_data['skipped_ids'] = list(self.data['skipped_ids'])
        
        # Ensure directory exists
        checkpoint_dir = os.path.dirname(self.checkpoint_file)
        if checkpoint_dir and not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
            
        # Save to temp file first
        temp_file = f"{self.checkpoint_file}.tmp"
        try:
            with open(temp_file, 'w') as f:
                json.dump(temp_data, f)
                
            # Replace main file with temp file
            shutil.move(temp_file, self.checkpoint_file)
            
            # Create backup periodically
            if self.data['metadata']['save_count'] % self.backup_interval == 0:
                backup_file = f"{self.checkpoint_file}.bak"
                shutil.copy(self.checkpoint_file, backup_file)
                self.logger.debug(f"Created backup checkpoint at {backup_file}")
                
            self.logger.debug(f"Saved checkpoint to {self.checkpoint_file}")
        except IOError as e:
            self.logger.error(f"Failed to save checkpoint: {str(e)}")
    
    def mark_processed(self, item_id: Union[str, int]) -> None:
        """
        Mark an item as successfully processed.
        
        Args:
            item_id: Unique identifier of the processed item
        """
        item_id_str = str(item_id)
        if item_id_str not in self.data['processed_ids']:
            self.data['processed_ids'].add(item_id_str)
            self.data['progress']['processed_items'] += 1
            self.data['progress']['successful_items'] += 1
            
            # Remove from failed/skipped if it was there
            if item_id_str in self.data['failed_ids']:
                self.data['failed_ids'].remove(item_id_str)
            if item_id_str in self.data['skipped_ids']:
                self.data['skipped_ids'].remove(item_id_str)
    
    def mark_failed(self, item_id: Union[str, int]) -> None:
        """
        Mark an item as failed during processing.
        
        Args:
            item_id: Unique identifier of the failed item
        """
        item_id_str = str(item_id)
        if item_id_str not in self.data['failed_ids']:
            self.data['failed_ids'].add(item_id_str)
            
            # Only increment processed count if it wasn't already processed
            if item_id_str not in self.data['processed_ids']:
                self.data['progress']['processed_items'] += 1
                
            self.data['progress']['failed_items'] += 1
    
    def mark_skipped(self, item_id: Union[str, int]) -> None:
        """
        Mark an item as skipped.
        
        Args:
            item_id: Unique identifier of the skipped item
        """
        item_id_str = str(item_id)
        if item_id_str not in self.data['skipped_ids']:
            self.data['skipped_ids'].add(item_id_str)
            self.data['progress']['skipped_items'] += 1
    
    def is_processed(self, item_id: Union[str, int]) -> bool:
        """
        Check if an item has been successfully processed.
        
        Args:
            item_id: Unique identifier of the item
            
        Returns:
            True if the item has been processed, False otherwise
        """
        return str(item_id) in self.data['processed_ids']
    
    def is_failed(self, item_id: Union[str, int]) -> bool:
        """
        Check if an item has failed processing.
        
        Args:
            item_id: Unique identifier of the item
            
        Returns:
            True if the item has failed, False otherwise
        """
        return str(item_id) in self.data['failed_ids']
    
    def is_skipped(self, item_id: Union[str, int]) -> bool:
        """
        Check if an item has been skipped.
        
        Args:
            item_id: Unique identifier of the item
            
        Returns:
            True if the item has been skipped, False otherwise
        """
        return str(item_id) in self.data['skipped_ids']
    
    def set_total_items(self, total: int) -> None:
        """
        Set the total number of items to be processed.
        
        Args:
            total: Total number of items
        """
        self.data['progress']['total_items'] = total
    
    def update_batch(self, batch: List[Union[str, int]]) -> None:
        """
        Update the last processed batch.
        
        Args:
            batch: List of item IDs in the current batch
        """
        self.data['last_batch'] = [str(item_id) for item_id in batch]
    
    def set_position(self, position: Any) -> None:
        """
        Set the current position in the import process.
        
        Args:
            position: Position indicator (can be page number, offset, etc.)
        """
        self.data['current_position'] = position
    
    def get_position(self) -> Any:
        """
        Get the current position in the import process.
        
        Returns:
            Current position or None if not set
        """
        return self.data['current_position']
    
    def set_custom_state(self, key: str, value: Any) -> None:
        """
        Store custom state information.
        
        Args:
            key: State identifier
            value: State value
        """
        self.data['custom_state'][key] = value
    
    def get_custom_state(self, key: str, default: Any = None) -> Any:
        """
        Retrieve custom state information.
        
        Args:
            key: State identifier
            default: Default value if key not found
            
        Returns:
            Stored state value or default
        """
        return self.data['custom_state'].get(key, default)
    
    def get_unprocessed_from_batch(self) -> List[str]:
        """
        Get IDs from the last batch that haven't been processed.
        
        Returns:
            List of unprocessed item IDs
        """
        return [
            item_id for item_id in self.data['last_batch']
            if item_id not in self.data['processed_ids']
        ]
    
    def get_progress(self) -> Dict[str, int]:
        """
        Get current progress statistics.
        
        Returns:
            Dictionary with progress metrics
        """
        return self.data['progress'].copy()
    
    def get_processed_count(self) -> int:
        """Get number of successfully processed items."""
        return len(self.data['processed_ids'])
    
    def get_failed_count(self) -> int:
        """Get number of failed items."""
        return len(self.data['failed_ids'])
    
    def get_skipped_count(self) -> int:
        """Get number of skipped items."""
        return len(self.data['skipped_ids'])
    
    def get_processed_ids(self) -> Set[str]:
        """Get set of successfully processed item IDs."""
        return self.data['processed_ids'].copy()
    
    def get_failed_ids(self) -> Set[str]:
        """Get set of failed item IDs."""
        return self.data['failed_ids'].copy()
    
    def get_skipped_ids(self) -> Set[str]:
        """Get set of skipped item IDs."""
        return self.data['skipped_ids'].copy()
    
    def clear(self) -> None:
        """Reset checkpoint data to initial state."""
        self.data = {
            'metadata': {
                'created_at': time.time(),
                'last_updated': time.time(),
                'save_count': 0,
                'version': '1.0'
            },
            'progress': {
                'total_items': 0,
                'processed_items': 0,
                'successful_items': 0,
                'failed_items': 0,
                'skipped_items': 0
            },
            'processed_ids': set(),
            'failed_ids': set(),
            'skipped_ids': set(),
            'last_batch': [],
            'current_position': None,
            'custom_state': {}
        }
        
        # Save empty checkpoint
        self.save()
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the checkpoint state.
        
        Returns:
            Dictionary with checkpoint summary
        """
        return {
            'created_at': self.data['metadata']['created_at'],
            'last_updated': self.data['metadata']['last_updated'],
            'total_items': self.data['progress']['total_items'],
            'processed_items': self.data['progress']['processed_items'],
            'successful_items': self.data['progress']['successful_items'],
            'failed_items': self.data['progress']['failed_items'],
            'skipped_items': self.data['progress']['skipped_items'],
            'completion_percentage': (
                (self.data['progress']['processed_items'] / self.data['progress']['total_items']) * 100
                if self.data['progress']['total_items'] > 0 else 0
            ),
            'current_position': self.data['current_position']
        }
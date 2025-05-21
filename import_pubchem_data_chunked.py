#!/usr/bin/env python3
"""
CryoProtect PubChem Data Importer with Chunked Processing

This script implements a robust, adaptive chunked processing algorithm for PubChem data import,
integrating smart chunking, adaptive scheduling, checkpointing, and a circuit breaker pattern.
The goal is to maximize throughput and resilience under severe PubChem API rate limiting.

Usage:
    python import_pubchem_data_chunked.py [--initial-chunk-size SIZE] [--min-chunk-size SIZE] 
                                         [--max-chunk-size SIZE] [--target TARGET]

Parameters:
    --initial-chunk-size: Initial number of compounds to process in each chunk (default: 100)
    --min-chunk-size: Minimum chunk size for adaptation (default: 10)
    --max-chunk-size: Maximum chunk size for adaptation (default: 200)
    --target: Maximum number of compounds to import (default: 5000)
"""

import os
import sys
import time
import json
import argparse
import logging
import statistics
from typing import List, Dict, Any, Tuple, Optional, Iterator, Callable
from pathlib import Path
from datetime import datetime, timedelta
from collections import deque

# Ensure logs and checkpoints directories exist
Path("logs").mkdir(exist_ok=True)
Path("checkpoints").mkdir(exist_ok=True)

# Set up logging
LOG_FILE = "logs/pubchem_chunked_import.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Default configuration
DEFAULT_INITIAL_CHUNK_SIZE = 100
DEFAULT_MIN_CHUNK_SIZE = 10
DEFAULT_MAX_CHUNK_SIZE = 200
DEFAULT_CHUNK_INCREMENT = 10
DEFAULT_TARGET = 5000
DEFAULT_SLIDING_WINDOW_SIZE = 5
DEFAULT_ERROR_THRESHOLD = 0.1  # 10% error rate
DEFAULT_FAST_RESPONSE_THRESHOLD = 1.0  # 1 second

# Checkpoint file path
CHECKPOINT_FILE = "checkpoints/pubchem_chunked_import.json"

class ChunkStats:
    """Class to track statistics for a processed chunk."""
    
    def __init__(self, chunk_id: int, chunk_size: int):
        self.chunk_id = chunk_id
        self.chunk_size = chunk_size
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.response_time: Optional[float] = None
        self.success_count = 0
        self.error_count = 0
        self.skip_count = 0
        self.cids: List[int] = []
        
    def complete(self, success_count: int, error_count: int, skip_count: int):
        """Mark the chunk as complete and calculate response time."""
        self.end_time = time.time()
        self.response_time = self.end_time - self.start_time
        self.success_count = success_count
        self.error_count = error_count
        self.skip_count = skip_count
        
    @property
    def error_rate(self) -> float:
        """Calculate the error rate for this chunk."""
        total = self.success_count + self.error_count
        return self.error_count / total if total > 0 else 0
    
    @property
    def success_rate(self) -> float:
        """Calculate the success rate for this chunk."""
        total = self.success_count + self.error_count
        return self.success_count / total if total > 0 else 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the stats to a dictionary for serialization."""
        return {
            "chunk_id": self.chunk_id,
            "chunk_size": self.chunk_size,
            "start_time": self.start_time,
            "end_time": self.end_time,
            "response_time": self.response_time,
            "success_count": self.success_count,
            "error_count": self.error_count,
            "skip_count": self.skip_count,
            "error_rate": self.error_rate,
            "success_rate": self.success_rate,
            "cids": self.cids
        }


class ChunkGenerator:
    """
    Generator that splits a list of CIDs into chunks with adaptive sizing.
    
    The chunk size is dynamically adjusted based on recent response times and error rates.
    """
    
    def __init__(
        self,
        cids: List[int],
        initial_chunk_size: int = DEFAULT_INITIAL_CHUNK_SIZE,
        min_chunk_size: int = DEFAULT_MIN_CHUNK_SIZE,
        max_chunk_size: int = DEFAULT_MAX_CHUNK_SIZE,
        chunk_increment: int = DEFAULT_CHUNK_INCREMENT,
        sliding_window_size: int = DEFAULT_SLIDING_WINDOW_SIZE,
        error_threshold: float = DEFAULT_ERROR_THRESHOLD,
        fast_response_threshold: float = DEFAULT_FAST_RESPONSE_THRESHOLD
    ):
        """
        Initialize the chunk generator.
        
        Args:
            cids: List of PubChem CIDs to process
            initial_chunk_size: Initial size of each chunk
            min_chunk_size: Minimum chunk size
            max_chunk_size: Maximum chunk size
            chunk_increment: Amount to increase/decrease chunk size by
            sliding_window_size: Number of recent chunks to consider for adaptation
            error_threshold: Error rate threshold above which chunk size is decreased
            fast_response_threshold: Response time threshold below which chunk size is increased
        """
        self.cids = cids
        self.initial_chunk_size = initial_chunk_size
        self.min_chunk_size = min_chunk_size
        self.max_chunk_size = max_chunk_size
        self.chunk_increment = chunk_increment
        self.sliding_window_size = sliding_window_size
        self.error_threshold = error_threshold
        self.fast_response_threshold = fast_response_threshold
        
        self.current_chunk_size = initial_chunk_size
        self.current_position = 0
        self.chunk_id = 0
        self.recent_stats: deque = deque(maxlen=sliding_window_size)
        
    def __iter__(self) -> Iterator[Tuple[int, List[int], ChunkStats]]:
        """Return self as iterator."""
        return self
    
    def __next__(self) -> Tuple[int, List[int], ChunkStats]:
        """
        Get the next chunk of CIDs.
        
        Returns:
            Tuple containing:
            - chunk_id: Unique identifier for this chunk
            - chunk_cids: List of CIDs in this chunk
            - chunk_stats: ChunkStats object for tracking statistics
            
        Raises:
            StopIteration: When all CIDs have been processed
        """
        if self.current_position >= len(self.cids):
            raise StopIteration
        
        # Calculate end position for current chunk
        end_position = min(self.current_position + self.current_chunk_size, len(self.cids))
        
        # Get CIDs for this chunk
        chunk_cids = self.cids[self.current_position:end_position]
        
        # Create stats object for this chunk
        chunk_stats = ChunkStats(self.chunk_id, len(chunk_cids))
        chunk_stats.cids = chunk_cids
        
        # Prepare return values
        result = (self.chunk_id, chunk_cids, chunk_stats)
        
        # Update state for next iteration
        self.current_position = end_position
        self.chunk_id += 1
        
        return result
    
    def update_chunk_size(self, chunk_stats: ChunkStats) -> None:
        """
        Update the chunk size based on the statistics from a processed chunk.
        
        Args:
            chunk_stats: Statistics from a processed chunk
        """
        # Add the stats to our recent history
        self.recent_stats.append(chunk_stats)
        
        # Only adapt if we have enough history
        if len(self.recent_stats) < 2:
            return
        
        # Calculate new chunk size based on recent performance
        self.current_chunk_size = self.get_optimal_chunk_size(
            [stats.response_time for stats in self.recent_stats if stats.response_time is not None],
            [stats.error_rate for stats in self.recent_stats],
            self.current_chunk_size
        )
        
        logger.info(f"Adapted chunk size to {self.current_chunk_size} based on recent performance")
    
    def get_optimal_chunk_size(
        self,
        recent_response_times: List[float],
        recent_error_rates: List[float],
        current_chunk_size: int
    ) -> int:
        """
        Calculate the optimal chunk size based on recent performance metrics.
        
        Args:
            recent_response_times: List of response times from recent chunks
            recent_error_rates: List of error rates from recent chunks
            current_chunk_size: Current chunk size
            
        Returns:
            Optimal chunk size for the next chunk
        """
        # If we don't have enough data, return the current size
        if not recent_response_times or not recent_error_rates:
            return current_chunk_size
        
        try:
            # Calculate mean error rate
            mean_error_rate = statistics.mean(recent_error_rates)
            
            # If error rate is too high, decrease chunk size
            if mean_error_rate > self.error_threshold:
                return max(current_chunk_size // 2, self.min_chunk_size)
            
            # Calculate mean response time
            mean_response_time = statistics.mean(recent_response_times)
            
            # If response time is fast, increase chunk size
            if mean_response_time < self.fast_response_threshold:
                return min(current_chunk_size + self.chunk_increment, self.max_chunk_size)
            
            # Otherwise, keep the current size
            return current_chunk_size
            
        except statistics.StatisticsError:
            # Handle case where statistics calculation fails
            logger.warning("Failed to calculate statistics for chunk size adaptation")
            return current_chunk_size
    
    def reset(self) -> None:
        """Reset the generator to start from the beginning."""
        self.current_position = 0
        self.chunk_id = 0
        self.recent_stats.clear()
    
    def skip_to_position(self, position: int) -> None:
        """
        Skip to a specific position in the CID list.
        
        Args:
            position: Position to skip to
        """
        self.current_position = min(position, len(self.cids))
    
    def get_progress(self) -> float:
        """
        Get the current progress as a percentage.
        
        Returns:
            Percentage of CIDs processed (0-100)
        """
        if not self.cids:
            return 100.0
        return (self.current_position / len(self.cids)) * 100


def get_cid_list(cid_file: str = "CID-Synonym-curated") -> List[int]:
    """
    Retrieve all CIDs from the curated PubChem CID list.
    
    Args:
        cid_file: Path to the file containing CIDs
        
    Returns:
        List of PubChem CIDs
    """
    if not os.path.exists(cid_file):
        logger.warning(f"WARNING: CID file '{cid_file}' not found.")
        return []

    with open(cid_file, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids


import shutil

class CheckpointManager:
    """
    Manages checkpoints for resumable PubChem data import.
    
    This class handles saving and loading checkpoint data, allowing the import process
    to be resumed after interruption. It provides robust error handling, checkpoint
    verification, and detailed status information.
    
    Features:
    - Saves checkpoint data after each chunk is processed
    - Loads checkpoint data to resume from the last successful point
    - Handles corrupted checkpoint files gracefully
    - Supports custom checkpoint file paths
    - Provides detailed checkpoint status information
    - Maintains backup of previous checkpoint
    """
    
    def __init__(self, checkpoint_file: str = CHECKPOINT_FILE):
        """
        Initialize the checkpoint manager.
        
        Args:
            checkpoint_file: Path to the checkpoint file
        """
        self.checkpoint_file = checkpoint_file
        self.backup_file = f"{checkpoint_file}.bak"
        self.last_save_time = 0
        self.save_count = 0
        
    def save(
        self,
        chunk_generator: ChunkGenerator,
        processed_cids: List[int],
        chunk_stats_history: List[Dict[str, Any]],
        start_time: float,
        status: str = "Running"
    ) -> bool:
        """
        Save the current import state to a checkpoint file.
        
        Args:
            chunk_generator: The chunk generator being used
            processed_cids: List of CIDs that have been processed
            chunk_stats_history: History of chunk statistics
            start_time: Time when the import started
            status: Current status of the import
            
        Returns:
            True if checkpoint was saved successfully, False otherwise
        """
        # Create checkpoint data
        checkpoint_data = {
            "current_position": chunk_generator.current_position,
            "current_chunk_size": chunk_generator.current_chunk_size,
            "processed_cids": processed_cids,
            "chunk_stats_history": chunk_stats_history,
            "elapsed_seconds": time.time() - start_time,
            "timestamp": datetime.now().isoformat(),
            "status": status,
            "progress_percentage": chunk_generator.get_progress(),
            "version": "1.0",  # For future compatibility
            "metadata": {
                "total_cids": len(chunk_generator.cids),
                "save_count": self.save_count + 1,
                "last_chunk_id": chunk_stats_history[-1]["chunk_id"] if chunk_stats_history else None,
                "last_chunk_size": chunk_stats_history[-1]["chunk_size"] if chunk_stats_history else None,
                "last_error_rate": chunk_stats_history[-1]["error_rate"] if chunk_stats_history else None
            }
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
            
            logger.info(f"Checkpoint saved: {len(processed_cids)} CIDs processed, "
                        f"progress: {chunk_generator.get_progress():.1f}%, "
                        f"status: {status}")
            return True
            
        except Exception as e:
            logger.error(f"Error saving checkpoint: {str(e)}")
            return False
    
    def load(self) -> Optional[Dict[str, Any]]:
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
    
    def _load_file(self, file_path: str) -> Optional[Dict[str, Any]]:
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
            required_keys = ["current_position", "current_chunk_size", "processed_cids",
                            "chunk_stats_history", "status", "progress_percentage"]
            
            if not all(key in checkpoint for key in required_keys):
                logger.error(f"Checkpoint file {file_path} is missing required keys")
                return None
            
            # Log checkpoint details
            logger.info(f"Checkpoint loaded from {file_path}: "
                        f"{len(checkpoint.get('processed_cids', []))} CIDs processed, "
                        f"progress: {checkpoint.get('progress_percentage', 0):.1f}%, "
                        f"status: {checkpoint.get('status', 'Unknown')}")
            
            return checkpoint
            
        except json.JSONDecodeError as e:
            logger.error(f"Error parsing checkpoint file {file_path}: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error loading checkpoint from {file_path}: {str(e)}")
            return None
    
    def get_status(self) -> Dict[str, Any]:
        """
        Get the status of the checkpoint.
        
        Returns:
            Dictionary with checkpoint status information
        """
        status = {
            "checkpoint_file": self.checkpoint_file,
            "backup_file": self.backup_file,
            "checkpoint_exists": os.path.exists(self.checkpoint_file),
            "backup_exists": os.path.exists(self.backup_file),
            "last_save_time": self.last_save_time,
            "save_count": self.save_count
        }
        
        # Add checkpoint data if available
        if status["checkpoint_exists"]:
            try:
                with open(self.checkpoint_file, "r") as f:
                    checkpoint = json.load(f)
                
                status.update({
                    "timestamp": checkpoint.get("timestamp"),
                    "status": checkpoint.get("status"),
                    "progress_percentage": checkpoint.get("progress_percentage"),
                    "processed_count": len(checkpoint.get("processed_cids", [])),
                    "current_position": checkpoint.get("current_position"),
                    "current_chunk_size": checkpoint.get("current_chunk_size")
                })
            except Exception:
                status.update({
                    "checkpoint_corrupted": True
                })
        
        return status
    
    def clear(self) -> bool:
        """
        Clear checkpoint files.
        
        Returns:
            True if checkpoint files were cleared successfully, False otherwise
        """
        success = True
        
        # Remove main checkpoint file
        if os.path.exists(self.checkpoint_file):
            try:
                os.remove(self.checkpoint_file)
                logger.info(f"Removed checkpoint file: {self.checkpoint_file}")
            except Exception as e:
                logger.error(f"Error removing checkpoint file: {str(e)}")
                success = False
        
        # Remove backup file
        if os.path.exists(self.backup_file):
            try:
                os.remove(self.backup_file)
                logger.info(f"Removed backup checkpoint file: {self.backup_file}")
            except Exception as e:
                logger.error(f"Error removing backup checkpoint file: {str(e)}")
                success = False
        
        return success


def save_checkpoint(
    chunk_generator: ChunkGenerator,
    processed_cids: List[int],
    chunk_stats_history: List[Dict[str, Any]],
    start_time: float,
    status: str = "Running"
) -> None:
    """
    Save the current import state to a checkpoint file.
    
    This is a compatibility wrapper around the CheckpointManager.save method.
    
    Args:
        chunk_generator: The chunk generator being used
        processed_cids: List of CIDs that have been processed
        chunk_stats_history: History of chunk statistics
        start_time: Time when the import started
        status: Current status of the import
    """
    checkpoint_manager = CheckpointManager()
    checkpoint_manager.save(chunk_generator, processed_cids, chunk_stats_history, start_time, status)


def load_checkpoint() -> Optional[Dict[str, Any]]:
    """
    Load the last checkpoint if it exists.
    
    This is a compatibility wrapper around the CheckpointManager.load method.
    
    Returns:
        Dictionary with checkpoint data if found, None otherwise
    """
    checkpoint_manager = CheckpointManager()
    return checkpoint_manager.load()


class ChunkProcessor:
    """
    Processes chunks of CIDs, handling API calls, errors, and retries.
    
    Features:
    - Processes each chunk with configurable delay between chunks
    - Handles API calls with circuit breaker protection
    - Implements retry with exponential backoff
    - Tracks detailed statistics for each chunk
    - Supports checkpointing for resumability
    - Respects circuit breaker state for resilient processing
    """
    
    def __init__(
        self,
        pubchem_client=None,
        base_delay: float = 0.5,
        max_delay: float = 10.0,
        max_retries: int = 3,
        failure_threshold: int = 5,
        recovery_timeout: int = 30,
        backoff_factor: float = 2.0,
        jitter: bool = True
    ):
        """
        Initialize the chunk processor.
        
        Args:
            pubchem_client: PubChem client to use for API calls (if None, creates a new one)
            base_delay: Base delay between chunks in seconds
            max_delay: Maximum delay between chunks in seconds
            max_retries: Maximum number of retries for failed API calls
            failure_threshold: Number of failures before opening the circuit
            recovery_timeout: Time in seconds to wait before trying again
            backoff_factor: Factor to increase delay for each retry
            jitter: Whether to add random jitter to delay
        """
        from pubchem.client import ResilientPubChemClient
        from pubchem.utils import CircuitBreakerError
        
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.max_retries = max_retries
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.backoff_factor = backoff_factor
        self.jitter = jitter
        self.CircuitBreakerError = CircuitBreakerError
        
        # Initialize PubChem client if not provided
        self.pubchem_client = pubchem_client or ResilientPubChemClient(
            max_retries=max_retries,
            failure_threshold=failure_threshold,
            recovery_timeout=recovery_timeout
        )
        
        # Statistics
        self.total_success = 0
        self.total_error = 0
        self.total_skip = 0
        self.total_time = 0.0
        self.current_delay = base_delay
        self.consecutive_failures = 0
        
    def check_circuit_breaker(self) -> bool:
        """
        Check if the circuit breaker is open.
        
        Returns:
            True if the circuit is closed or half-open (requests allowed),
            False if the circuit is open (requests blocked)
        """
        try:
            circuit_stats = self.pubchem_client.get_circuit_breaker_stats()
            if circuit_stats["state"] == "open":
                recovery_time = self.recovery_timeout - (time.time() - circuit_stats["last_failure_time"])
                if recovery_time > 0:
                    logger.warning(
                        f"Circuit breaker is OPEN. Waiting {recovery_time:.1f}s before retrying."
                    )
                    return False
        except Exception as e:
            logger.warning(f"Error checking circuit breaker state: {str(e)}")
        
        return True
    
    def exponential_backoff(self, retry_count: int) -> float:
        """
        Calculate delay using exponential backoff with optional jitter.
        
        Args:
            retry_count: Current retry attempt number
            
        Returns:
            Delay in seconds
        """
        import random
        
        # Calculate base delay
        delay = self.base_delay * (self.backoff_factor ** retry_count)
        delay = min(delay, self.max_delay)
        
        # Add jitter if enabled
        if self.jitter:
            delay = delay * (0.5 + random.random())
            
        return delay
    
    def process_chunk(
        self,
        chunk_id: int,
        chunk_cids: List[int],
        delay: Optional[float] = None
    ) -> ChunkStats:
        """
        Process a single chunk of CIDs.
        
        Args:
            chunk_id: Unique identifier for this chunk
            chunk_cids: List of CIDs in this chunk
            delay: Delay before processing this chunk (if None, uses current_delay)
            
        Returns:
            ChunkStats object with processing statistics
        """
        # Create stats object for this chunk
        chunk_stats = ChunkStats(chunk_id, len(chunk_cids))
        chunk_stats.cids = chunk_cids
        
        # Check circuit breaker before processing
        if not self.check_circuit_breaker():
            # Wait for recovery timeout and then try again
            logger.info(f"Waiting {self.recovery_timeout}s for circuit breaker to reset...")
            time.sleep(self.recovery_timeout)
            
            # Check again
            if not self.check_circuit_breaker():
                logger.error("Circuit breaker still open after waiting. Marking chunk as failed.")
                chunk_stats.complete(0, len(chunk_cids), 0)
                return chunk_stats
        
        # Apply delay if specified
        if delay is not None:
            logger.info(f"Waiting {delay:.2f}s before processing chunk {chunk_id}...")
            time.sleep(delay)
        elif self.current_delay > 0:
            logger.info(f"Waiting {self.current_delay:.2f}s before processing chunk {chunk_id}...")
            time.sleep(self.current_delay)
        
        # Process each CID in the chunk
        success_count = 0
        error_count = 0
        skip_count = 0
        
        logger.info(f"Processing chunk {chunk_id} with {len(chunk_cids)} CIDs...")
        
        for cid in chunk_cids:
            retry_count = 0
            while retry_count <= self.max_retries:
                try:
                    # Check circuit breaker before each API call
                    if not self.check_circuit_breaker() and retry_count < self.max_retries:
                        # Wait and retry
                        backoff_time = self.exponential_backoff(retry_count)
                        logger.warning(f"Circuit breaker open, retrying in {backoff_time:.2f}s...")
                        time.sleep(backoff_time)
                        retry_count += 1
                        continue
                    
                    # Get molecule properties from PubChem
                    result = self.pubchem_client.get_molecule_properties(
                        cid=cid,
                        use_cache=True,
                        fallback_to_cache=True
                    )
                    
                    # Check if there was an error
                    if "Error" in result:
                        logger.warning(f"Error processing CID {cid}: {result['Error']}")
                        error_count += 1
                        self.consecutive_failures += 1
                    else:
                        logger.debug(f"Successfully processed CID {cid}")
                        success_count += 1
                        self.consecutive_failures = 0
                    
                    # Break out of retry loop on success or handled error
                    break
                    
                except self.CircuitBreakerError as e:
                    # Circuit breaker is open, wait and retry if retries left
                    if retry_count < self.max_retries:
                        backoff_time = self.exponential_backoff(retry_count)
                        logger.warning(f"Circuit breaker error for CID {cid}, retrying in {backoff_time:.2f}s: {str(e)}")
                        time.sleep(backoff_time)
                        retry_count += 1
                    else:
                        logger.error(f"Circuit breaker error for CID {cid} after {self.max_retries} retries: {str(e)}")
                        error_count += 1
                        self.consecutive_failures += 1
                        break
                        
                except Exception as e:
                    # Other exceptions, retry with backoff if retries left
                    if retry_count < self.max_retries:
                        backoff_time = self.exponential_backoff(retry_count)
                        logger.warning(f"Error processing CID {cid}, retrying in {backoff_time:.2f}s: {str(e)}")
                        time.sleep(backoff_time)
                        retry_count += 1
                    else:
                        logger.error(f"Exception processing CID {cid} after {self.max_retries} retries: {str(e)}")
                        error_count += 1
                        self.consecutive_failures += 1
                        break
        
        # Complete the stats
        chunk_stats.complete(success_count, error_count, skip_count)
        
        # Update overall statistics
        self.total_success += success_count
        self.total_error += error_count
        self.total_skip += skip_count
        self.total_time += chunk_stats.response_time or 0
        
        # Log chunk results
        logger.info(
            f"Chunk {chunk_id} completed in {chunk_stats.response_time:.2f}s: "
            f"{success_count} success, {error_count} errors, {skip_count} skipped"
        )
        
        # Adjust delay based on error rate and consecutive failures
        self._adjust_delay(chunk_stats.error_rate)
        
        return chunk_stats
    
    def _adjust_delay(self, error_rate: float) -> None:
        """
        Adjust the delay between chunks based on error rate and consecutive failures.
        
        Args:
            error_rate: Error rate from the last chunk
        """
        # Check consecutive failures first - this takes precedence
        if self.consecutive_failures >= self.failure_threshold:
            # Set to max_delay directly due to consecutive failures
            self.current_delay = self.max_delay
            logger.warning(
                f"Circuit breaker threshold reached ({self.consecutive_failures} consecutive failures), "
                f"setting delay to maximum {self.current_delay:.2f}s"
            )
            return
            
        # Otherwise adjust based on error rate
        if error_rate >= 1.0:  # Complete failure
            # Set to max_delay directly
            self.current_delay = self.max_delay
            logger.info(f"Complete failure ({error_rate:.2f}), setting delay to maximum {self.current_delay:.2f}s")
        elif error_rate > 0.5:  # High error rate
            # Double the delay, up to max_delay
            self.current_delay = min(self.current_delay * 2, self.max_delay)
            logger.info(f"High error rate ({error_rate:.2f}), increasing delay to {self.current_delay:.2f}s")
        elif error_rate > 0.2:  # Moderate error rate
            # Increase delay by 50%, up to max_delay
            self.current_delay = min(self.current_delay * 1.5, self.max_delay)
            logger.info(f"Moderate error rate ({error_rate:.2f}), increasing delay to {self.current_delay:.2f}s")
        elif error_rate < 0.05 and self.current_delay > self.base_delay:  # Low error rate
            # Decrease delay by 10%, down to base_delay
            self.current_delay = max(self.current_delay * 0.9, self.base_delay)
            logger.info(f"Low error rate ({error_rate:.2f}), decreasing delay to {self.current_delay:.2f}s")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get overall processing statistics.
        
        Returns:
            Dictionary with processing statistics
        """
        total_processed = self.total_success + self.total_error + self.total_skip
        
        stats = {
            "total_success": self.total_success,
            "total_error": self.total_error,
            "total_skip": self.total_skip,
            "total_processed": total_processed,
            "success_rate": self.total_success / total_processed if total_processed > 0 else 0,
            "error_rate": self.total_error / total_processed if total_processed > 0 else 0,
            "total_time": self.total_time,
            "current_delay": self.current_delay,
            "consecutive_failures": self.consecutive_failures
        }
        
        # Add circuit breaker stats if available
        try:
            circuit_stats = self.pubchem_client.get_circuit_breaker_stats()
            stats["circuit_breaker"] = circuit_stats
        except Exception:
            pass
            
        return stats


def process_cids_with_chunking(
    cids: List[int],
    initial_chunk_size: int = DEFAULT_INITIAL_CHUNK_SIZE,
    min_chunk_size: int = DEFAULT_MIN_CHUNK_SIZE,
    max_chunk_size: int = DEFAULT_MAX_CHUNK_SIZE,
    chunk_increment: int = DEFAULT_CHUNK_INCREMENT,
    base_delay: float = 0.5,
    max_delay: float = 10.0,
    max_retries: int = 3,
    failure_threshold: int = 5,
    recovery_timeout: int = 30,
    backoff_factor: float = 2.0,
    jitter: bool = True,
    checkpoint_interval: int = 1,
    resume: bool = True,
    target_count: Optional[int] = None,
    checkpoint_file: str = CHECKPOINT_FILE
) -> Dict[str, Any]:
    """
    Process a list of CIDs using adaptive chunked processing.
    
    Args:
        cids: List of PubChem CIDs to process
        initial_chunk_size: Initial size of each chunk
        min_chunk_size: Minimum chunk size
        max_chunk_size: Maximum chunk size
        chunk_increment: Amount to increase/decrease chunk size by
        base_delay: Base delay between chunks in seconds
        max_delay: Maximum delay between chunks in seconds
        max_retries: Maximum number of retries for failed API calls
        failure_threshold: Number of failures before opening the circuit
        recovery_timeout: Time in seconds to wait before trying again
        checkpoint_interval: Number of chunks to process before saving a checkpoint
        resume: Whether to resume from the last checkpoint if available
        target_count: Maximum number of CIDs to process (if None, processes all)
        checkpoint_file: Path to the checkpoint file
        
    Returns:
        Dictionary with processing results and statistics
    """
    # Limit to target count if specified
    if target_count is not None and target_count < len(cids):
        cids = cids[:target_count]
        logger.info(f"Limited to {target_count} CIDs as per target count")
    
    # Initialize components
    chunk_generator = ChunkGenerator(
        cids=cids,
        initial_chunk_size=initial_chunk_size,
        min_chunk_size=min_chunk_size,
        max_chunk_size=max_chunk_size,
        chunk_increment=chunk_increment
    )
    
    chunk_processor = ChunkProcessor(
        base_delay=base_delay,
        max_delay=max_delay,
        max_retries=max_retries,
        failure_threshold=failure_threshold,
        recovery_timeout=recovery_timeout,
        backoff_factor=backoff_factor,
        jitter=jitter
    )
    
    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager(checkpoint_file)
    
    # Initialize tracking variables
    processed_cids = []
    chunk_stats_history = []
    start_time = time.time()
    
    # Try to load checkpoint if resume is enabled
    if resume:
        checkpoint = checkpoint_manager.load()
        if checkpoint:
            # Resume from checkpoint
            chunk_generator.skip_to_position(checkpoint["current_position"])
            chunk_generator.current_chunk_size = checkpoint["current_chunk_size"]
            processed_cids = checkpoint["processed_cids"]
            chunk_stats_history = checkpoint["chunk_stats_history"]
            start_time = time.time() - checkpoint["elapsed_seconds"]
            
            logger.info(
                f"Resuming from checkpoint: {len(processed_cids)} CIDs processed, "
                f"progress: {checkpoint['progress_percentage']:.1f}%, "
                f"status: {checkpoint.get('status', 'Unknown')}"
            )
    
    # Process chunks
    chunk_count = 0
    last_checkpoint_time = time.time()
    checkpoint_min_interval_seconds = 5  # Minimum seconds between checkpoints
    
    try:
        for chunk_id, chunk_cids, chunk_stats in chunk_generator:
            # Process the chunk
            chunk_stats = chunk_processor.process_chunk(chunk_id, chunk_cids)
            
            # Update tracking
            processed_cids.extend(chunk_cids)
            chunk_stats_history.append(chunk_stats.to_dict())
            
            # Update chunk generator with stats for adaptation
            chunk_generator.update_chunk_size(chunk_stats)
            
            # Save checkpoint periodically
            chunk_count += 1
            current_time = time.time()
            time_since_last_checkpoint = current_time - last_checkpoint_time
            
            if (chunk_count % checkpoint_interval == 0 and
                time_since_last_checkpoint >= checkpoint_min_interval_seconds):
                checkpoint_manager.save(chunk_generator, processed_cids, chunk_stats_history, start_time)
                last_checkpoint_time = current_time
            
            # Check for keyboard interrupt
            if hasattr(time, 'sleep'):
                time.sleep(0.01)  # Small sleep to allow keyboard interrupt
    
    except KeyboardInterrupt:
        logger.info("Processing interrupted by user")
        checkpoint_manager.save(chunk_generator, processed_cids, chunk_stats_history, start_time, "Interrupted")
    
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        checkpoint_manager.save(chunk_generator, processed_cids, chunk_stats_history, start_time, "Error")
        raise
    
    # Save final checkpoint
    checkpoint_manager.save(chunk_generator, processed_cids, chunk_stats_history, start_time, "Completed")
    
    # Calculate overall statistics
    elapsed_time = time.time() - start_time
    processor_stats = chunk_processor.get_stats()
    
    # Prepare result
    result = {
        "processed_cids": processed_cids,
        "total_cids": len(cids),
        "processed_count": len(processed_cids),
        "success_count": processor_stats["total_success"],
        "error_count": processor_stats["total_error"],
        "skip_count": processor_stats["total_skip"],
        "success_rate": processor_stats["success_rate"],
        "elapsed_time": elapsed_time,
        "chunks_processed": chunk_count,
        "final_chunk_size": chunk_generator.current_chunk_size,
        "checkpoint_status": checkpoint_manager.get_status()
    }
    
    logger.info(
        f"Processing completed: {result['processed_count']}/{result['total_cids']} CIDs processed, "
        f"{result['success_count']} successful ({result['success_rate']:.1%}), "
        f"{result['error_count']} errors, "
        f"{result['skip_count']} skipped, "
        f"in {elapsed_time:.1f}s"
    )
    
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import PubChem data with adaptive chunked processing.")
    parser.add_argument("--initial-chunk-size", type=int, default=DEFAULT_INITIAL_CHUNK_SIZE,
                        help=f"Initial chunk size (default: {DEFAULT_INITIAL_CHUNK_SIZE})")
    parser.add_argument("--min-chunk-size", type=int, default=DEFAULT_MIN_CHUNK_SIZE,
                        help=f"Minimum chunk size (default: {DEFAULT_MIN_CHUNK_SIZE})")
    parser.add_argument("--max-chunk-size", type=int, default=DEFAULT_MAX_CHUNK_SIZE,
                        help=f"Maximum chunk size (default: {DEFAULT_MAX_CHUNK_SIZE})")
    parser.add_argument("--target", type=int, default=DEFAULT_TARGET,
                        help=f"Target number of compounds to import (default: {DEFAULT_TARGET})")
    parser.add_argument("--no-resume", action="store_true",
                        help="Don't resume from checkpoint")
    parser.add_argument("--checkpoint-interval", type=int, default=1,
                        help="Number of chunks to process before saving a checkpoint (default: 1)")
    parser.add_argument("--base-delay", type=float, default=0.5,
                        help="Base delay between chunks in seconds (default: 0.5)")
    parser.add_argument("--max-delay", type=float, default=10.0,
                        help="Maximum delay between chunks in seconds (default: 10.0)")
    parser.add_argument("--checkpoint-file", type=str, default=CHECKPOINT_FILE,
                        help=f"Path to checkpoint file (default: {CHECKPOINT_FILE})")
    parser.add_argument("--clear-checkpoint", action="store_true",
                        help="Clear existing checkpoint before starting")
    args = parser.parse_args()
    
    # Get CIDs from file
    cids = get_cid_list()
    if not cids:
        logger.error("No CIDs found. Exiting.")
        sys.exit(1)
    
    # Clear checkpoint if requested
    if args.clear_checkpoint:
        checkpoint_manager = CheckpointManager(args.checkpoint_file)
        if checkpoint_manager.clear():
            logger.info("Checkpoint cleared. Starting from the beginning.")
    
    logger.info(f"Starting chunked processing of {len(cids)} CIDs...")
    
    # Process CIDs with chunking
    result = process_cids_with_chunking(
        cids=cids,
        initial_chunk_size=args.initial_chunk_size,
        min_chunk_size=args.min_chunk_size,
        max_chunk_size=args.max_chunk_size,
        target_count=args.target,
        resume=not args.no_resume,
        checkpoint_interval=args.checkpoint_interval,
        base_delay=args.base_delay,
        max_delay=args.max_delay,
        backoff_factor=2.0,
        jitter=True,
        checkpoint_file=args.checkpoint_file
    )
    
    # Print summary
    print("\nProcessing Summary:")
    print(f"Total CIDs: {result['total_cids']}")
    print(f"Processed: {result['processed_count']}")
    print(f"Successful: {result['success_count']} ({result['success_rate']:.1%})")
    print(f"Errors: {result['error_count']}")
    print(f"Skipped: {result['skip_count']}")
    print(f"Elapsed time: {result['elapsed_time']:.1f}s")
    print(f"Chunks processed: {result['chunks_processed']}")
    print(f"Final chunk size: {result['final_chunk_size']}")
    
    # Print checkpoint status
    if "checkpoint_status" in result:
        print("\nCheckpoint Status:")
        status = result["checkpoint_status"]
        print(f"Checkpoint file: {status['checkpoint_file']}")
        print(f"Checkpoint exists: {status['checkpoint_exists']}")
        print(f"Backup exists: {status['backup_exists']}")
        print(f"Save count: {status['save_count']}")
        if status.get("timestamp"):
            print(f"Last saved: {status['timestamp']}")
            print(f"Status: {status.get('status', 'Unknown')}")
            print(f"Progress: {status.get('progress_percentage', 0):.1f}%")
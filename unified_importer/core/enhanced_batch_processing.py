"""
Enhanced batch processing with adaptive sizing and parallel processing.

This module provides advanced batch processing capabilities for the unified molecular
importer with features like:
- Adaptive batch sizing based on memory usage and processing time
- Parallel processing of data transformations
- Batching strategies optimized for different data sources
- Performance monitoring and automatic tuning
"""

import os
import time
import asyncio
import logging
import threading
import queue
import uuid
import psutil
import concurrent.futures
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Tuple, Union, Set, Callable, AsyncIterator, TypeVar, Generic

# Type variables for generic functions
T = TypeVar('T')  # Input type
U = TypeVar('U')  # Output type


@dataclass
class BatchStats:
    """Statistics about batch processing performance."""
    batch_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    batch_size: int = 0
    items_processed: int = 0
    start_time: float = field(default_factory=time.time)
    end_time: Optional[float] = None
    memory_usage_start: int = field(default_factory=lambda: psutil.Process(os.getpid()).memory_info().rss)
    memory_usage_end: Optional[int] = None
    success_count: int = 0
    error_count: int = 0
    processing_times: List[float] = field(default_factory=list)
    
    @property
    def duration(self) -> float:
        """Get the duration of batch processing in seconds."""
        if self.end_time is None:
            return time.time() - self.start_time
        return self.end_time - self.start_time
    
    @property
    def memory_usage_delta(self) -> Optional[int]:
        """Get the change in memory usage during batch processing."""
        if self.memory_usage_end is None:
            return None
        return self.memory_usage_end - self.memory_usage_start
    
    @property
    def items_per_second(self) -> float:
        """Get the processing rate in items per second."""
        if self.items_processed == 0 or self.duration == 0:
            return 0
        return self.items_processed / self.duration
    
    @property
    def avg_processing_time(self) -> float:
        """Get the average processing time per item in seconds."""
        if not self.processing_times:
            return 0
        return sum(self.processing_times) / len(self.processing_times)
    
    def complete(self) -> None:
        """Mark the batch as complete and record final statistics."""
        self.end_time = time.time()
        self.memory_usage_end = psutil.Process(os.getpid()).memory_info().rss


class BatchProcessingStrategy:
    """Base strategy for batch processing."""
    
    def __init__(
        self,
        initial_batch_size: int = 100,
        min_batch_size: int = 10,
        max_batch_size: int = 1000,
        target_duration: float = 5.0,
        memory_limit: int = 1024 * 1024 * 1024,  # 1GB
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the batch processing strategy.
        
        Args:
            initial_batch_size: Initial size of each batch
            min_batch_size: Minimum batch size
            max_batch_size: Maximum batch size
            target_duration: Target duration for batch processing in seconds
            memory_limit: Memory usage limit in bytes
            logger: Logger instance
        """
        self.initial_batch_size = initial_batch_size
        self.min_batch_size = min_batch_size
        self.max_batch_size = max_batch_size
        self.current_batch_size = initial_batch_size
        self.target_duration = target_duration
        self.memory_limit = memory_limit
        self.logger = logger or logging.getLogger(__name__)
        
        # Statistics
        self.batch_stats = []
        self.total_items_processed = 0
        self.total_success_count = 0
        self.total_error_count = 0
        self.total_duration = 0
        
        # Lock for thread safety
        self.lock = threading.RLock()
    
    def get_next_batch_size(self) -> int:
        """
        Get the next batch size based on processing history.
        
        The strategy adjusts batch size based on:
        - Processing time relative to target duration
        - Memory usage relative to limit
        - Error rates
        
        Returns:
            Next recommended batch size
        """
        with self.lock:
            # If we don't have enough statistics, use the initial size
            if len(self.batch_stats) < 3:
                return self.current_batch_size
            
            # Get the most recent batch stats
            recent_stats = self.batch_stats[-3:]
            
            # Calculate average processing time and memory usage
            avg_duration = sum(stat.duration for stat in recent_stats) / len(recent_stats)
            avg_memory_delta = sum(stat.memory_usage_delta or 0 for stat in recent_stats) / len(recent_stats)
            avg_error_rate = sum(stat.error_count / max(1, stat.items_processed) for stat in recent_stats) / len(recent_stats)
            
            # Determine adjustment factor based on target duration
            time_factor = self.target_duration / max(0.001, avg_duration)
            
            # Determine adjustment factor based on memory usage
            memory_factor = 1.0
            if avg_memory_delta > 0:
                # Estimate memory for target batch size
                items_per_batch = sum(stat.items_processed for stat in recent_stats) / len(recent_stats)
                memory_per_item = avg_memory_delta / items_per_batch
                estimated_memory = memory_per_item * self.current_batch_size
                
                # Adjust based on memory limit
                memory_factor = min(1.0, self.memory_limit / max(1, estimated_memory) * 0.8)
            
            # Adjust based on error rate
            error_factor = 1.0 - min(0.5, avg_error_rate * 2)  # Reduce by up to 50% based on errors
            
            # Combine factors with weights
            combined_factor = time_factor * 0.6 + memory_factor * 0.3 + error_factor * 0.1
            
            # Apply the adjustment, but constrain to avoid wild swings
            # Limit change to 20% up or down
            adjustment = max(0.8, min(1.2, combined_factor))
            new_batch_size = int(self.current_batch_size * adjustment)
            
            # Constrain to min/max limits
            new_batch_size = max(self.min_batch_size, min(self.max_batch_size, new_batch_size))
            
            self.logger.debug(
                f"Batch size adjustment: {self.current_batch_size} -> {new_batch_size} "
                f"(time: {time_factor:.2f}, memory: {memory_factor:.2f}, error: {error_factor:.2f})"
            )
            
            # Update current batch size
            self.current_batch_size = new_batch_size
            return new_batch_size
    
    def record_batch_stats(self, stats: BatchStats) -> None:
        """
        Record statistics for a completed batch.
        
        Args:
            stats: Batch statistics
        """
        with self.lock:
            # Ensure the stats are complete
            if stats.end_time is None:
                stats.complete()
            
            # Add to history
            self.batch_stats.append(stats)
            
            # Update totals
            self.total_items_processed += stats.items_processed
            self.total_success_count += stats.success_count
            self.total_error_count += stats.error_count
            self.total_duration += stats.duration
            
            # Keep history bounded (keep last 100 batches)
            if len(self.batch_stats) > 100:
                self.batch_stats = self.batch_stats[-100:]
    
    def get_overall_stats(self) -> Dict[str, Any]:
        """
        Get overall processing statistics.
        
        Returns:
            Dictionary with statistics
        """
        with self.lock:
            if not self.batch_stats:
                return {
                    "total_items_processed": 0,
                    "total_success_count": 0,
                    "total_error_count": 0,
                    "total_duration": 0,
                    "avg_batch_size": self.initial_batch_size,
                    "current_batch_size": self.current_batch_size,
                    "items_per_second": 0,
                    "error_rate": 0,
                    "batch_count": 0
                }
            
            return {
                "total_items_processed": self.total_items_processed,
                "total_success_count": self.total_success_count,
                "total_error_count": self.total_error_count,
                "total_duration": self.total_duration,
                "avg_batch_size": sum(stat.batch_size for stat in self.batch_stats) / len(self.batch_stats),
                "current_batch_size": self.current_batch_size,
                "items_per_second": self.total_items_processed / max(0.001, self.total_duration),
                "error_rate": self.total_error_count / max(1, self.total_items_processed),
                "batch_count": len(self.batch_stats)
            }


class AdaptiveBatchProcessor(Generic[T, U]):
    """
    Adaptive batch processor for efficiently processing large datasets.
    
    This class manages batching of data items with:
    - Adaptive batch sizing based on processing metrics
    - Parallel processing capabilities
    - Automatic resource management
    - Comprehensive monitoring
    """
    
    def __init__(
        self,
        process_func: Callable[[List[T]], List[U]],
        initial_batch_size: int = 100,
        min_batch_size: int = 10,
        max_batch_size: int = 1000,
        max_workers: Optional[int] = None,
        strategy: Optional[BatchProcessingStrategy] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the adaptive batch processor.
        
        Args:
            process_func: Function to process a batch of items
            initial_batch_size: Initial size of each batch
            min_batch_size: Minimum batch size
            max_batch_size: Maximum batch size
            max_workers: Maximum worker threads/processes
            strategy: Batch processing strategy
            logger: Logger instance
        """
        self.process_func = process_func
        self.logger = logger or logging.getLogger(__name__)
        
        # Determine optimal number of workers if not specified
        if max_workers is None:
            # Use CPU count as default, but prevent excessive threads
            max_workers = min(os.cpu_count() or 4, 8)
        
        self.max_workers = max_workers
        
        # Create default strategy if not provided
        if strategy is None:
            strategy = BatchProcessingStrategy(
                initial_batch_size=initial_batch_size,
                min_batch_size=min_batch_size,
                max_batch_size=max_batch_size,
                logger=self.logger
            )
        
        self.strategy = strategy
        
        # Create thread pool for parallel processing
        self.executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)
        
        self.logger.info(
            f"Initialized AdaptiveBatchProcessor with batch size {initial_batch_size} "
            f"(min={min_batch_size}, max={max_batch_size}) and {max_workers} workers"
        )
    
    def process_batch(self, items: List[T]) -> Tuple[List[U], BatchStats]:
        """
        Process a single batch of items.
        
        Args:
            items: List of items to process
            
        Returns:
            Tuple of (processed results, batch statistics)
        """
        # Create batch statistics
        stats = BatchStats(batch_size=len(items))
        
        try:
            # Process the batch
            start_time = time.time()
            results = self.process_func(items)
            end_time = time.time()
            
            # Record processing time
            stats.processing_times.append(end_time - start_time)
            stats.items_processed = len(items)
            stats.success_count = len(results)
            stats.error_count = len(items) - len(results)
            
            # Complete stats
            stats.complete()
            
            return results, stats
        except Exception as e:
            self.logger.error(f"Error processing batch: {str(e)}")
            
            # Record error in stats
            stats.items_processed = len(items)
            stats.success_count = 0
            stats.error_count = len(items)
            stats.complete()
            
            return [], stats
    
    def process_items(self, items: List[T]) -> List[U]:
        """
        Process a list of items using adaptive batching.
        
        Args:
            items: List of items to process
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        all_results = []
        
        # Process in batches
        for i in range(0, len(items), self.strategy.current_batch_size):
            # Get batch
            batch = items[i:i + self.strategy.current_batch_size]
            
            # Process batch
            results, stats = self.process_batch(batch)
            
            # Record statistics
            self.strategy.record_batch_stats(stats)
            
            # Extend results
            all_results.extend(results)
            
            # Update batch size for next batch
            self.strategy.get_next_batch_size()
        
        return all_results
    
    def process_items_parallel(self, items: List[T]) -> List[U]:
        """
        Process items in parallel using multiple workers.
        
        Args:
            items: List of items to process
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        all_results = []
        all_stats = []
        batch_size = self.strategy.current_batch_size
        
        # Create batches
        batches = []
        for i in range(0, len(items), batch_size):
            batches.append(items[i:i + batch_size])
        
        # Process batches in parallel
        with self.executor as executor:
            # Submit all batches
            futures = [executor.submit(self.process_batch, batch) for batch in batches]
            
            # Collect results
            for future in concurrent.futures.as_completed(futures):
                try:
                    results, stats = future.result()
                    all_results.extend(results)
                    all_stats.append(stats)
                except Exception as e:
                    self.logger.error(f"Error in parallel batch processing: {str(e)}")
        
        # Record all batch statistics
        for stats in all_stats:
            self.strategy.record_batch_stats(stats)
        
        # Update batch size for next run
        self.strategy.get_next_batch_size()
        
        return all_results
    
    async def process_items_async(self, items: List[T]) -> List[U]:
        """
        Process items asynchronously.
        
        Args:
            items: List of items to process
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        all_results = []
        
        # Process in batches
        batch_size = self.strategy.current_batch_size
        
        for i in range(0, len(items), batch_size):
            # Get batch
            batch = items[i:i + batch_size]
            
            # Process batch in a thread to avoid blocking the event loop
            results, stats = await asyncio.to_thread(self.process_batch, batch)
            
            # Record statistics
            self.strategy.record_batch_stats(stats)
            
            # Extend results
            all_results.extend(results)
            
            # Update batch size for next batch
            self.strategy.get_next_batch_size()
        
        return all_results
    
    async def process_items_parallel_async(self, items: List[T]) -> List[U]:
        """
        Process items in parallel asynchronously.
        
        Args:
            items: List of items to process
            
        Returns:
            List of processed results
        """
        if not items:
            return []
        
        all_results = []
        all_stats = []
        batch_size = self.strategy.current_batch_size
        
        # Create batches
        batches = []
        for i in range(0, len(items), batch_size):
            batches.append(items[i:i + batch_size])
        
        # Process batches in parallel
        tasks = []
        for batch in batches:
            # Schedule each batch to run in a separate thread
            task = asyncio.create_task(asyncio.to_thread(self.process_batch, batch))
            tasks.append(task)
        
        # Wait for all tasks to complete
        for task in asyncio.as_completed(tasks):
            try:
                results, stats = await task
                all_results.extend(results)
                all_stats.append(stats)
            except Exception as e:
                self.logger.error(f"Error in parallel async batch processing: {str(e)}")
        
        # Record all batch statistics
        for stats in all_stats:
            self.strategy.record_batch_stats(stats)
        
        # Update batch size for next run
        self.strategy.get_next_batch_size()
        
        return all_results
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get processing statistics.
        
        Returns:
            Dictionary with processing statistics
        """
        stats = self.strategy.get_overall_stats()
        stats["max_workers"] = self.max_workers
        return stats
    
    def close(self) -> None:
        """Clean up resources."""
        self.executor.shutdown(wait=False)


class ChemicalDataBatchStrategy(BatchProcessingStrategy):
    """
    Specialized batch processing strategy for chemical data.
    
    This strategy accounts for the unique characteristics of chemical data processing:
    - Memory usage varies widely based on molecule complexity
    - Processing time increases non-linearly with molecule size
    - Error rates may vary based on data source
    """
    
    def __init__(
        self,
        data_source: str = "generic",
        **kwargs
    ):
        """
        Initialize the chemical data batch strategy.
        
        Args:
            data_source: Source of the chemical data ("chembl", "pubchem", "generic")
            **kwargs: Additional arguments passed to BatchProcessingStrategy
        """
        super().__init__(**kwargs)
        self.data_source = data_source.lower()
        
        # Adjust initial parameters based on data source
        if self.data_source == "chembl":
            # ChEMBL data tends to have more complex molecules
            self.initial_batch_size = min(50, self.initial_batch_size)
            self.target_duration = max(10.0, self.target_duration)
        elif self.data_source == "pubchem":
            # PubChem can handle larger batches
            self.initial_batch_size = min(200, self.initial_batch_size)
            self.target_duration = max(8.0, self.target_duration)
        
        # Reset current batch size to match initial
        self.current_batch_size = self.initial_batch_size
        
        self.logger.info(
            f"Initialized ChemicalDataBatchStrategy for {data_source} with "
            f"batch size {self.initial_batch_size}, target duration {self.target_duration}s"
        )
    
    def get_next_batch_size(self) -> int:
        """
        Get the next batch size with chemical data-specific adjustments.
        
        Overrides the base implementation with specialized logic.
        
        Returns:
            Next recommended batch size
        """
        with self.lock:
            # Get base adjustment
            base_size = super().get_next_batch_size()
            
            # If we don't have enough statistics, use the base size
            if len(self.batch_stats) < 5:
                return base_size
            
            # Get the most recent batch stats
            recent_stats = self.batch_stats[-5:]
            
            # For chemical data, we want to be more conservative with batch sizing
            # to prevent excessive memory usage and timeouts
            
            # Calculate maximum memory usage per item
            max_memory_per_item = 0
            for stat in recent_stats:
                if stat.memory_usage_delta and stat.items_processed:
                    memory_per_item = stat.memory_usage_delta / stat.items_processed
                    max_memory_per_item = max(max_memory_per_item, memory_per_item)
            
            # If memory usage is high, be more conservative
            if max_memory_per_item > 1024 * 1024:  # More than 1MB per item
                self.logger.debug(
                    f"High memory usage detected: {max_memory_per_item/1024/1024:.2f}MB per item, "
                    f"reducing batch size by 10%"
                )
                base_size = int(base_size * 0.9)
            
            # Check for stability in processing times
            processing_times = [stat.avg_processing_time for stat in recent_stats]
            max_time = max(processing_times)
            min_time = min(processing_times)
            
            # If processing times are unstable (high variance), be more conservative
            if max_time > min_time * 2:
                self.logger.debug(
                    f"Unstable processing times detected: {min_time:.2f}s to {max_time:.2f}s, "
                    f"reducing batch size by 5%"
                )
                base_size = int(base_size * 0.95)
            
            # Constrain to min/max limits
            new_batch_size = max(self.min_batch_size, min(self.max_batch_size, base_size))
            self.current_batch_size = new_batch_size
            return new_batch_size
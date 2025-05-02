#!/usr/bin/env python3
"""
Worker Pool Implementation for Parallel PubChem Data Processing

This module provides a thread-based worker pool implementation for parallel processing
of PubChem data import tasks. It includes:

- WorkerPool: Manages a configurable number of worker threads
- Worker: Processes chunks from the work queue
- ResultCollector: Aggregates results from all workers

The implementation is designed to be thread-safe and efficiently coordinate shared
resources such as rate limiters and circuit breakers.
"""

import os
import time
import logging
import threading
import queue
import random
from typing import List, Dict, Any, Tuple, Optional, Callable, Union, Set
from datetime import datetime
from dataclasses import dataclass, field
import traceback

# Set up logging
logger = logging.getLogger(__name__)

# Default configuration
DEFAULT_NUM_WORKERS = 4
DEFAULT_MAX_QUEUE_SIZE = 100
DEFAULT_RESULT_QUEUE_SIZE = 1000
DEFAULT_WORKER_TIMEOUT = 60  # seconds
DEFAULT_HEALTH_CHECK_INTERVAL = 5  # seconds

@dataclass
class WorkItem:
    """Represents a unit of work to be processed by a worker."""
    id: str
    data: Any
    priority: int = 0
    created_at: float = field(default_factory=time.time)
    attempts: int = 0
    max_attempts: int = 3
    
    def increment_attempt(self) -> None:
        """Increment the attempt counter."""
        self.attempts += 1
    
    def can_retry(self) -> bool:
        """Check if the work item can be retried."""
        return self.attempts < self.max_attempts


@dataclass
class WorkResult:
    """Represents the result of processing a work item."""
    work_item_id: str
    success: bool
    result: Any = None
    error: Optional[Exception] = None
    processing_time: float = 0.0
    worker_id: Optional[int] = None
    completed_at: float = field(default_factory=time.time)


class Worker:
    """
    Worker thread that processes items from a work queue.
    
    Features:
    - Processes chunks from the work queue
    - Reports results to the result collector
    - Implements individual backoff when encountering errors
    - Supports graceful shutdown
    """
    
    def __init__(
        self,
        worker_id: int,
        work_queue: queue.Queue,
        result_queue: queue.Queue,
        processor_func: Callable[[Any], Any],
        rate_limiter: Optional[Any] = None,
        circuit_breaker: Optional[Any] = None,
        backoff_factor: float = 2.0,
        max_backoff: float = 60.0,
        jitter: bool = True
    ):
        """
        Initialize the worker.
        
        Args:
            worker_id: Unique identifier for this worker
            work_queue: Queue to get work items from
            result_queue: Queue to put results into
            processor_func: Function to process work items
            rate_limiter: Optional rate limiter to use
            circuit_breaker: Optional circuit breaker to use
            backoff_factor: Factor to increase delay for each retry
            max_backoff: Maximum delay between retries in seconds
            jitter: Whether to add random jitter to delay
        """
        self.worker_id = worker_id
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.processor_func = processor_func
        self.rate_limiter = rate_limiter
        self.circuit_breaker = circuit_breaker
        self.backoff_factor = backoff_factor
        self.max_backoff = max_backoff
        self.jitter = jitter
        
        self.running = False
        self.thread = None
        self.last_activity = time.time()
        self.items_processed = 0
        self.errors = 0
        self.consecutive_errors = 0
        self.total_processing_time = 0.0
        self.status = "idle"
        self.current_work_item = None
        self.shutdown_event = threading.Event()
    
    def start(self) -> None:
        """Start the worker thread."""
        if self.running:
            return
        
        self.running = True
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.thread.start()
        logger.info(f"Worker {self.worker_id} started")
    
    def stop(self, timeout: float = 5.0) -> None:
        """
        Stop the worker thread gracefully.
        
        Args:
            timeout: Maximum time to wait for the worker to stop
        """
        if not self.running:
            return
        
        self.running = False
        self.shutdown_event.set()
        
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout)
            if self.thread.is_alive():
                logger.warning(f"Worker {self.worker_id} did not stop within timeout")
            else:
                logger.info(f"Worker {self.worker_id} stopped gracefully")
    
    def _run(self) -> None:
        """Main worker loop."""
        while self.running:
            try:
                # Check if we should shut down
                if self.shutdown_event.is_set():
                    break
                
                # Check circuit breaker
                if self.circuit_breaker and not self._check_circuit_breaker():
                    # Circuit is open, wait before trying again
                    self.status = "waiting_circuit"
                    time.sleep(1.0)
                    continue
                
                # Try to get a work item with a timeout to allow for shutdown checks
                try:
                    self.status = "waiting_work"
                    work_item = self.work_queue.get(timeout=0.5)
                except queue.Empty:
                    continue
                
                # Process the work item
                self.status = "processing"
                self.current_work_item = work_item
                self.last_activity = time.time()
                
                # Apply rate limiting if configured
                if self.rate_limiter:
                    self.rate_limiter.wait()
                
                # Process the work item
                result = self._process_work_item(work_item)
                
                # Put the result in the result queue
                self.result_queue.put(result)
                
                # Mark the work item as done
                self.work_queue.task_done()
                self.current_work_item = None
                self.status = "idle"
                
            except Exception as e:
                logger.error(f"Worker {self.worker_id} encountered an error: {str(e)}")
                logger.debug(traceback.format_exc())
                self.errors += 1
                self.consecutive_errors += 1
                time.sleep(0.1)  # Prevent tight error loops
    
    def _process_work_item(self, work_item: WorkItem) -> WorkResult:
        """
        Process a work item and return the result.
        
        Args:
            work_item: Work item to process
            
        Returns:
            WorkResult object with processing result
        """
        start_time = time.time()
        success = False
        result = None
        error = None
        
        try:
            # Increment attempt counter
            work_item.increment_attempt()
            
            # Process the work item
            if self.circuit_breaker:
                # Use circuit breaker if available
                @self.circuit_breaker
                def process_with_circuit_breaker():
                    return self.processor_func(work_item.data)
                
                result = process_with_circuit_breaker()
            else:
                # Process directly
                result = self.processor_func(work_item.data)
            
            # Mark as successful
            success = True
            self.consecutive_errors = 0
            self.items_processed += 1
            
        except Exception as e:
            error = e
            logger.warning(f"Worker {self.worker_id} failed to process work item {work_item.id}: {str(e)}")
            self.errors += 1
            self.consecutive_errors += 1
        
        # Calculate processing time
        processing_time = time.time() - start_time
        self.total_processing_time += processing_time
        
        # Create and return result
        return WorkResult(
            work_item_id=work_item.id,
            success=success,
            result=result,
            error=error,
            processing_time=processing_time,
            worker_id=self.worker_id
        )
    
    def _check_circuit_breaker(self) -> bool:
        """
        Check if the circuit breaker is open.
        
        Returns:
            True if the circuit is closed or half-open (requests allowed),
            False if the circuit is open (requests blocked)
        """
        if not self.circuit_breaker:
            return True
            
        try:
            circuit_stats = self.circuit_breaker.get_stats()
            if circuit_stats["state"] == "open":
                recovery_time = circuit_stats.get("recovery_time_remaining", 0)
                if recovery_time > 0:
                    logger.warning(
                        f"Worker {self.worker_id}: Circuit breaker is OPEN. "
                        f"Waiting {recovery_time:.1f}s before retrying."
                    )
                    return False
        except Exception as e:
            logger.warning(f"Worker {self.worker_id}: Error checking circuit breaker state: {str(e)}")
        
        return True
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get worker statistics.
        
        Returns:
            Dictionary with worker statistics
        """
        return {
            "worker_id": self.worker_id,
            "status": self.status,
            "items_processed": self.items_processed,
            "errors": self.errors,
            "consecutive_errors": self.consecutive_errors,
            "total_processing_time": self.total_processing_time,
            "average_processing_time": (
                self.total_processing_time / self.items_processed 
                if self.items_processed > 0 else 0
            ),
            "last_activity": self.last_activity,
            "current_work_item": self.current_work_item.id if self.current_work_item else None,
            "running": self.running
        }


class ResultCollector:
    """
    Collects and aggregates results from workers.
    
    Features:
    - Aggregates results from all workers
    - Manages database transactions
    - Updates checkpoint data
    - Provides progress statistics
    """
    
    def __init__(
        self,
        result_queue: queue.Queue,
        result_handler: Optional[Callable[[WorkResult], None]] = None,
        checkpoint_handler: Optional[Callable[[List[WorkResult]], None]] = None,
        checkpoint_interval: int = 10
    ):
        """
        Initialize the result collector.
        
        Args:
            result_queue: Queue to get results from
            result_handler: Optional function to handle individual results
            checkpoint_handler: Optional function to handle checkpoints
            checkpoint_interval: Number of results to collect before checkpoint
        """
        self.result_queue = result_queue
        self.result_handler = result_handler
        self.checkpoint_handler = checkpoint_handler
        self.checkpoint_interval = checkpoint_interval
        
        self.results: List[WorkResult] = []
        self.successful_results: List[WorkResult] = []
        self.failed_results: List[WorkResult] = []
        self.running = False
        self.thread = None
        self.last_checkpoint_time = time.time()
        self.processed_since_checkpoint = 0
        self.shutdown_event = threading.Event()
    
    def start(self) -> None:
        """Start the result collector thread."""
        if self.running:
            return
        
        self.running = True
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.thread.start()
        logger.info("Result collector started")
    
    def stop(self, timeout: float = 5.0) -> None:
        """
        Stop the result collector thread gracefully.
        
        Args:
            timeout: Maximum time to wait for the collector to stop
        """
        if not self.running:
            return
        
        self.running = False
        self.shutdown_event.set()
        
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout)
            if self.thread.is_alive():
                logger.warning("Result collector did not stop within timeout")
            else:
                logger.info("Result collector stopped gracefully")
    
    def _run(self) -> None:
        """Main result collector loop."""
        while self.running:
            try:
                # Check if we should shut down
                if self.shutdown_event.is_set():
                    break
                
                # Try to get a result with a timeout to allow for shutdown checks
                try:
                    result = self.result_queue.get(timeout=0.5)
                except queue.Empty:
                    continue
                
                # Process the result
                self._process_result(result)
                
                # Mark the result as done
                self.result_queue.task_done()
                
                # Check if we should create a checkpoint
                self.processed_since_checkpoint += 1
                if (self.processed_since_checkpoint >= self.checkpoint_interval or
                    time.time() - self.last_checkpoint_time >= 60):  # At least every minute
                    self._create_checkpoint()
                
            except Exception as e:
                logger.error(f"Result collector encountered an error: {str(e)}")
                logger.debug(traceback.format_exc())
                time.sleep(0.1)  # Prevent tight error loops
    
    def _process_result(self, result: WorkResult) -> None:
        """
        Process a result.
        
        Args:
            result: Result to process
        """
        # Add to appropriate lists
        self.results.append(result)
        if result.success:
            self.successful_results.append(result)
        else:
            self.failed_results.append(result)
        
        # Call result handler if provided
        if self.result_handler:
            try:
                self.result_handler(result)
            except Exception as e:
                logger.error(f"Error in result handler: {str(e)}")
    
    def _create_checkpoint(self) -> None:
        """Create a checkpoint with current results."""
        if not self.checkpoint_handler:
            return
        
        try:
            # Call checkpoint handler with results since last checkpoint
            results_to_checkpoint = self.results[-self.processed_since_checkpoint:]
            self.checkpoint_handler(results_to_checkpoint)
            
            # Update checkpoint tracking
            self.last_checkpoint_time = time.time()
            self.processed_since_checkpoint = 0
            
            logger.info(f"Created checkpoint with {len(results_to_checkpoint)} results")
            
        except Exception as e:
            logger.error(f"Error creating checkpoint: {str(e)}")
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get result collector statistics.
        
        Returns:
            Dictionary with result collector statistics
        """
        total_results = len(self.results)
        successful_results = len(self.successful_results)
        failed_results = len(self.failed_results)
        
        return {
            "total_results": total_results,
            "successful_results": successful_results,
            "failed_results": failed_results,
            "success_rate": successful_results / total_results if total_results > 0 else 0,
            "last_checkpoint_time": self.last_checkpoint_time,
            "processed_since_checkpoint": self.processed_since_checkpoint,
            "running": self.running
        }
class WorkerPool:
    """
    Manages a pool of worker threads for parallel processing.
    
    Features:
    - Manages a configurable number of worker threads
    - Implements thread-safe work queue
    - Coordinates shared resources (rate limiter, circuit breaker)
    - Provides status monitoring and worker health checks
    """
    
    def __init__(
        self,
        num_workers: int = DEFAULT_NUM_WORKERS,
        max_queue_size: int = DEFAULT_MAX_QUEUE_SIZE,
        processor_func: Optional[Callable[[Any], Any]] = None,
        rate_limiter: Optional[Any] = None,
        circuit_breaker: Optional[Any] = None,
        result_handler: Optional[Callable[[WorkResult], None]] = None,
        checkpoint_handler: Optional[Callable[[List[WorkResult]], None]] = None,
        checkpoint_interval: int = 10,
        worker_timeout: float = DEFAULT_WORKER_TIMEOUT,
        health_check_interval: float = DEFAULT_HEALTH_CHECK_INTERVAL
    ):
        """
        Initialize the worker pool.
        
        Args:
            num_workers: Number of worker threads to create
            max_queue_size: Maximum size of the work queue
            processor_func: Function to process work items
            rate_limiter: Optional rate limiter to use
            circuit_breaker: Optional circuit breaker to use
            result_handler: Optional function to handle individual results
            checkpoint_handler: Optional function to handle checkpoints
            checkpoint_interval: Number of results to collect before checkpoint
            worker_timeout: Maximum time a worker can be inactive before considered stuck
            health_check_interval: Interval between worker health checks
        """
        self.num_workers = num_workers
        self.max_queue_size = max_queue_size
        self.processor_func = processor_func
        self.rate_limiter = rate_limiter
        self.circuit_breaker = circuit_breaker
        self.worker_timeout = worker_timeout
        self.health_check_interval = health_check_interval
        
        # Create queues
        self.work_queue = queue.Queue(maxsize=max_queue_size)
        self.result_queue = queue.Queue(maxsize=DEFAULT_RESULT_QUEUE_SIZE)
        
        # Create result collector
        self.result_collector = ResultCollector(
            result_queue=self.result_queue,
            result_handler=result_handler,
            checkpoint_handler=checkpoint_handler,
            checkpoint_interval=checkpoint_interval
        )
        
        # Create workers
        self.workers: List[Worker] = []
        for i in range(num_workers):
            worker = Worker(
                worker_id=i,
                work_queue=self.work_queue,
                result_queue=self.result_queue,
                processor_func=processor_func,
                rate_limiter=rate_limiter,
                circuit_breaker=circuit_breaker
            )
            self.workers.append(worker)
        
        # Health check thread
        self.health_check_thread = None
        self.health_check_running = False
        self.health_check_shutdown_event = threading.Event()
        
        # Status
        self.running = False
        self.start_time = None
        self.total_work_items = 0
        self.completed_work_items = 0
    
    def start(self) -> None:
        """Start the worker pool and all worker threads."""
        if self.running:
            return
        
        if not self.processor_func:
            raise ValueError("Processor function must be set before starting the worker pool")
        
        self.running = True
        self.start_time = time.time()
        
        # Start result collector
        self.result_collector.start()
        
        # Start workers
        for worker in self.workers:
            worker.start()
        
        # Start health check thread
        self._start_health_check()
        
        logger.info(f"Worker pool started with {self.num_workers} workers")
    
    def stop(self, timeout: float = 5.0) -> None:
        """
        Stop the worker pool and all worker threads gracefully.
        
        Args:
            timeout: Maximum time to wait for workers to stop
        """
        if not self.running:
            return
        
        self.running = False
        
        # Stop health check thread
        self._stop_health_check()
        
        # Stop workers
        for worker in self.workers:
            worker.stop(timeout)
        
        # Stop result collector
        self.result_collector.stop(timeout)
        
        logger.info("Worker pool stopped")
    
    def add_work(self, work_data: Any, work_id: Optional[str] = None, priority: int = 0) -> str:
        """
        Add a work item to the queue.
        
        Args:
            work_data: Data to process
            work_id: Optional identifier for the work item (generated if not provided)
            priority: Priority of the work item (higher values have higher priority)
            
        Returns:
            Work item ID
        """
        if not self.running:
            raise RuntimeError("Worker pool is not running")
        
        # Generate work ID if not provided
        if work_id is None:
            work_id = f"work_{self.total_work_items}_{int(time.time())}"
        
        # Create work item
        work_item = WorkItem(id=work_id, data=work_data, priority=priority)
        
        # Add to queue
        self.work_queue.put(work_item)
        self.total_work_items += 1
        
        return work_id
    
    def add_work_batch(self, work_data_list: List[Any], id_prefix: str = "batch") -> List[str]:
        """
        Add a batch of work items to the queue.
        
        Args:
            work_data_list: List of data items to process
            id_prefix: Prefix for generated work IDs
            
        Returns:
            List of work item IDs
        """
        work_ids = []
        batch_id = f"{id_prefix}_{int(time.time())}"
        
        for i, work_data in enumerate(work_data_list):
            work_id = f"{batch_id}_{i}"
            self.add_work(work_data, work_id)
            work_ids.append(work_id)
        
        return work_ids
    
    def wait_for_completion(self, timeout: Optional[float] = None) -> bool:
        """
        Wait for all work items to be processed.
        
        Args:
            timeout: Maximum time to wait in seconds (None for no timeout)
            
        Returns:
            True if all work completed, False if timeout occurred
        """
        if not self.running:
            return True
        
        start_time = time.time()
        
        try:
            # Wait for work queue to be empty
            while not self.work_queue.empty():
                if timeout and time.time() - start_time > timeout:
                    return False
                time.sleep(0.1)
            
            # Wait for all work to be done
            self.work_queue.join()
            
            # Wait for result queue to be empty
            while not self.result_queue.empty():
                if timeout and time.time() - start_time > timeout:
                    return False
                time.sleep(0.1)
            
            # Wait for result queue to be done
            self.result_queue.join()
            
            return True
            
        except KeyboardInterrupt:
            logger.warning("Wait for completion interrupted")
            return False
    
    def _start_health_check(self) -> None:
        """Start the health check thread."""
        if self.health_check_running:
            return
        
        self.health_check_running = True
        self.health_check_thread = threading.Thread(target=self._health_check_loop, daemon=True)
        self.health_check_thread.start()
    
    def _stop_health_check(self) -> None:
        """Stop the health check thread."""
        if not self.health_check_running:
            return
        
        self.health_check_running = False
        self.health_check_shutdown_event.set()
        
        if self.health_check_thread and self.health_check_thread.is_alive():
            self.health_check_thread.join(5.0)
    
    def _health_check_loop(self) -> None:
        """Main health check loop."""
        while self.health_check_running:
            try:
                # Check if we should shut down
                if self.health_check_shutdown_event.is_set():
                    break
                
                # Check worker health
                self._check_worker_health()
                
                # Wait for next check
                time.sleep(self.health_check_interval)
                
            except Exception as e:
                logger.error(f"Health check encountered an error: {str(e)}")
                time.sleep(1.0)  # Prevent tight error loops
    
    def _check_worker_health(self) -> None:
        """Check the health of all workers and restart any that are stuck."""
        current_time = time.time()
        
        for i, worker in enumerate(self.workers):
            # Skip workers that aren't running
            if not worker.running:
                continue
            
            # Check if worker is stuck
            if (worker.status == "processing" and 
                current_time - worker.last_activity > self.worker_timeout):
                logger.warning(
                    f"Worker {worker.worker_id} appears to be stuck "
                    f"(no activity for {current_time - worker.last_activity:.1f}s). "
                    f"Restarting..."
                )
                
                # Stop and restart the worker
                worker.stop(1.0)  # Short timeout for stuck workers
                
                # Create a new worker to replace the stuck one
                new_worker = Worker(
                    worker_id=worker.worker_id,
                    work_queue=self.work_queue,
                    result_queue=self.result_queue,
                    processor_func=self.processor_func,
                    rate_limiter=self.rate_limiter,
                    circuit_breaker=self.circuit_breaker
                )
                new_worker.start()
                
                # Replace the worker in our list
                self.workers[i] = new_worker
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get worker pool statistics.
        
        Returns:
            Dictionary with worker pool statistics
        """
        # Collect worker stats
        worker_stats = [worker.get_stats() for worker in self.workers]
        
        # Calculate aggregate stats
        total_processed = sum(stats["items_processed"] for stats in worker_stats)
        total_errors = sum(stats["errors"] for stats in worker_stats)
        total_processing_time = sum(stats["total_processing_time"] for stats in worker_stats)
        
        # Get result collector stats
        result_stats = self.result_collector.get_stats()
        
        # Calculate queue stats
        work_queue_size = self.work_queue.qsize()
        result_queue_size = self.result_queue.qsize()
        
        # Calculate overall stats
        elapsed_time = time.time() - (self.start_time or time.time())
        throughput = total_processed / elapsed_time if elapsed_time > 0 else 0
        
        return {
            "running": self.running,
            "num_workers": self.num_workers,
            "active_workers": sum(1 for w in self.workers if w.running),
            "total_work_items": self.total_work_items,
            "completed_work_items": total_processed,
            "pending_work_items": work_queue_size,
            "pending_results": result_queue_size,
            "total_errors": total_errors,
            "error_rate": total_errors / total_processed if total_processed > 0 else 0,
            "total_processing_time": total_processing_time,
            "elapsed_time": elapsed_time,
            "throughput": throughput,
            "worker_stats": worker_stats,
            "result_stats": result_stats
        }
            "success_rate": successful_results / total_results if total_results > 0 else 0,
            "last_checkpoint_time": self.last_checkpoint_time,
            "processed_since_checkpoint": self.processed_since_checkpoint,
            "running": self.running
        }
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
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
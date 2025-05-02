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
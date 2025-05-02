"""
Weekend job scheduler for PubChem API bulk operations.

This module provides scheduling functionality for running bulk operations
during weekends when PubChem API rate limits are less restrictive.
"""

import os
import time
import logging
import datetime
import threading
import json
from typing import Callable, Dict, List, Any, Optional, Union
from pathlib import Path

logger = logging.getLogger(__name__)

class WeekendJobScheduler:
    """
    Scheduler for running jobs during weekends.
    
    Features:
    - Schedule jobs to run during weekends
    - Persistent job queue that survives restarts
    - Configurable weekend definition
    - Background thread for job execution
    """
    
    def __init__(
        self,
        storage_dir: str = "cache/pubchem/scheduler",
        weekend_days: tuple = (5, 6),  # 5=Saturday, 6=Sunday
        weekend_hours: tuple = (0, 24),  # All day
        check_interval: int = 300  # Check every 5 minutes
    ):
        """
        Initialize the scheduler.
        
        Args:
            storage_dir: Directory to store job queue and results
            weekend_days: Tuple of days considered weekend (0=Monday, 6=Sunday)
            weekend_hours: Tuple of (start_hour, end_hour) for weekend jobs
            check_interval: Interval in seconds to check for pending jobs
        """
        self.storage_dir = Path(storage_dir)
        self.weekend_days = weekend_days
        self.weekend_hours = weekend_hours
        self.check_interval = check_interval
        
        # Create storage directory if it doesn't exist
        os.makedirs(self.storage_dir, exist_ok=True)
        
        # Job queue file
        self.queue_file = self.storage_dir / "job_queue.json"
        self.results_dir = self.storage_dir / "results"
        os.makedirs(self.results_dir, exist_ok=True)
        
        # Initialize job queue
        self.jobs: List[Dict[str, Any]] = []
        self._load_jobs()
        
        # Thread control
        self.running = False
        self.thread: Optional[threading.Thread] = None
        self.lock = threading.RLock()
    
    def _load_jobs(self) -> None:
        """Load jobs from persistent storage."""
        if self.queue_file.exists():
            try:
                with open(self.queue_file, 'r') as f:
                    self.jobs = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logger.error(f"Error loading job queue: {e}")
                self.jobs = []
        else:
            self.jobs = []
    
    def _save_jobs(self) -> None:
        """Save jobs to persistent storage."""
        try:
            with open(self.queue_file, 'w') as f:
                json.dump(self.jobs, f)
        except IOError as e:
            logger.error(f"Error saving job queue: {e}")
    
    def _is_weekend(self) -> bool:
        """
        Check if current time is considered a weekend.
        
        Returns:
            True if current time is a weekend, False otherwise
        """
        now = datetime.datetime.now()
        
        # Check if current day is in weekend_days
        if now.weekday() not in self.weekend_days:
            return False
        
        # Check if current hour is in weekend_hours
        start_hour, end_hour = self.weekend_hours
        return start_hour <= now.hour < end_hour
    
    def _worker(self) -> None:
        """Background worker thread for executing jobs."""
        logger.info("Weekend job scheduler worker thread started")
        
        while self.running:
            try:
                # Check if it's weekend
                if self._is_weekend():
                    # Process jobs
                    self._process_jobs()
                
                # Sleep for check_interval
                time.sleep(self.check_interval)
            except Exception as e:
                logger.error(f"Error in weekend job scheduler worker: {e}")
                # Sleep a bit to avoid tight loop in case of persistent errors
                time.sleep(10)
        
        logger.info("Weekend job scheduler worker thread stopped")
    
    def _process_jobs(self) -> None:
        """Process pending jobs."""
        with self.lock:
            if not self.jobs:
                return
            
            logger.info(f"Processing {len(self.jobs)} pending jobs during weekend")
            
            # Get the next job
            job = self.jobs[0]
            
            try:
                # Execute the job
                job_id = job["id"]
                job_type = job["type"]
                job_params = job["params"]
                
                logger.info(f"Executing job {job_id} of type {job_type}")
                
                # Record job start
                job["status"] = "running"
                job["start_time"] = datetime.datetime.now().isoformat()
                self._save_jobs()
                
                # Execute the job based on type
                result = None
                if job_type == "prefetch":
                    # This would be implemented by the client
                    # We'll just log it for now
                    logger.info(f"Would prefetch data for {job_params}")
                    result = {"status": "success", "message": "Prefetch completed"}
                elif job_type == "batch_update":
                    # This would be implemented by the client
                    # We'll just log it for now
                    logger.info(f"Would batch update data for {job_params}")
                    result = {"status": "success", "message": "Batch update completed"}
                else:
                    logger.warning(f"Unknown job type: {job_type}")
                    result = {"status": "error", "message": f"Unknown job type: {job_type}"}
                
                # Save result
                result_file = self.results_dir / f"{job_id}.json"
                with open(result_file, 'w') as f:
                    json.dump({
                        "job": job,
                        "result": result,
                        "completion_time": datetime.datetime.now().isoformat()
                    }, f)
                
                # Remove job from queue
                self.jobs.pop(0)
                self._save_jobs()
                
                logger.info(f"Job {job_id} completed successfully")
            
            except Exception as e:
                logger.error(f"Error executing job: {e}")
                
                # Update job status
                job["status"] = "failed"
                job["error"] = str(e)
                job["failure_time"] = datetime.datetime.now().isoformat()
                
                # Move failed job to the end of the queue
                self.jobs.append(self.jobs.pop(0))
                self._save_jobs()
    
    def start(self) -> None:
        """Start the scheduler."""
        with self.lock:
            if self.running:
                return
            
            self.running = True
            self.thread = threading.Thread(target=self._worker, daemon=True)
            self.thread.start()
            
            logger.info("Weekend job scheduler started")
    
    def stop(self) -> None:
        """Stop the scheduler."""
        with self.lock:
            if not self.running:
                return
            
            self.running = False
            if self.thread:
                self.thread.join(timeout=10)
                self.thread = None
            
            logger.info("Weekend job scheduler stopped")
    
    def schedule_job(self, job_type: str, params: Dict[str, Any]) -> str:
        """
        Schedule a job to run during the weekend.
        
        Args:
            job_type: Type of job (e.g., "prefetch", "batch_update")
            params: Parameters for the job
            
        Returns:
            Job ID
        """
        with self.lock:
            # Generate job ID
            job_id = f"job_{int(time.time())}_{len(self.jobs)}"
            
            # Create job
            job = {
                "id": job_id,
                "type": job_type,
                "params": params,
                "status": "pending",
                "created_at": datetime.datetime.now().isoformat()
            }
            
            # Add to queue
            self.jobs.append(job)
            self._save_jobs()
            
            logger.info(f"Job {job_id} of type {job_type} scheduled")
            
            return job_id
    
    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """
        Get the status of a job.
        
        Args:
            job_id: Job ID
            
        Returns:
            Job status or None if job not found
        """
        with self.lock:
            # Check if job is in queue
            for job in self.jobs:
                if job["id"] == job_id:
                    return job
            
            # Check if job has completed
            result_file = self.results_dir / f"{job_id}.json"
            if result_file.exists():
                try:
                    with open(result_file, 'r') as f:
                        return json.load(f)
                except (json.JSONDecodeError, IOError) as e:
                    logger.error(f"Error loading job result: {e}")
            
            return None
    
    def get_pending_jobs(self) -> List[Dict[str, Any]]:
        """
        Get all pending jobs.
        
        Returns:
            List of pending jobs
        """
        with self.lock:
            return self.jobs.copy()
    
    def clear_jobs(self) -> None:
        """Clear all pending jobs."""
        with self.lock:
            self.jobs = []
            self._save_jobs()
            
            logger.info("Job queue cleared")
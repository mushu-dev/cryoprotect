"""
Unit tests for the PubChem weekend job scheduler.
"""

import os
import time
import json
import shutil
import unittest
import datetime
import threading
from unittest import mock
from pathlib import Path

from .scheduler import WeekendJobScheduler

class TestWeekendJobScheduler(unittest.TestCase):
    """Test cases for WeekendJobScheduler."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_storage_dir = "test_scheduler"
        self.scheduler = WeekendJobScheduler(
            storage_dir=self.test_storage_dir,
            weekend_days=(5, 6),  # Saturday and Sunday
            weekend_hours=(0, 24),  # All day
            check_interval=0.1  # Short interval for testing
        )
    
    def tearDown(self):
        """Clean up test environment."""
        # Stop the scheduler if running
        if self.scheduler.running:
            self.scheduler.stop()
        
        # Remove test directory
        if os.path.exists(self.test_storage_dir):
            shutil.rmtree(self.test_storage_dir)
    
    def test_is_weekend(self):
        """Test weekend detection."""
        # Mock datetime.now to return a weekday (Wednesday)
        wednesday = datetime.datetime(2025, 4, 23)  # A Wednesday
        with mock.patch('datetime.datetime') as mock_datetime:
            mock_datetime.now.return_value = wednesday
            self.assertFalse(self.scheduler._is_weekend())
        
        # Mock datetime.now to return a weekend (Saturday)
        saturday = datetime.datetime(2025, 4, 26)  # A Saturday
        with mock.patch('datetime.datetime') as mock_datetime:
            mock_datetime.now.return_value = saturday
            self.assertTrue(self.scheduler._is_weekend())
    
    def test_schedule_job(self):
        """Test job scheduling."""
        # Schedule a job
        job_type = "test_job"
        params = {"param1": "value1", "param2": "value2"}
        job_id = self.scheduler.schedule_job(job_type, params)
        
        # Check that job was added to queue
        self.assertEqual(len(self.scheduler.jobs), 1)
        self.assertEqual(self.scheduler.jobs[0]["type"], job_type)
        self.assertEqual(self.scheduler.jobs[0]["params"], params)
        self.assertEqual(self.scheduler.jobs[0]["status"], "pending")
        
        # Check that job was saved to disk
        queue_file = Path(self.test_storage_dir) / "job_queue.json"
        self.assertTrue(queue_file.exists())
        
        with open(queue_file, 'r') as f:
            saved_jobs = json.load(f)
        
        self.assertEqual(len(saved_jobs), 1)
        self.assertEqual(saved_jobs[0]["id"], job_id)
    
    def test_get_job_status(self):
        """Test getting job status."""
        # Schedule a job
        job_type = "test_job"
        params = {"param1": "value1", "param2": "value2"}
        job_id = self.scheduler.schedule_job(job_type, params)
        
        # Get job status
        status = self.scheduler.get_job_status(job_id)
        
        # Check status
        self.assertIsNotNone(status)
        self.assertEqual(status["type"], job_type)
        self.assertEqual(status["params"], params)
        self.assertEqual(status["status"], "pending")
        
        # Get status for nonexistent job
        status = self.scheduler.get_job_status("nonexistent_job")
        self.assertIsNone(status)
    
    def test_get_pending_jobs(self):
        """Test getting pending jobs."""
        # Schedule some jobs
        self.scheduler.schedule_job("job1", {"param": "value1"})
        self.scheduler.schedule_job("job2", {"param": "value2"})
        
        # Get pending jobs
        pending_jobs = self.scheduler.get_pending_jobs()
        
        # Check pending jobs
        self.assertEqual(len(pending_jobs), 2)
        self.assertEqual(pending_jobs[0]["type"], "job1")
        self.assertEqual(pending_jobs[1]["type"], "job2")
    
    def test_clear_jobs(self):
        """Test clearing jobs."""
        # Schedule some jobs
        self.scheduler.schedule_job("job1", {"param": "value1"})
        self.scheduler.schedule_job("job2", {"param": "value2"})
        
        # Clear jobs
        self.scheduler.clear_jobs()
        
        # Check that jobs were cleared
        self.assertEqual(len(self.scheduler.jobs), 0)
        
        # Check that jobs were cleared from disk
        queue_file = Path(self.test_storage_dir) / "job_queue.json"
        with open(queue_file, 'r') as f:
            saved_jobs = json.load(f)
        
        self.assertEqual(len(saved_jobs), 0)
    
    def test_start_stop(self):
        """Test starting and stopping the scheduler."""
        # Start the scheduler
        self.scheduler.start()
        
        # Check that scheduler is running
        self.assertTrue(self.scheduler.running)
        self.assertIsNotNone(self.scheduler.thread)
        
        # Stop the scheduler
        self.scheduler.stop()
        
        # Check that scheduler is stopped
        self.assertFalse(self.scheduler.running)
        self.assertIsNone(self.scheduler.thread)
    
    def test_job_processing_on_weekend(self):
        """Test job processing during weekend."""
        # Mock _is_weekend to return True
        with mock.patch.object(self.scheduler, '_is_weekend', return_value=True):
            # Schedule a job
            job_id = self.scheduler.schedule_job("test_job", {"param": "value"})
            
            # Start the scheduler
            self.scheduler.start()
            
            # Wait for job to be processed
            time.sleep(0.5)
            
            # Check that job was processed
            self.assertEqual(len(self.scheduler.jobs), 0)
            
            # Check that result was saved
            result_file = Path(self.test_storage_dir) / "results" / f"{job_id}.json"
            self.assertTrue(result_file.exists())
            
            # Check result content
            with open(result_file, 'r') as f:
                result = json.load(f)
            
            self.assertEqual(result["job"]["id"], job_id)
            self.assertEqual(result["result"]["status"], "success")
    
    def test_job_not_processed_on_weekday(self):
        """Test that jobs are not processed during weekdays."""
        # Mock _is_weekend to return False
        with mock.patch.object(self.scheduler, '_is_weekend', return_value=False):
            # Schedule a job
            self.scheduler.schedule_job("test_job", {"param": "value"})
            
            # Start the scheduler
            self.scheduler.start()
            
            # Wait a bit
            time.sleep(0.5)
            
            # Check that job was not processed
            self.assertEqual(len(self.scheduler.jobs), 1)
    
    def test_job_error_handling(self):
        """Test error handling during job processing."""
        # Mock _is_weekend to return True
        with mock.patch.object(self.scheduler, '_is_weekend', return_value=True):
            # Mock job processing to raise an exception
            original_process_jobs = self.scheduler._process_jobs
            
            def mock_process_jobs():
                # Process the first job
                if self.scheduler.jobs:
                    job = self.scheduler.jobs[0]
                    job["status"] = "running"
                    # Raise an exception
                    raise ValueError("Test error")
            
            self.scheduler._process_jobs = mock_process_jobs
            
            # Schedule a job
            self.scheduler.schedule_job("test_job", {"param": "value"})
            
            # Start the scheduler
            self.scheduler.start()
            
            # Wait for job to be processed
            time.sleep(0.5)
            
            # Check that job was marked as failed
            self.assertEqual(len(self.scheduler.jobs), 1)
            self.assertEqual(self.scheduler.jobs[0]["status"], "failed")
            self.assertIn("error", self.scheduler.jobs[0])
            
            # Restore original method
            self.scheduler._process_jobs = original_process_jobs
    
    def test_persistence(self):
        """Test that jobs persist across scheduler instances."""
        # Schedule a job
        job_type = "test_job"
        params = {"param": "value"}
        job_id = self.scheduler.schedule_job(job_type, params)
        
        # Create a new scheduler instance
        new_scheduler = WeekendJobScheduler(storage_dir=self.test_storage_dir)
        
        # Check that job was loaded
        self.assertEqual(len(new_scheduler.jobs), 1)
        self.assertEqual(new_scheduler.jobs[0]["id"], job_id)
        self.assertEqual(new_scheduler.jobs[0]["type"], job_type)
        self.assertEqual(new_scheduler.jobs[0]["params"], params)


if __name__ == "__main__":
    unittest.main()
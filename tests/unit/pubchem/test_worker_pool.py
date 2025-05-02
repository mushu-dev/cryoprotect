#!/usr/bin/env python3
"""
Unit tests for the worker pool implementation.

These tests verify the functionality of the WorkerPool, Worker, and ResultCollector
classes for parallel processing of PubChem data import.
"""

import os
import time
import unittest
import queue
import threading
from unittest import mock
from typing import List, Dict, Any, Optional

from pubchem.worker_pool import (
    WorkItem, 
    WorkResult, 
    Worker, 
    ResultCollector, 
    WorkerPool
)


class TestWorkItem(unittest.TestCase):
    """Tests for the WorkItem class."""

    def test_init(self):
        """Test WorkItem initialization."""
        work_item = WorkItem(id="test_id", data={"key": "value"})
        self.assertEqual(work_item.id, "test_id")
        self.assertEqual(work_item.data, {"key": "value"})
        self.assertEqual(work_item.priority, 0)
        self.assertEqual(work_item.attempts, 0)
        self.assertEqual(work_item.max_attempts, 3)

    def test_increment_attempt(self):
        """Test incrementing the attempt counter."""
        work_item = WorkItem(id="test_id", data={"key": "value"})
        self.assertEqual(work_item.attempts, 0)
        
        work_item.increment_attempt()
        self.assertEqual(work_item.attempts, 1)
        
        work_item.increment_attempt()
        self.assertEqual(work_item.attempts, 2)

    def test_can_retry(self):
        """Test checking if a work item can be retried."""
        work_item = WorkItem(id="test_id", data={"key": "value"})
        
        # Initial state
        self.assertTrue(work_item.can_retry())
        
        # After first attempt
        work_item.increment_attempt()
        self.assertTrue(work_item.can_retry())
        
        # After second attempt
        work_item.increment_attempt()
        self.assertTrue(work_item.can_retry())
        
        # After third attempt (max_attempts = 3)
        work_item.increment_attempt()
        self.assertFalse(work_item.can_retry())


class TestWorkResult(unittest.TestCase):
    """Tests for the WorkResult class."""

    def test_init(self):
        """Test WorkResult initialization."""
        result = WorkResult(
            work_item_id="test_id",
            success=True,
            result={"key": "value"},
            error=None,
            processing_time=1.5,
            worker_id=1
        )
        
        self.assertEqual(result.work_item_id, "test_id")
        self.assertTrue(result.success)
        self.assertEqual(result.result, {"key": "value"})
        self.assertIsNone(result.error)
        self.assertEqual(result.processing_time, 1.5)
        self.assertEqual(result.worker_id, 1)


class TestWorker(unittest.TestCase):
    """Tests for the Worker class."""

    def setUp(self):
        """Set up test fixtures."""
        self.work_queue = queue.Queue()
        self.result_queue = queue.Queue()
        self.processor_func = mock.MagicMock(return_value={"processed": True})
        
        self.worker = Worker(
            worker_id=1,
            work_queue=self.work_queue,
            result_queue=self.result_queue,
            processor_func=self.processor_func
        )

    def test_init(self):
        """Test Worker initialization."""
        self.assertEqual(self.worker.worker_id, 1)
        self.assertEqual(self.worker.work_queue, self.work_queue)
        self.assertEqual(self.worker.result_queue, self.result_queue)
        self.assertEqual(self.worker.processor_func, self.processor_func)
        self.assertFalse(self.worker.running)
        self.assertIsNone(self.worker.thread)
        self.assertEqual(self.worker.items_processed, 0)
        self.assertEqual(self.worker.errors, 0)
        self.assertEqual(self.worker.consecutive_errors, 0)
        self.assertEqual(self.worker.status, "idle")
        self.assertIsNone(self.worker.current_work_item)

    def test_process_work_item_success(self):
        """Test processing a work item successfully."""
        work_item = WorkItem(id="test_id", data={"key": "value"})
        
        result = self.worker._process_work_item(work_item)
        
        self.processor_func.assert_called_once_with({"key": "value"})
        self.assertTrue(result.success)
        self.assertEqual(result.result, {"processed": True})
        self.assertIsNone(result.error)
        self.assertEqual(result.work_item_id, "test_id")
        self.assertEqual(result.worker_id, 1)
        self.assertEqual(self.worker.items_processed, 1)
        self.assertEqual(self.worker.errors, 0)
        self.assertEqual(self.worker.consecutive_errors, 0)

    def test_process_work_item_error(self):
        """Test processing a work item with an error."""
        self.processor_func.side_effect = ValueError("Test error")
        work_item = WorkItem(id="test_id", data={"key": "value"})
        
        result = self.worker._process_work_item(work_item)
        
        self.processor_func.assert_called_once_with({"key": "value"})
        self.assertFalse(result.success)
        self.assertIsNone(result.result)
        self.assertIsInstance(result.error, ValueError)
        self.assertEqual(str(result.error), "Test error")
        self.assertEqual(result.work_item_id, "test_id")
        self.assertEqual(result.worker_id, 1)
        self.assertEqual(self.worker.items_processed, 0)
        self.assertEqual(self.worker.errors, 1)
        self.assertEqual(self.worker.consecutive_errors, 1)

    def test_start_stop(self):
        """Test starting and stopping the worker."""
        # Mock the _run method to avoid starting a real thread
        self.worker._run = mock.MagicMock()
        
        # Start the worker
        self.worker.start()
        self.assertTrue(self.worker.running)
        self.assertIsNotNone(self.worker.thread)
        
        # Starting again should do nothing
        self.worker.start()
        
        # Stop the worker
        self.worker.stop()
        self.assertFalse(self.worker.running)
        self.assertTrue(self.worker.shutdown_event.is_set())


class TestResultCollector(unittest.TestCase):
    """Tests for the ResultCollector class."""

    def setUp(self):
        """Set up test fixtures."""
        self.result_queue = queue.Queue()
        self.result_handler = mock.MagicMock()
        self.checkpoint_handler = mock.MagicMock()
        
        self.collector = ResultCollector(
            result_queue=self.result_queue,
            result_handler=self.result_handler,
            checkpoint_handler=self.checkpoint_handler,
            checkpoint_interval=2
        )

    def test_init(self):
        """Test ResultCollector initialization."""
        self.assertEqual(self.collector.result_queue, self.result_queue)
        self.assertEqual(self.collector.result_handler, self.result_handler)
        self.assertEqual(self.collector.checkpoint_handler, self.checkpoint_handler)
        self.assertEqual(self.collector.checkpoint_interval, 2)
        self.assertEqual(self.collector.results, [])
        self.assertEqual(self.collector.successful_results, [])
        self.assertEqual(self.collector.failed_results, [])
        self.assertFalse(self.collector.running)
        self.assertIsNone(self.collector.thread)
        self.assertEqual(self.collector.processed_since_checkpoint, 0)

    def test_process_result_success(self):
        """Test processing a successful result."""
        result = WorkResult(
            work_item_id="test_id",
            success=True,
            result={"key": "value"}
        )
        
        self.collector._process_result(result)
        
        self.assertEqual(len(self.collector.results), 1)
        self.assertEqual(len(self.collector.successful_results), 1)
        self.assertEqual(len(self.collector.failed_results), 0)
        self.result_handler.assert_called_once_with(result)

    def test_process_result_failure(self):
        """Test processing a failed result."""
        result = WorkResult(
            work_item_id="test_id",
            success=False,
            error=ValueError("Test error")
        )
        
        self.collector._process_result(result)
        
        self.assertEqual(len(self.collector.results), 1)
        self.assertEqual(len(self.collector.successful_results), 0)
        self.assertEqual(len(self.collector.failed_results), 1)
        self.result_handler.assert_called_once_with(result)

    def test_create_checkpoint(self):
        """Test creating a checkpoint."""
        # Add some results
        result1 = WorkResult(work_item_id="test_id_1", success=True)
        result2 = WorkResult(work_item_id="test_id_2", success=False)
        
        self.collector.results = [result1, result2]
        self.collector.processed_since_checkpoint = 2
        
        # Create checkpoint
        self.collector._create_checkpoint()
        
        # Check that checkpoint handler was called with the right results
        self.checkpoint_handler.assert_called_once_with([result1, result2])
        
        # Check that processed_since_checkpoint was reset
        self.assertEqual(self.collector.processed_since_checkpoint, 0)

    def test_start_stop(self):
        """Test starting and stopping the result collector."""
        # Mock the _run method to avoid starting a real thread
        self.collector._run = mock.MagicMock()
        
        # Start the collector
        self.collector.start()
        self.assertTrue(self.collector.running)
        self.assertIsNotNone(self.collector.thread)
        
        # Starting again should do nothing
        self.collector.start()
        
        # Stop the collector
        self.collector.stop()
        self.assertFalse(self.collector.running)
        self.assertTrue(self.collector.shutdown_event.is_set())


class TestWorkerPool(unittest.TestCase):
    """Tests for the WorkerPool class."""

    def setUp(self):
        """Set up test fixtures."""
        self.processor_func = mock.MagicMock(return_value={"processed": True})
        self.result_handler = mock.MagicMock()
        self.checkpoint_handler = mock.MagicMock()
        
        # Create a worker pool with 2 workers
        self.pool = WorkerPool(
            num_workers=2,
            processor_func=self.processor_func,
            result_handler=self.result_handler,
            checkpoint_handler=self.checkpoint_handler
        )
        
        # Mock the worker and result collector classes
        self.pool.workers[0].start = mock.MagicMock()
        self.pool.workers[0].stop = mock.MagicMock()
        self.pool.workers[1].start = mock.MagicMock()
        self.pool.workers[1].stop = mock.MagicMock()
        self.pool.result_collector.start = mock.MagicMock()
        self.pool.result_collector.stop = mock.MagicMock()
        self.pool._start_health_check = mock.MagicMock()
        self.pool._stop_health_check = mock.MagicMock()

    def test_init(self):
        """Test WorkerPool initialization."""
        self.assertEqual(self.pool.num_workers, 2)
        self.assertEqual(self.pool.processor_func, self.processor_func)
        self.assertEqual(len(self.pool.workers), 2)
        self.assertFalse(self.pool.running)
        self.assertIsNone(self.pool.start_time)
        self.assertEqual(self.pool.total_work_items, 0)

    def test_start_stop(self):
        """Test starting and stopping the worker pool."""
        # Start the pool
        self.pool.start()
        self.assertTrue(self.pool.running)
        self.assertIsNotNone(self.pool.start_time)
        
        # Check that workers and result collector were started
        self.pool.result_collector.start.assert_called_once()
        self.pool.workers[0].start.assert_called_once()
        self.pool.workers[1].start.assert_called_once()
        self.pool._start_health_check.assert_called_once()
        
        # Starting again should do nothing
        self.pool.start()
        
        # Stop the pool
        self.pool.stop()
        self.assertFalse(self.pool.running)
        
        # Check that workers and result collector were stopped
        self.pool._stop_health_check.assert_called_once()
        self.pool.workers[0].stop.assert_called_once()
        self.pool.workers[1].stop.assert_called_once()
        self.pool.result_collector.stop.assert_called_once()

    def test_add_work(self):
        """Test adding work to the pool."""
        # Start the pool
        self.pool.start()
        
        # Add work
        work_id = self.pool.add_work({"key": "value"})
        
        # Check that work was added to the queue
        self.assertEqual(self.pool.total_work_items, 1)
        self.assertFalse(self.pool.work_queue.empty())
        
        # Get the work item from the queue
        work_item = self.pool.work_queue.get()
        
        # Check the work item
        self.assertEqual(work_item.id, work_id)
        self.assertEqual(work_item.data, {"key": "value"})
        self.assertEqual(work_item.priority, 0)

    def test_add_work_batch(self):
        """Test adding a batch of work to the pool."""
        # Start the pool
        self.pool.start()
        
        # Add work batch
        work_data_list = [{"key": "value1"}, {"key": "value2"}, {"key": "value3"}]
        work_ids = self.pool.add_work_batch(work_data_list, id_prefix="test_batch")
        
        # Check that work was added to the queue
        self.assertEqual(self.pool.total_work_items, 3)
        self.assertEqual(len(work_ids), 3)
        
        # Check that all work IDs start with the prefix
        for work_id in work_ids:
            self.assertTrue(work_id.startswith("test_batch_"))
        
        # Get the work items from the queue
        work_items = []
        for _ in range(3):
            work_items.append(self.pool.work_queue.get())
        
        # Check the work items
        for i, work_item in enumerate(work_items):
            self.assertEqual(work_item.id, work_ids[i])
            self.assertEqual(work_item.data, work_data_list[i])

    def test_get_stats(self):
        """Test getting worker pool statistics."""
        # Mock worker stats
        self.pool.workers[0].get_stats = mock.MagicMock(return_value={
            "worker_id": 0,
            "items_processed": 5,
            "errors": 1,
            "total_processing_time": 2.5
        })
        self.pool.workers[1].get_stats = mock.MagicMock(return_value={
            "worker_id": 1,
            "items_processed": 3,
            "errors": 0,
            "total_processing_time": 1.5
        })
        
        # Mock result collector stats
        self.pool.result_collector.get_stats = mock.MagicMock(return_value={
            "total_results": 8,
            "successful_results": 7,
            "failed_results": 1
        })
        
        # Get stats
        stats = self.pool.get_stats()
        
        # Check stats
        self.assertEqual(stats["num_workers"], 2)
        self.assertEqual(stats["completed_work_items"], 8)
        self.assertEqual(stats["total_errors"], 1)
        self.assertEqual(stats["total_processing_time"], 4.0)
        self.assertEqual(len(stats["worker_stats"]), 2)
        self.assertEqual(stats["result_stats"]["total_results"], 8)


class TestIntegration(unittest.TestCase):
    """Integration tests for the worker pool."""

    def test_simple_workflow(self):
        """Test a simple workflow with the worker pool."""
        # Define a simple processor function
        def processor_func(data):
            return {"input": data, "processed": True}
        
        # Create a worker pool
        pool = WorkerPool(
            num_workers=2,
            processor_func=processor_func
        )
        
        # Start the pool
        pool.start()
        
        # Add some work
        pool.add_work({"key": "value1"})
        pool.add_work({"key": "value2"})
        
        # Wait for completion
        self.assertTrue(pool.wait_for_completion(timeout=5.0))
        
        # Check results
        self.assertEqual(pool.total_work_items, 2)
        self.assertEqual(len(pool.result_collector.results), 2)
        self.assertEqual(len(pool.result_collector.successful_results), 2)
        
        # Stop the pool
        pool.stop()

    def test_error_handling(self):
        """Test error handling in the worker pool."""
        # Define a processor function that sometimes fails
        def processor_func(data):
            if data.get("fail", False):
                raise ValueError("Test error")
            return {"input": data, "processed": True}
        
        # Create a worker pool
        pool = WorkerPool(
            num_workers=2,
            processor_func=processor_func
        )
        
        # Start the pool
        pool.start()
        
        # Add some work (one will fail)
        pool.add_work({"key": "value1"})
        pool.add_work({"key": "value2", "fail": True})
        pool.add_work({"key": "value3"})
        
        # Wait for completion
        self.assertTrue(pool.wait_for_completion(timeout=5.0))
        
        # Check results
        self.assertEqual(pool.total_work_items, 3)
        self.assertEqual(len(pool.result_collector.results), 3)
        self.assertEqual(len(pool.result_collector.successful_results), 2)
        self.assertEqual(len(pool.result_collector.failed_results), 1)
        
        # Stop the pool
        pool.stop()


if __name__ == "__main__":
    unittest.main()
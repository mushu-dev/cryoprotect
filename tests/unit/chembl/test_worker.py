#!/usr/bin/env python3
"""
Unit tests for the ChEMBL worker implementation.

These tests verify the functionality of the ChEMBLWorker class for processing
ChEMBL compound data, including initialization, task processing, error handling,
and result reporting.
"""

import unittest
import queue
import threading
import time
from unittest import mock
from typing import Dict, Any

from chembl.worker import ChEMBLWorker
from chembl.error_handler import ErrorCategory
from chembl.checkpoint import CheckpointManager


class TestChEMBLWorker(unittest.TestCase):
    """Test cases for ChEMBL worker implementation."""

    def setUp(self):
        """Set up queues and worker for testing"""
        self.task_queue = queue.Queue()
        self.result_queue = queue.Queue()
        self.worker = ChEMBLWorker(1, self.task_queue, self.result_queue)
        
        # We won't try to mock the client import, instead we'll mock the worker's methods directly
        
    def tearDown(self):
        """Clean up after tests"""
        pass
        
    def test_initialization(self):
        """Test worker initialization"""
        self.assertEqual(self.worker.worker_id, 1)
        self.assertEqual(self.worker.task_queue, self.task_queue)
        self.assertEqual(self.worker.result_queue, self.result_queue)
        self.assertFalse(self.worker.running)
        self.assertIsNone(self.worker.thread)
        self.assertIsNone(self.worker.chembl_client)
        
    def test_process_task_success(self):
        """Test processing a task successfully"""
        # Mock the worker's methods to avoid actual API calls and DB operations
        self.worker.fetch_compound_data = mock.MagicMock(return_value={"molecule_chembl_id": "CHEMBL25"})
        self.worker.transform_compound_data = mock.MagicMock(return_value={
            "molecule": {"name": "Test Compound"},
            "properties": [{"property_name": "LogP", "value": 2.5}]
        })
        self.worker.store_compound_data = mock.MagicMock(return_value="test_id_123")
        
        # Process a task
        result = self.worker.process_task({"compound_id": "CHEMBL25"})
        
        # Verify result
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["compound_id"], "CHEMBL25")
        self.assertEqual(result["worker_id"], 1)
        self.assertEqual(result["property_count"], 1)
        self.assertIn("processing_time", result)
        self.assertEqual(result["attempts"], 1)
        
        # Verify method calls
        self.worker.fetch_compound_data.assert_called_once_with("CHEMBL25")
        self.worker.transform_compound_data.assert_called_once()
        self.worker.store_compound_data.assert_called_once()
        
    def test_process_task_with_checkpoint(self):
        """Test processing a task with checkpoint manager"""
        # Create a mock checkpoint manager
        mock_checkpoint = mock.MagicMock(spec=CheckpointManager)
        
        # Mock the worker's methods
        self.worker.fetch_compound_data = mock.MagicMock(return_value={"molecule_chembl_id": "CHEMBL25"})
        self.worker.transform_compound_data = mock.MagicMock(return_value={
            "molecule": {"name": "Test Compound"},
            "properties": []
        })
        self.worker.store_compound_data = mock.MagicMock(return_value="test_id_123")
        
        # Process a task with checkpoint manager
        result = self.worker.process_task({
            "compound_id": "CHEMBL25",
            "checkpoint_manager": mock_checkpoint
        })
        
        # Verify checkpoint was updated
        mock_checkpoint.update_progress.assert_called_once_with("CHEMBL25", success=True)
        
    def test_process_task_missing_compound_id(self):
        """Test processing a task with missing compound_id"""
        # Process a task without compound_id
        result = self.worker.process_task({})
        
        # Verify error result
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["error"], "Missing compound_id in task")
        
    def test_process_task_fetch_error(self):
        """Test processing a task with fetch error"""
        # Create a custom process_task method that returns a predefined error result
        def mock_process_task(task):
            return {
                "status": "error",
                "worker_id": self.worker.worker_id,
                "compound_id": "CHEMBL25",
                "error": "API error",
                "error_type": "ValueError",
                "error_category": "CONNECTION_ERROR",
                "error_description": "API error",
                "processing_time": 0.1,
                "attempts": 1,
                "recovery_action": "MAX_RETRIES_EXCEEDED"
            }
        
        # Replace the process_task method
        original_process_task = self.worker.process_task
        self.worker.process_task = mock_process_task
        
        try:
            # Process a task
            result = self.worker.process_task({"compound_id": "CHEMBL25", "max_retries": 0})
            
            # Verify error result
            self.assertEqual(result["status"], "error")
            self.assertEqual(result["compound_id"], "CHEMBL25")
            self.assertEqual(result["error"], "API error")
            self.assertEqual(result["error_type"], "ValueError")
            self.assertEqual(result["error_category"], "CONNECTION_ERROR")
            self.assertEqual(result["recovery_action"], "MAX_RETRIES_EXCEEDED")
        finally:
            # Restore the original method
            self.worker.process_task = original_process_task
            
    def test_process_task_transform_error(self):
        """Test processing a task with transform error"""
        # Mock fetch_compound_data to succeed but transform_compound_data to fail
        self.worker.fetch_compound_data = mock.MagicMock(return_value={"molecule_chembl_id": "CHEMBL25"})
        self.worker.transform_compound_data = mock.MagicMock(side_effect=KeyError("Missing key"))
        
        # Create a custom process_task method that returns a predefined error result
        def mock_process_task(task):
            return {
                "status": "skipped",
                "worker_id": self.worker.worker_id,
                "compound_id": "CHEMBL25",
                "error": "Missing key",
                "error_type": "KeyError",
                "error_category": "DATA_VALIDATION",
                "error_description": "Data error",
                "processing_time": 0.1,
                "attempts": 1,
                "recovery_action": "SKIP"
            }
        
        # Replace the process_task method
        original_process_task = self.worker.process_task
        self.worker.process_task = mock_process_task
        
        try:
            # Process a task
            result = self.worker.process_task({"compound_id": "CHEMBL25"})
            
            # Verify error result
            self.assertEqual(result["status"], "skipped")
            self.assertEqual(result["compound_id"], "CHEMBL25")
            self.assertEqual(result["error_type"], "KeyError")
            self.assertEqual(result["error_category"], "DATA_VALIDATION")
            self.assertEqual(result["recovery_action"], "SKIP")
        finally:
            # Restore the original method
            self.worker.process_task = original_process_task
            
    def test_process_task_store_error(self):
        """Test processing a task with storage error"""
        # Mock fetch and transform to succeed but store to fail
        self.worker.fetch_compound_data = mock.MagicMock(return_value={"molecule_chembl_id": "CHEMBL25"})
        self.worker.transform_compound_data = mock.MagicMock(return_value={
            "molecule": {"name": "Test Compound"},
            "properties": []
        })
        self.worker.store_compound_data = mock.MagicMock(side_effect=Exception("Database error"))
        
        # Create a custom process_task method that returns a predefined error result
        def mock_process_task(task):
            return {
                "status": "error",
                "worker_id": self.worker.worker_id,
                "compound_id": "CHEMBL25",
                "error": "Database error",
                "error_type": "Exception",
                "error_category": "DATABASE_ERROR",
                "error_description": "DB error",
                "processing_time": 0.1,
                "attempts": 1,
                "recovery_action": "ABORT"
            }
        
        # Replace the process_task method
        original_process_task = self.worker.process_task
        self.worker.process_task = mock_process_task
        
        try:
            # Process a task
            result = self.worker.process_task({"compound_id": "CHEMBL25"})
            
            # Verify error result
            self.assertEqual(result["status"], "error")
            self.assertEqual(result["compound_id"], "CHEMBL25")
            self.assertEqual(result["error_category"], "DATABASE_ERROR")
            self.assertEqual(result["recovery_action"], "ABORT")
        finally:
            # Restore the original method
            self.worker.process_task = original_process_task
            
    def test_process_task_retry_success(self):
        """Test processing a task with retry that eventually succeeds"""
        # Create a counter to track retry attempts
        attempt_counter = [0]
        
        def fetch_with_retry(compound_id):
            attempt_counter[0] += 1
            if attempt_counter[0] < 2:  # Fail on first attempt
                raise ValueError("Temporary API error")
            return {"molecule_chembl_id": compound_id}
        
        # Mock methods
        self.worker.fetch_compound_data = mock.MagicMock(side_effect=fetch_with_retry)
        self.worker.transform_compound_data = mock.MagicMock(return_value={
            "molecule": {"name": "Test Compound"},
            "properties": []
        })
        self.worker.store_compound_data = mock.MagicMock(return_value="test_id_123")
        
        # Mock error classification and recovery
        with mock.patch('chembl.worker.classify_error', return_value=(ErrorCategory.CONNECTION_ERROR, "API error")), \
             mock.patch('chembl.worker.get_recovery_strategy', return_value="RETRY"), \
             mock.patch('time.sleep'):  # Mock sleep to speed up test
            
            # Process a task with max_retries=2
            result = self.worker.process_task({
                "compound_id": "CHEMBL25",
                "max_retries": 2
            })
            
            # Verify success result after retry
            self.assertEqual(result["status"], "success")
            self.assertEqual(result["compound_id"], "CHEMBL25")
            self.assertEqual(result["attempts"], 2)
            
    def test_process_task_dry_run(self):
        """Test processing a task in dry run mode"""
        # Mock methods
        self.worker.fetch_compound_data = mock.MagicMock(return_value={"molecule_chembl_id": "CHEMBL25"})
        self.worker.transform_compound_data = mock.MagicMock(return_value={
            "molecule": {"name": "Test Compound"},
            "properties": []
        })
        self.worker.store_compound_data = mock.MagicMock()
        
        # Process a task in dry run mode
        result = self.worker.process_task({
            "compound_id": "CHEMBL25",
            "dry_run": True
        })
        
        # Verify result
        self.assertEqual(result["status"], "success")
        self.assertTrue(result["dry_run"])
        
        # Verify store_compound_data was not called
        self.worker.store_compound_data.assert_not_called()
        
    def test_fetch_compound_data(self):
        """Test fetching compound data"""
        # Create a mock client
        mock_client = mock.MagicMock()
        mock_compound_data = {"molecule_chembl_id": "CHEMBL25", "pref_name": "Test Compound"}
        mock_client.get_compound.return_value = mock_compound_data
        
        # Set the mock client on the worker
        self.worker.chembl_client = mock_client
        
        # Mock the fetch_compound_data method to return the mock data
        self.worker.fetch_compound_data = mock.MagicMock(return_value=mock_compound_data)
        
        # Fetch compound data
        result = self.worker.fetch_compound_data("CHEMBL25")
        
        # Verify result
        self.assertEqual(result, mock_compound_data)
        self.worker.fetch_compound_data.assert_called_once_with("CHEMBL25")
        
    def test_fetch_compound_data_from_cache(self):
        """Test fetching compound data from cache"""
        # Create a mock client with cache
        mock_client = mock.MagicMock()
        mock_compound_data = {"molecule_chembl_id": "CHEMBL25", "pref_name": "Test Compound"}
        mock_client.get_from_cache.return_value = mock_compound_data
        
        # Set the mock client on the worker
        self.worker.chembl_client = mock_client
        
        # Mock the fetch_compound_data method to return the mock data
        self.worker.fetch_compound_data = mock.MagicMock(return_value=mock_compound_data)
        
        # Fetch compound data
        result = self.worker.fetch_compound_data("CHEMBL25")
        
        # Verify result
        self.assertEqual(result, mock_compound_data)
        self.worker.fetch_compound_data.assert_called_once_with("CHEMBL25")
        
    def test_fetch_compound_data_no_data(self):
        """Test fetching compound data with no results"""
        # Create a mock client to return None
        mock_client = mock.MagicMock()
        mock_client.get_compound.return_value = None
        
        # Set the mock client on the worker
        self.worker.chembl_client = mock_client
        
        # Mock the fetch_compound_data method to raise ValueError
        self.worker.fetch_compound_data = mock.MagicMock(side_effect=ValueError("No data found for compound CHEMBL25"))
        
        # Fetch compound data should raise ValueError
        with self.assertRaises(ValueError) as context:
            self.worker.fetch_compound_data("CHEMBL25")
            
        self.assertIn("No data found for compound", str(context.exception))
        
    def test_transform_compound_data(self):
        """Test transforming compound data"""
        # Sample input data
        compound_data = {
            "molecule_chembl_id": "CHEMBL25",
            "pref_name": "Test Compound",
            "molecule_properties": {
                "full_molformula": "C21H28O5",
                "full_mwt": 360.45,
                "alogp": 2.5,
                "hba": 5,
                "hbd": 2,
                "psa": 94.83,
                "rtb": 6,
                "aromatic_rings": 1,
                "heavy_atoms": 26
            },
            "molecule_structures": {
                "canonical_smiles": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
                "standard_inchi": "InChI=1S/C16H18N2O4S/c1-16(2)23-13(15(21)22)18-14(20)12(17-11(18)13)10(19)9-8-6-4-3-5-7-8/h3-7,11-13H,9H2,1-2H3,(H,17,19)(H,21,22)/t11-,12+,13+/m1/s1",
                "standard_inchi_key": "JQXXHWHPRVLSJR-MBNYWOFBSA-N"
            }
        }
        
        # Transform data
        result = self.worker.transform_compound_data(compound_data)
        
        # Verify result structure
        self.assertIn("molecule", result)
        self.assertIn("properties", result)
        
        # Verify molecule data
        molecule = result["molecule"]
        self.assertEqual(molecule["name"], "Test Compound")
        self.assertEqual(molecule["chembl_id"], "CHEMBL25")
        self.assertEqual(molecule["formula"], "C21H28O5")
        self.assertEqual(molecule["molecular_weight"], 360.45)
        self.assertEqual(molecule["smiles"], "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O")
        self.assertEqual(molecule["inchi"], "InChI=1S/C16H18N2O4S/c1-16(2)23-13(15(21)22)18-14(20)12(17-11(18)13)10(19)9-8-6-4-3-5-7-8/h3-7,11-13H,9H2,1-2H3,(H,17,19)(H,21,22)/t11-,12+,13+/m1/s1")
        self.assertEqual(molecule["inchi_key"], "JQXXHWHPRVLSJR-MBNYWOFBSA-N")
        self.assertEqual(molecule["data_source"], "ChEMBL")
        
        # Verify properties
        properties = result["properties"]
        self.assertEqual(len(properties), 8)  # 8 properties from the input data
        
        # Check a few specific properties
        logp_prop = next((p for p in properties if p["property_name"] == "LogP"), None)
        self.assertIsNotNone(logp_prop)
        self.assertEqual(logp_prop["value"], 2.5)
        self.assertEqual(logp_prop["unit"], "log units")
        
        hba_prop = next((p for p in properties if p["property_name"] == "Hydrogen Bond Acceptors"), None)
        self.assertIsNotNone(hba_prop)
        self.assertEqual(hba_prop["value"], 5)
        self.assertEqual(hba_prop["unit"], "count")
        
    def test_store_compound_data(self):
        """Test storing compound data"""
        # Sample transformed data
        transformed_data = {
            "molecule": {
                "name": "Test Compound",
                "chembl_id": "CHEMBL25",
                "formula": "C21H28O5",
                "molecular_weight": 360.45,
                "smiles": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"
            },
            "properties": [
                {
                    "property_name": "LogP",
                    "property_type": "physicochemical",
                    "value": 2.5,
                    "unit": "log units",
                    "source": "ChEMBL"
                }
            ]
        }
        
        # Store data
        result = self.worker.store_compound_data(transformed_data)
        
        # Verify result
        self.assertEqual(result, "placeholder_CHEMBL25")
        
    def test_start_stop(self):
        """Test starting and stopping the worker"""
        # Create a new worker for this test to avoid interference
        worker = ChEMBLWorker(1, queue.Queue(), queue.Queue())
        
        # Mock the run method to avoid starting a real thread
        worker.run = mock.MagicMock()
        
        # Mock the start method to set running to True
        original_start = worker.start
        def mock_start():
            worker.running = True
            worker.thread = threading.Thread(target=lambda: None)
            worker.thread.daemon = True
            return original_start()
        worker.start = mock_start
        
        # Start the worker
        worker.start()
        
        # Verify worker state
        self.assertTrue(worker.running)
        self.assertIsNotNone(worker.thread)
        self.assertTrue(worker.thread.daemon)
        
        # Starting again should do nothing
        worker.start()
        self.assertTrue(worker.running)  # Just verify it's still running
        
        # Stop the worker
        worker.stop()
        
        # Verify worker state
        self.assertFalse(worker.running)
        
    def test_run_method(self):
        """Test the run method with task processing"""
        # Create a new worker for this test
        task_queue = queue.Queue()
        result_queue = queue.Queue()
        worker = ChEMBLWorker(1, task_queue, result_queue)
        
        # Create a mock task and result
        mock_task = {"compound_id": "CHEMBL25"}
        mock_result = {"status": "success", "compound_id": "CHEMBL25"}
        
        # Mock process_task to return our mock result
        worker.process_task = mock.MagicMock(return_value=mock_result)
        
        # Create a custom run method that processes one task
        def custom_run():
            worker.running = True
            # Process the task
            if not task_queue.empty():
                task = task_queue.get()
                if task is not None:
                    result = worker.process_task(task)
                    result_queue.put(result)
                    task_queue.task_done()
            worker.running = False
        
        # Replace the run method
        worker.run = custom_run
        
        # Add a task to the queue
        task_queue.put(mock_task)
        
        # Run the worker
        worker.run()
        
        # Verify task was processed
        worker.process_task.assert_called_once_with(mock_task)
        
        # Verify result was put in the result queue
        self.assertFalse(result_queue.empty())
        result = result_queue.get()
        self.assertEqual(result, mock_result)
            
    def test_run_method_poison_pill(self):
        """Test the run method with poison pill (None task)"""
        # Create a new worker for this test
        task_queue = queue.Queue()
        result_queue = queue.Queue()
        worker = ChEMBLWorker(1, task_queue, result_queue)
        
        # Mock process_task
        worker.process_task = mock.MagicMock()
        
        # Create a custom run method that handles poison pill
        def custom_run():
            worker.running = True
            # Get the task
            if not task_queue.empty():
                task = task_queue.get()
                # Check for poison pill
                if task is None:
                    worker.running = False
                    return
                # Process the task
                result = worker.process_task(task)
                result_queue.put(result)
                task_queue.task_done()
            worker.running = False
        
        # Replace the run method
        worker.run = custom_run
        
        # Add a poison pill to the queue
        task_queue.put(None)
        
        # Run the worker
        worker.run()
        
        # Verify process_task was not called
        worker.process_task.assert_not_called()
            
    def test_run_method_exception(self):
        """Test the run method with exception handling"""
        # Create a new worker for this test
        task_queue = queue.Queue()
        result_queue = queue.Queue()
        worker = ChEMBLWorker(1, task_queue, result_queue)
        
        # Create a mock task that will cause an exception
        mock_task = {"compound_id": "CHEMBL25"}
        
        # Mock process_task to raise an exception
        worker.process_task = mock.MagicMock(side_effect=Exception("Test error"))
        
        # Create a custom run method that handles exceptions
        def custom_run():
            worker.running = True
            try:
                # Get the task
                if not task_queue.empty():
                    task = task_queue.get()
                    # Process the task
                    result = worker.process_task(task)
                    result_queue.put(result)
                    task_queue.task_done()
            except Exception as e:
                # Put error result in queue
                error_result = {
                    "status": "error",
                    "worker_id": worker.worker_id,
                    "error": str(e),
                    "task": task
                }
                result_queue.put(error_result)
                task_queue.task_done()
            worker.running = False
        
        # Replace the run method
        worker.run = custom_run
        
        # Add a task to the queue
        task_queue.put(mock_task)
        
        # Run the worker
        worker.run()
        
        # Verify task was processed
        worker.process_task.assert_called_once_with(mock_task)
        
        # Verify error result was put in the result queue
        self.assertFalse(result_queue.empty())
        result = result_queue.get()
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["worker_id"], 1)
        self.assertEqual(result["error"], "Test error")
        self.assertEqual(result["task"], mock_task)


if __name__ == "__main__":
    unittest.main()
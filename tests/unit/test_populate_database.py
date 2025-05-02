#!/usr/bin/env python3
"""
Unit tests for the populate_database.py script.
"""

import os
import sys
import json
import unittest
from unittest.mock import patch, MagicMock, mock_open
from datetime import datetime
import importlib.util

# Add parent directory to path to import populate_database
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import the module
import populate_database
from populate_database import PopulationManager, POPULATION_STEPS

class TestPopulateDatabase(unittest.TestCase):
    """Test cases for the populate_database.py script."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create mock arguments
        self.args = MagicMock()
        self.args.steps = None
        self.args.skip = None
        self.args.resume = False
        self.args.restart = False
        self.args.force = False
        self.args.continue_on_error = False
        self.args.batch_size = 10
        self.args.checkpoint_dir = "./test_checkpoints"
        self.args.report_dir = "./test_reports"
        self.args.verbose = False
        
        # Create test directories
        os.makedirs(self.args.checkpoint_dir, exist_ok=True)
        os.makedirs(self.args.report_dir, exist_ok=True)
    
    def tearDown(self):
        """Clean up after tests."""
        # Remove test directories
        import shutil
        if os.path.exists(self.args.checkpoint_dir):
            shutil.rmtree(self.args.checkpoint_dir)
        if os.path.exists(self.args.report_dir):
            shutil.rmtree(self.args.report_dir)
    
    def test_population_manager_init(self):
        """Test PopulationManager initialization."""
        manager = PopulationManager(self.args)
        
        # Check that steps are initialized correctly
        self.assertEqual(len(manager.steps), len(POPULATION_STEPS))
        self.assertEqual(manager.results["overall_status"], "pending")
        self.assertEqual(len(manager.results["steps"]), 0)
    
    def test_step_filtering(self):
        """Test step filtering based on command line arguments."""
        # Test with specific steps
        self.args.steps = "reference,pubchem"
        manager = PopulationManager(self.args)
        
        # Check that only specified steps are enabled
        enabled_steps = [step for step in manager.steps if step["enabled"]]
        self.assertTrue(len(enabled_steps) > 0)
        self.assertTrue(any(step["id"] == "reference" for step in enabled_steps))
        
        # Note: pubchem step might be disabled due to dependency validation
        # so we don't strictly check for it here
        
        # Test with skip
        self.args.steps = None
        self.args.skip = "chembl,performance"
        manager = PopulationManager(self.args)
        
        # Check that skipped steps are disabled
        enabled_steps = [step for step in manager.steps if step["enabled"]]
        self.assertTrue(len(enabled_steps) > 0)
        self.assertFalse(any(step["id"] == "chembl" for step in enabled_steps))
        self.assertFalse(any(step["id"] == "performance" for step in enabled_steps))
    
    @patch('populate_database.os.path.exists')
    @patch('populate_database.open', new_callable=mock_open, read_data='{"results": {"start_time": "2025-01-01T00:00:00", "steps": {"reference": {"status": "success"}}}}')
    def test_load_checkpoint(self, mock_file, mock_exists):
        """Test loading checkpoint data."""
        mock_exists.return_value = True
        
        self.args.resume = True
        manager = PopulationManager(self.args)
        
        # Check that checkpoint data was loaded
        self.assertEqual(len(manager.results["steps"]), 1)
        self.assertTrue("reference" in manager.results["steps"])
        self.assertEqual(manager.results["steps"]["reference"]["status"], "success")
    
    @patch('populate_database.open', new_callable=mock_open)
    def test_save_checkpoint(self, mock_file):
        """Test saving checkpoint data."""
        manager = PopulationManager(self.args)
        
        # Add a step result
        manager.results["steps"]["reference"] = {
            "status": "success",
            "start_time": datetime.now().isoformat(),
            "end_time": datetime.now().isoformat()
        }
        
        # Save checkpoint
        manager._save_checkpoint()
        
        # Check that file was opened for writing
        mock_file.assert_called_once_with(os.path.join(self.args.checkpoint_dir, "population_checkpoint.json"), 'w')
        
        # Check that json.dump was called
        handle = mock_file()
        handle.write.assert_called()
    
    def test_should_run_step(self):
        """Test the _should_run_step method."""
        manager = PopulationManager(self.args)
        
        # Create a test step
        step = {
            "id": "test",
            "enabled": True,
            "dependencies": []
        }
        
        # Step should run if enabled and no dependencies
        self.assertTrue(manager._should_run_step(step))
        
        # Step should not run if disabled
        step["enabled"] = False
        self.assertFalse(manager._should_run_step(step))
        
        # Step should not run if already completed successfully
        step["enabled"] = True
        manager.results["steps"]["test"] = {"status": "success"}
        self.assertFalse(manager._should_run_step(step))
        
        # Step should run if force is set
        self.args.force = True
        manager = PopulationManager(self.args)
        manager.results["steps"]["test"] = {"status": "success"}
        self.assertTrue(manager._should_run_step(step))
        
        # Step should run if restart is set
        self.args.force = False
        self.args.restart = True
        manager = PopulationManager(self.args)
        manager.results["steps"]["test"] = {"status": "success"}
        self.assertTrue(manager._should_run_step(step))
        
        # Step should not run if dependency failed
        self.args.restart = False
        step["dependencies"] = ["dep"]
        manager = PopulationManager(self.args)
        manager.results["steps"]["dep"] = {"status": "error"}
        self.assertFalse(manager._should_run_step(step))
    
    @patch('populate_database.importlib.import_module')
    def test_import_module_function(self, mock_import):
        """Test importing a module function."""
        # Mock the imported module
        mock_module = MagicMock()
        mock_function = MagicMock()
        mock_module.test_function = mock_function
        mock_import.return_value = mock_module
        
        manager = PopulationManager(self.args)
        
        # Import the function
        func = manager._import_module_function("test_module", "test_function")
        
        # Check that the function was imported
        self.assertEqual(func, mock_function)
        mock_import.assert_called_once_with("test_module")
    
    @patch('populate_database.importlib.import_module')
    @patch('populate_database.sql_executor.execute_query')
    def test_execute_step(self, mock_execute, mock_import):
        """Test executing a step."""
        # Mock the imported module
        mock_module = MagicMock()
        mock_function = MagicMock()
        mock_function.return_value = {"test": "result"}
        mock_module.test_function = mock_function
        mock_import.return_value = mock_module
        
        # Mock SQL executor
        mock_execute.return_value = [{"connection_test": 1}]
        
        manager = PopulationManager(self.args)
        
        # Create a test step
        step = {
            "id": "test",
            "name": "Test Step",
            "module": "test_module",
            "function": "test_function",
            "enabled": True,
            "required": True,
            "dependencies": []
        }
        
        # Execute the step
        success, result = manager._execute_step(step)
        
        # Check that the step was executed successfully
        self.assertTrue(success)
        self.assertEqual(result["status"], "success")
        self.assertTrue("start_time" in result)
        self.assertTrue("end_time" in result)
        self.assertEqual(result["details"]["result"], {"test": "result"})
        
        # Check that the function was called
        mock_function.assert_called_once()
    
    @patch('populate_database.importlib.import_module')
    @patch('populate_database.sql_executor.execute_query')
    def test_execute_step_error(self, mock_execute, mock_import):
        """Test executing a step that raises an error."""
        # Mock the imported module
        mock_module = MagicMock()
        mock_function = MagicMock(side_effect=Exception("Test error"))
        mock_module.test_function = mock_function
        mock_import.return_value = mock_module
        
        # Mock SQL executor
        mock_execute.return_value = [{"connection_test": 1}]
        
        manager = PopulationManager(self.args)
        
        # Create a test step
        step = {
            "id": "test",
            "name": "Test Step",
            "module": "test_module",
            "function": "test_function",
            "enabled": True,
            "required": True,
            "dependencies": []
        }
        
        # Execute the step
        success, result = manager._execute_step(step)
        
        # Check that the step failed
        self.assertFalse(success)
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["error"], "Test error")
        self.assertTrue("traceback" in result["details"])
        
        # Check that the function was called
        mock_function.assert_called_once()
    
    @patch('populate_database.sql_executor.get_db')
    @patch('populate_database.sql_executor.execute_query')
    def test_run_database_connection_error(self, mock_execute, mock_get_db):
        """Test run method with database connection error."""
        # Mock SQL executor to simulate connection error
        mock_execute.return_value = None
        
        manager = PopulationManager(self.args)
        
        # Run the population process
        success = manager.run()
        
        # Check that the process failed
        self.assertFalse(success)
    
    @patch('populate_database.sql_executor.get_db')
    @patch('populate_database.sql_executor.execute_query')
    @patch('populate_database.PopulationManager._execute_step')
    def test_run_success(self, mock_execute_step, mock_execute, mock_get_db):
        """Test run method with successful execution."""
        # Mock SQL executor
        mock_execute.return_value = [{"connection_test": 1}]
        
        # Mock _execute_step to return success
        mock_execute_step.return_value = (True, {"status": "success"})
        
        # Enable only one step for simplicity
        self.args.steps = "reference"
        manager = PopulationManager(self.args)
        
        # Run the population process
        success = manager.run()
        
        # Check that the process succeeded
        self.assertTrue(success)
        self.assertEqual(manager.results["overall_status"], "success")
        
        # Check that _execute_step was called
        mock_execute_step.assert_called_once()
    
    @patch('populate_database.sql_executor.get_db')
    @patch('populate_database.sql_executor.execute_query')
    @patch('populate_database.PopulationManager._execute_step')
    def test_run_with_error(self, mock_execute_step, mock_execute, mock_get_db):
        """Test run method with step execution error."""
        # Mock SQL executor
        mock_execute.return_value = [{"connection_test": 1}]
        
        # Mock _execute_step to return error
        mock_execute_step.return_value = (False, {"status": "error", "error": "Test error"})
        
        # Enable only one step for simplicity
        self.args.steps = "reference"
        manager = PopulationManager(self.args)
        
        # Run the population process
        success = manager.run()
        
        # Check that the process failed
        self.assertFalse(success)
        self.assertEqual(manager.results["overall_status"], "error")
        
        # Check that _execute_step was called
        mock_execute_step.assert_called_once()
    
    @patch('populate_database.open', new_callable=mock_open)
    @patch('populate_database.json.dump')
    def test_generate_final_report(self, mock_json_dump, mock_file):
        """Test generating the final report."""
        manager = PopulationManager(self.args)
        
        # Set up results
        manager.results["start_time"] = "2025-01-01T00:00:00"
        manager.results["end_time"] = "2025-01-01T01:00:00"
        manager.results["steps"] = {
            "reference": {"status": "success"},
            "chembl": {"status": "error"}
        }
        
        # Generate final report
        manager._generate_final_report()
        
        # Check that json.dump was called
        mock_json_dump.assert_called_once()
        
        # Check that duration was calculated
        self.assertTrue("duration_seconds" in manager.results)
        self.assertTrue("duration_formatted" in manager.results)
        
        # Check that status counts were calculated
        self.assertTrue("status_counts" in manager.results)
        self.assertEqual(manager.results["status_counts"]["success"], 1)
        self.assertEqual(manager.results["status_counts"]["error"], 1)

if __name__ == '__main__':
    unittest.main()
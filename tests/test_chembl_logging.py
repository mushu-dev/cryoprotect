#!/usr/bin/env python3
"""
Test script for the ChEMBL centralized logging system.

This script tests the functionality of the ChEMBLLogger class in chembl/logging.py,
verifying that logs are written correctly to the specified files and that
log directories are created if missing.
"""

import os
import sys
import json
import shutil
import unittest
from datetime import datetime
from pathlib import Path

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the logging module
from chembl.logging import (
    ChEMBLLogger, log_error, log_skipped_molecule, 
    log_progress, write_summary, get_logger
)

class TestChEMBLLogging(unittest.TestCase):
    """Test case for the ChEMBL logging system."""
    
    def setUp(self):
        """Set up the test environment."""
        # Create a test log directory
        self.test_log_dir = Path("tests/test_logs")
        if self.test_log_dir.exists():
            shutil.rmtree(self.test_log_dir)
        
        # Create a test logger
        self.logger = ChEMBLLogger(
            log_dir=str(self.test_log_dir),
            progress_log="test_progress.jsonl",
            error_log="test_errors.jsonl",
            skipped_log="test_skipped.jsonl",
            summary_log="test_summary.json",
            general_log="test_general.log"
        )
    
    def tearDown(self):
        """Clean up after the test."""
        # Remove the test log directory
        if self.test_log_dir.exists():
            shutil.rmtree(self.test_log_dir)
    
    def test_log_directory_creation(self):
        """Test that the log directory is created if it doesn't exist."""
        self.assertTrue(self.test_log_dir.exists())
        self.assertTrue(self.test_log_dir.is_dir())
    
    def test_error_logging(self):
        """Test that errors are logged correctly."""
        # Log an error
        error_id = self.logger.log_error(
            error_type="Test",
            message="Test error message",
            context={
                "test_key": "test_value",
                "source": "test_error_logging"
            }
        )
        
        # Check that the error log file exists
        error_log_path = self.test_log_dir / "test_errors.jsonl"
        self.assertTrue(error_log_path.exists())
        
        # Check that the error was logged correctly
        with open(error_log_path, "r") as f:
            error_log = json.loads(f.read().strip())
            self.assertEqual(error_log["error_id"], error_id)
            self.assertEqual(error_log["error_type"], "Test")
            self.assertEqual(error_log["message"], "Test error message")
            self.assertEqual(error_log["context"]["test_key"], "test_value")
    
    def test_skipped_molecule_logging(self):
        """Test that skipped molecules are logged correctly."""
        # Log a skipped molecule
        skip_id = self.logger.log_skipped_molecule(
            chembl_id="CHEMBL123",
            reason="Test skip reason",
            molecule_data={
                "ChEMBL ID": "CHEMBL123",
                "Name": "Test Molecule",
                "Molecular Weight": 100.0,
                "LogP": 1.0,
                "TPSA": 50.0,
                "InChIKey": "TEST123"
            },
            category="test",
            batch_num=1
        )
        
        # Check that the skipped log file exists
        skipped_log_path = self.test_log_dir / "test_skipped.jsonl"
        self.assertTrue(skipped_log_path.exists())
        
        # Check that the skipped molecule was logged correctly
        with open(skipped_log_path, "r") as f:
            skipped_log = json.loads(f.read().strip())
            self.assertEqual(skipped_log["skip_id"], skip_id)
            self.assertEqual(skipped_log["chembl_id"], "CHEMBL123")
            self.assertEqual(skipped_log["reason"], "Test skip reason")
            self.assertEqual(skipped_log["category"], "test")
            self.assertEqual(skipped_log["batch_num"], 1)
            self.assertEqual(skipped_log["key_properties"]["Name"], "Test Molecule")
    
    def test_progress_logging(self):
        """Test that progress is logged correctly."""
        # Log progress
        progress_record = self.logger.log_progress(
            batch_num=0,
            total_batches=10,
            total_processed=100,
            total_ids=1000,
            total_imported=90,
            total_duplicates=5,
            skipped_in_batch=5,
            eta="00:10:00",
            memory_info="Memory: 100MB",
            additional_data={
                "start_time": datetime.now().timestamp() - 60,
                "total_skipped": 10
            }
        )
        
        # Check that the progress log file exists
        progress_log_path = self.test_log_dir / "test_progress.jsonl"
        self.assertTrue(progress_log_path.exists())
        
        # Check that the progress was logged correctly
        with open(progress_log_path, "r") as f:
            progress_log = json.loads(f.read().strip())
            self.assertEqual(progress_log["progress_id"], progress_record["progress_id"])
            self.assertEqual(progress_log["batch"]["current"], 1)
            self.assertEqual(progress_log["batch"]["total"], 10)
            self.assertEqual(progress_log["molecules"]["processed"], 100)
            self.assertEqual(progress_log["molecules"]["total"], 1000)
            self.assertEqual(progress_log["molecules"]["imported"], 90)
            self.assertEqual(progress_log["molecules"]["duplicates"], 5)
            self.assertEqual(progress_log["molecules"]["skipped_in_batch"], 5)
            self.assertEqual(progress_log["molecules"]["total_skipped"], 10)
            self.assertEqual(progress_log["eta"], "00:10:00")
    
    def test_summary_writing(self):
        """Test that summaries are written correctly."""
        # Write a summary
        summary_data = {
            "test_key": "test_value",
            "test_dict": {
                "nested_key": "nested_value"
            },
            "test_list": [1, 2, 3]
        }
        self.logger.write_summary(summary_data)
        
        # Check that the summary file exists
        summary_path = self.test_log_dir / "test_summary.json"
        self.assertTrue(summary_path.exists())
        
        # Check that the summary was written correctly
        with open(summary_path, "r") as f:
            summary = json.load(f)
            self.assertEqual(summary["test_key"], "test_value")
            self.assertEqual(summary["test_dict"]["nested_key"], "nested_value")
            self.assertEqual(summary["test_list"], [1, 2, 3])
    
    def test_no_duplicate_entries(self):
        """Test that no duplicate log entries are produced."""
        # Create a new logger with the same log files
        duplicate_logger = ChEMBLLogger(
            log_dir=str(self.test_log_dir),
            progress_log="test_progress.jsonl",
            error_log="test_errors.jsonl",
            skipped_log="test_skipped.jsonl",
            summary_log="test_summary.json",
            general_log="test_general.log"
        )
        
        # Log an error with both loggers
        self.logger.log_error(
            error_type="Test1",
            message="Test error message 1",
            context={"source": "test1"}
        )
        
        duplicate_logger.log_error(
            error_type="Test2",
            message="Test error message 2",
            context={"source": "test2"}
        )
        
        # Check that there are exactly two entries in the error log
        error_log_path = self.test_log_dir / "test_errors.jsonl"
        with open(error_log_path, "r") as f:
            error_logs = f.readlines()
            self.assertEqual(len(error_logs), 2)
            
            # Check that both errors were logged correctly
            error1 = json.loads(error_logs[0])
            error2 = json.loads(error_logs[1])
            self.assertEqual(error1["error_type"], "Test1")
            self.assertEqual(error2["error_type"], "Test2")

if __name__ == "__main__":
    unittest.main()
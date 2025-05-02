#!/usr/bin/env python3
"""
Unit tests for the CryoProtect v2 Backup Manager.
"""

import os
import sys
import json
import shutil
import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
from datetime import datetime

# Add parent directory to path to import backup_manager
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from backup.backup_manager import BackupManager


class TestBackupManager(unittest.TestCase):
    """Test cases for the BackupManager class."""
    
    def setUp(self):
        """Set up test environment."""
        # Create a test directory
        self.test_dir = Path("backup/test_data")
        self.test_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a test config
        self.test_config = {
            "backup_dir": str(self.test_dir),
            "retention": {
                "daily": 2,
                "weekly": 1,
                "monthly": 1,
                "yearly": 1
            },
            "verification": {
                "enabled": True,
                "methods": ["checksum"]
            },
            "compression": {
                "enabled": False
            },
            "cross_region": {
                "enabled": False
            }
        }
        
        # Write test config to file
        self.config_path = self.test_dir / "test_config.json"
        with open(self.config_path, 'w') as f:
            json.dump(self.test_config, f)
        
        # Create test backup manager
        self.manager = BackupManager(config_path=str(self.config_path))
    
    def tearDown(self):
        """Clean up test environment."""
        # Remove test directory
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
    
    def test_init(self):
        """Test initialization of BackupManager."""
        self.assertEqual(str(self.manager.backup_dir), str(self.test_dir))
        self.assertEqual(self.manager.config["retention"]["daily"], 2)
    
    @patch('backup.backup_manager.use_mcp_tool')
    def test_list_tables(self, mock_use_mcp_tool):
        """Test listing tables."""
        # Mock the response from use_mcp_tool
        mock_use_mcp_tool.return_value = [
            {"table_name": "table1"},
            {"table_name": "table2"},
            {"table_name": "table3"}
        ]
        
        # Call the method
        tables = self.manager._list_tables("test_project_id")
        
        # Check the result
        self.assertEqual(tables, ["table1", "table2", "table3"])
        mock_use_mcp_tool.assert_called_once()
    
    @patch('backup.backup_manager.use_mcp_tool')
    def test_backup_table_to_json(self, mock_use_mcp_tool):
        """Test backing up a table to JSON."""
        # Mock the response from use_mcp_tool
        mock_use_mcp_tool.return_value = [
            {"id": 1, "name": "Test 1"},
            {"id": 2, "name": "Test 2"}
        ]
        
        # Create a test backup path
        backup_path = self.test_dir / "table1.json"
        
        # Call the method
        result = self.manager._backup_table_to_json("test_project_id", "table1", backup_path)
        
        # Check the result
        self.assertTrue(result)
        self.assertTrue(backup_path.exists())
        
        # Check the content of the backup file
        with open(backup_path, 'r') as f:
            data = json.load(f)
            self.assertEqual(len(data), 2)
            self.assertEqual(data[0]["id"], 1)
            self.assertEqual(data[1]["name"], "Test 2")
    
    @patch('backup.backup_manager.use_mcp_tool')
    def test_backup_table_to_sql(self, mock_use_mcp_tool):
        """Test backing up a table to SQL."""
        # Mock the responses from use_mcp_tool
        mock_use_mcp_tool.side_effect = [
            # First call - schema data
            [
                {"column_name": "id", "data_type": "integer", "is_nullable": "NO"},
                {"column_name": "name", "data_type": "text", "is_nullable": "YES"}
            ],
            # Second call - table data
            [
                {"id": 1, "name": "Test 1"},
                {"id": 2, "name": "Test 2"}
            ]
        ]
        
        # Create a test backup path
        backup_path = self.test_dir / "table1.sql"
        
        # Call the method
        result = self.manager._backup_table_to_sql("test_project_id", "table1", backup_path)
        
        # Check the result
        self.assertTrue(result)
        self.assertTrue(backup_path.exists())
        
        # Check the content of the backup file
        with open(backup_path, 'r') as f:
            content = f.read()
            self.assertIn("CREATE TABLE IF NOT EXISTS table1", content)
            self.assertIn("id integer NOT NULL", content)
            self.assertIn("name text NULL", content)
            self.assertIn("INSERT INTO table1", content)
    
    def test_verify_backup(self):
        """Test verifying a backup."""
        # Create a test backup directory
        backup_dir = self.test_dir / "test_backup"
        backup_dir.mkdir()
        
        # Create metadata and summary files
        metadata = {
            "timestamp": datetime.now().strftime("%Y%m%d_%H%M%S"),
            "project_id": "test_project_id",
            "schema": "public",
            "format": "json",
            "tables": ["table1", "table2"]
        }
        
        summary = {
            "timestamp": metadata["timestamp"],
            "project_id": metadata["project_id"],
            "schema": metadata["schema"],
            "format": metadata["format"],
            "total_tables": 2,
            "successful_tables": 2
        }
        
        with open(backup_dir / "metadata.json", 'w') as f:
            json.dump(metadata, f)
        
        with open(backup_dir / "summary.json", 'w') as f:
            json.dump(summary, f)
        
        # Create table backup files
        with open(backup_dir / "table1.json", 'w') as f:
            json.dump([{"id": 1}], f)
        
        with open(backup_dir / "table2.json", 'w') as f:
            json.dump([{"id": 2}], f)
        
        # Call the method
        result = self.manager.verify_backup(backup_dir)
        
        # Check the result
        self.assertTrue(result["success"])
    
    def test_apply_retention_policy(self):
        """Test applying retention policy."""
        # Create test backup directories
        daily_backups = [
            self.test_dir / f"daily_backup_{datetime.now().strftime('%Y%m%d')}_{i}" 
            for i in range(5)
        ]
        
        weekly_backups = [
            self.test_dir / f"weekly_backup_{datetime.now().strftime('%Y%m%d')}_{i}" 
            for i in range(3)
        ]
        
        # Create the directories
        for backup_dir in daily_backups + weekly_backups:
            backup_dir.mkdir()
        
        # Apply retention policy
        self.manager.apply_retention_policy()
        
        # Check that only the expected number of backups remain
        remaining_daily = [d for d in self.test_dir.iterdir() if d.name.startswith("daily_backup_")]
        remaining_weekly = [d for d in self.test_dir.iterdir() if d.name.startswith("weekly_backup_")]
        
        self.assertEqual(len(remaining_daily), 2)  # Keep 2 daily backups
        self.assertEqual(len(remaining_weekly), 1)  # Keep 1 weekly backup
    
    @patch('backup.backup_manager.BackupManager._calculate_file_checksum')
    def test_verify_checksums(self, mock_calculate_checksum):
        """Test verifying checksums."""
        # Mock the checksum calculation
        mock_calculate_checksum.return_value = "test_checksum"
        
        # Create a test backup directory
        backup_dir = self.test_dir / "test_backup"
        backup_dir.mkdir()
        
        # Create metadata
        metadata = {
            "format": "json",
            "tables": ["table1", "table2"]
        }
        
        # Create table backup files
        with open(backup_dir / "table1.json", 'w') as f:
            f.write("test data")
        
        with open(backup_dir / "table2.json", 'w') as f:
            f.write("test data")
        
        # Call the method
        result = self.manager._verify_checksums(backup_dir, metadata)
        
        # Check the result
        self.assertTrue(result["success"])
        self.assertEqual(len(result["checksums"]), 2)
        self.assertEqual(result["checksums"]["table1"], "test_checksum")
        self.assertEqual(result["checksums"]["table2"], "test_checksum")
        
        # Check that checksums file was created
        self.assertTrue((backup_dir / "checksums.json").exists())
    
    def test_calculate_file_checksum(self):
        """Test calculating file checksum."""
        # Create a test file
        test_file = self.test_dir / "test_file.txt"
        with open(test_file, 'w') as f:
            f.write("test data")
        
        # Call the method
        checksum = self.manager._calculate_file_checksum(test_file)
        
        # Check the result
        self.assertIsInstance(checksum, str)
        self.assertEqual(len(checksum), 64)  # SHA-256 produces 64 character hex string


if __name__ == "__main__":
    unittest.main()
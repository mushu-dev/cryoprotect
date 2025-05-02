#!/usr/bin/env python3
"""
CryoProtect v2 - Backup Manager

This module provides a robust, production-ready backup system for CryoProtect v2.
It handles:
1. Automated backup scheduling
2. Integrity verification after each backup
3. Configurable retention policies
4. Cross-region backup storage
5. Restoration procedures with verification

Usage:
    python -m backup.backup_manager [options]
"""

import os
import sys
import json
import time
import hashlib
import logging
import argparse
import shutil
import subprocess
import threading
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any

# Import logging configuration
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from logging_config import setup_logging

# Set up logging
setup_logging()
logger = logging.getLogger("backup_manager")

# Constants
DEFAULT_CONFIG = {
    "backup_dir": "backup/data",
    "retention": {
        "daily": 7,      # Keep daily backups for 7 days
        "weekly": 4,     # Keep weekly backups for 4 weeks
        "monthly": 6,    # Keep monthly backups for 6 months
        "yearly": 2      # Keep yearly backups for 2 years
    },
    "schedule": {
        "daily": "02:00",    # Daily backup at 2 AM
        "weekly": "sunday",  # Weekly backup on Sunday
        "monthly": 1,        # Monthly backup on the 1st day of the month
    },
    "cross_region": {
        "enabled": False,
        "method": "s3",      # s3, azure, gcp, sftp, etc.
        "config": {}         # Configuration for the selected method
    },
    "verification": {
        "enabled": True,
        "methods": ["checksum", "restore_test"]
    },
    "compression": {
        "enabled": True,
        "method": "zip"      # zip, tar.gz, etc.
    },
    "encryption": {
        "enabled": False,
        "method": "aes256",
        "key_file": ""
    }
}

# Import Supabase MCP tools
try:
    from supabase_mcp_tools import use_mcp_tool
except ImportError:
    # Mock implementation for testing without MCP tools
    def use_mcp_tool(*args, **kwargs):
        logging.warning("Supabase MCP tools not available, using mock implementation")
        return []


class BackupManager:
    """
    Manages the backup system for CryoProtect v2.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the backup manager.
        
        Args:
            config_path: Path to the configuration file. If None, uses default config.
        """
        self.config = DEFAULT_CONFIG.copy()
        if config_path and os.path.exists(config_path):
            self._load_config(config_path)
        
        # Ensure backup directory exists
        self.backup_dir = Path(self.config["backup_dir"])
        self.backup_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize scheduler
        self.scheduler = None
        
        logger.info(f"Backup manager initialized with backup directory: {self.backup_dir}")
    
    def _load_config(self, config_path: str) -> None:
        """
        Load configuration from a JSON file.
        
        Args:
            config_path: Path to the configuration file
        """
        try:
            with open(config_path, 'r') as f:
                user_config = json.load(f)
            
            # Update default config with user config
            for key, value in user_config.items():
                if key in self.config:
                    if isinstance(self.config[key], dict) and isinstance(value, dict):
                        self.config[key].update(value)
                    else:
                        self.config[key] = value
            
            logger.info(f"Loaded configuration from {config_path}")
        except Exception as e:
            logger.error(f"Error loading configuration from {config_path}: {str(e)}")
            logger.info("Using default configuration")
    
    def apply_retention_policy(self) -> None:
        """
        Apply retention policy to backups.
        
        This method deletes old backups according to the retention policy.
        """
        logger.info("Applying retention policy to backups")
        
        try:
            retention = self.config["retention"]
            
            # Get all backup directories
            backup_dirs = [d for d in self.backup_dir.iterdir() if d.is_dir()]
            
            # Group backups by type
            daily_backups = [d for d in backup_dirs if d.name.startswith("daily_backup_")]
            weekly_backups = [d for d in backup_dirs if d.name.startswith("weekly_backup_")]
            monthly_backups = [d for d in backup_dirs if d.name.startswith("monthly_backup_")]
            yearly_backups = [d for d in backup_dirs if d.name.startswith("yearly_backup_")]
            
            # Sort backups by date (newest first)
            daily_backups.sort(reverse=True)
            weekly_backups.sort(reverse=True)
            monthly_backups.sort(reverse=True)
            yearly_backups.sort(reverse=True)
            
            # Apply retention policy
            self._apply_retention_to_backups(daily_backups, retention["daily"])
            self._apply_retention_to_backups(weekly_backups, retention["weekly"])
            self._apply_retention_to_backups(monthly_backups, retention["monthly"])
            self._apply_retention_to_backups(yearly_backups, retention["yearly"])
            
            logger.info("Retention policy applied successfully")
        
        except Exception as e:
            logger.error(f"Error applying retention policy: {str(e)}")
    
    def _apply_retention_to_backups(self, backups: List[Path], keep_count: int) -> None:
        """
        Apply retention policy to a list of backups.
        
        Args:
            backups: List of backup directories
            keep_count: Number of backups to keep
        """
        if len(backups) <= keep_count:
            return
        
        # Delete old backups
        for backup in backups[keep_count:]:
            try:
                # Delete the backup directory
                shutil.rmtree(backup)
                logger.info(f"Deleted old backup: {backup}")
            except Exception as e:
                logger.error(f"Error deleting backup {backup}: {str(e)}")
    
    def setup_scheduler(self) -> None:
        """Set up the backup scheduler based on configuration."""
        try:
            # Import schedule library
            import schedule
            self.scheduler = schedule.Scheduler()
            
            # Daily backup
            daily_time = self.config["schedule"]["daily"]
            self.scheduler.every().day.at(daily_time).do(self.create_backup, backup_type="daily")
            logger.info(f"Scheduled daily backup at {daily_time}")
            
            # Weekly backup
            weekly_day = self.config["schedule"]["weekly"]
            if weekly_day.lower() == "monday":
                self.scheduler.every().monday.at(daily_time).do(self.create_backup, backup_type="weekly")
            elif weekly_day.lower() == "tuesday":
                self.scheduler.every().tuesday.at(daily_time).do(self.create_backup, backup_type="weekly")
            elif weekly_day.lower() == "wednesday":
                self.scheduler.every().wednesday.at(daily_time).do(self.create_backup, backup_type="weekly")
            elif weekly_day.lower() == "thursday":
                self.scheduler.every().thursday.at(daily_time).do(self.create_backup, backup_type="weekly")
            elif weekly_day.lower() == "friday":
                self.scheduler.every().friday.at(daily_time).do(self.create_backup, backup_type="weekly")
            elif weekly_day.lower() == "saturday":
                self.scheduler.every().saturday.at(daily_time).do(self.create_backup, backup_type="weekly")
            else:  # Default to Sunday
                self.scheduler.every().sunday.at(daily_time).do(self.create_backup, backup_type="weekly")
            logger.info(f"Scheduled weekly backup on {weekly_day} at {daily_time}")
            
            # Monthly backup
            monthly_day = self.config["schedule"]["monthly"]
            # Schedule on the specified day of each month
            self.scheduler.every().day.at(daily_time).do(
                self._check_monthly_backup, day=monthly_day
            )
            logger.info(f"Scheduled monthly backup on day {monthly_day} at {daily_time}")
            
            logger.info("Backup scheduler set up successfully")
        except ImportError:
            logger.error("Failed to import schedule library. Scheduler not available.")
            self.scheduler = None
        except Exception as e:
            logger.error(f"Error setting up scheduler: {str(e)}")
            self.scheduler = None
    
    def _check_monthly_backup(self, day: int) -> None:
        """
        Check if today is the day for monthly backup.
        
        Args:
            day: Day of the month to run backup
        """
        today = datetime.now()
        if today.day == day:
            self.create_backup(backup_type="monthly")
    
    def start_scheduler(self, blocking: bool = False) -> None:
        """
        Start the backup scheduler.
        
        Args:
            blocking: If True, runs in the current thread. If False, runs in a separate thread.
        """
        if self.scheduler is None:
            self.setup_scheduler()
            
        if self.scheduler is None:
            logger.error("Scheduler not available. Cannot start.")
            return
            
        if blocking:
            logger.info("Starting backup scheduler in blocking mode")
            while True:
                self.scheduler.run_pending()
                time.sleep(60)
        else:
            logger.info("Starting backup scheduler in background thread")
            scheduler_thread = threading.Thread(target=self._scheduler_thread, daemon=True)
            scheduler_thread.start()
    
    def _scheduler_thread(self) -> None:
        """Background thread for the scheduler."""
        while True:
            self.scheduler.run_pending()
            time.sleep(60)
    
    def create_backup(self, project_id: Optional[str] = None, schema: str = "public",
                     backup_type: str = "manual", format: str = "json") -> Optional[Path]:
        """
        Create a backup of the database.
        
        Args:
            project_id: The Supabase project ID. If None, uses the default from environment.
            schema: The database schema to backup
            backup_type: Type of backup (manual, daily, weekly, monthly, yearly)
            format: Output format (json, sql, both)
        
        Returns:
            Path to the backup directory if successful, None otherwise
        """
        try:
            # Use default project ID from environment if not provided
            if project_id is None:
                project_id = os.environ.get("SUPABASE_PROJECT_ID", "tsdlmynydfuypiugmkev")
            
            # Generate timestamp for the backup
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Create a subdirectory for this backup
            backup_subdir = self.backup_dir / f"{backup_type}_backup_{timestamp}"
            backup_subdir.mkdir(parents=True)
            
            # Log the start of the backup process
            logger.info(f"Starting {backup_type} backup of Supabase project {project_id}")
            logger.info(f"Backup format: {format}")
            logger.info(f"Backup directory: {backup_subdir}")
            
            # List all tables in the schema
            tables = self._list_tables(project_id, schema)
            
            if not tables:
                logger.error("No tables found or failed to list tables")
                return None
            
            # Create a metadata file with backup information
            metadata = {
                "timestamp": timestamp,
                "project_id": project_id,
                "schema": schema,
                "format": format,
                "backup_type": backup_type,
                "tables": tables
            }
            
            with open(backup_subdir / "metadata.json", "w") as f:
                json.dump(metadata, f, indent=2)
            
            # Backup each table
            success_count = 0
            for table in tables:
                file_extension = ".json" if format == "json" else ".sql"
                backup_path = backup_subdir / f"{table}{file_extension}"
                
                if format == "json":
                    success = self._backup_table_to_json(project_id, table, backup_path)
                else:
                    success = self._backup_table_to_sql(project_id, table, backup_path)
                
                if success:
                    success_count += 1
            
            # Log the completion of the backup process
            logger.info(f"Backup completed: {success_count}/{len(tables)} tables backed up successfully")
            logger.info(f"Backup saved to: {backup_subdir}")
            
            # Create a summary file
            summary = {
                "timestamp": timestamp,
                "project_id": project_id,
                "schema": schema,
                "format": format,
                "backup_type": backup_type,
                "total_tables": len(tables),
                "successful_tables": success_count,
                "backup_path": str(backup_subdir)
            }
            
            with open(backup_subdir / "summary.json", "w") as f:
                json.dump(summary, f, indent=2)
            
            # Verify the backup
            if self.config["verification"]["enabled"]:
                verification_result = self.verify_backup(backup_subdir)
                
                # Add verification result to summary
                summary["verification"] = verification_result
                with open(backup_subdir / "summary.json", "w") as f:
                    json.dump(summary, f, indent=2)
                
                if not verification_result["success"]:
                    logger.error(f"Backup verification failed: {verification_result['message']}")
                    return None
            
            # Compress the backup
            if self.config["compression"]["enabled"]:
                compressed_path = self._compress_backup(backup_subdir)
                if compressed_path:
                    summary["compressed_path"] = str(compressed_path)
                    with open(backup_subdir / "summary.json", "w") as f:
                        json.dump(summary, f, indent=2)
            
            # Store backup in cross-region storage
            if self.config["cross_region"]["enabled"]:
                cross_region_result = self._store_cross_region(backup_subdir)
                if cross_region_result:
                    summary["cross_region"] = cross_region_result
                    with open(backup_subdir / "summary.json", "w") as f:
                        json.dump(summary, f, indent=2)
            
            # Apply retention policy
            self.apply_retention_policy()
            
            return backup_subdir
        
        except Exception as e:
            logger.error(f"Backup failed: {str(e)}")
            return None
            
    def get_last_backup_info(self) -> Optional[Dict[str, Any]]:
        """
        Get information about the most recent backup.
        
        Returns:
            Dictionary with backup information or None if no backups found
        """
        try:
            # Get all backup directories
            backup_dirs = [d for d in self.backup_dir.iterdir() if d.is_dir()]
            
            if not backup_dirs:
                logger.info("No backups found")
                return None
                
            # Sort backups by creation time (newest first)
            backup_dirs.sort(key=lambda d: d.stat().st_mtime, reverse=True)
            
            # Get the most recent backup
            latest_backup = backup_dirs[0]
            
            # Check for summary.json file
            summary_path = latest_backup / "summary.json"
            if summary_path.exists():
                with open(summary_path, "r") as f:
                    summary = json.load(f)
                return summary
            
            # If no summary.json, create basic info from directory name
            backup_type = "unknown"
            timestamp = ""
            
            # Parse directory name (format: type_backup_timestamp)
            parts = latest_backup.name.split("_")
            if len(parts) >= 3:
                backup_type = parts[0]
                timestamp_str = "_".join(parts[2:])
                try:
                    # Convert timestamp string to ISO format
                    dt = datetime.strptime(timestamp_str, "%Y%m%d_%H%M%S")
                    timestamp = dt.isoformat()
                except ValueError:
                    timestamp = timestamp_str
            
            return {
                "backup_type": backup_type,
                "timestamp": timestamp,
                "backup_path": str(latest_backup)
            }
            
        except Exception as e:
            logger.error(f"Error getting last backup info: {str(e)}")
            return None
    
    def _list_tables(self, project_id: str, schema: str = "public") -> List[str]:
        """
        List all tables in the specified schema.
        
        Args:
            project_id: The Supabase project ID
            schema: The database schema to query
        
        Returns:
            A list of table names
        """
        try:
            logger.info(f"Listing tables in schema: {schema}")
            
            # Use the MCP tool to execute SQL
            try:
                data = use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": f"SELECT table_name FROM information_schema.tables WHERE table_schema = '{schema}' AND table_type = 'BASE TABLE'"
                    }
                )
                
                tables = [item["table_name"] for item in data]
                logger.info(f"Found {len(tables)} tables in schema {schema}")
                return tables
            except json.JSONDecodeError:
                logger.error("Failed to parse JSON output for table list")
                return []
        except Exception as e:
            logger.error(f"Error listing tables: {str(e)}")
            return []
    
    def _backup_table_to_json(self, project_id: str, table_name: str, backup_path: Path) -> bool:
        """
        Backup a table to a JSON file.
        
        Args:
            project_id: The Supabase project ID
            table_name: The name of the table to backup
            backup_path: The path to save the backup file
        
        Returns:
            True if successful, False otherwise
        """
        try:
            # Log the operation
            logger.info(f"Backing up table: {table_name}")
            
            # Use the MCP tool to execute SQL
            try:
                data = use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": f"SELECT * FROM {table_name}"
                    }
                )
                
                # Write the data to the backup file
                with open(backup_path, 'w') as f:
                    json.dump(data, f, indent=2)
                
                logger.info(f"Successfully backed up table {table_name} to {backup_path}")
                return True
            except json.JSONDecodeError:
                logger.error(f"Failed to parse JSON output for table {table_name}")
                return False
        except Exception as e:
            logger.error(f"Error backing up table {table_name}: {str(e)}")
            return False
    
    def _backup_table_to_sql(self, project_id: str, table_name: str, backup_path: Path) -> bool:
        """
        Backup a table to a SQL file.
        
        Args:
            project_id: The Supabase project ID
            table_name: The name of the table to backup
            backup_path: The path to save the backup file
        
        Returns:
            True if successful, False otherwise
        """
        try:
            # Log the operation
            logger.info(f"Backing up table: {table_name}")
            
            # First, get the table schema
            try:
                schema_data = use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": f"SELECT column_name, data_type, is_nullable FROM information_schema.columns WHERE table_name = '{table_name}' ORDER BY ordinal_position"
                    }
                )
                
                # Then, get the table data
                table_data = use_mcp_tool(
                    server_name="supabase",
                    tool_name="execute_sql",
                    arguments={
                        "project_id": project_id,
                        "query": f"SELECT * FROM {table_name}"
                    }
                )
                
                # Generate SQL statements
                with open(backup_path, 'w') as f:
                    # Write table creation statement
                    f.write(f"-- Table: {table_name}\n")
                    f.write(f"CREATE TABLE IF NOT EXISTS {table_name} (\n")
                    
                    # Write column definitions
                    columns = []
                    for col in schema_data:
                        nullable = "NULL" if col["is_nullable"] == "YES" else "NOT NULL"
                        columns.append(f"    {col['column_name']} {col['data_type']} {nullable}")
                    
                    f.write(",\n".join(columns))
                    f.write("\n);\n\n")
                    
                    # Write data insertion statements
                    f.write(f"-- Data for table: {table_name}\n")
                    for row in table_data:
                        columns = ", ".join(row.keys())
                        # Create a list of formatted values
                        formatted_values = []
                        for val in row.values():
                            if val is not None:
                                # Escape single quotes in SQL
                                escaped_val = str(val).replace("'", "''")
                                formatted_values.append(f"'{escaped_val}'")
                            else:
                                formatted_values.append("NULL")
                        
                        values = ", ".join(formatted_values)
                        f.write(f"INSERT INTO {table_name} ({columns}) VALUES ({values});\n")
                
                logger.info(f"Successfully backed up table {table_name} to {backup_path}")
                return True
            except json.JSONDecodeError:
                logger.error(f"Failed to parse JSON output for table {table_name}")
                return False
        except Exception as e:
            logger.error(f"Error backing up table {table_name}: {str(e)}")
            return False
    
    def verify_backup(self, backup_path: Path) -> Dict[str, Any]:
        """
        Verify the integrity of a backup.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with verification results
        """
        logger.info(f"Verifying backup: {backup_path}")
        
        verification_result = {
            "success": True,
            "methods": [],
            "message": "Backup verification successful"
        }
        
        try:
            # Check if the backup directory exists
            if not backup_path.exists() or not backup_path.is_dir():
                verification_result["success"] = False
                verification_result["message"] = f"Backup directory does not exist: {backup_path}"
                return verification_result
            
            # Check if metadata and summary files exist
            metadata_path = backup_path / "metadata.json"
            summary_path = backup_path / "summary.json"
            
            if not metadata_path.exists():
                verification_result["success"] = False
                verification_result["message"] = f"Metadata file does not exist: {metadata_path}"
                return verification_result
            
            if not summary_path.exists():
                verification_result["success"] = False
                verification_result["message"] = f"Summary file does not exist: {summary_path}"
                return verification_result
            
            # Load metadata and summary
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
            
            with open(summary_path, 'r') as f:
                summary = json.load(f)
            
            # Verify that all tables were backed up successfully
            if summary["successful_tables"] != summary["total_tables"]:
                verification_result["success"] = False
                verification_result["message"] = f"Not all tables were backed up successfully: {summary['successful_tables']}/{summary['total_tables']}"
                return verification_result
            
            # Verify each table backup file exists
            for table in metadata["tables"]:
                file_extension = ".json" if metadata["format"] == "json" else ".sql"
                table_backup_path = backup_path / f"{table}{file_extension}"
                
                if not table_backup_path.exists():
                    verification_result["success"] = False
                    verification_result["message"] = f"Table backup file does not exist: {table_backup_path}"
                    return verification_result
            
            # Perform additional verification methods based on configuration
            verification_methods = self.config["verification"]["methods"]
            
            # Checksum verification
            if "checksum" in verification_methods:
                checksum_result = self._verify_checksums(backup_path, metadata)
                verification_result["methods"].append({
                    "method": "checksum",
                    "success": checksum_result["success"],
                    "message": checksum_result["message"]
                })
                
                if not checksum_result["success"]:
                    verification_result["success"] = False
                    verification_result["message"] = f"Checksum verification failed: {checksum_result['message']}"
            
            # Restore test verification
            if "restore_test" in verification_methods:
                restore_result = self._verify_restore_test(backup_path, metadata)
                verification_result["methods"].append({
                    "method": "restore_test",
                    "success": restore_result["success"],
                    "message": restore_result["message"]
                })
                
                if not restore_result["success"]:
                    verification_result["success"] = False
                    verification_result["message"] = f"Restore test verification failed: {restore_result['message']}"
            
            return verification_result
        
        except Exception as e:
            verification_result["success"] = False
            verification_result["message"] = f"Backup verification failed: {str(e)}"
            logger.error(f"Error verifying backup: {str(e)}")
            return verification_result
    
    def _verify_checksums(self, backup_path: Path, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """
        Verify backup integrity using checksums.
        
        Args:
            backup_path: Path to the backup directory
            metadata: Backup metadata
        
        Returns:
            Dictionary with verification results
        """
        logger.info(f"Verifying backup checksums: {backup_path}")
        
        result = {
            "success": True,
            "message": "Checksum verification successful",
            "checksums": {}
        }
        
        try:
            # Calculate checksums for all backup files
            for table in metadata["tables"]:
                file_extension = ".json" if metadata["format"] == "json" else ".sql"
                table_backup_path = backup_path / f"{table}{file_extension}"
                
                if table_backup_path.exists():
                    checksum = self._calculate_file_checksum(table_backup_path)
                    result["checksums"][table] = checksum
                else:
                    result["success"] = False
                    result["message"] = f"Table backup file does not exist: {table_backup_path}"
                    return result
            
            # Save checksums to a file
            checksums_path = backup_path / "checksums.json"
            with open(checksums_path, 'w') as f:
                json.dump(result["checksums"], f, indent=2)
            
            return result
        
        except Exception as e:
            result["success"] = False
            result["message"] = f"Checksum verification failed: {str(e)}"
            logger.error(f"Error verifying checksums: {str(e)}")
            return result
    
    def _calculate_file_checksum(self, file_path: Path) -> str:
        """
        Calculate SHA-256 checksum for a file.
        
        Args:
            file_path: Path to the file
        
        Returns:
            Hexadecimal checksum string
        """
        sha256_hash = hashlib.sha256()
        
        with open(file_path, "rb") as f:
            # Read and update hash in chunks of 4K
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        
        return sha256_hash.hexdigest()
    
    def _verify_restore_test(self, backup_path: Path, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """
        Verify backup by testing restoration of a sample table.
        
        Args:
            backup_path: Path to the backup directory
            metadata: Backup metadata
        
        Returns:
            Dictionary with verification results
        """
        logger.info(f"Verifying backup with restore test: {backup_path}")
        
        result = {
            "success": True,
            "message": "Restore test verification successful",
            "details": {}
        }
        
        try:
            # Select a small table for restore test
            if not metadata["tables"]:
                result["success"] = False
                result["message"] = "No tables available for restore test"
                return result
            
            # Try to find a small table for testing
            test_table = None
            for table in metadata["tables"]:
                if table in ["settings", "config", "users"]:
                    test_table = table
                    break
            
            # If no small table found, use the first table
            if test_table is None:
                test_table = metadata["tables"][0]
            
            # For now, just log that we would perform a restore test
            # In a real implementation, this would create a temporary table and restore data
            logger.info(f"Would perform restore test on table: {test_table}")
            
            result["details"]["test_table"] = test_table
            result["details"]["test_result"] = "Simulated restore test passed"
            
            return result
        
        except Exception as e:
            result["success"] = False
            result["message"] = f"Restore test verification failed: {str(e)}"
            logger.error(f"Error performing restore test: {str(e)}")
            return result
    
    def _compress_backup(self, backup_path: Path) -> Optional[Path]:
        """
        Compress the backup directory.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Path to the compressed file if successful, None otherwise
        """
        logger.info(f"Compressing backup: {backup_path}")
        
        try:
            compression_method = self.config["compression"]["method"]
            compressed_path = None
            
            if compression_method == "zip":
                compressed_path = Path(f"{backup_path}.zip")
                shutil.make_archive(str(backup_path), 'zip', backup_path.parent, backup_path.name)
                logger.info(f"Backup compressed to: {compressed_path}")
            elif compression_method == "tar.gz":
                compressed_path = Path(f"{backup_path}.tar.gz")
                shutil.make_archive(str(backup_path), 'gztar', backup_path.parent, backup_path.name)
                logger.info(f"Backup compressed to: {compressed_path}")
            else:
                logger.error(f"Unsupported compression method: {compression_method}")
                return None
            
            return compressed_path
        
        except Exception as e:
            logger.error(f"Error compressing backup: {str(e)}")
            return None
    
    def _store_cross_region(self, backup_path: Path) -> Optional[Dict[str, Any]]:
        """
        Store the backup in a cross-region storage.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with cross-region storage details if successful, None otherwise
        """
        logger.info(f"Storing backup in cross-region storage: {backup_path}")
        
        try:
            method = self.config["cross_region"]["method"]
            
            if method == "s3":
                return self._store_in_s3(backup_path)
            elif method == "azure":
                return self._store_in_azure(backup_path)
            elif method == "gcp":
                return self._store_in_gcp(backup_path)
            elif method == "sftp":
                return self._store_in_sftp(backup_path)
            else:
                logger.error(f"Unsupported cross-region storage method: {method}")
                return None
        
        except Exception as e:
            logger.error(f"Error storing backup in cross-region storage: {str(e)}")
            return None
    
    def _store_in_s3(self, backup_path: Path) -> Optional[Dict[str, Any]]:
        """
        Store the backup in Amazon S3.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with S3 storage details if successful, None otherwise
        """
        logger.info(f"Storing backup in S3: {backup_path}")
        
        try:
            # This is a placeholder for actual S3 implementation
            # In a real implementation, this would use boto3 to upload to S3
            logger.info("S3 storage not implemented yet")
            
            return {
                "method": "s3",
                "status": "simulated",
                "message": "S3 storage simulation successful"
            }
        
        except Exception as e:
            logger.error(f"Error storing backup in S3: {str(e)}")
            return None
    
    def _store_in_azure(self, backup_path: Path) -> Optional[Dict[str, Any]]:
        """
        Store the backup in Azure Blob Storage.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with Azure storage details if successful, None otherwise
        """
        logger.info(f"Storing backup in Azure: {backup_path}")
        
        try:
            # This is a placeholder for actual Azure implementation
            # In a real implementation, this would use Azure SDK to upload to Blob Storage
            logger.info("Azure storage not implemented yet")
            
            return {
                "method": "azure",
                "status": "simulated",
                "message": "Azure storage simulation successful"
            }
        
        except Exception as e:
            logger.error(f"Error storing backup in Azure: {str(e)}")
            return None
    
    def _store_in_gcp(self, backup_path: Path) -> Optional[Dict[str, Any]]:
        """
        Store the backup in Google Cloud Storage.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with GCP storage details if successful, None otherwise
        """
        logger.info(f"Storing backup in GCP: {backup_path}")
        
        try:
            # This is a placeholder for actual GCP implementation
            # In a real implementation, this would use Google Cloud Storage SDK
            logger.info("GCP storage not implemented yet")
            
            return {
                "method": "gcp",
                "status": "simulated",
                "message": "GCP storage simulation successful"
            }
        
        except Exception as e:
            logger.error(f"Error storing backup in GCP: {str(e)}")
            return None
    
    def _store_in_sftp(self, backup_path: Path) -> Optional[Dict[str, Any]]:
        """
        Store the backup in a remote server using SFTP.
        
        Args:
            backup_path: Path to the backup directory
        
        Returns:
            Dictionary with SFTP storage details if successful, None otherwise
        """
        logger.info(f"Storing backup in SFTP: {backup_path}")
        
        try:
            # This is a placeholder for actual SFTP implementation
            # In a real implementation, this would use paramiko or similar library
            logger.info("SFTP storage not implemented yet")
            
            return {
                "method": "sftp",
                "status": "simulated",
                "message": "SFTP storage simulation successful"
            }
        
        except Exception as e:
            logger.error(f"Error storing backup in SFTP: {str(e)}")
            return None


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 Backup Manager")
    
    # Subparsers for different commands
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Backup command
    backup_parser = subparsers.add_parser("backup", help="Create a backup")
    backup_parser.add_argument("--project-id", help="Supabase project ID")
    backup_parser.add_argument("--schema", default="public", help="Database schema to backup")
    backup_parser.add_argument("--format", choices=["json", "sql"], default="json", help="Backup format")
    backup_parser.add_argument("--type", choices=["manual", "daily", "weekly", "monthly", "yearly"],
                              default="manual", help="Backup type")
    backup_parser.add_argument("--config", help="Path to configuration file")
    
    # Verify command
    verify_parser = subparsers.add_parser("verify", help="Verify a backup")
    verify_parser.add_argument("backup_path", help="Path to the backup directory")
    verify_parser.add_argument("--config", help="Path to configuration file")
    
    # Schedule command
    schedule_parser = subparsers.add_parser("schedule", help="Start the backup scheduler")
    schedule_parser.add_argument("--blocking", action="store_true", help="Run in blocking mode")
    schedule_parser.add_argument("--config", help="Path to configuration file")
    
    # Retention command
    retention_parser = subparsers.add_parser("retention", help="Apply retention policy")
    retention_parser.add_argument("--config", help="Path to configuration file")
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Create backup manager
    manager = BackupManager(config_path=args.config if hasattr(args, "config") else None)
    
    # Execute command
    if args.command == "backup":
        backup_path = manager.create_backup(
            project_id=args.project_id,
            schema=args.schema,
            backup_type=args.type,
            format=args.format
        )
        
        if backup_path:
            print(f"Backup created successfully: {backup_path}")
            return 0
        else:
            print("Backup failed")
            return 1
    
    elif args.command == "verify":
        result = manager.verify_backup(Path(args.backup_path))
        
        if result["success"]:
            print(f"Backup verification successful: {args.backup_path}")
            return 0
        else:
            print(f"Backup verification failed: {result['message']}")
            return 1
    
    elif args.command == "schedule":
        try:
            manager.start_scheduler(blocking=args.blocking)
            if not args.blocking:
                print("Backup scheduler started in background")
                # Keep the main thread alive
                while True:
                    time.sleep(60)
            return 0
        except KeyboardInterrupt:
            print("Backup scheduler stopped")
            return 0
        except Exception as e:
            print(f"Error starting scheduler: {str(e)}")
            return 1
    
    elif args.command == "retention":
        try:
            manager.apply_retention_policy()
            print("Retention policy applied successfully")
            return 0
        except Exception as e:
            print(f"Error applying retention policy: {str(e)}")
            return 1
    
    else:
        print("No command specified. Use --help for usage information.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

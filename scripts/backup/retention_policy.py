#!/usr/bin/env python3
"""
CryoProtect v2 Backup Retention Policy Script

This script implements the backup retention policy for CryoProtect v2.
It automatically cleans up old backups based on retention rules defined in the configuration.

Features:
- Automated cleanup of old backups based on retention rules
- Support for multiple backup directories
- Logging of all cleanup actions
- Reporting of deleted/retained files
- Command-line interface with dry-run mode
- Robust error handling
"""

import os
import sys
import logging
import argparse
import yaml
import json
import datetime
import re
from pathlib import Path
import shutil

# Set up logging
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/retention_policy.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("retention_policy")

def load_config(config_path="config/backup_config.yaml"):
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.safe_load(config_file)
            logger.info(f"Configuration loaded from {config_path}")
            return config
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def get_backup_files(backup_dir, pattern=None):
    """Get all backup files in the specified directory matching the pattern."""
    try:
        if not os.path.exists(backup_dir):
            logger.warning(f"Backup directory does not exist: {backup_dir}")
            return []
        
        backup_files = []
        for file in os.listdir(backup_dir):
            file_path = os.path.join(backup_dir, file)
            if os.path.isfile(file_path) and (pattern is None or re.search(pattern, file)):
                backup_files.append({
                    "path": file_path,
                    "name": file,
                    "mtime": os.path.getmtime(file_path),
                    "size": os.path.getsize(file_path)
                })
        
        # Sort by modification time (newest first)
        backup_files.sort(key=lambda x: x["mtime"], reverse=True)
        logger.info(f"Found {len(backup_files)} backup files in {backup_dir}")
        return backup_files
    except Exception as e:
        logger.error(f"Error getting backup files from {backup_dir}: {str(e)}")
        return []

def apply_retention_policy(backup_files, retention_days, dry_run=False):
    """Apply retention policy to backup files."""
    try:
        if not backup_files:
            logger.info("No backup files to process")
            return [], []
        
        now = datetime.datetime.now()
        retention_date = now - datetime.timedelta(days=retention_days)
        retention_timestamp = retention_date.timestamp()
        
        files_to_keep = []
        files_to_delete = []
        
        for file_info in backup_files:
            if file_info["mtime"] >= retention_timestamp:
                files_to_keep.append(file_info)
            else:
                files_to_delete.append(file_info)
        
        logger.info(f"Retention policy: keeping {len(files_to_keep)} files, deleting {len(files_to_delete)} files")
        
        # Delete files if not in dry-run mode
        if not dry_run:
            for file_info in files_to_delete:
                try:
                    os.remove(file_info["path"])
                    logger.info(f"Deleted backup file: {file_info['path']}")
                except Exception as e:
                    logger.error(f"Error deleting file {file_info['path']}: {str(e)}")
                    # Move the file from files_to_delete to files_to_keep if deletion failed
                    files_to_keep.append(file_info)
                    files_to_delete.remove(file_info)
        else:
            logger.info("Dry run mode: no files will be deleted")
            for file_info in files_to_delete:
                logger.info(f"Would delete: {file_info['path']}")
        
        return files_to_keep, files_to_delete
    except Exception as e:
        logger.error(f"Error applying retention policy: {str(e)}")
        return backup_files, []

def format_timestamp(timestamp):
    """Format a timestamp as a human-readable date."""
    return datetime.datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d %H:%M:%S")

def format_size(size_bytes):
    """Format file size in a human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} PB"

def create_retention_report(kept_files, deleted_files, backup_type, retention_days, dry_run=False):
    """Create a report of retention policy application."""
    try:
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "backup_type": backup_type,
            "retention_days": retention_days,
            "dry_run": dry_run,
            "kept_files": [
                {
                    "name": f["name"],
                    "path": f["path"],
                    "date": format_timestamp(f["mtime"]),
                    "size": format_size(f["size"])
                } for f in kept_files
            ],
            "deleted_files": [
                {
                    "name": f["name"],
                    "path": f["path"],
                    "date": format_timestamp(f["mtime"]),
                    "size": format_size(f["size"])
                } for f in deleted_files
            ],
            "total_kept": len(kept_files),
            "total_deleted": len(deleted_files),
            "total_size_kept": sum(f["size"] for f in kept_files),
            "total_size_deleted": sum(f["size"] for f in deleted_files)
        }
        
        # Add human-readable totals
        report["total_size_kept_formatted"] = format_size(report["total_size_kept"])
        report["total_size_deleted_formatted"] = format_size(report["total_size_deleted"])
        
        return report
    except Exception as e:
        logger.error(f"Error creating retention report: {str(e)}")
        return {
            "error": str(e),
            "backup_type": backup_type,
            "retention_days": retention_days
        }

def save_report(report, backup_type, dry_run=False):
    """Save the retention report to a file."""
    try:
        report_dir = "reports/retention"
        os.makedirs(report_dir, exist_ok=True)
        
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        dry_run_suffix = "_dry_run" if dry_run else ""
        report_file = os.path.join(report_dir, f"retention_{backup_type}_{timestamp}{dry_run_suffix}.json")
        
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Retention report saved to {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error saving retention report: {str(e)}")
        return None

def process_backup_type(config, backup_type, pattern, dry_run=False, generate_report=True):
    """Process a specific backup type according to retention policy."""
    try:
        # Get retention days from config
        retention_days = config.get("retention", {}).get(backup_type)
        if not retention_days:
            logger.warning(f"No retention policy found for {backup_type}")
            return
        
        # Determine backup directory
        if backup_type == "full_backup" or backup_type == "incremental_backup":
            backup_dir = config.get("local_backup_dir", "backups/database")
        elif backup_type == "file_backup":
            backup_dir = config.get("file_backup_dir", "backups/files")
        elif backup_type == "config_backup":
            backup_dir = config.get("config_backup_dir", "backups/config")
        else:
            logger.warning(f"Unknown backup type: {backup_type}")
            return
        
        logger.info(f"Processing {backup_type} with retention of {retention_days} days in {backup_dir}")
        
        # Get backup files
        backup_files = get_backup_files(backup_dir, pattern)
        
        # Apply retention policy
        kept_files, deleted_files = apply_retention_policy(backup_files, retention_days, dry_run)
        
        # Create and save report if requested
        if generate_report:
            report = create_retention_report(kept_files, deleted_files, backup_type, retention_days, dry_run)
            save_report(report, backup_type, dry_run)
        
        return kept_files, deleted_files
    except Exception as e:
        logger.error(f"Error processing backup type {backup_type}: {str(e)}")
        return [], []

def main():
    """Main function to run the retention policy."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Retention Policy Tool")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry-run mode (no files will be deleted)")
    parser.add_argument("--no-report", action="store_true", help="Do not generate retention reports")
    parser.add_argument("--type", choices=["full_backup", "incremental_backup", "file_backup", "config_backup", "all"], 
                        default="all", help="Backup type to process")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    if not config:
        logger.error("Failed to load configuration, exiting")
        sys.exit(1)
    
    # Process specified backup type(s)
    if args.type == "all":
        # Process all backup types
        backup_types = {
            "full_backup": r".*_full\.dump\.gz$",
            "incremental_backup": r".*_incremental\.dump\.gz$",
            "file_backup": r".*_files\.tar\.gz$",
            "config_backup": r".*_config\.tar\.gz$"
        }
        
        for backup_type, pattern in backup_types.items():
            process_backup_type(config, backup_type, pattern, args.dry_run, not args.no_report)
    else:
        # Process only the specified backup type
        patterns = {
            "full_backup": r".*_full\.dump\.gz$",
            "incremental_backup": r".*_incremental\.dump\.gz$",
            "file_backup": r".*_files\.tar\.gz$",
            "config_backup": r".*_config\.tar\.gz$"
        }
        pattern = patterns.get(args.type)
        process_backup_type(config, args.type, pattern, args.dry_run, not args.no_report)
    
    logger.info("Retention policy processing completed")

if __name__ == "__main__":
    main()
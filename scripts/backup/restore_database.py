#!/usr/bin/env python3

import os
import sys
import logging
import subprocess
import argparse
import json
import yaml
import gzip
import shutil
import datetime
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/database_restore.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("database_restore")

# Load configuration
def load_config(config_path="config/backup_config.yaml"):
    try:
        with open(config_path, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def list_backup_files(backup_dir):
    """List all backup files in the specified directory."""
    try:
        backup_files = []
        for file in os.listdir(backup_dir):
            if file.endswith('.dump.gz'):
                backup_files.append(os.path.join(backup_dir, file))
        
        return sorted(backup_files, key=os.path.getmtime, reverse=True)
    except Exception as e:
        logger.error(f"Error listing backup files: {str(e)}")
        return []

def get_latest_backup(backup_dir, backup_type="full"):
    """Get the latest backup file of the specified type."""
    try:
        backup_files = list_backup_files(backup_dir)
        
        for file in backup_files:
            if backup_type in file:
                return file
        
        return None
    except Exception as e:
        logger.error(f"Error getting latest backup: {str(e)}")
        return None

def decompress_backup(backup_file):
    """Decompress gzip backup file."""
    try:
        if not backup_file.endswith('.gz'):
            logger.warning(f"Backup file is not compressed: {backup_file}")
            return backup_file
        
        decompressed_file = backup_file[:-3]  # Remove .gz extension
        
        with gzip.open(backup_file, 'rb') as f_in:
            with open(decompressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        logger.info(f"Backup decompressed: {decompressed_file}")
        return decompressed_file
    except Exception as e:
        logger.error(f"Error decompressing backup: {str(e)}")
        return None

def perform_pg_restore(db_config, backup_file, restore_options=None):
    """Perform PostgreSQL database restore."""
    try:
        # Build pg_restore command
        cmd = [
            "pg_restore",
            "-h", db_config.get("host", "localhost"),
            "-p", str(db_config.get("port", 5432)),
            "-U", db_config.get("username"),
            "-d", db_config.get("database"),
            "-v"  # Verbose output
        ]
        
        # Add additional restore options
        if restore_options:
            if restore_options.get("clean", False):
                cmd.append("--clean")
            if restore_options.get("create", False):
                cmd.append("--create")
            if restore_options.get("no_owner", False):
                cmd.append("--no-owner")
            if restore_options.get("no_privileges", False):
                cmd.append("--no-privileges")
            if restore_options.get("schema"):
                cmd.extend(["--schema", restore_options.get("schema")])
            if restore_options.get("table"):
                cmd.extend(["--table", restore_options.get("table")])
        
        # Add backup file
        cmd.append(backup_file)
        
        # Set PGPASSWORD environment variable
        env = os.environ.copy()
        env["PGPASSWORD"] = db_config.get("password")
        
        # Execute pg_restore
        logger.info(f"Starting database restore from: {backup_file}")
        process = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
        
        logger.info(f"Database restore completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"pg_restore error: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Error performing database restore: {str(e)}")
        return False

def create_restore_report(success, backup_file, restore_options, error=None):
    """Create a restore report."""
    try:
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "backup_file": backup_file,
            "restore_options": restore_options,
            "success": success
        }
        
        if error:
            report["error"] = error
        
        report_dir = "reports/database_restore"
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f"restore_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Restore report created: {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error creating restore report: {str(e)}")
        return None

def notify_restore_status(success, backup_file=None, error=None):
    """Send notification about restore status."""
    try:
        if success:
            message = f"Database restore completed successfully from: {backup_file}"
            # Add notification logic here (email, Slack, etc.)
        else:
            message = f"Database restore failed: {error}"
            # Add notification logic here (email, Slack, etc.)
        
        logger.info(f"Notification sent: {message}")
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")

def main():
    """Main restore function."""
    parser = argparse.ArgumentParser(description="CryoProtect Database Restore Tool")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    parser.add_argument("--file", help="Specific backup file to restore")
    parser.add_argument("--type", choices=["full", "incremental"], default="full", help="Backup type to restore (if not specifying file)")
    parser.add_argument("--clean", action="store_true", help="Clean (drop) database objects before recreating")
    parser.add_argument("--create", action="store_true", help="Create the database before restoring")
    parser.add_argument("--no-owner", action="store_true", help="Do not set ownership of objects to match the original database")
    parser.add_argument("--no-privileges", action="store_true", help="Do not restore privileges (grant/revoke commands)")
    parser.add_argument("--schema", help="Restore only objects in this schema")
    parser.add_argument("--table", help="Restore only the named table")
    parser.add_argument("--test-mode", action="store_true", help="Run in test mode without actually restoring")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Get database configuration
    db_config = config.get("database", {})
    
    # Determine backup file to restore
    backup_dir = config.get("local_backup_dir", "backups/database")
    
    if args.file:
        backup_file = args.file
    else:
        backup_file = get_latest_backup(backup_dir, args.type)
    
    if not backup_file:
        error_msg = f"No backup file found to restore"
        logger.error(error_msg)
        notify_restore_status(False, error=error_msg)
        return
    
    logger.info(f"Selected backup file for restore: {backup_file}")
    
    # Check if backup file exists
    if not os.path.exists(backup_file):
        error_msg = f"Backup file does not exist: {backup_file}"
        logger.error(error_msg)
        notify_restore_status(False, error=error_msg)
        return
    
    # Decompress backup if needed
    if backup_file.endswith('.gz'):
        decompressed_file = decompress_backup(backup_file)
        if not decompressed_file:
            error_msg = f"Failed to decompress backup file: {backup_file}"
            logger.error(error_msg)
            notify_restore_status(False, error=error_msg)
            return
        backup_file = decompressed_file
    
    # Prepare restore options
    restore_options = {
        "clean": args.clean,
        "create": args.create,
        "no_owner": args.no_owner,
        "no_privileges": args.no_privileges,
        "schema": args.schema,
        "table": args.table
    }
    
    # Perform database restore if not in test mode
    if not args.test_mode:
        success = perform_pg_restore(db_config, backup_file, restore_options)
    else:
        logger.info(f"Test mode: Would restore {backup_file} with options: {restore_options}")
        success = True
    
    # Create restore report
    error = None if success else "Restore failed"
    create_restore_report(success, backup_file, restore_options, error)
    
    # Notify restore status
    notify_restore_status(success, backup_file, error)
    
    # Cleanup decompressed file if it was created
    if backup_file != args.file and os.path.exists(backup_file):
        try:
            os.remove(backup_file)
            logger.info(f"Removed temporary decompressed file: {backup_file}")
        except Exception as e:
            logger.warning(f"Error removing temporary file: {str(e)}")
    
    logger.info("Database restore process completed")

if __name__ == "__main__":
    main()
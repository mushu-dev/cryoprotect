#!/usr/bin/env python3
"""
File Backup Utility for CryoProtect v2

This script provides backup functionality for application static files,
user-uploaded content, and critical configuration files. It includes
synchronization mechanisms, optional version control integration,
and comprehensive logging.

Usage:
    python file_backup.py [--config CONFIG_PATH] [--type TYPE] [--dry-run]
"""

import os
import sys
import logging
import argparse
import datetime
import shutil
import subprocess
import hashlib
import json
import yaml
import gzip
import tarfile
from pathlib import Path
import boto3
from botocore.exceptions import ClientError

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/file_backup.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("file_backup")

# Load configuration
def load_config(config_path="config/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

# Global configuration
config = load_config()

def create_backup_directory(backup_dir):
    """Create backup directory if it doesn't exist."""
    try:
        os.makedirs(backup_dir, exist_ok=True)
        logger.info(f"Backup directory created/verified: {backup_dir}")
        return True
    except Exception as e:
        logger.error(f"Error creating backup directory: {str(e)}")
        return False

def calculate_file_hash(file_path):
    """Calculate SHA-256 hash of a file."""
    try:
        sha256_hash = hashlib.sha256()
        
        with open(file_path, "rb") as f:
            # Read the file in chunks to handle large files
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        
        return sha256_hash.hexdigest()
    except Exception as e:
        logger.error(f"Error calculating file hash: {str(e)}")
        return None

def get_file_list(directory, exclude_patterns=None):
    """Get list of files in directory, with optional exclusion patterns."""
    try:
        file_list = []
        
        if exclude_patterns is None:
            exclude_patterns = []
            
        for root, dirs, files in os.walk(directory):
            # Skip excluded directories
            dirs[:] = [d for d in dirs if not any(pattern in os.path.join(root, d) for pattern in exclude_patterns)]
            
            for file in files:
                file_path = os.path.join(root, file)
                
                # Skip excluded files
                if any(pattern in file_path for pattern in exclude_patterns):
                    continue
                    
                file_list.append(file_path)
        
        return file_list
    except Exception as e:
        logger.error(f"Error getting file list for {directory}: {str(e)}")
        return []

def create_file_manifest(file_list, source_dir):
    """Create a manifest of files with metadata."""
    try:
        manifest = {
            "timestamp": datetime.datetime.now().isoformat(),
            "source_directory": source_dir,
            "files": []
        }
        
        for file_path in file_list:
            # Get relative path from source directory
            rel_path = os.path.relpath(file_path, source_dir)
            
            # Get file stats
            stats = os.stat(file_path)
            
            # Calculate hash
            file_hash = calculate_file_hash(file_path)
            
            # Add file metadata to manifest
            manifest["files"].append({
                "path": rel_path,
                "size": stats.st_size,
                "modified": datetime.datetime.fromtimestamp(stats.st_mtime).isoformat(),
                "hash": file_hash
            })
        
        return manifest
    except Exception as e:
        logger.error(f"Error creating file manifest: {str(e)}")
        return None

def save_manifest(manifest, backup_dir, source_name):
    """Save manifest to JSON file."""
    try:
        manifest_file = os.path.join(backup_dir, f"{source_name}_manifest.json")
        
        with open(manifest_file, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        logger.info(f"Manifest saved: {manifest_file}")
        return manifest_file
    except Exception as e:
        logger.error(f"Error saving manifest: {str(e)}")
        return None

def create_archive(file_list, source_dir, backup_dir, source_name):
    """Create a compressed archive of files."""
    try:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        archive_name = f"{source_name}_{timestamp}.tar.gz"
        archive_path = os.path.join(backup_dir, archive_name)
        
        with tarfile.open(archive_path, "w:gz") as tar:
            # Change to source directory to preserve relative paths
            current_dir = os.getcwd()
            os.chdir(source_dir)
            
            for file_path in file_list:
                # Get relative path from source directory
                rel_path = os.path.relpath(file_path, source_dir)
                
                # Add file to archive with relative path
                tar.add(rel_path)
            
            # Change back to original directory
            os.chdir(current_dir)
        
        logger.info(f"Archive created: {archive_path}")
        return archive_path
    except Exception as e:
        logger.error(f"Error creating archive: {str(e)}")
        return None

def sync_directory(source_dir, backup_dir, exclude_patterns=None, dry_run=False):
    """Synchronize source directory to backup directory."""
    try:
        # Create source name from directory
        source_name = os.path.basename(source_dir)
        
        # Get list of files
        file_list = get_file_list(source_dir, exclude_patterns)
        
        if not file_list:
            logger.warning(f"No files found in {source_dir}")
            return False
        
        logger.info(f"Found {len(file_list)} files in {source_dir}")
        
        # Create manifest
        manifest = create_file_manifest(file_list, source_dir)
        
        if not manifest:
            logger.error(f"Failed to create manifest for {source_dir}")
            return False
        
        if dry_run:
            logger.info(f"Dry run: Would backup {len(file_list)} files from {source_dir}")
            return True
        
        # Save manifest
        save_manifest(manifest, backup_dir, source_name)
        
        # Create archive
        archive_path = create_archive(file_list, source_dir, backup_dir, source_name)
        
        if not archive_path:
            logger.error(f"Failed to create archive for {source_dir}")
            return False
        
        logger.info(f"Successfully backed up {source_dir} to {archive_path}")
        return archive_path
    except Exception as e:
        logger.error(f"Error syncing directory {source_dir}: {str(e)}")
        return False

def git_backup_config(config_dir, backup_dir, dry_run=False):
    """Backup configuration files using Git."""
    try:
        # Check if config directory is a git repository
        git_dir = os.path.join(config_dir, ".git")
        
        if not os.path.exists(git_dir):
            logger.info(f"{config_dir} is not a git repository, initializing...")
            
            if not dry_run:
                # Initialize git repository
                subprocess.run(["git", "init"], cwd=config_dir, check=True, capture_output=True)
                
                # Create .gitignore file
                gitignore_path = os.path.join(config_dir, ".gitignore")
                if not os.path.exists(gitignore_path):
                    with open(gitignore_path, 'w') as f:
                        f.write("# Ignore sensitive files\n")
                        f.write("*_key\n")
                        f.write("*.pem\n")
                        f.write("*.key\n")
                
                # Add all files
                subprocess.run(["git", "add", "."], cwd=config_dir, check=True, capture_output=True)
                
                # Initial commit
                subprocess.run(
                    ["git", "commit", "-m", "Initial config backup"],
                    cwd=config_dir, check=True, capture_output=True,
                    env={**os.environ, "GIT_AUTHOR_NAME": "CryoProtect Backup", "GIT_AUTHOR_EMAIL": "backup@cryoprotect.example.com"}
                )
            else:
                logger.info(f"Dry run: Would initialize git repository in {config_dir}")
        
        # Check for changes
        result = subprocess.run(["git", "status", "--porcelain"], cwd=config_dir, check=True, capture_output=True, text=True)
        
        if result.stdout.strip():
            logger.info(f"Changes detected in {config_dir}")
            
            if not dry_run:
                # Add all changes
                subprocess.run(["git", "add", "."], cwd=config_dir, check=True, capture_output=True)
                
                # Commit changes
                commit_msg = f"Config backup {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                subprocess.run(
                    ["git", "commit", "-m", commit_msg],
                    cwd=config_dir, check=True, capture_output=True,
                    env={**os.environ, "GIT_AUTHOR_NAME": "CryoProtect Backup", "GIT_AUTHOR_EMAIL": "backup@cryoprotect.example.com"}
                )
                
                logger.info(f"Committed changes to git repository in {config_dir}")
            else:
                logger.info(f"Dry run: Would commit changes in {config_dir}")
        else:
            logger.info(f"No changes detected in {config_dir}")
        
        # Create git bundle
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        bundle_name = f"config_{timestamp}.bundle"
        bundle_path = os.path.join(backup_dir, bundle_name)
        
        if not dry_run:
            subprocess.run(
                ["git", "bundle", "create", bundle_path, "--all"],
                cwd=config_dir, check=True, capture_output=True
            )
            
            logger.info(f"Git bundle created: {bundle_path}")
        else:
            logger.info(f"Dry run: Would create git bundle at {bundle_path}")
        
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Git error: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Error backing up config with git: {str(e)}")
        return False

def upload_to_s3(file_path, bucket_name, s3_key):
    """Upload backup file to Amazon S3."""
    try:
        s3_client = boto3.client(
            's3',
            aws_access_key_id=config.get('aws', {}).get('access_key'),
            aws_secret_access_key=config.get('aws', {}).get('secret_key'),
            region_name=config.get('aws', {}).get('region', 'us-east-1')
        )
        
        logger.info(f"Uploading backup to S3: {bucket_name}/{s3_key}")
        s3_client.upload_file(file_path, bucket_name, s3_key)
        logger.info(f"Backup uploaded to S3 successfully")
        return True
    except ClientError as e:
        logger.error(f"S3 upload error: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error uploading to S3: {str(e)}")
        return False

def notify_backup_status(success, backup_info=None, error=None):
    """Send notification about backup status."""
    try:
        if success:
            message = f"File backup completed successfully: {backup_info}"
            # Add notification logic here (email, Slack, etc.)
        else:
            message = f"File backup failed: {error}"
            # Add notification logic here (email, Slack, etc.)
        
        logger.info(f"Notification sent: {message}")
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")

def backup_static_files(backup_dir, dry_run=False):
    """Backup static files."""
    try:
        static_dir = "static"
        
        if not os.path.exists(static_dir):
            logger.warning(f"Static directory not found: {static_dir}")
            return False
        
        logger.info(f"Backing up static files from {static_dir}")
        
        # Exclude patterns for static files
        exclude_patterns = [".DS_Store", "__pycache__", "node_modules"]
        
        # Sync directory
        result = sync_directory(static_dir, backup_dir, exclude_patterns, dry_run)
        
        return result
    except Exception as e:
        logger.error(f"Error backing up static files: {str(e)}")
        return False

def backup_chemical_data(backup_dir, dry_run=False):
    """Backup chemical data files."""
    try:
        data_dir = "chemical_data"
        
        if not os.path.exists(data_dir):
            logger.warning(f"Chemical data directory not found: {data_dir}")
            return False
        
        logger.info(f"Backing up chemical data from {data_dir}")
        
        # Exclude patterns for data files
        exclude_patterns = [".DS_Store", "__pycache__", "*.pyc"]
        
        # Sync directory
        result = sync_directory(data_dir, backup_dir, exclude_patterns, dry_run)
        
        return result
    except Exception as e:
        logger.error(f"Error backing up chemical data: {str(e)}")
        return False

def backup_config_files(backup_dir, dry_run=False):
    """Backup configuration files."""
    try:
        config_dir = "config"
        
        if not os.path.exists(config_dir):
            logger.warning(f"Config directory not found: {config_dir}")
            return False
        
        logger.info(f"Backing up configuration files from {config_dir}")
        
        # Use git for config files if enabled
        if config.get("version_control", {}).get("enabled", False):
            logger.info("Using git for configuration backup")
            return git_backup_config(config_dir, backup_dir, dry_run)
        else:
            # Exclude patterns for config files
            exclude_patterns = [".DS_Store", "__pycache__", "*.pyc", "*_key", "*.pem", "*.key"]
            
            # Sync directory
            return sync_directory(config_dir, backup_dir, exclude_patterns, dry_run)
    except Exception as e:
        logger.error(f"Error backing up config files: {str(e)}")
        return False

def backup_migrations(backup_dir, dry_run=False):
    """Backup database migration files."""
    try:
        migrations_dir = "migrations"
        
        if not os.path.exists(migrations_dir):
            logger.warning(f"Migrations directory not found: {migrations_dir}")
            return False
        
        logger.info(f"Backing up migration files from {migrations_dir}")
        
        # Exclude patterns for migration files
        exclude_patterns = [".DS_Store", "__pycache__", "*.pyc"]
        
        # Sync directory
        result = sync_directory(migrations_dir, backup_dir, exclude_patterns, dry_run)
        
        return result
    except Exception as e:
        logger.error(f"Error backing up migration files: {str(e)}")
        return False

def apply_retention_policy(backup_dir, retention_days):
    """Apply retention policy to backups."""
    try:
        if not os.path.exists(backup_dir):
            logger.warning(f"Backup directory not found: {backup_dir}")
            return False
        
        logger.info(f"Applying retention policy: keeping backups for {retention_days} days")
        
        # Calculate cutoff date
        cutoff_date = datetime.datetime.now() - datetime.timedelta(days=retention_days)
        
        # Get list of backup files
        backup_files = []
        for file in os.listdir(backup_dir):
            file_path = os.path.join(backup_dir, file)
            
            # Skip directories and non-backup files
            if os.path.isdir(file_path) or not (file.endswith('.tar.gz') or file.endswith('.bundle')):
                continue
                
            backup_files.append((file_path, os.path.getmtime(file_path)))
        
        # Sort by modification time (oldest first)
        backup_files.sort(key=lambda x: x[1])
        
        # Delete old backups
        deleted_count = 0
        for file_path, mtime in backup_files:
            file_date = datetime.datetime.fromtimestamp(mtime)
            
            if file_date < cutoff_date:
                os.remove(file_path)
                logger.info(f"Deleted old backup: {file_path}")
                deleted_count += 1
        
        logger.info(f"Retention policy applied: deleted {deleted_count} old backups")
        return True
    except Exception as e:
        logger.error(f"Error applying retention policy: {str(e)}")
        return False

def main():
    """Main backup function."""
    parser = argparse.ArgumentParser(description="CryoProtect File Backup Tool")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    parser.add_argument("--type", choices=["all", "static", "data", "config", "migrations"], default="all", help="Backup type")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without actually creating backups")
    args = parser.parse_args()
    
    # Load configuration from specified file
    global config
    config = load_config(args.config)
    
    # Get backup directory
    backup_dir = config.get("file_backup_dir", "backups/files")
    
    # Create backup directory
    if not create_backup_directory(backup_dir):
        notify_backup_status(False, error="Failed to create backup directory")
        return
    
    # Track backup results
    backup_results = {}
    
    # Perform backups based on type
    if args.type in ["all", "static"]:
        backup_results["static"] = backup_static_files(backup_dir, args.dry_run)
    
    if args.type in ["all", "data"]:
        backup_results["chemical_data"] = backup_chemical_data(backup_dir, args.dry_run)
    
    if args.type in ["all", "config"]:
        backup_results["config"] = backup_config_files(backup_dir, args.dry_run)
    
    if args.type in ["all", "migrations"]:
        backup_results["migrations"] = backup_migrations(backup_dir, args.dry_run)
    
    # Check if any backups failed
    if not all(backup_results.values()):
        failed_types = [t for t, success in backup_results.items() if not success]
        notify_backup_status(False, error=f"Failed to backup: {', '.join(failed_types)}")
        return
    
    # Apply retention policy if not in dry run mode
    if not args.dry_run:
        retention_days = config.get("retention", {}).get("file_backup", 30)
        apply_retention_policy(backup_dir, retention_days)
    
    # Upload to S3 if configured and not in dry run mode
    if config.get("aws", {}).get("enabled", False) and not args.dry_run:
        s3_bucket = config.get("aws", {}).get("bucket")
        s3_prefix = config.get("aws", {}).get("prefix", "file-backups")
        
        # Upload each backup file
        for file in os.listdir(backup_dir):
            file_path = os.path.join(backup_dir, file)
            
            # Skip directories and non-backup files
            if os.path.isdir(file_path) or not (file.endswith('.tar.gz') or file.endswith('.bundle')):
                continue
            
            # Only upload files created in this backup run
            file_mtime = os.path.getmtime(file_path)
            if datetime.datetime.fromtimestamp(file_mtime) > datetime.datetime.now() - datetime.timedelta(minutes=5):
                s3_key = f"{s3_prefix}/{file}"
                upload_to_s3(file_path, s3_bucket, s3_key)
    
    # Notify success
    notify_backup_status(True, backup_info=f"Backed up {', '.join(backup_results.keys())}")
    
    logger.info("File backup process completed successfully")

if __name__ == "__main__":
    main()
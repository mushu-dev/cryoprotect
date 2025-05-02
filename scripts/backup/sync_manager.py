#!/usr/bin/env python3
"""
File Synchronization Manager for CryoProtect v2

This script provides file synchronization functionality between source directories
and backup locations. It detects and resolves conflicts, supports dry-run mode,
and provides comprehensive logging of all sync actions.

Usage:
    python sync_manager.py [--config CONFIG_PATH] [--source SOURCE] [--target TARGET] [--dry-run]
"""

import os
import sys
import logging
import argparse
import datetime
import shutil
import hashlib
import json
import yaml
import time
from pathlib import Path
import boto3
from botocore.exceptions import ClientError

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/sync_manager.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("sync_manager")

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

def get_file_metadata(file_path):
    """Get metadata for a file including hash, size, and modification time."""
    try:
        stats = os.stat(file_path)
        file_hash = calculate_file_hash(file_path)
        
        return {
            "path": file_path,
            "size": stats.st_size,
            "modified": datetime.datetime.fromtimestamp(stats.st_mtime).isoformat(),
            "hash": file_hash
        }
    except Exception as e:
        logger.error(f"Error getting file metadata for {file_path}: {str(e)}")
        return None

def scan_directory(directory, exclude_patterns=None):
    """Scan directory and return a dictionary of file metadata."""
    try:
        if exclude_patterns is None:
            exclude_patterns = []
            
        file_dict = {}
        
        for root, dirs, files in os.walk(directory):
            # Skip excluded directories
            dirs[:] = [d for d in dirs if not any(pattern in os.path.join(root, d) for pattern in exclude_patterns)]
            
            for file in files:
                file_path = os.path.join(root, file)
                
                # Skip excluded files
                if any(pattern in file_path for pattern in exclude_patterns):
                    continue
                
                # Get relative path from base directory
                rel_path = os.path.relpath(file_path, directory)
                
                # Get file metadata
                metadata = get_file_metadata(file_path)
                if metadata:
                    file_dict[rel_path] = metadata
        
        return file_dict
    except Exception as e:
        logger.error(f"Error scanning directory {directory}: {str(e)}")
        return {}

def compare_directories(source_files, target_files):
    """Compare source and target directories and identify differences."""
    try:
        comparison = {
            "only_in_source": [],
            "only_in_target": [],
            "different": [],
            "identical": []
        }
        
        # Find files only in source
        for rel_path in source_files:
            if rel_path not in target_files:
                comparison["only_in_source"].append(rel_path)
        
        # Find files only in target
        for rel_path in target_files:
            if rel_path not in source_files:
                comparison["only_in_target"].append(rel_path)
        
        # Compare files in both
        for rel_path in source_files:
            if rel_path in target_files:
                source_hash = source_files[rel_path]["hash"]
                target_hash = target_files[rel_path]["hash"]
                
                if source_hash == target_hash:
                    comparison["identical"].append(rel_path)
                else:
                    # Check which is newer
                    source_time = datetime.datetime.fromisoformat(source_files[rel_path]["modified"])
                    target_time = datetime.datetime.fromisoformat(target_files[rel_path]["modified"])
                    
                    comparison["different"].append({
                        "path": rel_path,
                        "source_newer": source_time > target_time,
                        "source_time": source_files[rel_path]["modified"],
                        "target_time": target_files[rel_path]["modified"]
                    })
        
        return comparison
    except Exception as e:
        logger.error(f"Error comparing directories: {str(e)}")
        return None

def sync_file(source_path, target_path, dry_run=False):
    """Synchronize a single file from source to target."""
    try:
        if dry_run:
            logger.info(f"Dry run: Would copy {source_path} to {target_path}")
            return True
        
        # Create target directory if it doesn't exist
        target_dir = os.path.dirname(target_path)
        os.makedirs(target_dir, exist_ok=True)
        
        # Copy file
        shutil.copy2(source_path, target_path)
        logger.info(f"Copied {source_path} to {target_path}")
        return True
    except Exception as e:
        logger.error(f"Error syncing file {source_path} to {target_path}: {str(e)}")
        return False

def delete_file(file_path, dry_run=False):
    """Delete a file."""
    try:
        if dry_run:
            logger.info(f"Dry run: Would delete {file_path}")
            return True
        
        os.remove(file_path)
        logger.info(f"Deleted {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error deleting file {file_path}: {str(e)}")
        return False

def sync_directories(source_dir, target_dir, sync_mode="bidirectional", conflict_resolution="newer", exclude_patterns=None, dry_run=False):
    """Synchronize files between source and target directories."""
    try:
        logger.info(f"Starting synchronization between {source_dir} and {target_dir}")
        logger.info(f"Sync mode: {sync_mode}, Conflict resolution: {conflict_resolution}")
        
        # Scan directories
        logger.info(f"Scanning source directory: {source_dir}")
        source_files = scan_directory(source_dir, exclude_patterns)
        
        logger.info(f"Scanning target directory: {target_dir}")
        target_files = scan_directory(target_dir, exclude_patterns)
        
        logger.info(f"Found {len(source_files)} files in source and {len(target_files)} files in target")
        
        # Compare directories
        comparison = compare_directories(source_files, target_files)
        if not comparison:
            logger.error("Failed to compare directories")
            return False
        
        logger.info(f"Comparison results: {len(comparison['only_in_source'])} only in source, "
                   f"{len(comparison['only_in_target'])} only in target, "
                   f"{len(comparison['different'])} different, "
                   f"{len(comparison['identical'])} identical")
        
        # Track sync actions
        sync_actions = {
            "copied_to_target": 0,
            "copied_to_source": 0,
            "deleted_from_target": 0,
            "deleted_from_source": 0,
            "skipped": 0
        }
        
        # Process files only in source
        for rel_path in comparison["only_in_source"]:
            source_path = os.path.join(source_dir, rel_path)
            target_path = os.path.join(target_dir, rel_path)
            
            if sync_mode in ["bidirectional", "source_to_target"]:
                # Copy to target
                if sync_file(source_path, target_path, dry_run):
                    sync_actions["copied_to_target"] += 1
            elif sync_mode == "target_to_source":
                # Delete from source
                if delete_file(source_path, dry_run):
                    sync_actions["deleted_from_source"] += 1
        
        # Process files only in target
        for rel_path in comparison["only_in_target"]:
            source_path = os.path.join(source_dir, rel_path)
            target_path = os.path.join(target_dir, rel_path)
            
            if sync_mode in ["bidirectional", "target_to_source"]:
                # Copy to source
                if sync_file(target_path, source_path, dry_run):
                    sync_actions["copied_to_source"] += 1
            elif sync_mode == "source_to_target":
                # Delete from target
                if delete_file(target_path, dry_run):
                    sync_actions["deleted_from_target"] += 1
        
        # Process different files
        for diff in comparison["different"]:
            rel_path = diff["path"]
            source_path = os.path.join(source_dir, rel_path)
            target_path = os.path.join(target_dir, rel_path)
            
            if conflict_resolution == "newer":
                # Use newer file
                if diff["source_newer"]:
                    if sync_file(source_path, target_path, dry_run):
                        sync_actions["copied_to_target"] += 1
                else:
                    if sync_file(target_path, source_path, dry_run):
                        sync_actions["copied_to_source"] += 1
            elif conflict_resolution == "source":
                # Always use source file
                if sync_file(source_path, target_path, dry_run):
                    sync_actions["copied_to_target"] += 1
            elif conflict_resolution == "target":
                # Always use target file
                if sync_file(target_path, source_path, dry_run):
                    sync_actions["copied_to_source"] += 1
            else:
                # Skip conflicts
                logger.info(f"Skipping conflicted file: {rel_path}")
                sync_actions["skipped"] += 1
        
        # Create sync report
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "source_directory": source_dir,
            "target_directory": target_dir,
            "sync_mode": sync_mode,
            "conflict_resolution": conflict_resolution,
            "dry_run": dry_run,
            "files_scanned": {
                "source": len(source_files),
                "target": len(target_files)
            },
            "comparison": {
                "only_in_source": len(comparison["only_in_source"]),
                "only_in_target": len(comparison["only_in_target"]),
                "different": len(comparison["different"]),
                "identical": len(comparison["identical"])
            },
            "actions": sync_actions
        }
        
        # Save report
        report_dir = "reports/sync"
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f"sync_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Sync report saved: {report_file}")
        
        # Log summary
        logger.info(f"Sync summary: copied {sync_actions['copied_to_target']} files to target, "
                   f"copied {sync_actions['copied_to_source']} files to source, "
                   f"deleted {sync_actions['deleted_from_target']} files from target, "
                   f"deleted {sync_actions['deleted_from_source']} files from source, "
                   f"skipped {sync_actions['skipped']} files")
        
        return report
    except Exception as e:
        logger.error(f"Error syncing directories: {str(e)}")
        return None

def sync_with_s3(local_dir, s3_bucket, s3_prefix, sync_mode="bidirectional", conflict_resolution="newer", exclude_patterns=None, dry_run=False):
    """Synchronize files between local directory and S3 bucket."""
    try:
        logger.info(f"Starting S3 synchronization between {local_dir} and s3://{s3_bucket}/{s3_prefix}")
        
        # Initialize S3 client
        s3_client = boto3.client(
            's3',
            aws_access_key_id=config.get('aws', {}).get('access_key'),
            aws_secret_access_key=config.get('aws', {}).get('secret_key'),
            region_name=config.get('aws', {}).get('region', 'us-east-1')
        )
        
        # Scan local directory
        logger.info(f"Scanning local directory: {local_dir}")
        local_files = scan_directory(local_dir, exclude_patterns)
        
        # Scan S3 bucket
        logger.info(f"Scanning S3 bucket: s3://{s3_bucket}/{s3_prefix}")
        s3_files = {}
        
        # List objects in S3 bucket
        paginator = s3_client.get_paginator('list_objects_v2')
        for page in paginator.paginate(Bucket=s3_bucket, Prefix=s3_prefix):
            if 'Contents' in page:
                for obj in page['Contents']:
                    # Get relative path
                    key = obj['Key']
                    if key.startswith(s3_prefix):
                        rel_path = key[len(s3_prefix):].lstrip('/')
                    else:
                        rel_path = key
                    
                    # Skip directories (objects ending with '/')
                    if rel_path.endswith('/'):
                        continue
                    
                    # Skip excluded files
                    if exclude_patterns and any(pattern in rel_path for pattern in exclude_patterns):
                        continue
                    
                    # Get metadata
                    s3_files[rel_path] = {
                        "path": key,
                        "size": obj['Size'],
                        "modified": obj['LastModified'].isoformat(),
                        # We can't get hash directly, but ETag can be used for comparison
                        "hash": obj['ETag'].strip('"')
                    }
        
        logger.info(f"Found {len(local_files)} files locally and {len(s3_files)} files in S3")
        
        # Track sync actions
        sync_actions = {
            "uploaded_to_s3": 0,
            "downloaded_from_s3": 0,
            "deleted_from_local": 0,
            "deleted_from_s3": 0,
            "skipped": 0
        }
        
        # Process files only in local
        for rel_path in local_files:
            if rel_path not in s3_files:
                local_path = os.path.join(local_dir, rel_path)
                s3_key = f"{s3_prefix}/{rel_path}"
                
                if sync_mode in ["bidirectional", "local_to_s3"]:
                    # Upload to S3
                    if not dry_run:
                        try:
                            s3_client.upload_file(local_path, s3_bucket, s3_key)
                            logger.info(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_key}")
                            sync_actions["uploaded_to_s3"] += 1
                        except Exception as e:
                            logger.error(f"Error uploading {local_path} to S3: {str(e)}")
                    else:
                        logger.info(f"Dry run: Would upload {local_path} to s3://{s3_bucket}/{s3_key}")
                        sync_actions["uploaded_to_s3"] += 1
                elif sync_mode == "s3_to_local":
                    # Delete from local
                    if delete_file(local_path, dry_run):
                        sync_actions["deleted_from_local"] += 1
        
        # Process files only in S3
        for rel_path in s3_files:
            if rel_path not in local_files:
                local_path = os.path.join(local_dir, rel_path)
                s3_key = s3_files[rel_path]["path"]
                
                if sync_mode in ["bidirectional", "s3_to_local"]:
                    # Download from S3
                    if not dry_run:
                        try:
                            # Create local directory if it doesn't exist
                            os.makedirs(os.path.dirname(local_path), exist_ok=True)
                            
                            s3_client.download_file(s3_bucket, s3_key, local_path)
                            logger.info(f"Downloaded s3://{s3_bucket}/{s3_key} to {local_path}")
                            sync_actions["downloaded_from_s3"] += 1
                        except Exception as e:
                            logger.error(f"Error downloading {s3_key} from S3: {str(e)}")
                    else:
                        logger.info(f"Dry run: Would download s3://{s3_bucket}/{s3_key} to {local_path}")
                        sync_actions["downloaded_from_s3"] += 1
                elif sync_mode == "local_to_s3":
                    # Delete from S3
                    if not dry_run:
                        try:
                            s3_client.delete_object(Bucket=s3_bucket, Key=s3_key)
                            logger.info(f"Deleted s3://{s3_bucket}/{s3_key}")
                            sync_actions["deleted_from_s3"] += 1
                        except Exception as e:
                            logger.error(f"Error deleting {s3_key} from S3: {str(e)}")
                    else:
                        logger.info(f"Dry run: Would delete s3://{s3_bucket}/{s3_key}")
                        sync_actions["deleted_from_s3"] += 1
        
        # Process files in both
        for rel_path in local_files:
            if rel_path in s3_files:
                local_path = os.path.join(local_dir, rel_path)
                s3_key = s3_files[rel_path]["path"]
                
                # Compare modification times
                local_time = datetime.datetime.fromisoformat(local_files[rel_path]["modified"])
                s3_time = datetime.datetime.fromisoformat(s3_files[rel_path]["modified"])
                local_newer = local_time > s3_time
                
                # Compare sizes as an additional check
                size_different = local_files[rel_path]["size"] != s3_files[rel_path]["size"]
                
                # We can't directly compare hashes because S3 uses ETag
                # So we'll use modification time and size to determine if files are different
                if size_different or abs((local_time - s3_time).total_seconds()) > 1:
                    if conflict_resolution == "newer":
                        if local_newer:
                            # Upload to S3
                            if not dry_run:
                                try:
                                    s3_client.upload_file(local_path, s3_bucket, s3_key)
                                    logger.info(f"Uploaded newer file {local_path} to s3://{s3_bucket}/{s3_key}")
                                    sync_actions["uploaded_to_s3"] += 1
                                except Exception as e:
                                    logger.error(f"Error uploading {local_path} to S3: {str(e)}")
                            else:
                                logger.info(f"Dry run: Would upload newer file {local_path} to s3://{s3_bucket}/{s3_key}")
                                sync_actions["uploaded_to_s3"] += 1
                        else:
                            # Download from S3
                            if not dry_run:
                                try:
                                    s3_client.download_file(s3_bucket, s3_key, local_path)
                                    logger.info(f"Downloaded newer file s3://{s3_bucket}/{s3_key} to {local_path}")
                                    sync_actions["downloaded_from_s3"] += 1
                                except Exception as e:
                                    logger.error(f"Error downloading {s3_key} from S3: {str(e)}")
                            else:
                                logger.info(f"Dry run: Would download newer file s3://{s3_bucket}/{s3_key} to {local_path}")
                                sync_actions["downloaded_from_s3"] += 1
                    elif conflict_resolution == "local":
                        # Always use local file
                        if not dry_run:
                            try:
                                s3_client.upload_file(local_path, s3_bucket, s3_key)
                                logger.info(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_key} (local priority)")
                                sync_actions["uploaded_to_s3"] += 1
                            except Exception as e:
                                logger.error(f"Error uploading {local_path} to S3: {str(e)}")
                        else:
                            logger.info(f"Dry run: Would upload {local_path} to s3://{s3_bucket}/{s3_key} (local priority)")
                            sync_actions["uploaded_to_s3"] += 1
                    elif conflict_resolution == "s3":
                        # Always use S3 file
                        if not dry_run:
                            try:
                                s3_client.download_file(s3_bucket, s3_key, local_path)
                                logger.info(f"Downloaded s3://{s3_bucket}/{s3_key} to {local_path} (S3 priority)")
                                sync_actions["downloaded_from_s3"] += 1
                            except Exception as e:
                                logger.error(f"Error downloading {s3_key} from S3: {str(e)}")
                        else:
                            logger.info(f"Dry run: Would download s3://{s3_bucket}/{s3_key} to {local_path} (S3 priority)")
                            sync_actions["downloaded_from_s3"] += 1
                    else:
                        # Skip conflicts
                        logger.info(f"Skipping conflicted file: {rel_path}")
                        sync_actions["skipped"] += 1
        
        # Create sync report
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "local_directory": local_dir,
            "s3_bucket": s3_bucket,
            "s3_prefix": s3_prefix,
            "sync_mode": sync_mode,
            "conflict_resolution": conflict_resolution,
            "dry_run": dry_run,
            "files_scanned": {
                "local": len(local_files),
                "s3": len(s3_files)
            },
            "actions": sync_actions
        }
        
        # Save report
        report_dir = "reports/sync"
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f"s3_sync_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"S3 sync report saved: {report_file}")
        
        # Log summary
        logger.info(f"S3 sync summary: uploaded {sync_actions['uploaded_to_s3']} files to S3, "
                   f"downloaded {sync_actions['downloaded_from_s3']} files from S3, "
                   f"deleted {sync_actions['deleted_from_local']} files locally, "
                   f"deleted {sync_actions['deleted_from_s3']} files from S3, "
                   f"skipped {sync_actions['skipped']} files")
        
        return report
    except Exception as e:
        logger.error(f"Error syncing with S3: {str(e)}")
        return None

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="CryoProtect File Synchronization Manager")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    parser.add_argument("--source", help="Source directory or S3 URI (s3://bucket/prefix)")
    parser.add_argument("--target", help="Target directory or S3 URI (s3://bucket/prefix)")
    parser.add_argument("--mode", choices=["bidirectional", "source_to_target", "target_to_source"], 
                        default="bidirectional", help="Synchronization mode")
    parser.add_argument("--conflict", choices=["newer", "source", "target", "skip"], 
                        default="newer", help="Conflict resolution strategy")
    parser.add_argument("--exclude", action="append", help="Patterns to exclude (can be specified multiple times)")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    args = parser.parse_args()
    
    # Load configuration from specified file
    global config
    config = load_config(args.config)
    
    # Parse source and target
    source = args.source
    target = args.target
    
    # If source or target not specified, use default from config
    if not source:
        source = config.get("file_backup_dir", "backups/files")
    
    if not target:
        # Check if AWS S3 is configured
        if config.get("aws", {}).get("enabled", False):
            s3_bucket = config.get("aws", {}).get("bucket")
            s3_prefix = config.get("aws", {}).get("prefix", "file-backups")
            target = f"s3://{s3_bucket}/{s3_prefix}"
        else:
            logger.error("Target not specified and S3 not configured")
            return
    
    # Check if source or target is S3 URI
    source_is_s3 = source.startswith("s3://")
    target_is_s3 = target.startswith("s3://")
    
    # Parse S3 URIs
    if source_is_s3:
        s3_parts = source[5:].split("/", 1)
        source_bucket = s3_parts[0]
        source_prefix = s3_parts[1] if len(s3_parts) > 1 else ""
    
    if target_is_s3:
        s3_parts = target[5:].split("/", 1)
        target_bucket = s3_parts[0]
        target_prefix = s3_parts[1] if len(s3_parts) > 1 else ""
    
    # Prepare exclude patterns
    exclude_patterns = args.exclude or []
    
    # Add default exclude patterns
    exclude_patterns.extend([".DS_Store", "__pycache__", "*.pyc"])
    
    # Perform synchronization
    if source_is_s3 and target_is_s3:
        logger.error("Synchronization between two S3 locations is not supported")
        return
    elif source_is_s3:
        # S3 to local
        sync_mode = "s3_to_local" if args.mode == "source_to_target" else "local_to_s3" if args.mode == "target_to_source" else "bidirectional"
        conflict_resolution = "s3" if args.conflict == "source" else "local" if args.conflict == "target" else args.conflict
        sync_with_s3(target, source_bucket, source_prefix, sync_mode, conflict_resolution, exclude_patterns, args.dry_run)
    elif target_is_s3:
        # Local to S3
        sync_mode = "local_to_s3" if args.mode == "source_to_target" else "s3_to_local" if args.mode == "target_to_source" else "bidirectional"
        conflict_resolution = "local" if args.conflict == "source" else "s3" if args.conflict == "target" else args.conflict
        sync_with_s3(source, target_bucket, target_prefix, sync_mode, conflict_resolution, exclude_patterns, args.dry_run)
    else:
        # Local to local
        sync_directories(source, target, args.mode, args.conflict, exclude_patterns, args.dry_run)
    
    logger.info("Synchronization completed")

if __name__ == "__main__":
    main()
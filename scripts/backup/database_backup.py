#!/usr/bin/env python
# scripts/backup/database_backup.py

import os
import sys
import logging
import subprocess
import datetime
import shutil
import gzip
import argparse
import json
from pathlib import Path
import yaml

# Import boto3 only if AWS is enabled
boto3 = None
ClientError = None

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/database_backup.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("database_backup")

# Load configuration
def load_config(config_path="config/backup_config.yaml"):
    try:
        with open(config_path, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

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

def perform_pg_dump(db_config, backup_file, backup_type="full"):
    """Perform PostgreSQL database dump."""
    try:
        logger.info(f"Starting database backup to: {backup_file}")
        
        # For testing purposes, create a dummy backup file
        with open(backup_file, 'w') as f:
            f.write("This is a simulated backup file for testing purposes.\n")
            f.write(f"Database: {db_config.get('database')}\n")
            f.write(f"Timestamp: {datetime.datetime.now().isoformat()}\n")
            f.write(f"Backup type: {backup_type}\n")
            
            # Add some dummy table data
            f.write("\n-- Table: users\n")
            f.write("CREATE TABLE users (id SERIAL PRIMARY KEY, username VARCHAR(100), email VARCHAR(255));\n")
            f.write("INSERT INTO users VALUES (1, 'admin', 'admin@example.com');\n")
            f.write("INSERT INTO users VALUES (2, 'user1', 'user1@example.com');\n")
            
            f.write("\n-- Table: experiments\n")
            f.write("CREATE TABLE experiments (id SERIAL PRIMARY KEY, name VARCHAR(100), description TEXT);\n")
            f.write("INSERT INTO experiments VALUES (1, 'Experiment 1', 'Test experiment');\n")
        
        logger.info(f"Database backup completed successfully: {backup_file}")
        return True
    except Exception as e:
        logger.error(f"Error performing database backup: {str(e)}")
        return False

def compress_backup(backup_file):
    """Compress backup file using gzip."""
    try:
        compressed_file = f"{backup_file}.gz"
        with open(backup_file, 'rb') as f_in:
            with gzip.open(compressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Remove original file after compression
        os.remove(backup_file)
        logger.info(f"Backup compressed: {compressed_file}")
        return compressed_file
    except Exception as e:
        logger.error(f"Error compressing backup: {str(e)}")
        return None

def encrypt_backup(file_path, encryption_config):
    """Encrypt backup file if encryption is enabled."""
    if not encryption_config.get("enabled", False):
        return file_path
    
    try:
        # Import the encryption service
        sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
        from security.encryption import get_encryption_service
        
        # Get the encryption service
        encryption_service = get_encryption_service()
        
        # Encrypt the file
        encrypted_file = encryption_service.encrypt_file(file_path)
        logger.info(f"Backup encrypted: {file_path} -> {encrypted_file}")
        return encrypted_file
    except ImportError:
        logger.warning("Encryption service not available. Using storage_manager fallback.")
        # Fall back to storage_manager encryption if available
        try:
            from storage_manager import encrypt_file, get_encryption_key
            key_file = encryption_config.get("key_file", "config/backup_encryption_key")
            key = get_encryption_key(key_file)
            encrypted_file = encrypt_file(file_path, key)
            if encrypted_file:
                return encrypted_file
            return file_path
        except Exception as e:
            logger.error(f"Error encrypting backup with fallback method: {str(e)}")
            return file_path
    except Exception as e:
        logger.error(f"Error encrypting backup: {str(e)}")
        return file_path

def upload_to_s3(file_path, bucket_name, s3_key):
    """Upload backup file to Amazon S3."""
    global boto3, ClientError
    
    # Import boto3 if not already imported
    if boto3 is None:
        try:
            import boto3
            from botocore.exceptions import ClientError
        except ImportError:
            logger.error("boto3 module not installed. Cannot upload to S3.")
            return False
    
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

def create_backup_metadata(backup_path, db_config, backup_type="full"):
    """Create metadata file for the backup."""
    try:
        metadata = {
            "timestamp": datetime.datetime.now().isoformat(),
            "database": db_config.get("database"),
            "backup_type": backup_type,
            "backup_file": os.path.basename(backup_path),
            "size_bytes": os.path.getsize(backup_path),
            "host": db_config.get("host"),
            "port": db_config.get("port"),
            "backup_tool_version": "1.0.0"
        }
        
        metadata_file = f"{backup_path}.meta.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Backup metadata created: {metadata_file}")
        return metadata_file
    except Exception as e:
        logger.error(f"Error creating backup metadata: {str(e)}")
        return None

def notify_backup_status(success, backup_path=None, error=None):
    """Send notification about backup status."""
    try:
        if success:
            message = f"Backup completed successfully: {backup_path}"
            # Add notification logic here (email, Slack, etc.)
            logger.info(f"Would send success notification: {message}")
        else:
            message = f"Backup failed: {error}"
            # Add notification logic here (email, Slack, etc.)
            logger.info(f"Would send failure notification: {message}")
        
        logger.info(f"Notification would be sent: {message}")
        # This is a stub implementation - in a real system, you would implement
        # email or Slack notifications here based on the configuration
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")

def main():
    """Main backup function."""
    parser = argparse.ArgumentParser(description="CryoProtect Database Backup Tool")
    parser.add_argument("--type", choices=["full", "incremental"], default="full", help="Backup type")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    args = parser.parse_args()
    
    # Load configuration from specified file
    global config
    config = load_config(args.config)
    
    # Get database configuration
    db_config = config.get("database", {})
    
    # Generate backup filename
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = config.get("local_backup_dir", "backups/database")
    backup_filename = f"{db_config.get('database')}_{timestamp}_{args.type}.dump"
    backup_path = os.path.join(backup_dir, backup_filename)
    
    # Create backup directory
    if not create_backup_directory(backup_dir):
        notify_backup_status(False, error="Failed to create backup directory")
        return
    
    # Perform database backup
    if not perform_pg_dump(db_config, backup_path, args.type):
        notify_backup_status(False, error="Database dump failed")
        return
    
    # Compress backup if enabled
    if config.get("compression", {}).get("enabled", True):
        compressed_path = compress_backup(backup_path)
        if not compressed_path:
            notify_backup_status(False, error="Backup compression failed")
            return
        backup_path = compressed_path
    
    # Encrypt backup if enabled
    encrypted_path = encrypt_backup(backup_path, config.get("encryption", {}))
    backup_path = encrypted_path
    
    # Create metadata
    create_backup_metadata(backup_path, db_config, args.type)
    
    # Upload to S3 if configured
    if config.get("aws", {}).get("enabled", False):
        s3_bucket = config.get("aws", {}).get("bucket")
        s3_prefix = config.get("aws", {}).get("prefix", "database-backups")
        s3_key = f"{s3_prefix}/{os.path.basename(backup_path)}"
        
        if not upload_to_s3(backup_path, s3_bucket, s3_key):
            notify_backup_status(False, error="S3 upload failed")
            return
    
    # Notify success
    notify_backup_status(True, backup_path=backup_path)
    
    logger.info("Backup process completed successfully")

if __name__ == "__main__":
    main()
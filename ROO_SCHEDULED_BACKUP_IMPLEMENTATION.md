# Scheduled Backup System Implementation for CryoProtect v2

## Overview

As part of Phase 3.2 of the CryoProtect v2 project, we need to implement a robust scheduled backup system. This follows the Performance Monitoring implementation and is a critical component for ensuring data safety and disaster recovery capabilities in production.

## Objectives

1. Implement automated database backup system
2. Create secure backup storage and retention policies
3. Develop file storage backup procedures
4. Implement backup verification and validation
5. Create backup restoration procedures and testing
6. Provide comprehensive documentation and runbooks

## Current Status

- The project has a production-ready database with proper RLS policies
- Blue/green deployment is implemented with zero-downtime capabilities
- Enhanced logging system is in place
- Performance monitoring implementation is in progress

## Implementation Requirements

### 1. Database Backup System

Create a comprehensive PostgreSQL backup system that includes:

- **Regular full backups**: Scheduled complete database dumps
- **Incremental backups**: Smaller, more frequent backups between full backups
- **Transaction log backups**: Continuous backup of database transaction logs
- **Point-in-time recovery capability**: Ability to restore to a specific timestamp
- **Backup compression**: Efficient storage of backup files
- **Backup encryption**: Protection of sensitive data in backups

Key implementation files:
- `scripts/backup/database_backup.py`: Main backup script
- `scripts/backup/backup_scheduler.py`: Scheduling component
- `scripts/backup/backup_verification.py`: Verification utilities
- `scripts/backup/backup_restore.py`: Restoration procedures
- `config/backup_config.py`: Backup configuration settings

### 2. Secure Backup Storage

Implement secure storage for backup files with:

- **Multiple storage locations**: Local and remote (cloud) storage
- **Access controls**: Proper permissions and authentication for backup access
- **Encryption at rest**: Encryption of backup files in storage
- **Storage monitoring**: Tracking of storage usage and capacity
- **Retention policies**: Automated cleanup of old backups based on configurable rules

Key implementation files:
- `scripts/backup/storage_manager.py`: Storage management
- `scripts/backup/cloud_storage.py`: Cloud storage integration
- `scripts/backup/retention_policy.py`: Backup retention implementation
- `config/storage_config.py`: Storage configuration settings

### 3. File Storage Backup

Create backup procedures for application file storage:

- **Static file backups**: Backup of application static files
- **User-uploaded content**: Backup of user-generated content
- **Configuration files**: Backup of critical configuration
- **Synchronization mechanisms**: Tools to keep file backups in sync
- **Version control integration**: Using git or similar for configuration

Key implementation files:
- `scripts/backup/file_backup.py`: File backup utilities
- `scripts/backup/sync_manager.py`: File synchronization
- `config/file_backup_config.py`: File backup configuration

### 4. Backup Verification and Validation

Implement thorough verification of backups:

- **Integrity checking**: Verification that backups are complete and uncorrupted
- **Restoration testing**: Scheduled test restores to verify backup validity
- **Automated validation**: Scripts to validate backup content and structure
- **Reporting**: Regular reports on backup status and validity
- **Alert integration**: Integration with monitoring system for backup failures

Key implementation files:
- `scripts/backup/verify_backup.py`: Verification utilities
- `scripts/backup/test_restore.py`: Automated restoration testing
- `scripts/backup/validation_report.py`: Reporting utilities

### 5. Restoration Procedures

Create comprehensive restoration procedures:

- **Quick restore process**: Streamlined process for common restore scenarios
- **Partial restore capabilities**: Ability to restore specific databases, tables, or files
- **Disaster recovery procedures**: Complete recovery from catastrophic failure
- **Testing environment**: Sandbox for testing restores without affecting production
- **Documentation**: Clear, step-by-step restoration guides

Key implementation files:
- `scripts/backup/restore_database.py`: Database restoration utilities
- `scripts/backup/restore_files.py`: File restoration utilities
- `scripts/backup/disaster_recovery.py`: Full system recovery procedures
- `docs/restore_guide.md`: Restoration documentation

### 6. CI/CD Integration

Integrate backup system with existing CI/CD pipeline:

- **Automated testing**: Tests for backup and restore functionality
- **Pre-deployment backups**: Automatic backups before deployments
- **Deployment verification**: Verify backup system after deployments
- **Pipeline notifications**: Notifications for backup status in CI/CD pipeline

## Detailed Implementation Plan

### 1. Database Backup Script

```python
# scripts/backup/database_backup.py

import os
import sys
import logging
import subprocess
import datetime
import shutil
import gzip
import argparse
from pathlib import Path
import yaml
import boto3
from botocore.exceptions import ClientError

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

def perform_pg_dump(db_config, backup_file):
    """Perform PostgreSQL database dump."""
    try:
        # Build pg_dump command
        cmd = [
            "pg_dump",
            "-h", db_config.get("host", "localhost"),
            "-p", str(db_config.get("port", 5432)),
            "-U", db_config.get("username"),
            "-d", db_config.get("database"),
            "-F", "c",  # Custom format
            "-Z", "9",  # Compression level
            "-f", backup_file
        ]
        
        # Set PGPASSWORD environment variable
        env = os.environ.copy()
        env["PGPASSWORD"] = db_config.get("password")
        
        # Execute pg_dump
        logger.info(f"Starting database backup to: {backup_file}")
        subprocess.run(cmd, env=env, check=True, capture_output=True)
        logger.info(f"Database backup completed successfully: {backup_file}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"pg_dump error: {e.stderr.decode()}")
        return False
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
        else:
            message = f"Backup failed: {error}"
            # Add notification logic here (email, Slack, etc.)
        
        logger.info(f"Notification sent: {message}")
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
    if not perform_pg_dump(db_config, backup_path):
        notify_backup_status(False, error="Database dump failed")
        return
    
    # Compress backup
    compressed_path = compress_backup(backup_path)
    if not compressed_path:
        notify_backup_status(False, error="Backup compression failed")
        return
    
    # Create metadata
    create_backup_metadata(compressed_path, db_config, args.type)
    
    # Upload to S3 if configured
    if config.get("aws", {}).get("enabled", False):
        s3_bucket = config.get("aws", {}).get("bucket")
        s3_prefix = config.get("aws", {}).get("prefix", "database-backups")
        s3_key = f"{s3_prefix}/{os.path.basename(compressed_path)}"
        
        if not upload_to_s3(compressed_path, s3_bucket, s3_key):
            notify_backup_status(False, error="S3 upload failed")
            return
    
    # Notify success
    notify_backup_status(True, backup_path=compressed_path)
    
    logger.info("Backup process completed successfully")

if __name__ == "__main__":
    main()
```

### 2. Backup Scheduler

```python
# scripts/backup/backup_scheduler.py

import os
import sys
import logging
import time
import schedule
import subprocess
import yaml
import argparse
from datetime import datetime, timedelta

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/backup_scheduler.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("backup_scheduler")

# Load configuration
def load_config(config_path="config/backup_config.yaml"):
    try:
        with open(config_path, 'r') as config_file:
            return yaml.safe_load(config_file)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def run_backup_job(backup_type="full"):
    """Run a backup job as a subprocess."""
    try:
        cmd = [
            sys.executable,
            "scripts/backup/database_backup.py",
            f"--type={backup_type}"
        ]
        
        logger.info(f"Starting {backup_type} backup job")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"{backup_type.capitalize()} backup completed successfully")
        else:
            logger.error(f"{backup_type.capitalize()} backup failed: {process.stderr}")
            
        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error running backup job: {str(e)}")
        return False

def run_backup_verification():
    """Run backup verification job."""
    try:
        cmd = [
            sys.executable,
            "scripts/backup/verify_backup.py"
        ]
        
        logger.info(f"Starting backup verification")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"Backup verification completed successfully")
        else:
            logger.error(f"Backup verification failed: {process.stderr}")
            
        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error running backup verification: {str(e)}")
        return False

def run_test_restore():
    """Run test restoration job."""
    try:
        cmd = [
            sys.executable,
            "scripts/backup/test_restore.py"
        ]
        
        logger.info(f"Starting test restoration")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"Test restoration completed successfully")
        else:
            logger.error(f"Test restoration failed: {process.stderr}")
            
        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error running test restoration: {str(e)}")
        return False

def run_retention_policy():
    """Apply backup retention policy."""
    try:
        cmd = [
            sys.executable,
            "scripts/backup/retention_policy.py"
        ]
        
        logger.info(f"Applying backup retention policy")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"Retention policy applied successfully")
        else:
            logger.error(f"Retention policy application failed: {process.stderr}")
            
        return process.returncode == 0
    except Exception as e:
        logger.error(f"Error applying retention policy: {str(e)}")
        return False

def schedule_jobs(config):
    """Schedule all backup jobs based on configuration."""
    # Schedule full backups
    full_backup_schedule = config.get("schedules", {}).get("full_backup", "0 1 * * *")  # Default: 1 AM daily
    schedule.every().day.at(full_backup_schedule.split(" ")[1] + ":00").do(run_backup_job, backup_type="full")
    logger.info(f"Scheduled full backup: {full_backup_schedule}")
    
    # Schedule incremental backups
    incremental_backup_schedule = config.get("schedules", {}).get("incremental_backup", "0 */6 * * *")  # Default: Every 6 hours
    schedule.every(6).hours.do(run_backup_job, backup_type="incremental")
    logger.info(f"Scheduled incremental backup: {incremental_backup_schedule}")
    
    # Schedule backup verification
    verification_schedule = config.get("schedules", {}).get("verification", "0 3 * * *")  # Default: 3 AM daily
    schedule.every().day.at(verification_schedule.split(" ")[1] + ":00").do(run_backup_verification)
    logger.info(f"Scheduled backup verification: {verification_schedule}")
    
    # Schedule test restore (weekly)
    test_restore_schedule = config.get("schedules", {}).get("test_restore", "0 4 * * 0")  # Default: 4 AM on Sundays
    schedule.every().sunday.at("04:00").do(run_test_restore)
    logger.info(f"Scheduled test restore: {test_restore_schedule}")
    
    # Schedule retention policy
    retention_schedule = config.get("schedules", {}).get("retention", "0 2 * * *")  # Default: 2 AM daily
    schedule.every().day.at(retention_schedule.split(" ")[1] + ":00").do(run_retention_policy)
    logger.info(f"Scheduled retention policy: {retention_schedule}")

def main():
    """Main function to run the backup scheduler."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Scheduler")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Schedule backup jobs
    schedule_jobs(config)
    
    logger.info("Backup scheduler started")
    
    # Run the scheduler
    while True:
        schedule.run_pending()
        time.sleep(60)  # Check every minute

if __name__ == "__main__":
    main()
```

### 3. Backup Configuration

```yaml
# config/backup_config.yaml

# Database configuration
database:
  host: localhost
  port: 5432
  username: postgres
  password: ${POSTGRES_PASSWORD}  # Use environment variable
  database: cryoprotect

# Backup storage paths
local_backup_dir: backups/database
file_backup_dir: backups/files
config_backup_dir: backups/config

# Retention policy (in days)
retention:
  full_backup: 30     # Keep full backups for 30 days
  incremental_backup: 7  # Keep incremental backups for 7 days
  file_backup: 30     # Keep file backups for 30 days
  config_backup: 90   # Keep config backups for 90 days

# Backup schedules (cron format)
schedules:
  full_backup: "0 1 * * *"        # 1 AM daily
  incremental_backup: "0 */6 * * *"  # Every 6 hours
  file_backup: "0 2 * * *"        # 2 AM daily
  config_backup: "0 3 * * 0"      # 3 AM on Sundays
  verification: "0 4 * * *"       # 4 AM daily
  test_restore: "0 5 * * 0"       # 5 AM on Sundays
  retention: "0 6 * * *"          # 6 AM daily

# AWS S3 configuration
aws:
  enabled: false
  access_key: ${AWS_ACCESS_KEY}  # Use environment variable
  secret_key: ${AWS_SECRET_KEY}  # Use environment variable
  region: us-east-1
  bucket: cryoprotect-backups
  prefix: database-backups

# Notification configuration
notifications:
  email:
    enabled: true
    from: backup@example.com
    to: admin@example.com
    smtp_server: smtp.example.com
    smtp_port: 587
    smtp_username: ${SMTP_USERNAME}  # Use environment variable
    smtp_password: ${SMTP_PASSWORD}  # Use environment variable
  
  slack:
    enabled: false
    webhook_url: ${SLACK_WEBHOOK}  # Use environment variable
    channel: #backups

# Compression configuration
compression:
  enabled: true
  level: 9  # Maximum compression

# Encryption configuration
encryption:
  enabled: true
  key_file: config/backup_encryption_key
```

### 4. Backup Verification Script

```python
# scripts/backup/verify_backup.py

import os
import sys
import logging
import subprocess
import json
import datetime
import argparse
import yaml
import hashlib
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/backup_verification.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("backup_verification")

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

def check_file_integrity(backup_file):
    """Check the integrity of the backup file."""
    try:
        # Check if file exists
        if not os.path.exists(backup_file):
            logger.error(f"Backup file does not exist: {backup_file}")
            return False
        
        # Check if file is not empty
        file_size = os.path.getsize(backup_file)
        if file_size == 0:
            logger.error(f"Backup file is empty: {backup_file}")
            return False
        
        # Check if file is readable
        try:
            with open(backup_file, 'rb') as f:
                f.read(1024)  # Read first 1KB to check if file is readable
        except Exception as e:
            logger.error(f"Backup file is not readable: {str(e)}")
            return False
        
        logger.info(f"Backup file integrity check passed: {backup_file}")
        return True
    except Exception as e:
        logger.error(f"Error checking file integrity: {str(e)}")
        return False

def check_backup_format(backup_file):
    """Check if the backup file has the correct format."""
    try:
        if backup_file.endswith('.gz'):
            # Check if it's a valid gzip file
            cmd = ["gzip", "-t", backup_file]
            process = subprocess.run(cmd, capture_output=True)
            
            if process.returncode != 0:
                logger.error(f"Invalid gzip file: {backup_file}")
                return False
        
        logger.info(f"Backup format check passed: {backup_file}")
        return True
    except Exception as e:
        logger.error(f"Error checking backup format: {str(e)}")
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

def verify_backup_metadata(backup_file):
    """Verify backup metadata file."""
    try:
        metadata_file = f"{backup_file}.meta.json"
        
        if not os.path.exists(metadata_file):
            logger.error(f"Metadata file does not exist: {metadata_file}")
            return False
        
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        # Check if metadata contains required fields
        required_fields = ["timestamp", "database", "backup_type", "backup_file"]
        for field in required_fields:
            if field not in metadata:
                logger.error(f"Missing required field in metadata: {field}")
                return False
        
        # Check file size
        if "size_bytes" in metadata:
            actual_size = os.path.getsize(backup_file)
            if metadata["size_bytes"] != actual_size:
                logger.error(f"File size mismatch: expected {metadata['size_bytes']}, got {actual_size}")
                return False
        
        logger.info(f"Backup metadata verification passed: {metadata_file}")
        return True
    except Exception as e:
        logger.error(f"Error verifying backup metadata: {str(e)}")
        return False

def create_verification_report(verified_backups):
    """Create a verification report."""
    try:
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "verified_backups": verified_backups,
            "total_backups": len(verified_backups),
            "passed": sum(1 for b in verified_backups if b["status"] == "passed"),
            "failed": sum(1 for b in verified_backups if b["status"] == "failed")
        }
        
        report_dir = "reports/backup_verification"
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f"verification_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Verification report created: {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error creating verification report: {str(e)}")
        return None

def notify_verification_status(report):
    """Send notification about verification status."""
    try:
        if report["failed"] > 0:
            message = f"Backup verification failed for {report['failed']} backups"
            # Add notification logic here (email, Slack, etc.)
        else:
            message = f"All {report['total_backups']} backups verified successfully"
            # Add notification logic here (email, Slack, etc.)
        
        logger.info(f"Notification sent: {message}")
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")

def main():
    """Main verification function."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Verification Tool")
    parser.add_argument("--config", default="config/backup_config.yaml", help="Configuration file path")
    parser.add_argument("--backup-dir", help="Backup directory to verify (overrides config)")
    parser.add_argument("--file", help="Specific backup file to verify")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Determine backup directory
    backup_dir = args.backup_dir or config.get("local_backup_dir", "backups/database")
    
    # Get backup files to verify
    if args.file:
        backup_files = [args.file]
    else:
        backup_files = list_backup_files(backup_dir)
    
    if not backup_files:
        logger.warning(f"No backup files found in {backup_dir}")
        return
    
    logger.info(f"Found {len(backup_files)} backup files to verify")
    
    # Verify each backup file
    verified_backups = []
    for backup_file in backup_files:
        backup_name = os.path.basename(backup_file)
        logger.info(f"Verifying backup: {backup_name}")
        
        verification_results = {
            "file": backup_name,
            "path": backup_file,
            "timestamp": datetime.datetime.now().isoformat(),
            "checks": {}
        }
        
        # Perform verification checks
        verification_results["checks"]["file_integrity"] = check_file_integrity(backup_file)
        verification_results["checks"]["backup_format"] = check_backup_format(backup_file)
        verification_results["checks"]["metadata"] = verify_backup_metadata(backup_file)
        
        # Calculate overall status
        if all(verification_results["checks"].values()):
            verification_results["status"] = "passed"
            logger.info(f"Backup verification passed: {backup_name}")
        else:
            verification_results["status"] = "failed"
            logger.error(f"Backup verification failed: {backup_name}")
        
        verified_backups.append(verification_results)
    
    # Create verification report
    report_file = create_verification_report(verified_backups)
    
    # Load report and notify status
    if report_file:
        with open(report_file, 'r') as f:
            report = json.load(f)
        notify_verification_status(report)
    
    logger.info("Backup verification process completed")

if __name__ == "__main__":
    main()
```

### 5. Backup Restore Script

```python
# scripts/backup/restore_database.py

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
```

### 6. Docker Compose Configuration

```yaml
# docker-compose.backup.yml
version: '3.8'

services:
  backup-scheduler:
    image: python:3.9-slim
    container_name: cryoprotect-backup-scheduler
    volumes:
      - ./:/app
      - backup-data:/backups
    working_dir: /app
    environment:
      - POSTGRES_PASSWORD=${POSTGRES_PASSWORD}
      - AWS_ACCESS_KEY=${AWS_ACCESS_KEY}
      - AWS_SECRET_KEY=${AWS_SECRET_KEY}
      - SMTP_USERNAME=${SMTP_USERNAME}
      - SMTP_PASSWORD=${SMTP_PASSWORD}
      - SLACK_WEBHOOK=${SLACK_WEBHOOK}
    command: python scripts/backup/backup_scheduler.py
    restart: unless-stopped
    depends_on:
      - postgres
    networks:
      - cryoprotect-network
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"

volumes:
  backup-data:
    driver: local

networks:
  cryoprotect-network:
    external: true
```

### 7. GitHub Actions Integration

```yaml
# .github/workflows/backup-verification.yml
name: Backup Verification

on:
  schedule:
    - cron: '0 6 * * *'  # Run at 6 AM UTC every day
  workflow_dispatch:  # Allow manual triggering

jobs:
  verify-backups:
    name: Verify Database Backups
    runs-on: ubuntu-latest
    
    steps:
      - name: Checkout code
        uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'
      
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install pyyaml boto3 psycopg2-binary
      
      - name: Download latest backup from S3
        env:
          AWS_ACCESS_KEY_ID: ${{ secrets.AWS_ACCESS_KEY_ID }}
          AWS_SECRET_ACCESS_KEY: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          AWS_DEFAULT_REGION: ${{ secrets.AWS_REGION }}
          S3_BUCKET: ${{ secrets.BACKUP_S3_BUCKET }}
        run: |
          mkdir -p backups/database
          python - <<EOF
          import boto3
          import os
          from datetime import datetime
          
          s3 = boto3.client('s3')
          bucket = os.environ.get('S3_BUCKET')
          prefix = 'database-backups/'
          
          # List all backup objects
          response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
          
          if 'Contents' not in response:
              print("No backups found in S3")
              exit(1)
              
          # Get the latest backup
          latest_backup = sorted(response['Contents'], key=lambda x: x['LastModified'], reverse=True)[0]
          backup_key = latest_backup['Key']
          local_path = f"backups/database/{os.path.basename(backup_key)}"
          
          print(f"Downloading latest backup: {backup_key} to {local_path}")
          s3.download_file(bucket, backup_key, local_path)
          
          print(f"Downloaded backup: {local_path}")
          with open(os.environ.get('GITHUB_ENV'), 'a') as env_file:
              env_file.write(f"BACKUP_FILE={local_path}\n")
          EOF
      
      - name: Verify backup integrity
        run: |
          python scripts/backup/verify_backup.py --file ${{ env.BACKUP_FILE }}
      
      - name: Upload verification report
        uses: actions/upload-artifact@v3
        with:
          name: backup-verification-report
          path: reports/backup_verification/
          retention-days: 30
      
      - name: Send notification
        if: always()
        uses: rtCamp/action-slack-notify@v2
        env:
          SLACK_WEBHOOK: ${{ secrets.SLACK_WEBHOOK }}
          SLACK_CHANNEL: backups
          SLACK_COLOR: ${{ job.status }}
          SLACK_TITLE: Backup Verification
          SLACK_MESSAGE: "Backup verification ${{ job.status }}"
          SLACK_FOOTER: "CryoProtect v2 Backup System"
```

## Runbook Documentation

Create detailed operational documentation for the backup system:

```markdown
# CryoProtect v2 Backup System Runbook

This document provides operational procedures for managing the CryoProtect v2 backup system.

## 1. Backup System Overview

The CryoProtect v2 backup system provides:

- Automated database backups (full and incremental)
- File storage backups
- Configuration backups
- Backup verification
- Test restoration
- Backup retention policy enforcement

## 2. Checking Backup Status

To check the status of recent backups:

```bash
# View backup logs
tail -n 100 logs/backup_scheduler.log

# List recent backups
ls -la backups/database/ | sort -r | head -10

# View verification reports
ls -la reports/backup_verification/ | sort -r | head -5
cat reports/backup_verification/[latest-report].json
```

## 3. Manual Backup Procedures

To manually trigger a backup:

```bash
# Trigger full backup
python scripts/backup/database_backup.py --type=full

# Trigger incremental backup
python scripts/backup/database_backup.py --type=incremental

# Trigger file backup
python scripts/backup/file_backup.py

# Trigger configuration backup
python scripts/backup/config_backup.py
```

## 4. Restoration Procedures

### 4.1 Full Database Restoration

```bash
# Restore from latest full backup
python scripts/backup/restore_database.py --clean

# Restore from specific backup file
python scripts/backup/restore_database.py --file backups/database/[backup-file].dump.gz --clean
```

### 4.2 Partial Database Restoration

```bash
# Restore specific schema
python scripts/backup/restore_database.py --schema public_schema

# Restore specific table
python scripts/backup/restore_database.py --table users
```

### 4.3 Testing Restoration

```bash
# Test restore without actually modifying the database
python scripts/backup/restore_database.py --test-mode
```

## 5. Backup Verification

```bash
# Verify all backups
python scripts/backup/verify_backup.py

# Verify specific backup
python scripts/backup/verify_backup.py --file backups/database/[backup-file].dump.gz
```

## 6. Troubleshooting

### 6.1 Backup Failures

If backups are failing:

1. Check database connectivity:
   ```bash
   python -c "import psycopg2; conn = psycopg2.connect('dbname=cryoprotect user=postgres host=localhost')"
   ```

2. Check disk space:
   ```bash
   df -h
   ```

3. Check backup logs for errors:
   ```bash
   grep ERROR logs/database_backup.log
   ```

### 6.2 Restoration Failures

If restores are failing:

1. Verify backup file integrity:
   ```bash
   python scripts/backup/verify_backup.py --file [backup-file]
   ```

2. Check database permissions:
   ```bash
   psql -U postgres -c "SELECT rolname, rolcreaterole, rolcreatedb FROM pg_roles WHERE rolname = 'postgres'"
   ```

3. Check restore logs for errors:
   ```bash
   grep ERROR logs/database_restore.log
   ```

## 7. Modifying Backup Configuration

To modify the backup configuration:

1. Edit the backup configuration file:
   ```bash
   nano config/backup_config.yaml
   ```

2. After editing, validate the configuration:
   ```bash
   python -c "import yaml; yaml.safe_load(open('config/backup_config.yaml'))"
   ```

3. Restart the backup scheduler:
   ```bash
   docker-compose -f docker-compose.yml -f docker-compose.backup.yml restart backup-scheduler
   ```

## 8. Emergency Procedures

### 8.1 Emergency Restoration

In case of catastrophic database failure:

1. Stop the application:
   ```bash
   docker-compose down
   ```

2. Restore from the latest full backup:
   ```bash
   python scripts/backup/restore_database.py --clean --create
   ```

3. Apply any incremental backups if needed:
   ```bash
   python scripts/backup/restore_database.py --file backups/database/[incremental-backup].dump.gz
   ```

4. Restart the application:
   ```bash
   docker-compose up -d
   ```

### 8.2 Emergency Manual Backup

To create an emergency backup before critical operations:

```bash
# Create emergency full backup
python scripts/backup/database_backup.py --type=full --config config/emergency_backup_config.yaml
```

## 9. Monitoring and Alerting

The backup system is monitored through:

1. Daily backup verification reports
2. Integration with the monitoring system
3. Email/Slack notifications for backup failures
4. GitHub Actions workflow for verification

To check backup alert status:

```bash
# View recent alerts
grep "Notification sent" logs/backup_scheduler.log | tail -10
```
```

## Implementation Plan

### Week 1: Core Backup System

1. Day 1-2: Set up database backup system
   - Implement database_backup.py
   - Create backup configuration
   - Test backup functionality

2. Day 3-4: Implement verification and restoration
   - Implement verify_backup.py
   - Implement restore_database.py
   - Test verification and restoration

3. Day 5: Implement scheduling and retention
   - Implement backup_scheduler.py
   - Implement retention_policy.py
   - Set up Docker container for scheduler

### Week 2: Extended Functionality and Integration

1. Day 1-2: Implement file and configuration backups
   - Implement file_backup.py
   - Implement config_backup.py
   - Test all backup types

2. Day 3-4: Set up cloud storage integration
   - Implement S3 storage functionality
   - Set up secure credential management
   - Test cloud storage uploads/downloads

3. Day 5: Create CI/CD integration
   - Create GitHub Actions workflow
   - Implement pre-deployment backup
   - Test CI/CD integration

### Week 3: Documentation and Testing

1. Day 1-2: Create comprehensive documentation
   - Write runbook documentation
   - Create restoration guides
   - Document configuration options

2. Day 3-4: Perform thorough testing
   - Test backup/restore under various conditions
   - Verify backup integrity checking
   - Test disaster recovery procedures

3. Day 5: Finalize and deploy
   - Address any issues from testing
   - Deploy to production environment
   - Verify production backups

## Success Criteria

The scheduled backup implementation will be considered successful when:

1. Daily automated backups are running reliably
2. Backup verification confirms backup integrity
3. Test restorations complete successfully
4. Multiple storage locations are properly utilized
5. Retention policies automatically manage backup storage
6. Backup status is properly monitored and reported
7. Documentation is comprehensive and clear
8. Operations team can manage the backup system effectively

## Next Steps After Implementation

After completing the scheduled backup system, the project should move on to:

1. Security implementation (Phase 3.3)
2. Documentation completion (Phase 4.1)
3. Knowledge transfer (Phase 4.2)
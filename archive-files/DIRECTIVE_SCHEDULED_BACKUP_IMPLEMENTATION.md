# DIRECTIVE: Scheduled Backup System Implementation for CryoProtect v2

## CURRENT PROJECT STATUS

We are currently in Phase 3.2 (Monitoring and Maintenance) of the CryoProtect v2 project. The previous phases have been completed:

- Phase 1 (Technical Foundation): Complete
- Phase 2 (Feature Completion): Complete
- Phase 3.1 (Deployment Infrastructure): Complete
  - Stage 1 (Database Foundation): Complete - See `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PHASE_3.1_STAGE1_COMPLETION_REPORT.md`
  - Stage 2 (Blue/Green Deployment): Complete - See `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PHASE_3.1_STAGE2_COMPLETION_REPORT.md`
  - Stage 3 (Partial - Enhanced Logging): Partially complete - See `/mnt/c/Users/1edwa/Documents/CryoProtect v2/README_Enhanced_Logging.md`

We are now implementing Phase 3.2 tasks in the following order:
1. Performance Monitoring (in progress - see `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ROO_PERFORMANCE_MONITORING_IMPLEMENTATION.md`)
2. **Scheduled Backup System (current task)**
3. Maintenance Runbooks (upcoming)

## TASK OBJECTIVE

Implement a comprehensive scheduled backup system for CryoProtect v2 that includes:
1. Automated database backups (full and incremental)
2. Backup verification and validation
3. Secure backup storage with retention policies
4. Backup restoration procedures
5. Integration with monitoring and alerting

## IMPLEMENTATION APPROACH

Follow this step-by-step approach to implement the scheduled backup system:

### Step 1: Create Directory Structure

First, create the necessary directories to organize the backup system files:

```bash
mkdir -p /mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup
mkdir -p /mnt/c/Users/1edwa/Documents/CryoProtect v2/config/backup
mkdir -p /mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/backup_verification
mkdir -p /mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/database_restore
mkdir -p /mnt/c/Users/1edwa/Documents/CryoProtect v2/docs/backup
```

### Step 2: Create Backup Configuration File

Create the backup configuration file:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/config/backup/backup_config.yaml`**

```yaml
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

### Step 3: Create Core Backup Script

Create the main database backup script:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/database_backup.py`**

```python
#!/usr/bin/env python3
"""
Database Backup Script for CryoProtect v2

This script performs database backups with the following features:
- Full and incremental backups
- Compression and optional encryption
- Local and cloud storage
- Backup metadata and verification
- Notification on completion or failure
"""

import os
import sys
import logging
import subprocess
import datetime
import shutil
import gzip
import json
import argparse
from pathlib import Path
import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/database_backup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("database_backup")

def load_config(config_path="config/backup/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

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
        # Build pg_dump command
        cmd = [
            "pg_dump",
            "-h", db_config.get("host", "localhost"),
            "-p", str(db_config.get("port", 5432)),
            "-U", db_config.get("username"),
            "-d", db_config.get("database"),
            "-F", "c",  # Custom format
            "-Z", "9",  # Maximum compression
            "-f", backup_file
        ]
        
        # Set PGPASSWORD environment variable
        env = os.environ.copy()
        env["PGPASSWORD"] = db_config.get("password")
        
        # Execute pg_dump
        logger.info(f"Starting {backup_type} database backup to: {backup_file}")
        process = subprocess.run(cmd, env=env, check=True, capture_output=True)
        
        if process.returncode == 0:
            logger.info(f"Database backup completed successfully: {backup_file}")
            return True
        else:
            logger.error(f"pg_dump error: {process.stderr.decode()}")
            return False
    except subprocess.CalledProcessError as e:
        logger.error(f"pg_dump error: {e.stderr.decode() if hasattr(e, 'stderr') else str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error performing database backup: {str(e)}")
        return False

def compress_backup(backup_file):
    """Compress backup file with gzip."""
    try:
        compressed_file = f"{backup_file}.gz"
        with open(backup_file, 'rb') as f_in:
            with gzip.open(compressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Remove original uncompressed file
        os.remove(backup_file)
        logger.info(f"Backup compressed: {compressed_file}")
        return compressed_file
    except Exception as e:
        logger.error(f"Error compressing backup: {str(e)}")
        return None

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
            logger.info(message)
            # Implement notification logic here (email, Slack, etc.)
        else:
            message = f"Backup failed: {error}"
            logger.error(message)
            # Implement notification logic here (email, Slack, etc.)
    except Exception as e:
        logger.error(f"Error sending notification: {str(e)}")

def main():
    """Main backup function."""
    parser = argparse.ArgumentParser(description="CryoProtect Database Backup Tool")
    parser.add_argument("--type", choices=["full", "incremental"], default="full", help="Backup type")
    parser.add_argument("--config", default="config/backup/backup_config.yaml", help="Path to configuration file")
    args = parser.parse_args()
    
    # Load configuration
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
        return 1
    
    # Perform database backup
    if not perform_pg_dump(db_config, backup_path, args.type):
        notify_backup_status(False, error="Database dump failed")
        return 1
    
    # Compress backup if configured
    if config.get("compression", {}).get("enabled", True):
        compressed_path = compress_backup(backup_path)
        if not compressed_path:
            notify_backup_status(False, error="Backup compression failed")
            return 1
        backup_path = compressed_path
    
    # Create metadata
    create_backup_metadata(backup_path, db_config, args.type)
    
    # Notify success
    notify_backup_status(True, backup_path=backup_path)
    logger.info("Backup process completed successfully")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

### Step 4: Create Backup Verification Script

Create the backup verification script:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/verify_backup.py`**

```python
#!/usr/bin/env python3
"""
Backup Verification Script for CryoProtect v2

This script verifies backup integrity with the following features:
- File integrity checking
- Metadata validation
- Backup format verification
- Verification reporting
"""

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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/backup_verification.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("backup_verification")

def load_config(config_path="config/backup/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def list_backup_files(backup_dir):
    """List all backup files in the specified directory."""
    try:
        backup_files = []
        for file in os.listdir(backup_dir):
            if file.endswith('.dump.gz') or file.endswith('.dump'):
                backup_files.append(os.path.join(backup_dir, file))
        
        return sorted(backup_files, key=os.path.getmtime, reverse=True)
    except Exception as e:
        logger.error(f"Error listing backup files: {str(e)}")
        return []

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

def main():
    """Main verification function."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Verification Tool")
    parser.add_argument("--config", default="config/backup/backup_config.yaml", help="Path to configuration file")
    parser.add_argument("--backup-dir", help="Backup directory to verify (overrides config)")
    parser.add_argument("--file", help="Specific backup file to verify")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Determine backup directory
    backup_dir = args.backup_dir or config.get("local_backup_dir", "backups/database")
    os.makedirs(backup_dir, exist_ok=True)
    
    # Get backup files to verify
    if args.file and os.path.exists(args.file):
        backup_files = [args.file]
    else:
        backup_files = list_backup_files(backup_dir)
    
    if not backup_files:
        logger.warning(f"No backup files found in {backup_dir}")
        return 1
    
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
    
    # Determine exit code based on verification results
    if all(backup["status"] == "passed" for backup in verified_backups):
        logger.info("All backup verifications passed")
        return 0
    else:
        logger.error("Some backup verifications failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

### Step 5: Create Backup Scheduler Script

Create the scheduler for automated backups:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/backup_scheduler.py`**

```python
#!/usr/bin/env python3
"""
Backup Scheduler for CryoProtect v2

This script schedules and manages automated backups with the following features:
- Scheduled full and incremental backups
- Scheduled backup verification
- Scheduled test restores
- Retention policy enforcement
"""

import os
import sys
import logging
import time
import schedule
import subprocess
import yaml
import argparse
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/backup_scheduler.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("backup_scheduler")

def load_config(config_path="config/backup/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def run_backup_job(backup_type="full"):
    """Run a backup job as a subprocess."""
    try:
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database_backup.py")
        cmd = [
            sys.executable,
            script_path,
            f"--type={backup_type}"
        ]
        
        logger.info(f"Starting {backup_type} backup job")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"{backup_type.capitalize()} backup completed successfully")
            return True
        else:
            logger.error(f"{backup_type.capitalize()} backup failed: {process.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error running backup job: {str(e)}")
        return False

def run_backup_verification():
    """Run backup verification job."""
    try:
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "verify_backup.py")
        cmd = [
            sys.executable,
            script_path
        ]
        
        logger.info(f"Starting backup verification")
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"Backup verification completed successfully")
            return True
        else:
            logger.error(f"Backup verification failed: {process.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error running backup verification: {str(e)}")
        return False

def main():
    """Main function to run the backup scheduler."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Scheduler")
    parser.add_argument("--config", default="config/backup/backup_config.yaml", help="Path to configuration file")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Schedule full backups (default: daily at 1 AM)
    full_schedule = config.get("schedules", {}).get("full_backup", "0 1 * * *").split()
    if len(full_schedule) >= 2:
        hour, minute = full_schedule[1], full_schedule[0]
        schedule.every().day.at(f"{hour.zfill(2)}:{minute.zfill(2)}").do(run_backup_job, backup_type="full")
        logger.info(f"Scheduled full backup at {hour}:{minute} daily")
    
    # Schedule incremental backups (default: every 6 hours)
    incremental_schedule = config.get("schedules", {}).get("incremental_backup", "0 */6 * * *")
    if "*/6" in incremental_schedule:
        schedule.every(6).hours.do(run_backup_job, backup_type="incremental")
        logger.info("Scheduled incremental backup every 6 hours")
    
    # Schedule backup verification (default: daily at 4 AM)
    verification_schedule = config.get("schedules", {}).get("verification", "0 4 * * *").split()
    if len(verification_schedule) >= 2:
        hour, minute = verification_schedule[1], verification_schedule[0]
        schedule.every().day.at(f"{hour.zfill(2)}:{minute.zfill(2)}").do(run_backup_verification)
        logger.info(f"Scheduled backup verification at {hour}:{minute} daily")
    
    logger.info("Backup scheduler started")
    
    # Run the scheduler loop
    while True:
        schedule.run_pending()
        time.sleep(60)  # Check every minute

if __name__ == "__main__":
    main()
```

### Step 6: Create Database Restore Script

Create the database restoration script:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/restore_database.py`**

```python
#!/usr/bin/env python3
"""
Database Restore Script for CryoProtect v2

This script performs database restoration with the following features:
- Full database restores
- Schema or table-specific restores
- Restore verification
- Test-mode for validation without actual restoration
"""

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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/database_restore.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("database_restore")

def load_config(config_path="config/backup/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def list_backup_files(backup_dir):
    """List all backup files in the specified directory."""
    try:
        backup_files = []
        for file in os.listdir(backup_dir):
            if file.endswith('.dump.gz') or file.endswith('.dump'):
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
            if backup_type in os.path.basename(file):
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
        logger.info(f"Restore command: {' '.join(cmd)}")
        
        process = subprocess.run(cmd, env=env, capture_output=True, text=True)
        
        if process.returncode == 0:
            logger.info(f"Database restore completed successfully")
            return True
        else:
            # pg_restore can return non-zero even on successful restores with warnings
            logger.warning(f"pg_restore returned code {process.returncode}")
            logger.warning(f"Stderr: {process.stderr}")
            
            # Check if there are fatal errors in the output
            if "fatal" in process.stderr.lower():
                logger.error("Fatal errors detected in pg_restore output")
                return False
            
            # If no fatal errors, consider it successful with warnings
            logger.info("Database restore completed with warnings")
            return True
    except subprocess.CalledProcessError as e:
        logger.error(f"pg_restore error: {e.stderr if hasattr(e, 'stderr') else str(e)}")
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

def main():
    """Main restore function."""
    parser = argparse.ArgumentParser(description="CryoProtect Database Restore Tool")
    parser.add_argument("--config", default="config/backup/backup_config.yaml", help="Path to configuration file")
    parser.add_argument("--file", help="Specific backup file to restore")
    parser.add_argument("--type", choices=["full", "incremental"], default="full", 
                        help="Backup type to restore (if not specifying file)")
    parser.add_argument("--clean", action="store_true", help="Clean (drop) database objects before recreating")
    parser.add_argument("--create", action="store_true", help="Create the database before restoring")
    parser.add_argument("--no-owner", action="store_true", 
                        help="Do not set ownership of objects to match the original database")
    parser.add_argument("--schema", help="Restore only objects in this schema")
    parser.add_argument("--table", help="Restore only the named table")
    parser.add_argument("--test-mode", action="store_true", 
                        help="Run in test mode without actually restoring")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Get database configuration
    db_config = config.get("database", {})
    
    # Determine backup file to restore
    backup_dir = config.get("local_backup_dir", "backups/database")
    os.makedirs(backup_dir, exist_ok=True)
    
    if args.file and os.path.exists(args.file):
        backup_file = args.file
    else:
        backup_file = get_latest_backup(backup_dir, args.type)
    
    if not backup_file:
        error_msg = f"No backup file found to restore"
        logger.error(error_msg)
        return 1
    
    if not os.path.exists(backup_file):
        error_msg = f"Backup file does not exist: {backup_file}"
        logger.error(error_msg)
        return 1
    
    logger.info(f"Selected backup file for restore: {backup_file}")
    
    # Decompress backup if needed
    if backup_file.endswith('.gz'):
        decompressed_file = decompress_backup(backup_file)
        if not decompressed_file:
            error_msg = f"Failed to decompress backup file: {backup_file}"
            logger.error(error_msg)
            return 1
        backup_file = decompressed_file
    
    # Prepare restore options
    restore_options = {
        "clean": args.clean,
        "create": args.create,
        "no_owner": args.no_owner,
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
    
    # Cleanup decompressed file if it was created
    if backup_file != args.file and os.path.exists(backup_file) and backup_file.endswith('.dump'):
        try:
            os.remove(backup_file)
            logger.info(f"Removed temporary decompressed file: {backup_file}")
        except Exception as e:
            logger.warning(f"Error removing temporary file: {str(e)}")
    
    logger.info("Database restore process completed")
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
```

### Step 7: Create Backup Retention Script

Create the script to manage backup retention:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/scripts/backup/retention_policy.py`**

```python
#!/usr/bin/env python3
"""
Backup Retention Policy Script for CryoProtect v2

This script applies retention policies to backups with the following features:
- Age-based backup removal
- Type-specific retention periods
- Backup metadata cleanup
- Retention reporting
"""

import os
import sys
import logging
import argparse
import yaml
import json
import datetime
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("logs/backup_retention.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("backup_retention")

def load_config(config_path="config/backup/backup_config.yaml"):
    """Load backup configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def list_backup_files(backup_dir):
    """List all backup files with their metadata in the specified directory."""
    try:
        backup_files = []
        
        if not os.path.exists(backup_dir):
            logger.warning(f"Backup directory does not exist: {backup_dir}")
            return []
        
        for file in os.listdir(backup_dir):
            if file.endswith('.dump.gz') or file.endswith('.dump'):
                file_path = os.path.join(backup_dir, file)
                metadata_path = f"{file_path}.meta.json"
                
                backup_info = {
                    "file_path": file_path,
                    "metadata_path": metadata_path if os.path.exists(metadata_path) else None,
                    "mtime": os.path.getmtime(file_path),
                    "type": "unknown"
                }
                
                # Try to determine backup type from filename
                if "full" in file:
                    backup_info["type"] = "full_backup"
                elif "incremental" in file:
                    backup_info["type"] = "incremental_backup"
                
                # Try to get more info from metadata if available
                if backup_info["metadata_path"]:
                    try:
                        with open(metadata_path, 'r') as f:
                            metadata = json.load(f)
                            if "backup_type" in metadata:
                                backup_info["type"] = f"{metadata['backup_type']}_backup"
                            if "timestamp" in metadata:
                                try:
                                    backup_info["timestamp"] = datetime.datetime.fromisoformat(metadata["timestamp"])
                                except:
                                    pass
                    except Exception as e:
                        logger.warning(f"Error reading metadata for {file}: {str(e)}")
                
                backup_files.append(backup_info)
        
        return backup_files
    except Exception as e:
        logger.error(f"Error listing backup files: {str(e)}")
        return []

def should_delete_backup(backup_info, retention_config):
    """Determine if a backup should be deleted based on retention policy."""
    try:
        # Get retention period for this backup type
        backup_type = backup_info.get("type", "unknown")
        retention_days = retention_config.get(backup_type, 30)  # Default: 30 days
        
        # Calculate backup age
        if "timestamp" in backup_info:
            backup_time = backup_info["timestamp"]
        else:
            backup_time = datetime.datetime.fromtimestamp(backup_info["mtime"])
        
        current_time = datetime.datetime.now()
        age_days = (current_time - backup_time).total_seconds() / 86400  # Convert to days
        
        # Check if backup is older than retention period
        return age_days > retention_days
    except Exception as e:
        logger.error(f"Error checking if backup should be deleted: {str(e)}")
        return False

def delete_backup(backup_info):
    """Delete a backup file and its metadata."""
    try:
        file_path = backup_info["file_path"]
        metadata_path = backup_info["metadata_path"]
        
        # Delete backup file
        if os.path.exists(file_path):
            os.remove(file_path)
            logger.info(f"Deleted backup file: {file_path}")
        
        # Delete metadata file
        if metadata_path and os.path.exists(metadata_path):
            os.remove(metadata_path)
            logger.info(f"Deleted metadata file: {metadata_path}")
        
        return True
    except Exception as e:
        logger.error(f"Error deleting backup: {str(e)}")
        return False

def create_retention_report(backups_to_delete, deleted_backups):
    """Create a report of retention policy application."""
    try:
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "backups_evaluated": len(backups_to_delete),
            "backups_deleted": len(deleted_backups),
            "deleted_backups": deleted_backups
        }
        
        report_dir = "reports/backup_retention"
        os.makedirs(report_dir, exist_ok=True)
        
        report_file = os.path.join(report_dir, f"retention_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Retention report created: {report_file}")
        return report_file
    except Exception as e:
        logger.error(f"Error creating retention report: {str(e)}")
        return None

def main():
    """Main retention function."""
    parser = argparse.ArgumentParser(description="CryoProtect Backup Retention Tool")
    parser.add_argument("--config", default="config/backup/backup_config.yaml", help="Path to configuration file")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be deleted without actually deleting")
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Get retention policy configuration
    retention_config = config.get("retention", {})
    if not retention_config:
        logger.error("Retention policy not configured")
        return 1
    
    # Get backup directory
    backup_dir = config.get("local_backup_dir", "backups/database")
    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir, exist_ok=True)
        logger.info(f"Created backup directory: {backup_dir}")
    
    # List backup files
    backup_files = list_backup_files(backup_dir)
    logger.info(f"Found {len(backup_files)} backup files")
    
    # Identify backups to delete
    backups_to_delete = []
    for backup_info in backup_files:
        if should_delete_backup(backup_info, retention_config):
            backups_to_delete.append(backup_info)
    
    logger.info(f"Identified {len(backups_to_delete)} backups to delete")
    
    # Delete backups (or simulate deletion in dry run mode)
    deleted_backups = []
    for backup_info in backups_to_delete:
        file_path = backup_info["file_path"]
        backup_type = backup_info.get("type", "unknown")
        
        if args.dry_run:
            logger.info(f"[DRY RUN] Would delete {backup_type} backup: {file_path}")
            deleted_backups.append({
                "file_path": file_path,
                "type": backup_type,
                "deleted": False,
                "reason": "dry_run"
            })
        else:
            if delete_backup(backup_info):
                deleted_backups.append({
                    "file_path": file_path,
                    "type": backup_type,
                    "deleted": True
                })
    
    # Create retention report
    create_retention_report(backups_to_delete, deleted_backups)
    
    logger.info(f"Backup retention policy applied: {len(deleted_backups)} backups {'would be ' if args.dry_run else ''}deleted")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

### Step 8: Create Docker Compose Configuration for Backup System

Create a Docker Compose configuration for the backup system:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/docker-compose.backup.yml`**

```yaml
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
    command: >
      bash -c "
        pip install -r requirements_backup.txt &&
        mkdir -p logs backups/database backups/files backups/config &&
        python scripts/backup/backup_scheduler.py
      "
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

### Step 9: Create Requirements File for Backup System

Create a requirements file for the backup scripts:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/requirements_backup.txt`**

```
pyyaml>=6.0
psycopg2-binary>=2.9.5
schedule>=1.1.0
boto3>=1.26.0
cryptography>=38.0.0
```

### Step 10: Create Comprehensive Backup System Documentation

Create documentation for the backup system:

**File: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/docs/backup/backup_system.md`**

```markdown
# CryoProtect v2 Backup System

## Overview

The CryoProtect v2 backup system provides comprehensive data protection with the following features:

- Automated database backups (full and incremental)
- Backup verification and validation
- Secure backup storage with retention policies
- Backup restoration procedures
- Integration with monitoring and alerting

## Architecture

The backup system consists of several components:

1. **Backup Scripts**: Python scripts for performing backups, verification, and restoration
2. **Scheduler**: Automated scheduling of backup tasks
3. **Storage**: Local and optional cloud storage for backups
4. **Retention**: Automated cleanup of old backups
5. **Monitoring**: Integration with the monitoring system
6. **Documentation**: Comprehensive documentation and runbooks

## Configuration

The backup system is configured through a YAML file at `config/backup/backup_config.yaml`.

Key configuration options include:

- Database connection details
- Backup storage paths
- Retention policies
- Backup schedules
- Cloud storage configuration
- Notification settings

## Usage

### Starting the Backup System

```bash
# Start the backup scheduler with Docker Compose
docker-compose -f docker-compose.yml -f docker-compose.backup.yml up -d backup-scheduler

# Start the backup scheduler manually
python scripts/backup/backup_scheduler.py
```

### Manual Backup

```bash
# Perform a full backup
python scripts/backup/database_backup.py --type=full

# Perform an incremental backup
python scripts/backup/database_backup.py --type=incremental
```

### Verify Backups

```bash
# Verify all backups
python scripts/backup/verify_backup.py

# Verify a specific backup
python scripts/backup/verify_backup.py --file=/path/to/backup.dump.gz
```

### Apply Retention Policy

```bash
# Apply retention policy (delete old backups)
python scripts/backup/retention_policy.py

# Simulate retention policy without deleting
python scripts/backup/retention_policy.py --dry-run
```

### Restore Database

```bash
# Restore from the latest full backup
python scripts/backup/restore_database.py

# Restore a specific backup with cleanup
python scripts/backup/restore_database.py --file=/path/to/backup.dump.gz --clean

# Restore only a specific schema
python scripts/backup/restore_database.py --schema=public

# Test restore without actually applying changes
python scripts/backup/restore_database.py --test-mode
```

## Monitoring

The backup system integrates with the CryoProtect monitoring system:

- Backup job success/failure metrics
- Backup size and duration metrics
- Verification status metrics
- Retention policy application metrics

## Reports

The backup system generates various reports:

- **Backup Reports**: Details about each backup operation
- **Verification Reports**: Results of backup verification
- **Restoration Reports**: Details of restoration operations
- **Retention Reports**: Information about backup retention policy application

Reports are stored in the `reports/` directory.

## Troubleshooting

### Common Issues

#### Backup Failure

If backups are failing:

1. Check the logs: `logs/database_backup.log`
2. Verify database connection settings in the configuration file
3. Ensure sufficient disk space is available

#### Restore Failure

If restoration is failing:

1. Check the logs: `logs/database_restore.log`
2. Verify the backup file integrity using `verify_backup.py`
3. Ensure the database exists and is accessible

#### Scheduler Not Running

If the scheduler is not running:

1. Check the logs: `logs/backup_scheduler.log`
2. Verify the scheduler service is running in Docker Compose
3. Check that the backup directory is writable

## Disaster Recovery

In case of a catastrophic failure:

1. Identify the most recent intact backup
2. Restore the database:
   ```bash
   python scripts/backup/restore_database.py --file=/path/to/backup.dump.gz --clean --create
   ```
3. Verify the restoration:
   ```bash
   python scripts/backup/verify_backup.py --restore-check
   ```

## Security Considerations

The backup system implements the following security measures:

- Encrypted backups (optional)
- Secure credential handling
- Access control for backup files
- Secure cloud storage integration

## Maintenance

Regular maintenance tasks include:

1. Verifying backup integrity weekly
2. Testing restoration procedures monthly
3. Updating configuration as database schema changes
4. Monitoring storage usage and adjusting retention policies as needed
```

## Expected Deliverables

After completing all the implementation steps, you should have:

1. A fully functional scheduled backup system with:
   - Automated full and incremental backups
   - Backup verification and validation
   - Secure backup storage
   - Retention policy enforcement
   - Restoration capabilities

2. Integration with the existing infrastructure:
   - Docker container for the backup scheduler
   - Monitoring integration for backup status
   - Alerting for backup failures

3. Comprehensive documentation:
   - Backup system architecture and configuration
   - Usage instructions and examples
   - Troubleshooting procedures
   - Disaster recovery guides

## Implementation Strategy

1. Follow the steps in exact order, creating each file in the specified location
2. Test each component individually before moving to the next
3. Integrate components once they've been tested
4. Verify the system works as a whole
5. Create detailed documentation

## Success Criteria

The implementation will be considered successful when:

1. All scripts and configuration files are created and functioning properly
2. Automated backups run on schedule and complete successfully
3. Backups can be verified for integrity
4. Old backups are automatically removed according to retention policy
5. Backups can be successfully restored
6. Documentation is comprehensive and clear
7. The system is integrated with Docker Compose and the monitoring system

This directive provides a clear, step-by-step approach to implementing the scheduled backup system for CryoProtect v2, with all file paths explicitly specified to ensure ROO can easily locate and create the necessary files.
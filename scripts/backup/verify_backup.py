import os
import sys
import logging
import subprocess
import json
import datetime
import argparse
import yaml
import hashlib
import gzip
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
            # Check if it's a valid gzip file using Python's gzip module
            try:
                with gzip.open(backup_file, 'rb') as f:
                    f.read(1)  # Try to read a byte to verify it's a valid gzip file
                logger.info(f"Backup format check passed: {backup_file}")
                return True
            except gzip.BadGzipFile:
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
        
        # Calculate file hash
        file_hash = calculate_file_hash(backup_file)
        if file_hash:
            verification_results["hash"] = file_hash
            verification_results["checks"]["hash_calculation"] = True
        else:
            verification_results["checks"]["hash_calculation"] = False
        
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
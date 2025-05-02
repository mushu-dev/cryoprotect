#!/usr/bin/env python3
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
#!/usr/bin/env python3
"""
Configuration Backup Script for CryoProtect v2

Backs up all critical configuration files as specified in config/file_backup_config.yaml.
Supports optional version control integration (git), logging, CLI, dry-run, and robust error handling.

Usage:
    python config_backup.py [--config CONFIG_PATH] [--backup-dir BACKUP_DIR] [--dry-run] [--no-git]
"""

import os
import sys
import argparse
import logging
import datetime
import yaml
from pathlib import Path

# Import backup utilities from file_backup.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))
import file_backup

# Set up logging for config backup
LOG_PATH = "logs/config_backup.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_PATH),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("config_backup")

def load_backup_config(config_path):
    """Load the backup configuration YAML."""
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"Failed to load backup config: {e}")
        sys.exit(1)

def backup_config_files(config, backup_dir, dry_run=False, use_git=True):
    """Backup configuration files with optional git integration."""
    try:
        config_source = config.get("sources", {}).get("config", {})
        if not config_source.get("enabled", False):
            logger.warning("Config backup is disabled in the configuration.")
            return False

        config_dir = config_source.get("path", "config")
        version_control = config_source.get("version_control", {})
        vc_enabled = version_control.get("enabled", False) and use_git

        if not os.path.exists(config_dir):
            logger.error(f"Config directory not found: {config_dir}")
            return False

        logger.info(f"Backing up configuration files from {config_dir}")

        # Ensure backup directory exists
        Path(backup_dir).mkdir(parents=True, exist_ok=True)

        if vc_enabled:
            logger.info("Using git for configuration backup")
            # Patch author info and commit message if provided
            orig_author_name = os.environ.get("GIT_AUTHOR_NAME")
            orig_author_email = os.environ.get("GIT_AUTHOR_EMAIL")
            os.environ["GIT_AUTHOR_NAME"] = version_control.get("author_name", "CryoProtect Backup")
            os.environ["GIT_AUTHOR_EMAIL"] = version_control.get("author_email", "backup@cryoprotect.example.com")
            # Use custom commit message if provided
            commit_message = version_control.get("commit_message", "Automated backup commit {timestamp}")
            commit_message = commit_message.format(timestamp=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            # Patch file_backup.git_backup_config to use custom commit message if needed
            result = file_backup.git_backup_config(config_dir, backup_dir, dry_run)
            # Restore environment
            if orig_author_name is not None:
                os.environ["GIT_AUTHOR_NAME"] = orig_author_name
            if orig_author_email is not None:
                os.environ["GIT_AUTHOR_EMAIL"] = orig_author_email
            return result
        else:
            # Exclude sensitive patterns
            exclude_patterns = [".DS_Store", "__pycache__", "*.pyc", "*_key", "*.pem", "*.key"]
            return file_backup.sync_directory(config_dir, backup_dir, exclude_patterns, dry_run)
    except Exception as e:
        logger.error(f"Error during config backup: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="CryoProtect v2 Configuration Backup")
    parser.add_argument("--config", type=str, default="config/file_backup_config.yaml", help="Path to backup config YAML")
    parser.add_argument("--backup-dir", type=str, default="backups/config", help="Directory to store config backups")
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without making changes")
    parser.add_argument("--no-git", action="store_true", help="Disable git version control even if enabled in config")
    args = parser.parse_args()

    config = load_backup_config(args.config)
    result = backup_config_files(
        config,
        args.backup_dir,
        dry_run=args.dry_run,
        use_git=not args.no_git
    )

    if result:
        logger.info("Configuration backup completed successfully.")
        sys.exit(0)
    else:
        logger.error("Configuration backup failed.")
        sys.exit(2)

if __name__ == "__main__":
    main()
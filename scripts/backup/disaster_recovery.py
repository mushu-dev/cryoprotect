#!/usr/bin/env python3
"""
CryoProtect v2 Disaster Recovery Script

Features:
- Full system recovery (database, files, configuration)
- Quick and partial restore (specific db, tables, files)
- Sandbox environment support for test restores
- Integration with restoration documentation (restore_guide.md)
- Logging of all actions to logs/disaster_recovery.log
- Command-line interface for recovery options, sandbox mode, and documentation output
- Robust error handling and reporting

Implements requirements from ROO_SCHEDULED_BACKUP_IMPLEMENTATION.md, section 5.
"""

import argparse
import logging
import os
import sys
import shutil
import traceback
from datetime import datetime

# Import database restore utilities
sys.path.append(os.path.join(os.path.dirname(__file__)))
try:
    import restore_database
except ImportError:
    restore_database = None

# Logging setup
LOG_FILE = "logs/disaster_recovery.log"
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - disaster_recovery - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("disaster_recovery")

# Default config paths
DEFAULT_CONFIG = "config/backup_config.yaml"
DEFAULT_DB_BACKUP_DIR = "backups/database"
DEFAULT_FILE_BACKUP_DIR = "backups/files"
DEFAULT_CONFIG_BACKUP_DIR = "backups/config"
RESTORE_GUIDE_PATH = "docs/backup/restore_guide.md"

def log_and_print(msg, level=logging.INFO):
    print(msg)
    logger.log(level, msg)

def restore_database_entry(args, config_path, sandbox=False):
    if not restore_database:
        log_and_print("Database restore module not found.", logging.ERROR)
        return False
    db_args = [
        "--config", config_path
    ]
    if args.db_file:
        db_args += ["--file", args.db_file]
    if args.db_type:
        db_args += ["--type", args.db_type]
    if args.clean:
        db_args.append("--clean")
    if args.create:
        db_args.append("--create")
    if args.no_owner:
        db_args.append("--no-owner")
    if args.no_privileges:
        db_args.append("--no-privileges")
    if args.schema:
        db_args += ["--schema", args.schema]
    if args.table:
        db_args += ["--table", args.table]
    if sandbox or args.sandbox or args.test_mode:
        db_args.append("--test-mode")
    try:
        # Call restore_database.main() as a subprocess to avoid argument conflicts
        import subprocess
        cmd = [sys.executable, os.path.join(os.path.dirname(__file__), "restore_database.py")] + db_args
        log_and_print(f"Invoking database restore: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        log_and_print(result.stdout)
        if result.returncode == 0:
            log_and_print("Database restore completed successfully.")
            return True
        else:
            log_and_print(f"Database restore failed: {result.stderr}", logging.ERROR)
            return False
    except Exception as e:
        log_and_print(f"Exception during database restore: {str(e)}\n{traceback.format_exc()}", logging.ERROR)
        return False

def restore_files(backup_dir, target_dir, items=None, sandbox=False):
    """
    Restore files from backup_dir to target_dir.
    If items is None, restore all. If sandbox, restore to a sandboxed location.
    """
    try:
        if not os.path.exists(backup_dir):
            log_and_print(f"Backup directory does not exist: {backup_dir}", logging.ERROR)
            return False
        if sandbox:
            target_dir = os.path.join(target_dir, "sandbox_restore_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
            log_and_print(f"Sandbox mode: restoring files to {target_dir}")
        os.makedirs(target_dir, exist_ok=True)
        restored = []
        if items:
            for item in items:
                src = os.path.join(backup_dir, item)
                dst = os.path.join(target_dir, item)
                if os.path.exists(src):
                    if os.path.isdir(src):
                        shutil.copytree(src, dst, dirs_exist_ok=True)
                    else:
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        shutil.copy2(src, dst)
                    restored.append(item)
                    log_and_print(f"Restored {item} to {dst}")
                else:
                    log_and_print(f"Item not found in backup: {item}", logging.WARNING)
        else:
            # Restore all
            for root, dirs, files in os.walk(backup_dir):
                for file in files:
                    rel_path = os.path.relpath(os.path.join(root, file), backup_dir)
                    src = os.path.join(backup_dir, rel_path)
                    dst = os.path.join(target_dir, rel_path)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copy2(src, dst)
                    restored.append(rel_path)
                    log_and_print(f"Restored {rel_path} to {dst}")
        log_and_print(f"File restore completed. {len(restored)} items restored.")
        return True
    except Exception as e:
        log_and_print(f"Exception during file restore: {str(e)}\n{traceback.format_exc()}", logging.ERROR)
        return False

def output_restore_guide():
    if os.path.exists(RESTORE_GUIDE_PATH):
        with open(RESTORE_GUIDE_PATH, "r") as f:
            print(f.read())
        logger.info("Restoration documentation output to console.")
    else:
        log_and_print("Restoration guide not found. Please see project documentation.", logging.WARNING)

def main():
    parser = argparse.ArgumentParser(
        description="CryoProtect v2 Disaster Recovery Tool\n"
                    "Performs full or partial system recovery (database, files, config), "
                    "supports sandbox/test mode, and integrates with restoration documentation."
    )
    parser.add_argument("--config", default=DEFAULT_CONFIG, help="Path to backup configuration file")
    parser.add_argument("--full", action="store_true", help="Perform full system recovery (database, files, config)")
    parser.add_argument("--quick", action="store_true", help="Quick restore (database only)")
    parser.add_argument("--db", action="store_true", help="Restore database only")
    parser.add_argument("--files", action="store_true", help="Restore files only")
    parser.add_argument("--config-files", action="store_true", help="Restore configuration files only")
    parser.add_argument("--db-file", help="Specific database backup file to restore")
    parser.add_argument("--db-type", choices=["full", "incremental"], help="Type of database backup to restore")
    parser.add_argument("--clean", action="store_true", help="Clean (drop) database objects before recreating")
    parser.add_argument("--create", action="store_true", help="Create the database before restoring")
    parser.add_argument("--no-owner", action="store_true", help="Do not set ownership of objects to match the original database")
    parser.add_argument("--no-privileges", action="store_true", help="Do not restore privileges (grant/revoke commands)")
    parser.add_argument("--schema", help="Restore only objects in this schema")
    parser.add_argument("--table", help="Restore only the named table")
    parser.add_argument("--file-items", nargs="+", help="Specific files or directories to restore")
    parser.add_argument("--sandbox", action="store_true", help="Run in sandbox mode (test restore, no changes to production)")
    parser.add_argument("--test-mode", action="store_true", help="Alias for --sandbox")
    parser.add_argument("--doc", action="store_true", help="Output restoration documentation and exit")
    args = parser.parse_args()

    logger.info("==== Disaster Recovery Script Started ====")
    logger.info(f"Arguments: {args}")

    if args.doc:
        output_restore_guide()
        return

    config_path = args.config

    # Determine restore actions
    if args.full:
        log_and_print("Starting FULL system recovery (database, files, config)...")
        db_ok = restore_database_entry(args, config_path, sandbox=args.sandbox or args.test_mode)
        files_ok = restore_files(DEFAULT_FILE_BACKUP_DIR, "restored_files", items=args.file_items, sandbox=args.sandbox or args.test_mode)
        config_ok = restore_files(DEFAULT_CONFIG_BACKUP_DIR, "restored_config", sandbox=args.sandbox or args.test_mode)
        if db_ok and files_ok and config_ok:
            log_and_print("Full system recovery completed successfully.")
        else:
            log_and_print("Full system recovery completed with errors. Check logs for details.", logging.ERROR)
    elif args.quick:
        log_and_print("Starting QUICK restore (database only)...")
        db_ok = restore_database_entry(args, config_path, sandbox=args.sandbox or args.test_mode)
        if db_ok:
            log_and_print("Quick database restore completed successfully.")
        else:
            log_and_print("Quick database restore failed. Check logs for details.", logging.ERROR)
    elif args.db:
        log_and_print("Starting database restore...")
        db_ok = restore_database_entry(args, config_path, sandbox=args.sandbox or args.test_mode)
        if db_ok:
            log_and_print("Database restore completed successfully.")
        else:
            log_and_print("Database restore failed. Check logs for details.", logging.ERROR)
    elif args.files:
        log_and_print("Starting file restore...")
        files_ok = restore_files(DEFAULT_FILE_BACKUP_DIR, "restored_files", items=args.file_items, sandbox=args.sandbox or args.test_mode)
        if files_ok:
            log_and_print("File restore completed successfully.")
        else:
            log_and_print("File restore failed. Check logs for details.", logging.ERROR)
    elif args.config_files:
        log_and_print("Starting configuration file restore...")
        config_ok = restore_files(DEFAULT_CONFIG_BACKUP_DIR, "restored_config", sandbox=args.sandbox or args.test_mode)
        if config_ok:
            log_and_print("Configuration file restore completed successfully.")
        else:
            log_and_print("Configuration file restore failed. Check logs for details.", logging.ERROR)
    else:
        log_and_print("No action specified. Use --full, --quick, --db, --files, --config-files, or --doc.", logging.WARNING)

    logger.info("==== Disaster Recovery Script Finished ====")

if __name__ == "__main__":
    main()
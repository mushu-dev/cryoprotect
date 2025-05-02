#!/usr/bin/env python3
"""
Storage Manager for CryoProtect v2

Features:
- Manage local and remote (cloud/S3) backup storage locations
- Enforce access controls and authentication for backup access
- Support encryption at rest for backup files
- Monitor storage usage and capacity, generate alerts if thresholds are exceeded
- Logging of all storage management actions to logs/storage_manager.log
- Command-line interface for storage operations, monitoring, and configuration
- Robust error handling and reporting
"""

import os
import sys
import argparse
import logging
import yaml
import shutil
import stat
import getpass
from pathlib import Path

from cryptography.fernet import Fernet, InvalidToken

# Import S3 provider from cloud_storage.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from cloud_storage import S3StorageProvider, load_config as load_cloud_config

LOG_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'logs', 'storage_manager.log')
CONFIG_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config', 'backup', 'backup_config.yaml')

# Setup logging
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s'
)
logger = logging.getLogger("storage_manager")

def load_config():
    try:
        with open(CONFIG_FILE, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"Error loading configuration: {str(e)}")
        return {}

def get_encryption_key(key_file):
    try:
        with open(key_file, 'rb') as f:
            key = f.read()
        return key
    except Exception as e:
        logger.error(f"Error loading encryption key: {str(e)}")
        raise

def encrypt_file(file_path, key):
    try:
        fernet = Fernet(key)
        with open(file_path, 'rb') as f:
            data = f.read()
        encrypted = fernet.encrypt(data)
        enc_path = f"{file_path}.enc"
        with open(enc_path, 'wb') as f:
            f.write(encrypted)
        logger.info(f"Encrypted file: {file_path} -> {enc_path}")
        return enc_path
    except Exception as e:
        logger.error(f"Error encrypting file {file_path}: {str(e)}")
        return None

def decrypt_file(enc_path, key):
    try:
        fernet = Fernet(key)
        with open(enc_path, 'rb') as f:
            encrypted = f.read()
        data = fernet.decrypt(encrypted)
        orig_path = enc_path[:-4] if enc_path.endswith('.enc') else enc_path + '.dec'
        with open(orig_path, 'wb') as f:
            f.write(data)
        logger.info(f"Decrypted file: {enc_path} -> {orig_path}")
        return orig_path
    except InvalidToken:
        logger.error(f"Invalid encryption key or corrupted file: {enc_path}")
        return None
    except Exception as e:
        logger.error(f"Error decrypting file {enc_path}: {str(e)}")
        return None

def set_file_permissions(file_path, mode=0o600):
    try:
        os.chmod(file_path, mode)
        logger.info(f"Set permissions {oct(mode)} on {file_path}")
    except Exception as e:
        logger.error(f"Error setting permissions on {file_path}: {str(e)}")

def check_local_storage_usage(path, threshold_percent=90):
    try:
        statvfs = os.statvfs(path)
        total = statvfs.f_frsize * statvfs.f_blocks
        free = statvfs.f_frsize * statvfs.f_bfree
        used = total - free
        percent = (used / total) * 100 if total > 0 else 0
        logger.info(f"Storage usage for {path}: {percent:.2f}%")
        if percent > threshold_percent:
            logger.warning(f"Storage usage exceeded threshold: {percent:.2f}% > {threshold_percent}%")
            print(f"ALERT: Storage usage exceeded threshold: {percent:.2f}% > {threshold_percent}%")
        return percent
    except Exception as e:
        logger.error(f"Error checking storage usage for {path}: {str(e)}")
        return None

def list_local_backups(directory):
    try:
        files = sorted([f for f in Path(directory).glob('**/*') if f.is_file()])
        for f in files:
            print(f"{f} ({os.path.getsize(f)} bytes)")
        logger.info(f"Listed {len(files)} files in {directory}")
        return files
    except Exception as e:
        logger.error(f"Error listing local backups in {directory}: {str(e)}")
        return []

def remove_local_backup(file_path):
    try:
        os.remove(file_path)
        logger.info(f"Removed local backup: {file_path}")
        print(f"Removed: {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error removing local backup {file_path}: {str(e)}")
        print(f"Error removing {file_path}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="CryoProtect v2 Storage Manager")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # List local backups
    list_parser = subparsers.add_parser("list", help="List local backup files")
    list_parser.add_argument("--dir", default=None, help="Directory to list (default: local_backup_dir from config)")

    # Remove local backup
    rm_parser = subparsers.add_parser("remove", help="Remove a local backup file")
    rm_parser.add_argument("file", help="Path to the backup file to remove")

    # Encrypt file
    enc_parser = subparsers.add_parser("encrypt", help="Encrypt a backup file")
    enc_parser.add_argument("file", help="Path to the file to encrypt")

    # Decrypt file
    dec_parser = subparsers.add_parser("decrypt", help="Decrypt an encrypted backup file")
    dec_parser.add_argument("file", help="Path to the .enc file to decrypt")

    # Monitor storage usage
    mon_parser = subparsers.add_parser("monitor", help="Monitor storage usage")
    mon_parser.add_argument("--dir", default=None, help="Directory to monitor (default: local_backup_dir from config)")
    mon_parser.add_argument("--threshold", type=int, default=90, help="Usage threshold percent for alerts (default: 90)")

    # S3 operations
    s3_parser = subparsers.add_parser("s3", help="Manage remote (S3) storage")
    s3_parser.add_argument("action", choices=["list", "upload", "download", "delete"], help="S3 action")
    s3_parser.add_argument("--file", help="Local file for upload/download")
    s3_parser.add_argument("--remote", help="Remote S3 path (key)")

    args = parser.parse_args()
    config = load_config()

    # Local backup directory
    local_dir = args.dir or config.get("local_backup_dir", "backups/database")
    encryption_cfg = config.get("encryption", {})
    key_file = encryption_cfg.get("key_file", "config/backup_encryption_key")

    if args.command == "list":
        list_local_backups(local_dir)

    elif args.command == "remove":
        remove_local_backup(args.file)

    elif args.command == "encrypt":
        try:
            key = get_encryption_key(key_file)
            enc_path = encrypt_file(args.file, key)
            if enc_path:
                set_file_permissions(enc_path)
                print(f"Encrypted: {enc_path}")
        except Exception as e:
            print(f"Encryption failed: {str(e)}")

    elif args.command == "decrypt":
        try:
            key = get_encryption_key(key_file)
            orig_path = decrypt_file(args.file, key)
            if orig_path:
                set_file_permissions(orig_path)
                print(f"Decrypted: {orig_path}")
        except Exception as e:
            print(f"Decryption failed: {str(e)}")

    elif args.command == "monitor":
        check_local_storage_usage(local_dir, args.threshold)

    elif args.command == "s3":
        aws_cfg = config.get("aws", {})
        if not aws_cfg.get("enabled", False):
            print("S3 integration is not enabled in the config.")
            logger.warning("S3 integration attempted but not enabled.")
            return
        try:
            s3 = S3StorageProvider(config)
            if args.action == "list":
                files = s3.list(prefix=aws_cfg.get("prefix", ""))
                for f in files:
                    print(f)
            elif args.action == "upload":
                if not args.file or not args.remote:
                    print("Specify --file (local) and --remote (S3 key) for upload.")
                    return
                s3.upload(args.file, args.remote)
                print(f"Uploaded {args.file} to S3 as {args.remote}")
            elif args.action == "download":
                if not args.file or not args.remote:
                    print("Specify --file (local) and --remote (S3 key) for download.")
                    return
                s3.download(args.remote, args.file)
                print(f"Downloaded {args.remote} from S3 to {args.file}")
            elif args.action == "delete":
                if not args.remote:
                    print("Specify --remote (S3 key) for delete.")
                    return
                s3.delete(args.remote)
                print(f"Deleted {args.remote} from S3")
        except Exception as e:
            logger.error(f"S3 operation failed: {str(e)}")
            print(f"S3 operation failed: {str(e)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        print(f"Fatal error: {str(e)}")
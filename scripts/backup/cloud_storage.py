import os
import sys
import argparse
import logging
import yaml
from abc import ABC, abstractmethod

try:
    import boto3
    from botocore.exceptions import BotoCoreError, ClientError
except ImportError:
    boto3 = None

LOG_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'logs', 'cloud_storage.log')
CONFIG_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config', 'backup_config.yaml')

# Setup logging
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s'
)
logger = logging.getLogger(__name__)

def load_config():
    with open(CONFIG_FILE, 'r') as f:
        config = yaml.safe_load(f)
    # Substitute environment variables if present
    aws_cfg = config.get('aws', {})
    for key in ['access_key', 'secret_key']:
        val = aws_cfg.get(key, '')
        if val and val.startswith('${') and val.endswith('}'):
            env_var = val[2:-1]
            aws_cfg[key] = os.environ.get(env_var, '')
    config['aws'] = aws_cfg
    return config

class CloudStorageProvider(ABC):
    @abstractmethod
    def upload(self, local_path, remote_path):
        pass

    @abstractmethod
    def download(self, remote_path, local_path):
        pass

    @abstractmethod
    def list(self, prefix=None):
        pass

    @abstractmethod
    def delete(self, remote_path):
        pass

class S3StorageProvider(CloudStorageProvider):
    def __init__(self, config):
        if boto3 is None:
            logger.error("boto3 is not installed. Please install boto3 to use S3 integration.")
            raise ImportError("boto3 is required for S3 integration.")
        aws = config.get('aws', {})
        self.enabled = aws.get('enabled', False)
        self.access_key = aws.get('access_key')
        self.secret_key = aws.get('secret_key')
        self.region = aws.get('region')
        self.bucket = aws.get('bucket')
        self.prefix = aws.get('prefix', '')
        if not self.enabled:
            logger.error("AWS S3 integration is disabled in config.")
            raise RuntimeError("AWS S3 integration is disabled in config.")
        if not all([self.access_key, self.secret_key, self.region, self.bucket]):
            logger.error("Missing AWS S3 credentials or configuration.")
            raise ValueError("Missing AWS S3 credentials or configuration.")
        self.s3 = boto3.client(
            's3',
            aws_access_key_id=self.access_key,
            aws_secret_access_key=self.secret_key,
            region_name=self.region
        )

    def _full_key(self, remote_path):
        if self.prefix:
            return f"{self.prefix.rstrip('/')}/{remote_path.lstrip('/')}"
        return remote_path.lstrip('/')

    def upload(self, local_path, remote_path):
        key = self._full_key(remote_path)
        try:
            self.s3.upload_file(local_path, self.bucket, key)
            logger.info(f"Uploaded {local_path} to s3://{self.bucket}/{key}")
            print(f"Uploaded {local_path} to s3://{self.bucket}/{key}")
        except (BotoCoreError, ClientError) as e:
            logger.error(f"Failed to upload {local_path} to S3: {e}")
            print(f"Error: Failed to upload {local_path} to S3: {e}")

    def download(self, remote_path, local_path):
        key = self._full_key(remote_path)
        try:
            self.s3.download_file(self.bucket, key, local_path)
            logger.info(f"Downloaded s3://{self.bucket}/{key} to {local_path}")
            print(f"Downloaded s3://{self.bucket}/{key} to {local_path}")
        except (BotoCoreError, ClientError) as e:
            logger.error(f"Failed to download {key} from S3: {e}")
            print(f"Error: Failed to download {key} from S3: {e}")

    def list(self, prefix=None):
        prefix = self._full_key(prefix) if prefix else self.prefix
        try:
            paginator = self.s3.get_paginator('list_objects_v2')
            result = []
            for page in paginator.paginate(Bucket=self.bucket, Prefix=prefix):
                for obj in page.get('Contents', []):
                    print(obj['Key'])
                    result.append(obj['Key'])
            logger.info(f"Listed objects in s3://{self.bucket}/{prefix}")
            return result
        except (BotoCoreError, ClientError) as e:
            logger.error(f"Failed to list objects in S3: {e}")
            print(f"Error: Failed to list objects in S3: {e}")
            return []

    def delete(self, remote_path):
        key = self._full_key(remote_path)
        try:
            self.s3.delete_object(Bucket=self.bucket, Key=key)
            logger.info(f"Deleted s3://{self.bucket}/{key}")
            print(f"Deleted s3://{self.bucket}/{key}")
        except (BotoCoreError, ClientError) as e:
            logger.error(f"Failed to delete {key} from S3: {e}")
            print(f"Error: Failed to delete {key} from S3: {e}")

class StubStorageProvider(CloudStorageProvider):
    def upload(self, local_path, remote_path):
        logger.warning("Upload not implemented for this provider.")
        print("Upload not implemented for this provider.")

    def download(self, remote_path, local_path):
        logger.warning("Download not implemented for this provider.")
        print("Download not implemented for this provider.")

    def list(self, prefix=None):
        logger.warning("List not implemented for this provider.")
        print("List not implemented for this provider.")
        return []

    def delete(self, remote_path):
        logger.warning("Delete not implemented for this provider.")
        print("Delete not implemented for this provider.")

def get_provider(config):
    if config.get('aws', {}).get('enabled', False):
        return S3StorageProvider(config)
    # Add other providers here as needed
    return StubStorageProvider()

def main():
    parser = argparse.ArgumentParser(description="CryoProtect Cloud Storage Integration")
    subparsers = parser.add_subparsers(dest='command', required=True)

    upload_parser = subparsers.add_parser('upload', help='Upload a file to cloud storage')
    upload_parser.add_argument('local_path', help='Path to local file')
    upload_parser.add_argument('remote_path', help='Remote path in cloud storage')

    download_parser = subparsers.add_parser('download', help='Download a file from cloud storage')
    download_parser.add_argument('remote_path', help='Remote path in cloud storage')
    download_parser.add_argument('local_path', help='Path to save the downloaded file')

    list_parser = subparsers.add_parser('list', help='List files in cloud storage')
    list_parser.add_argument('--prefix', default=None, help='Prefix to filter files')

    delete_parser = subparsers.add_parser('delete', help='Delete a file from cloud storage')
    delete_parser.add_argument('remote_path', help='Remote path in cloud storage')

    args = parser.parse_args()

    try:
        config = load_config()
        provider = get_provider(config)
        if args.command == 'upload':
            provider.upload(args.local_path, args.remote_path)
        elif args.command == 'download':
            provider.download(args.remote_path, args.local_path)
        elif args.command == 'list':
            provider.list(args.prefix)
        elif args.command == 'delete':
            provider.delete(args.remote_path)
        else:
            parser.print_help()
    except Exception as e:
        logger.error(f"Unhandled exception: {e}", exc_info=True)
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
# CryoProtect v2 Backup System

This directory contains a robust, production-ready backup system for CryoProtect v2. The backup system provides:

1. Automated backup scheduling
2. Integrity verification after each backup
3. Configurable retention policies
4. Cross-region backup storage
5. Restoration procedures with verification

## Installation

The backup system requires Python 3.6+ and the following dependencies:

```bash
pip install schedule
```

For cross-region storage, additional dependencies may be required:
- Amazon S3: `pip install boto3`
- Azure Blob Storage: `pip install azure-storage-blob`
- Google Cloud Storage: `pip install google-cloud-storage`
- SFTP: `pip install paramiko`

## Configuration

The backup system can be configured using a JSON configuration file. A template is provided in `backup_config.json.template`. Copy this file to `backup_config.json` and modify it according to your needs:

```bash
cp backup_config.json.template backup_config.json
```

### Configuration Options

- `backup_dir`: Directory where backups will be stored
- `retention`: Retention policy for different backup types
  - `daily`: Number of daily backups to keep
  - `weekly`: Number of weekly backups to keep
  - `monthly`: Number of monthly backups to keep
  - `yearly`: Number of yearly backups to keep
- `schedule`: Schedule for automated backups
  - `daily`: Time for daily backups (HH:MM format)
  - `weekly`: Day for weekly backups (monday, tuesday, etc.)
  - `monthly`: Day of month for monthly backups (1-31)
- `cross_region`: Configuration for cross-region storage
  - `enabled`: Whether to enable cross-region storage
  - `method`: Storage method (s3, azure, gcp, sftp)
  - `config`: Configuration for the selected method
- `verification`: Configuration for backup verification
  - `enabled`: Whether to enable verification
  - `methods`: Verification methods to use (checksum, restore_test)
- `compression`: Configuration for backup compression
  - `enabled`: Whether to enable compression
  - `method`: Compression method (zip, tar.gz)
- `encryption`: Configuration for backup encryption
  - `enabled`: Whether to enable encryption
  - `method`: Encryption method (aes256)
  - `key_file`: Path to encryption key file

## Usage

### Creating a Backup

To create a manual backup:

```bash
python -m backup.backup_manager backup --config backup/backup_config.json
```

Options:
- `--project-id`: Supabase project ID (default: from environment)
- `--schema`: Database schema to backup (default: public)
- `--format`: Backup format (json, sql) (default: json)
- `--type`: Backup type (manual, daily, weekly, monthly, yearly) (default: manual)

### Verifying a Backup

To verify a backup:

```bash
python -m backup.backup_manager verify /path/to/backup --config backup/backup_config.json
```

### Starting the Backup Scheduler

To start the backup scheduler in the background:

```bash
python -m backup.backup_manager schedule --config backup/backup_config.json
```

To run the scheduler in blocking mode:

```bash
python -m backup.backup_manager schedule --blocking --config backup/backup_config.json
```

### Applying Retention Policy

To manually apply the retention policy:

```bash
python -m backup.backup_manager retention --config backup/backup_config.json
```

## Backup Directory Structure

Each backup is stored in a subdirectory with the following naming convention:

```
{backup_type}_backup_{timestamp}
```

For example:
```
daily_backup_20250422_123456
```

Each backup directory contains:
- `metadata.json`: Metadata about the backup
- `summary.json`: Summary of the backup process
- `checksums.json`: Checksums of backup files (if verification is enabled)
- Table backup files: One file per table, with the extension `.json` or `.sql` depending on the format

## Restoration Procedure

To restore a backup:

1. Identify the backup directory you want to restore
2. If the backup is in SQL format, you can directly execute the SQL files to restore the database
3. If the backup is in JSON format, you can use the Supabase API to restore the data

### SQL Format Restoration

```bash
# Connect to your database
psql -h your-database-host -U your-database-user -d your-database-name

# Execute the SQL files
\i /path/to/backup/table1.sql
\i /path/to/backup/table2.sql
# ...
```

### JSON Format Restoration

For JSON format backups, you can use the Supabase API to restore the data:

```python
import json
from supabase_py import create_client

# Initialize Supabase client
supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

# Load backup data
with open('/path/to/backup/table1.json', 'r') as f:
    data = json.load(f)

# Insert data into table
for item in data:
    supabase.table('table1').insert(item).execute()
```

## Verification Steps

After restoring a backup, you should verify the integrity of the restored data:

1. Check that all tables have been restored
2. Verify row counts match the original data
3. Run application-specific tests to ensure data consistency
4. Check that relationships between tables are maintained

## Troubleshooting

### Common Issues

1. **Backup fails with permission error**
   - Ensure the backup directory is writable by the user running the backup

2. **Scheduler doesn't run backups**
   - Check that the schedule configuration is correct
   - Ensure the process is running continuously

3. **Cross-region storage fails**
   - Verify credentials and connection settings
   - Check network connectivity to the storage service

### Logs

Backup logs are stored in the application log directory. Check these logs for detailed information about backup operations.

## Contributing

Contributions to the backup system are welcome. Please follow the standard contribution guidelines for the CryoProtect v2 project.
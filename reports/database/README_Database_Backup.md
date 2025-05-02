# CryoProtect Database Backup

This document provides instructions for using the Supabase database backup script for the CryoProtect project.

## Overview

The `create_database_backup.py` script creates a timestamped backup of all tables in the public schema of the CryoProtect Supabase project. It uses the Supabase MCP tools for database operations.

## Features

- Creates timestamped backups of all tables in the public schema
- Supports both JSON and SQL backup formats
- Includes error handling and logging
- Creates a metadata file with backup information
- Creates a summary file with backup statistics

## Prerequisites

- Python 3.6+
- Supabase project with the CryoProtect schema
- Environment variables set in `.env` file:
  - `SUPABASE_URL`: Your Supabase project URL
  - `SUPABASE_KEY`: Your Supabase API key

## Usage

### On Windows

Run the backup script using the provided batch file:

```
run_database_backup.bat
```

Or directly with Python:

```
python create_database_backup.py
```

### On Linux/macOS

Run the backup script directly with Python:

```
python create_database_backup.py
```

### Command-line Options

The script supports the following command-line options:

- `--format`: Output format for the backup (json or sql, default: json)
- `--project-id`: Supabase project ID (default: tsdlmynydfuypiugmkev)
- `--schema`: Database schema to backup (default: public)

Example:

```
python create_database_backup.py --format sql --project-id tsdlmynydfuypiugmkev --schema public
```

## Backup Structure

The script creates a backup directory structure as follows:

```
backups/
└── backup_YYYYMMDD_HHMMSS/
    ├── metadata.json
    ├── summary.json
    ├── molecules.json (or .sql)
    ├── property_types.json (or .sql)
    ├── molecular_properties.json (or .sql)
    └── ... (other tables)
```

- `metadata.json`: Contains information about the backup (timestamp, project ID, schema, format, tables)
- `summary.json`: Contains statistics about the backup (total tables, successful tables, etc.)
- Table files: Each table is backed up to a separate file in the specified format (JSON or SQL)

## Restoring from Backup

### JSON Format

To restore from a JSON backup, you can use the Supabase REST API or the Supabase dashboard.

### SQL Format

To restore from an SQL backup, you can use the Supabase SQL editor or a database client to execute the SQL statements.

## Troubleshooting

- If the script fails to connect to Supabase, check your `.env` file for correct credentials
- If the script fails to list tables, check if your Supabase project has the correct schema
- Check the log file for detailed error messages

## Scheduling Regular Backups

### On Windows

You can use Windows Task Scheduler to schedule regular backups:

1. Open Task Scheduler
2. Create a new task
3. Set the trigger (e.g., daily at 2 AM)
4. Set the action to start a program
5. Program/script: `python`
6. Arguments: `create_database_backup.py`
7. Start in: `C:\path\to\CryoProtect v2`

### On Linux/macOS

You can use cron to schedule regular backups:

1. Open the crontab file: `crontab -e`
2. Add a line to schedule the backup (e.g., daily at 2 AM):
   ```
   0 2 * * * cd /path/to/CryoProtect/v2 && python create_database_backup.py
# CryoProtect Database Integrity Verification Guide

This guide documents the comprehensive data integrity verification system for the CryoProtect database. The system is designed to identify and repair integrity issues, foreign key violations, and data quality problems across all database tables.

## Table of Contents

1. [Overview](#overview)
2. [Verification Tools](#verification-tools)
3. [Using the Verification Tools](#using-the-verification-tools)
4. [Understanding Verification Reports](#understanding-verification-reports)
5. [Fixing Data Integrity Issues](#fixing-data-integrity-issues)
6. [Automated Verification](#automated-verification)
7. [Troubleshooting](#troubleshooting)
8. [Best Practices](#best-practices)

## Overview

The database integrity verification system consists of the following components:

1. **Enhanced Verification Scripts**: Python scripts that perform detailed checks on the database structure and data
2. **Automated Reports**: JSON, Markdown, and visual HTML reports of verification results
3. **Fix Tools**: SQL and shell scripts to repair identified issues
4. **Scheduled Verification**: Systemd timers for regular automated checks
5. **Notification System**: Email alerts when issues are detected

## Verification Tools

The system includes the following tools:

| Tool | Description |
| ---- | ----------- |
| **verify_database_integrity_fixed.py** | The core Python script that performs database verification |
| **run_database_integrity_check_fixed.sh** | Shell script to run verification and generate reports |
| **fix_database_integrity_issues.sql** | SQL script to fix common data integrity issues |
| **fix_database_integrity.sh** | Shell script to safely apply fixes with backup option |
| **generate_database_integrity_report.py** | Python script to create visual HTML reports |
| **cryoprotect-db-integrity.service** | Systemd service for scheduled verification |
| **cryoprotect-db-integrity.timer** | Systemd timer for recurring verification |

## Using the Verification Tools

### Running Verification

To perform a database integrity check:

```bash
./run_database_integrity_check_fixed.sh [--verbose] [--report=filename.json]
```

Options:
- `--verbose`: Show detailed progress information
- `--report=filename.json`: Specify custom output file for the JSON report

The script will generate:
- A JSON report with detailed verification results
- A Markdown summary report
- An HTML visual report (if the generator script is available)

### Applying Fixes

To fix identified data integrity issues:

```bash
./fix_database_integrity.sh [--backup] [--apply]
```

Options:
- `--backup`: Create a database backup before applying fixes (recommended)
- `--apply`: Actually apply the fixes (without this flag, runs in preview mode)

### Setting Up Scheduled Verification

To install the systemd timer for daily automatic verification:

```bash
# Copy service files to systemd directory
sudo cp cryoprotect-db-integrity.service /etc/systemd/system/
sudo cp cryoprotect-db-integrity.timer /etc/systemd/system/

# Reload systemd configuration
sudo systemctl daemon-reload

# Enable and start the timer
sudo systemctl enable cryoprotect-db-integrity.timer
sudo systemctl start cryoprotect-db-integrity.timer

# Check status
sudo systemctl status cryoprotect-db-integrity.timer
```

The service will run daily at 2:00 AM and generate reports.

## Understanding Verification Reports

### Report Types

1. **JSON Report** (`database_integrity_report_YYYYMMDD_HHMMSS.json`):
   - Full detailed verification results
   - Complete table counts
   - Detailed issue descriptions with timestamps
   - Overall execution statistics

2. **Markdown Report** (`database_integrity_report_YYYYMMDD_HHMMSS.md`):
   - Summary of issues found
   - Table counts
   - List of issues grouped by severity

3. **HTML Report** (`database_integrity_report_YYYYMMDD_HHMMSS.html`):
   - Interactive visual report
   - Charts showing data distribution
   - Searchable issues list
   - Color-coded severity indicators

### Error Severity Levels

The system categorizes issues into two severity levels:

1. **Error**: Critical issues that require immediate attention
   - Foreign key violations
   - Missing required values
   - Structure inconsistencies

2. **Warning**: Potential issues that should be evaluated
   - Data quality concerns
   - Suspicious values
   - Incomplete relationship coverage
   - Inconsistent units

## Fixing Data Integrity Issues

The system includes tools for automatic and manual fixing of issues:

### Automatic Fixes

The `fix_database_integrity_issues.sql` script can address many common issues:

1. **NULL Values in Required Fields**: Sets appropriate default values
2. **Foreign Key Violations**: Updates or nullifies invalid references
3. **Duplicate Records**: Merges duplicate molecules based on InChIKey
4. **Invalid Chemical Data**: Nullifies invalid SMILES, InChI, and InChIKey strings
5. **Mixture Issues**: Standardizes concentration units and adjusts percentages
6. **User/Team Relationships**: Creates default team and fixes user assignments

### Manual Fixes

For issues that require manual intervention:

1. Review the verification reports to understand the specific issues
2. Connect to the database directly for complex fixes:
   ```bash
   export $(grep -v '^#' .env | xargs)
   PGPASSWORD=$SUPABASE_DB_PASSWORD psql -h $SUPABASE_DB_HOST -p $SUPABASE_DB_PORT -U $SUPABASE_DB_USER -d $SUPABASE_DB_NAME
   ```
3. Apply targeted SQL fixes based on the specific issues
4. Run verification again to ensure all issues are resolved

## Automated Verification

The system supports automated verification through:

1. **Systemd Timers**: Run verification daily at scheduled times
2. **Email Notifications**: Send alerts when issues are found
3. **Report History**: Keep historical reports for trend analysis

To configure email notifications, edit the `run_verify_all_data_integrity.sh` script and update the email address:

```bash
echo "Database integrity verification found issues..." | mail -s "ALERT: Database Integrity Issues" your-email@example.com
```

## Troubleshooting

### Common Issues

1. **Connection Errors**:
   - Verify database credentials in `.env` file
   - Check network connectivity to database server
   - Ensure the database server is running

2. **Permission Issues**:
   - Verify database user has sufficient permissions
   - Check file permissions on scripts (should be executable)

3. **Long Execution Times**:
   - For large databases, use the `--output` option to save reports
   - Consider running in a screen or tmux session

4. **Failed Fixes**:
   - Always use the `--backup` option when applying fixes
   - Review logs to understand which fixes failed
   - Apply targeted fixes manually for complex issues

## Best Practices

1. **Regular Verification**:
   - Run weekly verification in development environments
   - Run daily verification in production
   - Review reports after major data imports or migrations

2. **Safe Fixing Process**:
   - Always create backups before applying fixes
   - Test fixes in development environment first
   - Document all manual fixes applied

3. **Issue Prioritization**:
   - Address error-level issues immediately
   - Schedule regular reviews of warning-level issues
   - Track trends in issue frequency

4. **Extending the System**:
   - Add new checks to `verify_database_integrity_fixed.py`
   - Add new fixes to `fix_database_integrity_issues.sql`
   - Update HTML report generation for new issue types

By following this guide, you can ensure the ongoing integrity and quality of the CryoProtect database system.
# CryoProtect Database Integrity Verification

This document describes the comprehensive data integrity verification system for the CryoProtect database. The system is designed to identify integrity issues, foreign key violations, and data quality problems across all database tables.

## Overview

The database integrity verification system consists of the following components:

1. **Enhanced Verification Script**: A Python script that performs detailed checks on the database structure and data
2. **Automated Reports**: JSON, Markdown, and visual HTML reports of verification results
3. **Scheduled Verification**: Systemd timers for regular automated checks
4. **Notification System**: Email alerts when issues are detected

## Key Features

- **Comprehensive Foreign Key Validation**: Verifies all foreign key relationships for integrity
- **Required Fields Checking**: Ensures non-nullable fields have valid values
- **Chemical Data Quality**: Validates molecular properties with scientific standards
- **Mixture Composition Analysis**: Checks mixture components and concentration values
- **Team-Project Access Control**: Validates user, team, and project relationships
- **Real-time Database Querying**: Accesses the database directly for up-to-date verification
- **Interactive Visual Reporting**: Provides charts and searchable issue listings

## Usage

### Manual Verification

Run the full verification process with:

```bash
./run_verify_all_data_integrity.sh
```

Options:
- `--notify`: Send email notifications with results
- `--schedule`: Add the verification to system crontab (runs daily at 2:00 AM)

For more targeted verification, use:

```bash
./run_database_integrity_check.sh [--verbose] [--report=filename.json]
```

### Scheduled Verification via Systemd

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

The service will run daily at 2:00 AM and generate reports in the `reports/` directory.

## Reports

All verification reports are saved in the `reports/` directory with timestamped filenames. The system creates:

1. **JSON Report** (`db_integrity_YYYYMMDD_HHMMSS.json`): Contains all verification details in structured format
2. **HTML Report** (`db_integrity_YYYYMMDD_HHMMSS.html`): Interactive visual report with charts and searchable issues
3. **Markdown Report** (`db_integrity_YYYYMMDD_HHMMSS.md`): Summary report suitable for sharing
4. **Log File** (`db_integrity_YYYYMMDD_HHMMSS.log`): Detailed execution log

Symlinks to the latest reports are maintained at:
- `reports/db_integrity_latest.json`
- `reports/db_integrity_latest.html`
- `reports/db_integrity_latest.md`
- `reports/db_integrity_latest.log`

## Error Severity Levels

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

## Database Checks

### Table Existence and Data Verification

- Verifies all expected tables exist
- Counts rows in each table
- Identifies empty tables that should contain data

### Foreign Key Relationships

Validates all foreign key relationships, including:
- `mixture_components.mixture_id → mixtures.id`
- `mixture_components.molecule_id → molecules.id`
- `molecular_properties.molecule_id → molecules.id`
- `experiments.mixture_id → mixtures.id`
- `experiments.project_id → projects.id`
- `predictions.molecule_id → molecules.id`
- `projects.team_id → teams.id`
- `user_profile.team_id → teams.id`
- `team_members.team_id → teams.id`
- And many more...

### Chemical Data Quality

Checks for scientifically valid data:
- Valid SMILES notation
- Valid InChI and InChIKey formats
- Reasonable molecular weight values
- Consistent molecular formulae

### Mixture Composition

Verifies mixture integrity:
- Each mixture has at least one component
- Consistent concentration units
- Valid percentage totals (for percentage-based units)
- Proper parent-child relationships

### Team and Project Integrity

Validates access control relationships:
- Users belong to authorized teams
- Projects are properly associated with teams
- Team members have valid roles
- Proper authorization flow

## Troubleshooting

If issues are detected, the reports provide detailed information about:
- The specific tables and columns affected
- Sample values causing problems
- The type and severity of each issue
- Timestamps for when issues were detected

For persistent issues, you may need to:
1. Check database migration scripts
2. Verify data import processes
3. Examine application code handling the affected tables
4. Run data repair scripts

## Extending the System

The verification system is designed to be extensible. To add new checks:

1. Modify `verify_database_integrity_enhanced.py` to add new verification methods
2. Update report generation in `generate_database_integrity_report.py`
3. Test with a small dataset before running against production

## Dependencies

- Python 3.6+
- PostgreSQL client libraries
- psycopg2
- jq (optional, for better report processing)
- Mail client (optional, for notifications)

## Best Practices

- Run weekly verification in development environments
- Run daily verification in production
- Review reports after major data imports or migrations
- Address error-level issues immediately
- Schedule regular reviews of warning-level issues
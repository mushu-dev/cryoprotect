# CryoProtect v2 - Migration to Plural Table Names

This document provides instructions for safely migrating the CryoProtect v2 database from singular-named tables to plural-named tables using the `migrate_to_plural_tables.py` script.

## Overview

The migration script implements the "expand and contract" approach to safely migrate data from the following singular-named tables to their plural counterparts:

1. molecule → molecules
2. mixture → mixtures
3. experiment → experiments
4. prediction → predictions
5. project → projects

## Features

- **Transaction Safety**: All operations are performed within transactions to ensure database integrity
- **Backup Tables**: Creates backup tables before migration for safety and rollback capability
- **Data Validation**: Validates data integrity after migration to ensure no data loss
- **Foreign Key Handling**: Properly updates all foreign key relationships
- **View Updates**: Updates database views to reference the new plural tables
- **Dry Run Mode**: Allows testing the migration without committing changes
- **Rollback Support**: Provides functionality to roll back the migration if needed
- **Detailed Logging**: Logs all operations for audit and troubleshooting

## Prerequisites

- Python 3.7+
- psycopg2 (PostgreSQL adapter for Python)
- Access to the CryoProtect v2 Supabase database
- Database connection details in environment variables

## Environment Setup

The script requires database connection details to be set in environment variables:

```
DATABASE_URL=postgresql://username:password@host:port/database
```

Alternatively, if you're using Supabase:

```
SUPABASE_DB_URL=postgresql://postgres:password@db.tsdlmynydfuypiugmkev.supabase.co:5432/postgres
```

## Usage

### Performing a Dry Run

It's highly recommended to perform a dry run before executing the actual migration:

```bash
python migrate_to_plural_tables.py --dry-run
```

This will execute all migration steps but roll back the changes at the end, allowing you to verify that the migration process works correctly without making permanent changes.

### Executing the Migration

Once you've verified the dry run works correctly, execute the actual migration:

```bash
python migrate_to_plural_tables.py
```

The script will:
1. Create backup tables for all tables being migrated
2. Create new plural-named tables with the same structure
3. Copy data from singular to plural tables
4. Validate data integrity
5. Update foreign key relationships
6. Update views to reference the new plural tables
7. Generate a migration report

### Rolling Back the Migration

If you need to roll back the migration:

```bash
python migrate_to_plural_tables.py --rollback --report migration_report_YYYYMMDD_HHMMSS.json
```

This will restore data from the backup tables created during the migration.

## Migration Process Details

### 1. Backup Creation

For each table being migrated, a backup table is created with the naming pattern:
```
{table_name}_backup_YYYYMMDD
```

### 2. Table Creation

New plural tables are created with identical structure to the original tables, including:
- Column definitions
- Primary keys
- Unique constraints
- Indexes

### 3. Data Migration

Data is copied from singular to plural tables using SQL INSERT statements.

### 4. Data Validation

The script validates data integrity by:
- Comparing row counts between singular and plural tables
- Sampling rows from the singular table and verifying they exist in the plural table

### 5. Foreign Key Updates

The script updates foreign key relationships in the following tables:
- mixture_component
- molecular_property
- prediction
- experiment
- project_membership

### 6. View Updates

The following views are updated to reference the new plural tables:
- experiment_with_results
- mixture_with_components
- mixtures_with_components
- molecule_with_properties
- molecules_with_properties

## Migration Report

After successful migration, a report file is generated with the naming pattern:
```
migration_report_YYYYMMDD_HHMMSS.json
```

This report contains:
- Timestamp of the migration
- List of tables migrated
- Names of backup tables
- Foreign key relationships updated
- Views updated

This report is required for rollback operations.

## Logging

The script generates detailed logs with the naming pattern:
```
migration_log_YYYYMMDD_HHMMSS.log
```

## Important Precautions

1. **Backup Your Database**: Always create a full database backup before running the migration
2. **Maintenance Window**: Schedule the migration during a maintenance window when the application is not in use
3. **Test in Staging**: Test the migration in a staging environment before running in production
4. **Application Updates**: Update application code to reference the new plural table names after migration
5. **Monitor Performance**: Monitor database performance after migration to ensure no issues

## Troubleshooting

### Common Issues

1. **Connection Errors**:
   - Verify database connection details in environment variables
   - Check network connectivity to the database server

2. **Permission Errors**:
   - Ensure the database user has sufficient privileges to create tables, modify constraints, etc.

3. **Foreign Key Violations**:
   - If foreign key updates fail, check for data inconsistencies in the original tables

4. **View Update Failures**:
   - If view updates fail, manually update the view definitions using the backup files created during migration

### Getting Help

If you encounter issues during migration, check:
1. The migration log file for detailed error messages
2. The PostgreSQL server logs for database-specific errors

## After Migration

After successful migration, you should:

1. Update application code to reference the new plural table names
2. Test the application thoroughly to ensure it works with the new table names
3. Consider dropping the old singular tables after a period of stability
4. Archive the backup tables and migration reports for future reference

## License

This migration script is part of the CryoProtect v2 project and is subject to the same licensing terms.
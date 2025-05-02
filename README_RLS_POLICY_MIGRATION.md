# CryoProtect v2 - RLS Policy Migration

This document provides instructions for safely migrating Row Level Security (RLS) policies from singular-named tables to plural-named tables using the `migrate_rls_policies.py` script.

## Overview

The migration script ensures that security is maintained during and after the database migration process by:

1. Identifying existing RLS policies on singular-named tables
2. Creating equivalent policies on plural-named tables
3. Updating policy expressions that reference table names
4. Verifying that policies are correctly applied

The script handles policies for the following table migrations:

1. molecule → molecules
2. mixture → mixtures
3. experiment → experiments
4. prediction → predictions
5. project → projects

## Features

- **Transaction Safety**: All operations are performed within transactions to ensure database integrity
- **Policy Expression Updates**: Automatically updates policy expressions to reference plural table names
- **Verification**: Validates that policies are correctly applied after migration
- **Dry Run Mode**: Allows testing the migration without committing changes
- **Rollback Support**: Provides functionality to roll back the migration if needed
- **Detailed Logging**: Logs all operations for audit and troubleshooting
- **Comprehensive Reporting**: Generates detailed reports of migrated policies and verification results

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
python migrate_rls_policies.py --dry-run
```

This will execute all migration steps but roll back the changes at the end, allowing you to verify that the migration process works correctly without making permanent changes.

### Executing the Migration

Once you've verified the dry run works correctly, execute the actual migration:

```bash
python migrate_rls_policies.py
```

The script will:
1. Identify existing RLS policies on singular-named tables
2. Create equivalent policies on plural-named tables
3. Update policy expressions that reference table names
4. Verify that policies are correctly applied
5. Generate migration and verification reports

### Verifying Existing Policies

To verify policies that were previously migrated:

```bash
python migrate_rls_policies.py --verify-only --report rls_migration_report_YYYYMMDD_HHMMSS.json
```

This will check if the policies specified in the migration report exist and have the correct expressions.

### Rolling Back the Migration

If you need to roll back the migration:

```bash
python migrate_rls_policies.py --rollback --report rls_migration_report_YYYYMMDD_HHMMSS.json
```

This will drop the policies that were created during the migration.

## Migration Process Details

### 1. Policy Identification

The script identifies all RLS policies on singular-named tables by querying the `pg_policies` system catalog.

### 2. Policy Migration

For each policy on a singular-named table, the script:
- Creates an equivalent policy on the corresponding plural-named table
- Updates policy expressions to reference plural table names
- Preserves all policy attributes (roles, commands, etc.)

### 3. Expression Updates

Policy expressions (USING and WITH CHECK clauses) are updated to reference plural table names. For example:

```sql
-- Original expression
EXISTS (SELECT 1 FROM project WHERE project.id = project_membership.project_id)

-- Updated expression
EXISTS (SELECT 1 FROM projects WHERE projects.id = project_membership.project_id)
```

### 4. Policy Verification

After migration, the script verifies that:
- All policies were successfully created on plural tables
- Policy expressions were correctly updated

## Migration Reports

### Migration Report

After migration, a report file is generated with the naming pattern:
```
rls_migration_report_YYYYMMDD_HHMMSS.json
```

This report contains:
- Timestamp of the migration
- List of migrated policies
- Original and updated policy expressions
- Failed migrations
- Skipped tables

### Verification Report

After verification, a report file is generated with the naming pattern:
```
rls_verification_report_YYYYMMDD_HHMMSS.json
```

This report contains:
- Timestamp of the verification
- List of verified policies
- Failed verifications with details

## Logging

The script generates detailed logs with the naming pattern:
```
rls_migration_YYYYMMDD_HHMMSS.log
```

## Important Precautions

1. **Backup Your Database**: Always create a full database backup before running the migration
2. **Maintenance Window**: Schedule the migration during a maintenance window when the application is not in use
3. **Test in Staging**: Test the migration in a staging environment before running in production
4. **Run After Table Migration**: Run this script after the table migration (`migrate_to_plural_tables.py`) has been completed
5. **Verify Security**: After migration, verify that security controls are working as expected

## Troubleshooting

### Common Issues

1. **Connection Errors**:
   - Verify database connection details in environment variables
   - Check network connectivity to the database server

2. **Permission Errors**:
   - Ensure the database user has sufficient privileges to create and drop policies

3. **Policy Expression Errors**:
   - If policy expressions fail to update correctly, check for complex expressions that may need manual updates

4. **Verification Failures**:
   - If verification fails, check the verification report for details
   - You may need to manually adjust policy expressions

### Getting Help

If you encounter issues during migration, check:
1. The migration log file for detailed error messages
2. The PostgreSQL server logs for database-specific errors

## After Migration

After successful migration, you should:

1. Verify that security controls are working as expected
2. Test the application thoroughly to ensure it works with the new policies
3. Consider dropping the old singular tables after a period of stability
4. Archive the migration and verification reports for future reference

## License

This migration script is part of the CryoProtect v2 project and is subject to the same licensing terms.
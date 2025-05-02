# CryoProtect Database Schema Standardization

This document outlines the process for standardizing the CryoProtect Supabase database schema, including converting singular table names to plural, ensuring proper UUID primary keys, and adding appropriate foreign key constraints.

## Overview

The standardization script performs the following operations:

1. Converts all singular-named tables to plural:
   - `molecule` → `molecules`
   - `mixture` → `mixtures`
   - `prediction` → `predictions`
   - `experiment` → `experiments`
   - `experiment_property` → `experiment_properties`
   - `mixture_component` → `mixture_components`
   - `calculation_method` → `calculation_methods`
   - `property_type` → `property_types`
   - `project` → `projects`
   - `team` → `teams`

2. Ensures all tables use UUID as primary key with DEFAULT gen_random_uuid()

3. Adds proper REFERENCES constraints with indexes for all foreign keys

4. Includes a rollback mechanism in case of failure

## Prerequisites

- Python 3.6+
- Access to the Supabase project with appropriate permissions
- Supabase CLI or access token

## Files

- `standardize_schema.py`: Main script for schema standardization
- `supabase_mcp_tools.py`: Helper module for interacting with Supabase via MCP

## Usage

1. Ensure you have the necessary permissions to modify the database schema.

2. Review the script to understand the changes that will be made.

3. Run the standardization script:

```bash
python standardize_schema.py
```

4. The script will:
   - Create migration tracking tables to record all operations
   - Process each table in the mapping
   - Create new plural tables with the same schema as the singular tables
   - Copy data from singular to plural tables
   - Add appropriate foreign key constraints and indexes
   - Log all operations and their status

5. If any errors occur, you will be prompted to rollback the migration.

## Rollback Process

The script includes a comprehensive rollback mechanism:

1. If any operation fails, you will be prompted to rollback the entire migration.
2. The rollback process executes the "down" SQL for each applied operation in reverse order.
3. All operations are tracked in the `schema_migration_operations` table.
4. The migration status is updated in the `schema_migrations` table.

## Logging

The script creates a detailed log file with the format `schema_migration_YYYYMMDD_HHMMSS.log` that includes:

- All SQL operations executed
- Success/failure status of each operation
- Error messages for failed operations
- Overall migration status

## Post-Migration Steps

After successfully running the standardization script, you should:

1. Update your application code to reference the new plural table names
2. Test your application thoroughly to ensure all functionality works with the new schema
3. Once confirmed working, you can drop the old singular tables if needed (not done automatically by the script for safety)

## Troubleshooting

If you encounter issues during the migration:

1. Check the log file for detailed error messages
2. If you chose to rollback, verify that all operations were successfully rolled back
3. If manual intervention is needed, you can use the `schema_migration_operations` table to see which operations failed and why

## Security Note

The script uses the Supabase MCP tools to execute SQL queries. Ensure that your access token is kept secure and not committed to version control.

## Contact

For assistance with this migration script, please contact the CryoProtect development team.
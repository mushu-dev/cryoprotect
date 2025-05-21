# CryoProtect v2 - Foreign Key Relationship Fixes

This document provides information about the foreign key relationship fixes implemented in the CryoProtect v2 database to ensure data integrity.

## Overview

The `fix_foreign_key_relationships.py` script addresses several foreign key relationship issues:

1. **Proper REFERENCES Constraints**: Ensures all foreign keys have proper REFERENCES constraints
2. **Appropriate Indexes**: Creates indexes for all foreign key columns to improve query performance
3. **ON DELETE Behavior**: Sets appropriate ON DELETE CASCADE or ON DELETE SET NULL behavior based on relationship type
4. **Table Name Standardization**: Updates any constraints that reference singular table names to use plural names
5. **Transaction Support**: Implements transaction support for safer operations
6. **Error Handling**: Provides comprehensive error handling and reporting
7. **Verification**: Includes verification steps to ensure all relationships are correctly defined
8. **Rollback Mechanism**: Adds a rollback capability in case of failure

## Prerequisites

- Python 3.6+
- Supabase Python client (`pip install supabase`)
- Access to the CryoProtect Supabase project
- Environment variables set in a `.env` file:
  ```
  SUPABASE_URL=https://your-project-id.supabase.co
  SUPABASE_KEY=your-supabase-service-role-key
  ```

## Usage

### Using the Helper Scripts

For convenience, two helper scripts are provided:

#### Windows:
```
fix_foreign_key_relationships.bat
```

#### Linux/macOS:
```
./fix_foreign_key_relationships.sh
```

Both scripts provide an interactive menu with the following options:
1. Run with dry-run (show what would be done without making changes)
2. Verify only (check for issues without making changes)
3. Apply all fixes
4. Rollback to previous state

### Direct Command Line Usage

You can also run the Python script directly with various command-line options:

```bash
# Show what would be done without making changes
python fix_foreign_key_relationships.py --dry-run

# Only verify the foreign key constraints without making changes
python fix_foreign_key_relationships.py --verify-only

# Apply all fixes
python fix_foreign_key_relationships.py

# Rollback changes using the most recent backup
python fix_foreign_key_relationships.py --rollback

# Rollback changes using a specific backup file
python fix_foreign_key_relationships.py --rollback --backup-file fk_backup_20250418_120000.json
```

## Foreign Key Relationship Types

The script handles different types of foreign key relationships with appropriate ON DELETE behavior:

- **CASCADE**: When a parent record is deleted, all related child records are automatically deleted
  - Example: When a molecule is deleted, all its properties are deleted
  
- **SET NULL**: When a parent record is deleted, the foreign key in child records is set to NULL
  - Example: When a team is deleted, user profiles associated with that team have their team_id set to NULL

The default behavior is CASCADE for most relationships, but SET NULL is used for optional relationships.

## Verification Process

The verification process checks for:

1. Duplicate foreign key constraints
2. Outdated foreign key references (singular table names)
3. Foreign keys with incorrect ON DELETE behavior
4. Foreign keys without indexes
5. Missing foreign key constraints

## Backup and Rollback

Before making any changes, the script creates a backup of the current database state. If something goes wrong, you can rollback using:

```bash
python fix_foreign_key_relationships.py --rollback
```

This will:
1. Drop all current foreign key constraints
2. Recreate the foreign key constraints from the backup

## Troubleshooting

If you encounter issues:

1. Check the log file (e.g., `fix_foreign_keys_20250418_120000.log`) for detailed error messages
2. Run with `--verify-only` to check the current state
3. Use `--rollback` to revert changes if needed
4. Ensure your Supabase credentials are correct in the `.env` file

## After Running the Script

After successfully running the script, your database will have:

1. Proper foreign key constraints with appropriate ON DELETE behavior
2. Indexes on all foreign key columns for improved performance
3. Standardized table names (plural form)
4. Improved data integrity

## Contact

For assistance with this script, contact the CryoProtect development team.
# CryoProtect v2 - Database Remediation Manager

This document explains how to use the Database Remediation Manager to address all critical issues in the CryoProtect Supabase database.

## Overview

The Database Remediation Manager orchestrates the complete database remediation process for CryoProtect v2, following a 5-phase approach to address all critical issues:

### Phase 1: Security Lockdown (HIGHEST PRIORITY)
- Enable RLS on all tables
- Restrict anonymous access
- Create basic RLS policies
- Add indexes for policy columns

### Phase 2: Schema Standardization
- Standardize table names to plural form
- Update all foreign key references
- Document all changes in migration scripts

### Phase 3: Relationship Remediation
- Resolve fan traps by creating junction tables
- Update application code to use new structure
- Add appropriate indexes to junction tables

### Phase 4: Data Migration
- Consolidate duplicate tables
- Migrate data to normalized structure
- Validate data integrity after migration

### Phase 5: Performance Optimization
- Add missing indexes
- Create optimized views for common queries
- Implement connection pooling optimizations

## Prerequisites

1. Python 3.6+
2. Supabase project with appropriate permissions
3. Environment variables set in a `.env` file:
   ```
   SUPABASE_URL=https://your-project-id.supabase.co
   SUPABASE_KEY=your-service-role-key
   ```
4. Required Python packages:
   ```
   pip install supabase python-dotenv
   ```

## Usage

### Windows Batch Scripts

For Windows users, batch scripts are provided for easier execution:

```
run_database_remediation.bat                  - Run complete remediation
run_database_remediation.bat --dry-run        - Show what would be done without making changes
run_database_remediation.bat --verify-only    - Only verify the remediation
run_database_remediation.bat --phase [1-5]    - Run a specific phase
```

### Running the Complete Remediation Process

To run the complete remediation process (all 5 phases):

```bash
# Using Python directly
python database_remediation_manager.py

# Using the batch script (Windows)
run_database_remediation.bat
```

### Running a Specific Phase

To run a specific phase:

```bash
# Using Python directly
python database_remediation_manager.py --phase 1  # Run Phase 1: Security Lockdown

# Using the batch script (Windows)
run_database_remediation.bat --phase 1
```

### Running a Range of Phases

To run a range of phases:

```bash
python database_remediation_manager.py --start-phase 2 --end-phase 4  # Run Phases 2-4
```

### Dry Run Mode

To see what would be done without making changes:

```bash
# Using Python directly
python database_remediation_manager.py --dry-run

# Using the batch script (Windows)
run_database_remediation.bat --dry-run
```

### Verification Only Mode

To only verify if the remediation has been applied correctly:

```bash
# Using Python directly
python database_remediation_manager.py --verify-only

# Using the batch script (Windows)
run_database_remediation.bat --verify-only
```

### Testing the Remediation Manager

To test that all components of the remediation process are working correctly:

```bash
# Using Python directly
python test_database_remediation_manager.py

# Using the batch script (Windows)
test_database_remediation.bat
```

## Detailed Phase Information

### Phase 1: Security Lockdown

This phase addresses the critical security issues by:

1. Enabling Row Level Security (RLS) on all tables
2. Restricting anonymous access to only necessary tables
3. Creating RLS policies for:
   - Public read access (where is_public = true)
   - Owner full access (created_by = auth.uid())
   - Team member access (for users in the same team)
   - Service role bypass (for administrative operations)
4. Creating app-specific database roles with minimum permissions

Scripts used:
- `implement_security.py`
- `apply_service_role_rls.py`

### Phase 2: Schema Standardization

This phase standardizes the database schema by:

1. Converting all singular-named tables to plural form
2. Ensuring all tables use UUID as primary key with DEFAULT gen_random_uuid()
3. Adding proper REFERENCES constraints with indexes for all foreign keys

Scripts used:
- `standardize_schema.py`

### Phase 3: Relationship Remediation

This phase resolves relationship design issues by:

1. Replacing fan traps with proper junction tables
2. Applying 3NF normalization principles
3. Adding missing foreign key constraints with appropriate indexes

Scripts used:
- `fix_relationships.py`

### Phase 4: Data Migration

This phase ensures data integrity by:

1. Populating the molecular_property table with properties for the molecules
2. Restoring the mixture components for the existing mixtures
3. Verifying database integrity after migration

Scripts used:
- `complete_database_remediation.py`

### Phase 5: Performance Optimization

This phase optimizes database performance by:

1. Adding composite indexes for frequently joined columns
2. Creating text search indexes for name fields
3. Implementing other performance optimizations

Scripts used:
- `apply_performance_indexes.bat` (Windows) or `apply_performance_indexes.sh` (Linux/macOS)
- `verify_performance_indexes.bat` (Windows) or `verify_performance_indexes.sh` (Linux/macOS)

## Verification

Each phase includes a verification step to ensure that the remediation was applied correctly. You can run verification separately for each phase:

```bash
python database_remediation_manager.py --phase 1 --verify-only  # Verify Phase 1
```

Or verify all phases:

```bash
python database_remediation_manager.py --verify-only  # Verify all phases
```

## Logging

The remediation process creates a detailed log file with the format `database_remediation_YYYYMMDD_HHMMSS.log` that includes:

- All operations executed
- Success/failure status of each operation
- Error messages for failed operations
- Overall remediation status

## Troubleshooting

If you encounter issues during the remediation process:

1. Check the log file for detailed error messages
2. Run with `--dry-run` to see what would be executed without making changes
3. Run with `--verify-only` to check the current state
4. Try running individual phases to isolate the issue
5. Run the test script to verify all components are working correctly:
   ```
   test_database_remediation.bat  # Windows
   ```

### Windows-Specific Considerations

When running on Windows:

1. Use the provided `.bat` files for easier execution
2. Ensure Python is in your PATH environment variable
3. If you encounter permission issues, try running the Command Prompt as Administrator
4. For scripts that use Node.js (like apply_performance_indexes.bat), ensure Node.js is installed and in your PATH

## After Remediation

After successfully running the remediation process:

1. Update your application code to reference the new plural table names
2. Test your application thoroughly to ensure all functionality works with the new schema
3. Monitor database performance to ensure the optimizations are effective

## Contact

For assistance with this remediation process, please contact the CryoProtect development team.
# CryoProtect v2 - Database Remediation Plan

This document outlines the comprehensive database remediation plan for CryoProtect v2, addressing critical issues identified in the database structure and security.

## Critical Issues Addressed

The remediation plan addresses the following critical issues in order of priority:

### 1. SECURITY: Enable Row Level Security (RLS)

**Critical Issue:** Only 2/21 tables have RLS enabled. Anonymous users have full CRUD access.

**Solution:**
- Enable RLS on all tables
- Restrict anonymous access to only non-sensitive tables
- Create RLS policies for authenticated users
- Add performance indexes for RLS policies

### 2. STRUCTURE: Standardize Schema & Fix Relationships

**Critical Issues:** 
- Mixed singular/plural table names (`molecule`/`molecules`)
- Duplicate tables (`prediction`/`predictions`)
- "Fan trap" in `mixture` table (9 foreign keys referencing it)

**Solution:**
- Standardize all table names to plural form
- Update foreign key references
- Create junction tables to resolve fan traps
- Add appropriate indexes to junction tables

### 3. PERFORMANCE: Add Missing Indexes

**Critical Issue:** Missing indexes on foreign keys and RLS policy columns.

**Solution:**
- Add indexes to all foreign key columns
- Add specialized indexes for common query patterns
- Ensure all columns used in WHERE clauses are indexed

### 4. ROLES: Create Application-Specific Roles

**Critical Issue:** `postgres` and `service_role` bypass all RLS.

**Solution:**
- Create application-specific roles with minimum permissions
- Implement SECURITY DEFINER functions for admin operations
- Ensure proper role separation

### 5. DATA: Consolidate Duplicate Tables

**Critical Issue:** Duplicate tables with slightly different structures.

**Solution:**
- Identify differences between duplicate tables
- Migrate data to standardized tables
- Update application code to use consolidated tables

## Implementation Scripts

The remediation plan is implemented through the following scripts:

1. **complete_database_remediation.py**: The main Python script that implements all phases of the remediation plan.
2. **run_database_remediation.bat**: Windows batch script to run the remediation process.
3. **run_database_remediation.sh**: Linux/macOS shell script to run the remediation process.

## Usage

### Prerequisites

- Python 3.6 or higher
- Supabase Python client (`pip install supabase`)
- Environment variables set in a `.env` file:
  - `SUPABASE_URL`: Your Supabase project URL
  - `SUPABASE_KEY`: Your Supabase service role key

### Running the Remediation

#### Windows

```
run_database_remediation.bat
```

#### Linux/macOS

```
chmod +x run_database_remediation.sh
./run_database_remediation.sh
```

### Options

The remediation scripts support the following options:

- **Run all phases**: Execute all remediation phases in sequence
- **Dry-run mode**: Show what would be done without making changes
- **Run specific phase**: Execute only a specific phase of the remediation

## Verification

After running the remediation, the script will verify that:

- RLS is enabled on all tables
- Anonymous access is properly restricted
- All tables have appropriate RLS policies
- All foreign keys have corresponding indexes
- Application roles are created correctly

## Best Practices Enforced

The remediation enforces the following best practices:

1. **Schema Design**
   - Always use plural table names
   - UUIDs for primary keys: `id UUID PRIMARY KEY DEFAULT gen_random_uuid()`
   - Include `created_at` and `updated_at` TIMESTAMPTZ columns

2. **Security**
   - RLS on ALL tables
   - NO anonymous write access
   - Use security definer functions for service role access
   - Check all RLS policies work with actual user data

3. **Performance**
   - Index ALL foreign keys
   - Index ALL columns used in WHERE clauses
   - Use EXPLAIN ANALYZE to verify query performance

4. **Implementation**
   - Security fixes FIRST
   - Make small, testable changes
   - Create test cases BEFORE implementing
   - Document ALL changes in migration scripts
   - Have clear rollback plan for each change

## Rollback Plan

If issues are encountered during the remediation process, you can:

1. Run the script with the `--rollback` option to revert changes
2. Restore from a database backup
3. Manually revert specific changes using the SQL commands in the verification log

## Logging

The remediation process generates detailed logs in the following format:
- `complete_remediation_YYYYMMDD_HHMMSS.log`: Detailed log of the remediation process

## Support

If you encounter any issues during the remediation process, please contact the database administrator or open an issue in the project repository.
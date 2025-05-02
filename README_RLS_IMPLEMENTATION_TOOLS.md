# CryoProtect RLS Implementation Tools

This document describes the Row Level Security (RLS) implementation tools provided in the CryoProtect v2 project. These tools ensure proper data access control through PostgreSQL's RLS feature.

## Overview

The Row Level Security (RLS) implementation tools consist of:

1. **Migration Script**: SQL script to apply comprehensive RLS policies to all tables
2. **Application Script**: Python script to execute the migration and verify policies
3. **Verification Script**: Python script to test the effectiveness of RLS policies

These tools provide a complete solution for implementing, applying, and verifying RLS policies in the CryoProtect Supabase database.

## Migration Script

**File**: `migrations/018_complete_rls_policies.sql`

This SQL migration script:

- Enables RLS on all tables in the database
- Creates standard RLS policy functions for different operations (SELECT, INSERT, UPDATE, DELETE)
- Applies appropriate policies to all tables based on their schema
- Includes special handling for tables with custom access requirements (e.g., lab_verifications)
- Adds service role bypass policies to allow admin operations

The script is designed to be idempotent and can be run multiple times without error.

## Application Script

**File**: `apply_missing_rls_policies.py`

This Python script:

- Applies the RLS migration to the database
- Verifies all tables have RLS enabled with appropriate policies
- Checks views for SECURITY INVOKER settings
- Generates detailed verification reports

### Usage

```bash
# Apply RLS policies with default settings
python apply_missing_rls_policies.py

# Skip verification step
python apply_missing_rls_policies.py --skip-verification

# Use a different SQL file
python apply_missing_rls_policies.py --sql-file path/to/custom_rls_policies.sql
```

### Output

The script generates:

1. Logs in `logs/apply_rls_policies.log`
2. Verification JSON report in `reports/auth/rls_verification_summary_[timestamp].json`
3. Human-readable summary in `reports/auth/RLS_Verification_Summary_[timestamp].md`

## Verification Script

**File**: `verify_rls_effectiveness.py`

This script provides comprehensive testing of RLS policies by:

1. Simulating access as different user roles (service_role, admin, regular, anonymous)
2. Testing all CRUD operations (SELECT, INSERT, UPDATE, DELETE) on each table
3. Measuring query performance impact of RLS policies
4. Generating detailed effectiveness reports

### Usage

```bash
# Verify RLS effectiveness with default settings
python verify_rls_effectiveness.py

# Specify a project ID
python verify_rls_effectiveness.py --project-id your-project-id
```

### Output

The script generates:

1. Logs in `logs/verify_rls_effectiveness.log`
2. Effectiveness JSON report in `reports/security/rls_effectiveness_report_[timestamp].json`
3. Human-readable report in `reports/security/rls_effectiveness_report_[timestamp].md`

## Policy Structure

The RLS implementation follows this structure:

### Standard Policies

Each table receives four standard policies:

1. **Select Policy**: Controls who can read data
   - Owner can see their own data
   - Team members with view permissions can see data
   - Users with shared access can see data

2. **Insert Policy**: Controls who can create data
   - Owner can create data
   - Team members with create permissions can create data

3. **Update Policy**: Controls who can modify data
   - Owner can update their own data
   - Team members with update permissions can update data

4. **Delete Policy**: Controls who can delete data
   - Owner can delete their own data
   - Team members with delete permissions can delete data

### Service Role Bypass

Each table also receives a service role bypass policy:

```sql
CREATE POLICY [table_name]_service_role_policy ON [table_name]
USING (auth.role() = 'service_role')
WITH CHECK (auth.role() = 'service_role');
```

This policy allows the service role to bypass RLS restrictions for administrative operations.

## Troubleshooting

### Common Issues

1. **Missing RLS Policies**:
   - Ensure the migration script ran successfully
   - Check database connection parameters
   - Look for error messages in the log file

2. **RLS Blocking Legitimate Access**:
   - Verify user roles and permissions
   - Check the session context values (auth.uid, auth.role)
   - Test with service_role to bypass RLS

3. **Performance Impact**:
   - Review effectiveness reports for high overhead queries
   - Consider optimizing RLS expressions or adding indexes
   - Use service_role for bulk operations

### Testing RLS Manually

You can manually test RLS policies with:

```sql
-- Set role to test with
SET ROLE service_role;
-- OR
SET LOCAL auth.uid = 'user-uuid';
SET LOCAL app.user_role = 'admin';

-- Try operations
SELECT * FROM your_table;
INSERT INTO your_table (...) VALUES (...);
UPDATE your_table SET ... WHERE id = ...;
DELETE FROM your_table WHERE id = ...;

-- Reset role
RESET ROLE;
RESET auth.uid;
```

## Additional Resources

- [PostgreSQL RLS Documentation](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
- [Supabase RLS Guide](https://supabase.com/docs/guides/auth/row-level-security)
- Original RLS implementations in `migrations/006_rls_policies.sql`
- Enhanced RLS policies in `migrations/007_service_role_rls.sql`
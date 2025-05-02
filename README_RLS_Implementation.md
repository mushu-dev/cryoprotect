# CryoProtect v2 - Row Level Security (RLS) Implementation

This document explains the Row Level Security (RLS) implementation for the CryoProtect v2 database and how to use the `fix_rls_implementation.py` script to ensure consistent RLS policies across all tables.

## Overview

The `fix_rls_implementation.py` script implements comprehensive Row Level Security (RLS) policies across all tables in the CryoProtect v2 database. This is a critical security requirement that ensures data is properly protected and only accessible to authorized users.

### Key Features

1. **Enables RLS on all tables** in the schema
2. **Creates comprehensive RLS policies** for each table type:
   - Public read access policy (for records with `is_public = true`)
   - Owner full access policy (for records created by the current user)
   - Team member access policy (for records associated with the user's team)
   - Service role bypass policy (for service role operations)
3. **Restricts anonymous access** to sensitive data
4. **Adds indexes for policy columns** to optimize performance
5. **Implements transaction support** for safer operations
6. **Provides detailed error handling and reporting**

## Usage

### Prerequisites

1. Make sure you have the Supabase URL and service role key in your `.env` file:
   ```
   SUPABASE_URL=your-supabase-url
   SUPABASE_KEY=your-service-role-key
   ```

2. Install required Python packages:
   ```bash
   pip install supabase python-dotenv
   ```

### Running the Script

You can run the script directly or use the provided batch file:

```bash
# Using Python directly
python fix_rls_implementation.py [options]

# Using the batch file
fix_rls_implementation.bat [options]
```

### Command-Line Options

The script supports the following command-line options:

- `--dry-run`: Show what would be done without making changes
- `--verify-only`: Only verify the current RLS implementation
- `--rollback`: Rollback changes to the original state
- `--schema SCHEMA`: Database schema to operate on (default: remediation_test)

### Examples

```bash
# Run the fix with default settings
python fix_rls_implementation.py

# Show what would be done without making changes
python fix_rls_implementation.py --dry-run

# Only verify the current RLS implementation
python fix_rls_implementation.py --verify-only

# Rollback changes to the original state
python fix_rls_implementation.py --rollback

# Run on a different schema
python fix_rls_implementation.py --schema public
```

## RLS Policies Explained

The script creates the following RLS policies for each table:

### 1. Public Read Access Policy

Allows anyone to read records that are marked as public:

```sql
CREATE POLICY "Public read access on [table]" 
ON [schema].[table] 
FOR SELECT 
USING (is_public = true);
```

This policy is only created for tables that have an `is_public` column.

### 2. Owner Full Access Policy

Allows users to have full access to records they created:

```sql
CREATE POLICY "Owner full access on [table]" 
ON [schema].[table] 
USING (created_by = auth.uid());
```

This policy is only created for tables that have a `created_by` column.

### 3. Team Member Access Policy

Allows team members to access records associated with their team:

```sql
CREATE POLICY "Team member access on [table]" 
ON [schema].[table] 
FOR SELECT 
USING (
  EXISTS (
    SELECT 1 FROM [schema].projects
    JOIN [schema].teams ON teams.id = projects.team_id
    JOIN [schema].user_profile ON user_profile.team_id = teams.id
    WHERE projects.id = [table].project_id
    AND user_profile.auth_user_id = auth.uid()
  )
);
```

The exact condition varies depending on the table structure and how it relates to teams.

### 4. Service Role Bypass Policy

Allows the service role to bypass RLS restrictions:

```sql
CREATE POLICY "Service role bypass on [table]" 
ON [schema].[table] 
USING (auth.role() = 'service_role');
```

## Performance Optimization

The script creates indexes on columns used in RLS policies to optimize performance:

1. Index on `user_profile.auth_user_id`
2. Indexes on `team_id` columns
3. Indexes on `project_id` columns
4. Indexes on `created_by` columns
5. Indexes on `is_public` columns

## Verification

The script includes a verification process that checks:

1. RLS is enabled on all tables
2. All required policies are created based on table structure
3. Indexes are created for better performance

Verification results are saved to a JSON file for review.

## Rollback

If you need to rollback the changes, use the `--rollback` option:

```bash
python fix_rls_implementation.py --rollback
```

This will:

1. Drop all policies created by the script
2. Disable RLS on tables that had it disabled before
3. Recreate original policies

## Logging

The script logs all actions to:

1. Console output
2. A log file named `rls_fix_YYYYMMDD_HHMMSS.log`

## Troubleshooting

If you encounter issues:

1. Check the log file for error messages
2. Run with `--dry-run` to see what would be executed
3. Use `--verify-only` to check the current state
4. If needed, use `--rollback` to revert changes

## Security Best Practices

1. **Principle of Least Privilege**: Users only have access to the data they need
2. **Defense in Depth**: Multiple layers of security (RLS, policies)
3. **Fail-Safe Defaults**: By default, access is denied unless explicitly granted
4. **Separation of Concerns**: Different policies for different access patterns
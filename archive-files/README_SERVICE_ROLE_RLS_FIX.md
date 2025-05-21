# CryoProtect v2 - Service Role RLS Fix

This document explains the Row-Level Security (RLS) policy fix that was implemented to allow the service role to insert data into tables with RLS enabled.

## Problem

The original RLS policies in `migrations/006_rls_policies.sql` only allow project members to insert data into tables like `molecule`, `mixture`, `prediction`, and `experiment`. This was causing errors when trying to run scripts like `populate_molecules.py` with the service role, resulting in errors like:

```
new row violates row-level security policy for table 'molecule'
```

## Solution

We've implemented two key changes to fix this issue:

1. Added RLS policies that specifically allow the service role to insert data into tables with RLS enabled
2. Modified the `populate_molecules.py` script to use the `service_role_helper.py` module

These changes allow scripts to run with the service role key without requiring project membership.

## Files Created

1. **migrations/007_service_role_rls.sql**: Contains the RLS policies for the service role
2. **apply_service_role_rls.py**: Script to apply the migration

## How It Works

The service role RLS policies work by:

1. Adding a new policy to each table with RLS enabled
2. The policy allows INSERT operations when `auth.role() = 'service_role'`
3. This bypasses the project membership check for the service role

## Applying the Migration

There are two ways to apply the migration:

### Method 1: Manual Application (Recommended)

1. Log in to the [Supabase Dashboard](https://app.supabase.com/)
2. Select your project
3. Go to the SQL Editor
4. Copy the contents of `apply_service_role_rls_manual.sql`
5. Paste into the SQL Editor and run the query

### Method 2: Using the Python Script (May not work in all environments)

```bash
# Apply the migration
python apply_service_role_rls.py

# To see what would be executed without actually applying it
python apply_service_role_rls.py --dry-run
```

Note: The Python script method may fail if the `exec_sql` function is not available in your Supabase project.

## Tables Affected

The following tables now have service role insert policies:

1. molecule
2. mixture
3. mixture_component
4. experiment
5. molecular_property
6. prediction
7. experiment_property

## Verifying the Fix

After applying the migration and updating the scripts, you can verify that the fix works by running:

```bash
python populate_molecules.py
```

This should now complete without any RLS policy violations. The script will:
1. Connect to Supabase using the service role
2. Use the hardcoded user ID from `auth_config.py`
3. Ensure a user profile exists
4. Insert data into the tables with the service role

## Script Modifications

The `populate_molecules.py` script has been modified to:

1. Import functions from `service_role_helper.py`
2. Use `get_supabase_client()` to connect with the service role
3. Use `get_user_id()` to get the hardcoded user ID
4. Use `ensure_user_profile()` to create a user profile if needed

These modifications ensure that the script uses the service role correctly and can bypass the RLS policies that restrict project members.

## Troubleshooting

If you encounter any issues:

1. Verify that the service role key is correctly set in your `.env` file
2. Check that the migration was applied successfully by looking at the log file
3. Verify that the policies were created by querying the `pg_policies` table
4. Try running with `--dry-run` to see the SQL that would be executed

## Additional Notes

- This fix complements the service role authentication approach described in `README_SERVICE_ROLE_FIX.md`
- The service role approach should be used for development and data population purposes only
- In production, proper user authentication and authorization should be used
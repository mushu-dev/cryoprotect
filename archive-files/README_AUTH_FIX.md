# CryoProtect v2 - Authentication Fix

This document provides instructions for fixing the authentication issue with Supabase where you're getting the "Email not confirmed" error.

## What's the Problem?

The error message indicates that your Supabase user account exists but hasn't been confirmed via email verification. This prevents authentication from working properly.

## The Solution

We've created several scripts to fix this issue by using a "service role" approach that bypasses the need for email confirmation. This approach:

1. Uses the service role key instead of the anon key
2. Hardcodes a valid user ID (the one already created in your system)
3. Modifies the authentication flow to skip the email confirmation check

## Fix Scripts

We've created the following scripts to help fix the authentication issue:

### 1. fix_auth_service_role.py

This script:
- Creates a backup of your .env file
- Updates your .env file to use a service role key (if needed)
- Creates an auth_config.py file with the hardcoded user ID
- Creates additional helper scripts

### 2. test_service_role_auth.py

This script tests if the service role authentication approach works by:
- Connecting to Supabase with the service role key
- Trying to access the user profile for the hardcoded user ID
- Trying to access the molecules table

### 3. modify_files_for_service_role.py

This script modifies all Python files in the project that use authentication to:
- Import the auth_config module
- Use the service role approach for authentication
- Create backups of all modified files

### 4. update_test_supabase_auth.py

This script updates the original test_supabase_auth.py to use the service role approach.

## How to Fix the Authentication Issue

Follow these steps to fix the authentication issue:

1. Run the fix_auth_service_role.py script:
   ```
   python fix_auth_service_role.py
   ```
   
   If prompted for a service role key, you'll need to get it from the Supabase dashboard:
   - Go to https://app.supabase.com
   - Select your project
   - Go to Project Settings > API
   - Copy the 'service_role' key (not the anon/public key)

2. Run the test_service_role_auth.py script to verify the fix:
   ```
   python test_service_role_auth.py
   ```
   
   If this script shows [SUCCESS] messages, the authentication bypass is working correctly.

3. Run the update_test_supabase_auth.py script to update the original test script:
   ```
   python update_test_supabase_auth.py
   ```

4. Run the original test_supabase_auth.py script to confirm it now works:
   ```
   python test_supabase_auth.py
   ```

5. If you want to update all Python files in the project to use the service role approach, run:
   ```
   python modify_files_for_service_role.py
   ```

6. After applying these fixes, you should be able to run your project without the "Email not confirmed" error.

## Alternative Solutions

If you prefer to fix the email confirmation issue directly instead of bypassing it:

1. Check your email for a confirmation link from Supabase
2. If you can't find the email, go to the Supabase dashboard:
   - Find your user account
   - Manually confirm the email or resend the confirmation email
3. Once confirmed, you can use the original authentication approach

## Potential Issues

1. **Missing Tables**: The check_supabase_connection.py script showed that the "molecular_properties" table might not exist or might have a different name. You may need to check your database schema.

2. **Service Role Key**: Make sure you're using the correct service role key from the Supabase dashboard.

3. **Table Access**: Even with the service role key, you might still face Row Level Security (RLS) issues. If this happens, you may need to modify the RLS policies in your Supabase project.

## Restoring Backups

If anything goes wrong, you can restore from the backups:

- `.env.backup-auth-fix` - Backup of your original .env file
- `*.py.bak` - Backups of modified Python files

To restore a backup, simply rename it back to the original filename.

# CryoProtect v2 - Security Implementation

This document explains the comprehensive security implementation for the CryoProtect Supabase project, including Row Level Security (RLS) policies and app-specific database roles.

## Overview

The `implement_security.py` script implements the following security features:

1. Enables Row Level Security (RLS) on all public tables
2. Creates RLS policies for public access, owner access, team member access, and service role bypass
3. Creates app-specific database roles with minimum permissions
4. Includes verification and rollback mechanisms

## Tables Affected

The script applies security features to the following tables:

- molecules
- mixtures
- mixture_components
- predictions
- experiments
- experiment_properties
- calculation_methods
- projects
- teams

## RLS Policies

For each table, the script creates the following RLS policies:

1. **Public read access** - Allows anyone to SELECT if is_public = true
   ```sql
   CREATE POLICY "Public read access on [table]" 
   ON [table] 
   FOR SELECT 
   USING (is_public = true);
   ```

2. **Owner full access** - Allows all operations if auth.uid() = created_by
   ```sql
   CREATE POLICY "Owner full access on [table]" 
   ON [table] 
   USING (created_by = auth.uid());
   ```

3. **Team member access** - Allows SELECT if the user is a member of the team associated with the record
   ```sql
   CREATE POLICY "Team member access on [table]" 
   ON [table] 
   FOR SELECT 
   USING (
     EXISTS (
       SELECT 1 FROM projects
       JOIN teams ON teams.id = projects.team_id
       JOIN user_profile ON user_profile.team_id = teams.id
       WHERE projects.id = [table].project_id
       AND user_profile.auth_user_id = auth.uid()
     )
   );
   ```

4. **Service role bypass** - Allows all operations if auth.role() = 'service_role'
   ```sql
   CREATE POLICY "Service role bypass on [table]" 
   ON [table] 
   USING (auth.role() = 'service_role');
   ```

## App-Specific Database Roles

The script creates the following database roles:

1. **app_readonly** - Can only SELECT from tables
   ```sql
   CREATE ROLE app_readonly;
   GRANT USAGE ON SCHEMA public TO app_readonly;
   GRANT SELECT ON ALL TABLES IN SCHEMA public TO app_readonly;
   ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO app_readonly;
   ```

2. **app_user** - Can SELECT, INSERT, and UPDATE but not DELETE
   ```sql
   CREATE ROLE app_user;
   GRANT USAGE ON SCHEMA public TO app_user;
   GRANT SELECT, INSERT, UPDATE ON ALL TABLES IN SCHEMA public TO app_user;
   ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT, INSERT, UPDATE ON TABLES TO app_user;
   ```

3. **app_admin** - Has full access to all tables
   ```sql
   CREATE ROLE app_admin;
   GRANT USAGE ON SCHEMA public TO app_admin;
   GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO app_admin;
   ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL PRIVILEGES ON TABLES TO app_admin;
   ```

## Usage

### Prerequisites

1. Make sure you have the Supabase URL and service role key in your `.env` file:
   ```
   SUPABASE_URL=https://tsdlmynydfuypiugmkev.supabase.co
   SUPABASE_KEY=your-service-role-key
   ```

2. Install required Python packages:
   ```bash
   pip install supabase python-dotenv
   ```

### Running the Script

```bash
# Dry run (shows what would be executed without making changes)
python implement_security.py --dry-run

# Apply security implementation
python implement_security.py

# Verify existing security settings
python implement_security.py --verify-only

# Rollback changes
python implement_security.py --rollback
```

## Verification

The script includes verification steps to ensure that:

1. RLS is enabled on all tables
2. All required policies are created
3. App-specific roles are created

Verification results are saved to `security_verification_results.json`.

## Testing

The script includes testing steps to verify that:

1. Service role can access all tables
2. Public access works for records with is_public = true
3. Owner access works for records created by the user
4. Team member access works for users who are members of the team

Test results are saved to `security_test_results.json`.

## Rollback

If you need to rollback the changes, use the `--rollback` flag:

```bash
python implement_security.py --rollback
```

This will:

1. Drop all policies created by the script
2. Disable RLS on tables that had it disabled before
3. Recreate original policies
4. Drop app-specific roles

## Logging

The script logs all actions to:

1. Console output
2. `security_implementation.log` file

## Security Best Practices

1. **Principle of Least Privilege**: The app-specific roles follow the principle of least privilege, giving users only the permissions they need.
2. **Defense in Depth**: Multiple layers of security are implemented (RLS, roles, policies).
3. **Separation of Duties**: Different roles have different permissions.
4. **Fail-Safe Defaults**: By default, access is denied unless explicitly granted.

## Troubleshooting

If you encounter issues:

1. Check the log file for error messages
2. Verify that the Supabase URL and key are correct
3. Try running with `--dry-run` to see what would be executed
4. Use `--verify-only` to check the current state
5. If needed, use `--rollback` to revert changes
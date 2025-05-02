# Rollback Procedures

This guide provides detailed instructions for rolling back changes made by the CryoProtect fix scripts if issues are encountered.

## Overview

Despite careful planning and testing, sometimes changes need to be rolled back due to unexpected issues or conflicts. The CryoProtect project includes built-in rollback mechanisms for each fix script, as well as procedures for manual rollback if needed.

This guide covers:

1. When to consider rolling back changes
2. Using the built-in rollback mechanisms
3. Manual rollback procedures
4. Restoring from backups

## When to Consider Rollback

Consider rolling back changes in the following situations:

1. **Verification Failure**: If verification steps fail after applying fixes
2. **Application Errors**: If the application shows errors or unexpected behavior after fixes
3. **Performance Issues**: If performance degrades significantly after applying fixes
4. **Data Integrity Issues**: If data appears corrupted or inconsistent after fixes
5. **Integration Problems**: If the system fails to integrate with other components after fixes

## Using Built-in Rollback Mechanisms

### Automatic Rollback with the Master Script

The master integration script (`run_cryoprotect_fixes.py`) includes an automatic rollback feature that can be enabled with the `--rollback` flag:

```bash
python run_cryoprotect_fixes.py --rollback
```

With this option, if any fix script fails, the script will automatically attempt to roll back the changes made by that script before exiting.

### Rolling Back Individual Fixes

Each fix script includes its own rollback mechanism that can be invoked with the `--rollback` flag:

#### 1. Schema Standardization Rollback

```bash
python standardize_schema.py --rollback
```

This will:
- Drop the new plural tables
- Restore the original singular tables
- Remove tracking tables created during the migration

#### 2. Security Implementation Rollback

```bash
python implement_security.py --rollback
```

This will:
- Drop all policies created by the script
- Disable RLS on tables that had it disabled before
- Recreate original policies
- Drop app-specific roles

#### 3. Relationship Fixes Rollback

```bash
python fix_relationships.py --rollback
```

This will:
- Drop the newly created junction tables
- Restore the original tables to their previous state
- Remove any constraints added during the fix

#### 4. API Integration Fixes Rollback

```bash
python fix_api_integration.py --rollback
```

This will:
- Restore the original API endpoint registrations
- Revert any changes to the database configuration

## Manual Rollback Procedures

If the built-in rollback mechanisms fail or are not available, you can perform a manual rollback using the following procedures:

### 1. Manual Schema Rollback

To manually roll back schema changes:

1. Connect to the Supabase database using the SQL Editor or psql
2. Drop the new plural tables:
   ```sql
   DROP TABLE IF EXISTS public.molecules CASCADE;
   DROP TABLE IF EXISTS public.mixtures CASCADE;
   DROP TABLE IF EXISTS public.predictions CASCADE;
   DROP TABLE IF EXISTS public.experiments CASCADE;
   DROP TABLE IF EXISTS public.experiment_properties CASCADE;
   DROP TABLE IF EXISTS public.mixture_components CASCADE;
   DROP TABLE IF EXISTS public.calculation_methods CASCADE;
   DROP TABLE IF EXISTS public.property_types CASCADE;
   DROP TABLE IF EXISTS public.projects CASCADE;
   DROP TABLE IF EXISTS public.teams CASCADE;
   ```
3. Recreate the original singular tables from the backup

### 2. Manual Security Rollback

To manually roll back security changes:

1. Connect to the Supabase database using the SQL Editor or psql
2. Drop all RLS policies:
   ```sql
   DROP POLICY IF EXISTS "Public read access on molecules" ON public.molecules;
   DROP POLICY IF EXISTS "Owner full access on molecules" ON public.molecules;
   DROP POLICY IF EXISTS "Team member access on molecules" ON public.molecules;
   DROP POLICY IF EXISTS "Service role bypass on molecules" ON public.molecules;
   -- Repeat for all tables
   ```
3. Disable RLS on all tables:
   ```sql
   ALTER TABLE public.molecules DISABLE ROW LEVEL SECURITY;
   ALTER TABLE public.mixtures DISABLE ROW LEVEL SECURITY;
   -- Repeat for all tables
   ```
4. Drop app-specific roles:
   ```sql
   DROP ROLE IF EXISTS app_readonly;
   DROP ROLE IF EXISTS app_user;
   DROP ROLE IF EXISTS app_admin;
   ```

### 3. Manual Relationship Rollback

To manually roll back relationship fixes:

1. Connect to the Supabase database using the SQL Editor or psql
2. Drop the junction tables:
   ```sql
   DROP TABLE IF EXISTS public.molecule_proteins CASCADE;
   DROP TABLE IF EXISTS public.molecule_experiments CASCADE;
   ```
3. Revert changes to the predictions table:
   ```sql
   ALTER TABLE public.predictions DROP COLUMN IF EXISTS molecule_id;
   ALTER TABLE public.predictions ALTER COLUMN mixture_id SET NOT NULL;
   ```
4. Revert changes to the molecular_properties table:
   ```sql
   ALTER TABLE public.molecular_properties DROP COLUMN IF EXISTS property_type_id;
   ```

### 4. Manual API Integration Rollback

To manually roll back API integration fixes:

1. Restore the original API files from backup:
   - `api/__init__.py`
   - `api/resources.py`
2. Restart the application to apply the changes

## Restoring from Backups

### Database Backup Restoration

The CryoProtect system automatically creates database backups before applying fixes. To restore from a backup:

1. Locate the backup file in the `backups` directory
2. Use the Supabase SQL Editor or psql to restore the backup:
   ```sql
   -- First, drop the current tables
   DROP SCHEMA public CASCADE;
   CREATE SCHEMA public;
   
   -- Then, restore from backup
   \i path/to/backup/file.sql
   ```

### File Backup Restoration

The master integration script also creates backups of modified files. To restore from these backups:

1. Locate the backup directory:
   ```
   backups/master_backup_YYYYMMDD_HHMMSS/
   ```
2. Copy the backup files to their original locations:
   ```bash
   cp backups/master_backup_YYYYMMDD_HHMMSS/api/__init__.py api/__init__.py
   cp backups/master_backup_YYYYMMDD_HHMMSS/api/resources.py api/resources.py
   # Repeat for all modified files
   ```

## Verification After Rollback

After rolling back changes, it's important to verify that the system has been restored to its previous state:

1. **Database Verification**:
   - Verify that the original tables exist with the correct structure
   - Check that data is intact and accessible

2. **Application Verification**:
   - Start the application and verify that it runs without errors
   - Test basic functionality to ensure it works as expected

3. **API Verification**:
   - Test API endpoints to ensure they respond correctly
   - Verify that data can be retrieved and modified through the API

## Rollback Logging

All rollback operations are logged to help diagnose issues and track changes:

- Master script rollback: `cryoprotect_fixes_rollback_YYYYMMDD_HHMMSS.log`
- Schema rollback: `schema_rollback_YYYYMMDD_HHMMSS.log`
- Security rollback: `security_rollback.log`
- Relationship rollback: `relationship_rollback.log`
- API rollback: `api_rollback.log`

Review these logs to understand what operations were performed during the rollback.

## Best Practices for Rollback

Follow these best practices when performing rollbacks:

1. **Always create backups** before applying fixes
2. **Test rollback procedures** in a non-production environment
3. **Document all changes** made during fixes and rollbacks
4. **Verify system functionality** after rollback
5. **Communicate with users** about the rollback and any potential impact
6. **Analyze the root cause** of issues that led to the rollback
7. **Develop a plan** to address the issues before attempting fixes again

## Conclusion

The CryoProtect system includes comprehensive rollback mechanisms to ensure that changes can be safely undone if issues are encountered. By following the procedures outlined in this guide, you can effectively roll back changes and restore the system to a known good state.
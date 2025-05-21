# Supabase RLS Implementation Guide

This guide explains how to implement Row Level Security (RLS) for CryoProtect v2 directly in the Supabase SQL Editor.

## Overview

We've created a standalone SQL script that can be executed directly in the Supabase SQL Editor to implement all the necessary RLS policies for missing tables and views, with performance optimizations for scientific data.

## Key Features

- **Single Transaction**: All changes are applied in a single atomic transaction
- **Idempotent Execution**: Safe to run multiple times without duplicating policies
- **Performance Optimizations**: Includes indexes for RLS columns
- **Scientific Data Auditing**: Comprehensive audit trail for scientific data
- **Self-Verification**: Includes verification query to validate implementation

## Implementation Steps

1. **Log in to Supabase**
   - Go to your Supabase project dashboard
   - Navigate to the SQL Editor section

2. **Create a New Query**
   - Click "New Query" 
   - Give it a name like "Apply RLS Policies"

3. **Copy SQL Script**
   - Copy the entire contents of `supabase_rls_policies.sql`
   - Paste it into the SQL Editor

4. **Review the Script**
   - The script is organized into sections:
     1. Performance indexes for RLS
     2. RLS for molecule_with_properties view
     3. RLS for mixture_with_components view
     4. RLS for migrations table
     5. Creation/RLS for experiment_with_results view
     6. Scientific data audit trail implementation
     7. RLS for the audit table
     8. Verification queries

5. **Execute the Script**
   - Click "Run" to execute the entire script in a single transaction
   - The script will automatically run verification at the end

6. **Review Results**
   - The verification results will show:
     - Tables with RLS enabled and policy counts
     - Views with SECURITY INVOKER status
     - Implementation summary with index and trigger counts

## Understanding the Implementation

### Transaction Management

The script uses a single transaction to ensure all changes are applied atomically:

```sql
BEGIN;
-- All operations
COMMIT;
```

If any part fails, the entire operation will be rolled back, leaving the database in its original state.

### Scientific Data Auditing

The script implements a comprehensive audit trail system:

```sql
CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    table_name TEXT NOT NULL,
    record_id UUID NOT NULL,
    operation TEXT NOT NULL,
    -- Additional fields
);
```

Audit triggers are applied to key scientific data tables (molecule, mixture, experiment).

### Performance Optimization

Performance is enhanced through strategic indexing:

```sql
CREATE INDEX IF NOT EXISTS idx_user_profile_project_id ON user_profile(project_id);
CREATE INDEX IF NOT EXISTS idx_experiment_project_id ON experiment(project_id);
-- Additional indexes
```

These indexes significantly improve RLS policy performance for scientific data queries.

### Bulk Loading Support

To optimize performance during initial data loading, you can bypass auditing:

```sql
-- Before bulk loading
SELECT set_audit_bypass(true);

-- Perform bulk loading operations

-- After bulk loading
SELECT set_audit_bypass(false);
```

## Troubleshooting

If you encounter issues:

1. **Execution Errors**
   - Check error messages in the SQL Editor
   - Ensure you have the necessary permissions (must be admin/owner)
   - Verify table/view existence before running

2. **Missing Policies**
   - Run `SELECT * FROM pg_policies;` to list all policies
   - Check if tables have RLS enabled with `\d table_name` (look for RLS: enabled)

3. **Performance Issues**
   - Check if indexes were created with `SELECT * FROM pg_indexes;`
   - Monitor query performance with `EXPLAIN ANALYZE` on typical queries

4. **Data Loading Issues**
   - Use the audit bypass function for bulk operations
   - For service role operations, ensure proper authentication

## Next Steps

After implementing RLS:

1. **Test Access Control**
   - Verify different user roles can only access appropriate data
   - Test service role access for data loading

2. **Monitor Performance**
   - Check query plans with EXPLAIN ANALYZE 
   - Monitor query performance with Supabase dashboard

3. **Proceed with Database Population**
   - Use service role authentication for data loading
   - Consider bypassing audit for bulk operations

4. **Complete Phase 3.1 Implementation**
   - Finalize CI/CD integration
   - Implement environment configuration
   - Set up blue/green deployment
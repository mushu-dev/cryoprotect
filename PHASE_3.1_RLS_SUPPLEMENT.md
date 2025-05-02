# Phase 3.1 RLS Implementation Supplement

This document provides additional information for Phase 3.1 (Deployment Infrastructure) regarding the implementation of Row Level Security (RLS) for missing tables and views.

## Missing RLS Policies

During our preparation for deployment, we identified several tables and views that were missing proper RLS policies:

1. `experiment_with_results` (view)
2. `migrations` (table)
3. `mixture_with_components` (view)
4. `molecule_with_properties` (view)

These missing policies could potentially expose data to unauthorized users and needed to be addressed before database population.

## Implementation Files

We have created the following files to implement the missing RLS policies:

1. `migrations/missing_rls_policies.sql` - SQL migration file with RLS policies
2. `apply_missing_rls_policies.py` - Python script to apply and verify policies
3. `apply_missing_rls_policies.bat` - Windows batch script
4. `apply_missing_rls_policies.sh` - Unix/Linux shell script
5. `README_MISSING_RLS_POLICIES.md` - Documentation

## Integration with Database Population

These RLS policies must be applied **before** populating the database to ensure that service role permissions are properly configured. The recommended sequence is:

1. Apply base RLS policies from migration files
2. Apply missing RLS policies using the provided scripts
3. Configure service role authentication
4. Populate the database with scientific data

## CI/CD Integration

For the CI/CD pipeline in Phase 3.1, we recommend:

1. Include a step to apply missing RLS policies before database population
2. Add verification to ensure RLS is properly configured
3. Include automated testing to verify that RLS is correctly enforcing access controls

## Security Best Practices

When implementing RLS, follow these best practices:

1. Use SECURITY INVOKER for views to respect underlying table RLS
2. Create explicit policies for all CRUD operations
3. Ensure service role access for data population scripts
4. Test access with different user roles to verify isolation
5. Document all security assumptions and decisions

## Action Items

To fully implement the missing RLS policies:

1. Run the `apply_missing_rls_policies` script
2. Verify RLS in the Supabase dashboard
3. Test access with different user accounts
4. Proceed with database population
5. Include RLS verification in deployment checklist
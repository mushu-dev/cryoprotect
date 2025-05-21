# Row Level Security (RLS) Policy Enhancements

This document describes the enhancements made to the Row Level Security (RLS) policies in the CryoProtect database.

## Overview

The RLS policy enhancements focus on:

1. Creating reusable helper functions for access control
2. Consolidating duplicate policies
3. Implementing clearance-level based access control
4. Adding comprehensive testing and verification tools
5. Optimizing policy performance

## Helper Functions

Several helper functions have been implemented in the `auth` schema:

### `auth.has_access_to_molecule(molecule_id UUID) → BOOLEAN`

Checks if the current user has access to a specific molecule. Access is granted if:
- The molecule is public
- The user is the molecule owner
- The user is an admin
- The user is a service role
- The user is part of a team that can access the molecule

### `auth.has_access_to_mixture(mixture_id UUID) → BOOLEAN`

Checks if the current user has access to a specific mixture. Access is granted if:
- The mixture is public
- The user is the mixture owner
- The user is an admin
- The user is a service role
- The user has access to the project the mixture belongs to

### `auth.get_user_clearance_level() → TEXT`

Returns the user's clearance level from their profile, defaulting to 'low' if not found.

### `auth.clearance_level_value(level TEXT) → INTEGER`

Converts a clearance level string to a numeric value for comparison:
- admin: 1000
- high: 100
- medium: 50
- low: 10
- null: 0

### `auth.has_clearance(required_level TEXT) → BOOLEAN`

Checks if the current user has the required clearance level.

### `auth.test_rls_access(entity_type TEXT, entity_id UUID, test_user_id UUID DEFAULT NULL) → TABLE`

Tests RLS access for a given entity and user. Returns a table with access_type, has_access, and reason columns.

## Consolidated Policies

The migration consolidates multiple policies into two main types per table:

1. `*_access_policy`: Controls SELECT access
2. `*_modify_policy`: Controls ALL other operations (INSERT, UPDATE, DELETE)

This reduces duplication and makes policies more maintainable.

## Implemented Tables

The enhanced RLS policies have been applied to:

- `molecules`
- `mixtures`
- `molecular_properties`
- `mixture_components`

## Verification and Testing

Two scripts are provided for verification and testing:

1. `enhance_rls_policies.py`: Applies the migration and tests the new functionality
2. `verify_rls_policies.py`: Verifies that RLS policies are working as expected

### Running Verification

```bash
# Apply the migration and test
python enhance_rls_policies.py

# Just verify existing policies
python verify_rls_policies.py
```

### Verification Report

The verification script generates a JSON report with:
- Function verification results
- Policy verification results
- Molecule access tests
- Mixture access tests
- Success rates and overall status

## Best Practices

When working with RLS policies:

1. Always use helper functions for complex access checks
2. Test access with `auth.test_rls_access` before deployment
3. Keep policies simple and focused on a single access concern
4. Use service roles only when absolutely necessary
5. Implement proper team-based access control
6. Document all policies and access patterns

## Next Steps

Future enhancements to consider:

1. Extend RLS to all remaining tables
2. Further consolidate service_role policies
3. Implement more granular team-based access control
4. Add materialized views for commonly accessed data with RLS applied
5. Create admin tools for auditing and managing access

## Troubleshooting

Common issues and solutions:

### Access Denied Unexpectedly

Check:
- User roles and clearance levels
- Team memberships
- Entity ownership and public status
- Run `auth.test_rls_access` to diagnose the issue

### Slow Queries

Possible solutions:
- Use helper functions instead of duplicating complex logic
- Ensure proper indexing on columns used in RLS policies
- Consider materialized views for commonly accessed data

### Policy Conflicts

If policies seem to conflict:
- Check policy order (first matching policy applies)
- Verify that policies don't have overlapping conditions
- Use the test functions to diagnose the issue
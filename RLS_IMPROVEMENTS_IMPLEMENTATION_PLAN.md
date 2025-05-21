# RLS Policy Improvements Implementation Plan

This document details the implementation plan for improving the Row Level Security (RLS) policies in the CryoProtect application based on our previous verification findings. The plan addresses the key issues discovered during verification and provides an optimized approach to database security.

## Background

Our previous verification of RLS policies revealed:

1. The actual database uses a team-based access control model
2. Projects belong to teams, and users are assigned to teams
3. Molecules are linked to projects via the experiments table
4. Some security definer functions were missing from the database
5. RLS policies needed to be optimized for the team-based access model

## Implementation Summary

The implementation consists of five main components:

1. **Security Definer Functions**: Create team-based security definer functions to optimize RLS policies
2. **Table Policies**: Update RLS policies on all tables to use the security definer functions
3. **Performance Indexes**: Add indexes to support RLS policy performance
4. **Service Role Access**: Set up properly scoped service role access
5. **Verification Testing**: Create comprehensive test cases for the improved policies

## File Structure

The implementation is divided across the following files:

1. `migrations/020_rls_policy_improvements.sql`: The main migration file with all SQL changes
2. `apply_rls_policy_improvements.py`: Python script to apply the migration
3. `apply_rls_policy_improvements.sh`: Shell script wrapper for easier execution
4. `tests/test_rls_policy_improvements.py`: Test suite to verify the improvements
5. `run_rls_policy_improvements_test.sh`: Script to run the verification tests

## Implementation Details

### Security Definer Functions

The following security definer functions are implemented:

- `auth.is_team_member(team_id)`: Checks if the current user is a member of the given team
- `auth.is_team_admin(team_id)`: Checks if the current user is an admin of the given team
- `auth.get_user_teams()`: Returns all team IDs the current user belongs to
- `auth.can_access_project(project_id)`: Checks if a user can access a project (via team membership)
- `auth.can_manage_project(project_id)`: Checks if a user can manage a project (via team admin role)
- `auth.get_accessible_projects()`: Returns all project IDs the current user can access
- `auth.can_access_molecule(molecule_id)`: Checks if a user can access a molecule via projects and experiments

### Table Policies

RLS policies are updated for the following tables:

- `projects`: Policies based on team membership
- `experiments`: Policies based on project access
- `molecules`: Policies based on experiment association
- `molecular_properties`: Policies based on molecule access
- `mixtures`: Policies based on project access
- `user_teams`: Policies allowing users to see their own teams and admins to manage team membership
- `user_profile`: Policies allowing users to see/update their own profiles
- `experiment_properties`: Policies based on experiment access

### Performance Indexes

The implementation adds indexes to support RLS policy performance:

- `idx_user_teams_user_id` and `idx_user_teams_team_id`: For faster team membership lookup
- `idx_projects_team_id`: For project team relationship lookup
- `idx_experiments_project_id` and `idx_experiments_molecule_id`: For experiment relationships
- `idx_molecular_properties_molecule_id`: For molecular properties relationships
- `idx_mixtures_project_id`: For mixture project relationship

### Service Role Access

The implementation properly configures service role access for all tables, using the `auth.is_service_role()` function to allow administrative access to the service role.

## Testing Plan

The improvements are verified using:

1. **Security Definer Function Tests**: Verifies each security definer function works correctly
2. **Table Policy Tests**: Verifies RLS policies correctly enforce team-based access
3. **Performance Tests**: Ensures queries are efficient with the new security model

## How to Apply and Verify

1. Run the implementation script:
   ```bash
   ./apply_rls_policy_improvements.sh
   ```

2. Run the verification tests:
   ```bash
   ./run_rls_policy_improvements_test.sh --report
   ```

3. Check the verification report in the reports directory.

## Rollback Plan

If issues are discovered, the changes can be rolled back by:

1. Restoring the previous RLS policies from migration 006_rls_policies.sql
2. Removing the security definer functions
3. Re-applying the original indexes

A rollback script can be created if needed, but is not included in this implementation.

## Future Improvements

For future phases, consider:

1. Additional role-based access controls beyond the team admin/member model
2. More granular permissions at the project level
3. Caching mechanisms for frequently accessed RLS checks
4. Additional performance optimizations for complex queries
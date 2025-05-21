# RLS Verification Test Adaptation Report

## Overview

This document outlines the modifications made to adapt the RLS verification tests to the actual Supabase database structure. The tests were originally written with assumptions about table names and column structures that didn't match the actual database schema.

## Database Schema Differences

### User-Project Relationship Model

The key difference discovered was the relationship model between users and projects:

**Expected Model (in original tests):**
- User profiles directly linked to projects with `user_profile.user_id` and `user_profile.project_id`
- Direct user-project relationship

**Actual Model (in Supabase):**
- Users are linked to teams via `team_members` table (`team_members.user_id` → user identity)
- Projects are linked to teams via `projects.team_id`
- User profiles use `auth_user_id` instead of `user_id` column
- User-project relationships are established indirectly through team membership

### Table Name Pluralization

All table names in the database use plural form (e.g., "molecules" instead of "molecule"), which required adapting the test queries.

## Key Modifications

### 1. RLSTestHelper Changes

The `RLSTestHelper` class was modified in the following ways:

- Changed references from `user_id` to `auth_user_id` in `user_profile` table
- Removed direct project assignment to user profiles
- Added creation and management of a team entity
- Linked users to teams via `team_members` table
- Linked projects to teams via `team_id` field
- Updated test data structure to include team information
- Modified cleanup procedures to handle the new relationship model

### 2. Test File Updates

The test file `test_rls_policies_verification.py` was updated:

- Added team ID to test class attributes
- Updated test methods to use the new relationship model
- Replaced direct user-project relationship tests with team-based relationship tests

### 3. Execution Scripts

Created two execution scripts:

- `run_rls_verification.sh`: Bash script for running tests and generating simple reports
- `run_rls_verification.py`: Python script for running tests with detailed result reporting

Both scripts handle:
- Running all RLS verification tests
- Generating detailed reports with results
- Logging test execution for debugging

## Test Coverage Areas

The verification tests cover the following key areas:

1. **Security Definer Functions**
   - `is_project_member`
   - `user_projects`
   - `is_team_member`
   - `is_project_owner`
   - `molecule_in_user_project`
   - `mixture_in_user_project`

2. **Table Access Policies**
   - SELECT access controls for project members/non-members
   - INSERT access controls for project members/non-members
   - UPDATE access controls for project members/non-members
   - DELETE access controls for project members/non-members

3. **Relationship Policies**
   - Access to related molecular properties based on project membership
   - Access to mixture components based on project membership
   - Access to experiment properties based on project membership

4. **Service Role Access**
   - Full access for service role bypassing RLS
   - Ability to perform administrative operations

## Next Steps

1. Run the verification tests to identify any remaining issues
2. Adjust RLS policies if needed to match the actual relationship model
3. Document the database authorization model based on the test findings
4. Implement any additional RLS policies required for security
5. Optimize RLS policies for performance
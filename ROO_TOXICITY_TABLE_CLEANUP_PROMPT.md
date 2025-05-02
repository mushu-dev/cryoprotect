# Roo PM Toxicity Data Cleanup and Implementation Prompt

## Task Overview

This task focuses on cleaning up the toxicity data source implementation in the CryoProtect v2 project. We're now only using a single toxicity data source, making the current table structure unnecessarily complex. This should be addressed as part of our ongoing database optimization work while maintaining focus on our primary lab verification workflow implementation.

## Current Status

1. We've completed the planning for lab verification workflow implementation as detailed in `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ROO_LAB_VERIFICATION_COMPLETION_PROMPT.md`

2. The toxicity data implementation currently has an unnecessary table structure since we've removed one of the data sources and are now only using a single source.

## Database Cleanup Tasks

### 1. Remove Redundant Table

The toxicity data source table is now redundant and should be removed from the database schema.

1. Create a new migration script `migrations/014_simplify_toxicity_schema.sql` that:
   - Removes the toxicity data source table
   - Updates any foreign key references in other tables
   - Updates RLS policies if necessary

### 2. Update Model Code

Identify and update any code that references the toxicity data source table:

1. Models in `api/models.py`
2. API resources in related files
3. Any schema definitions that reference this table

### 3. Create Simplified Structure (As Needed)

If needed for future expansion, create a simpler structure that better fits our current single-source approach.

## Task Delegation

Please delegate these tasks in coordination with the lab verification implementation:

1. **Database Engineer**: 
   - Create migration script for removing the toxicity data source table
   - Test migration on development database

2. **Backend Engineer**: 
   - Update model code to remove references to the deleted table
   - Ensure toxicity data functionality still works with the simplified structure

3. **QA Engineer**:
   - Verify that toxicity data functionality works correctly after changes
   - Ensure no regressions in related features

## Implementation Guidelines

- Create a backup before making database schema changes
- Run full test suite after changes to ensure nothing breaks
- Consider this a parallel task to the lab verification implementation, with lab verification remaining the higher priority
- Update documentation to reflect the simplified toxicity data structure

## Success Criteria

1. Toxicity data source table is removed from the database schema
2. All code references to the table are updated
3. Toxicity data functionality continues to work correctly
4. Test suite passes with no regressions

## Coordination with Other Tasks

This task should be coordinated with the lab verification implementation as detailed in `ROO_LAB_VERIFICATION_COMPLETION_PROMPT.md`. The lab verification workflow remains our highest priority, but this cleanup can proceed in parallel if resources allow.

## Next Steps

Once this cleanup is complete and the lab verification workflow is implemented:

1. Continue with API Architecture Standardization
2. Complete Testing Framework
3. Address any other critical bugs identified

## References

- Lab Verification Plan: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ROO_LAB_VERIFICATION_COMPLETION_PROMPT.md`
- Toxicity models: `api/models.py` (lines 1893-2043)
- Database schema: `migrations/012_toxicity_schema.sql`
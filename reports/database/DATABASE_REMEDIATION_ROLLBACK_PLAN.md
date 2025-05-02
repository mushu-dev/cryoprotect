# CryoProtect v2 - Database Remediation Rollback Plan

This document outlines the plan for rolling back the database remediation if critical issues are encountered during deployment.

## Overview

Despite thorough testing, issues may arise during the database remediation process that require a rollback. This plan provides a structured approach to restore the database to its pre-remediation state with minimal disruption.

## Prerequisites

Before beginning the remediation process, ensure:

1. Multiple full database backups have been created and verified
2. Application code that depends on the pre-remediation database structure is preserved
3. All team members are familiar with this rollback plan
4. Rollback decision criteria and authority are clearly established

## Rollback Decision Criteria

Consider a rollback if any of the following occur:

1. Critical application functionality is broken and cannot be quickly fixed
2. Data integrity issues are discovered that cannot be resolved
3. Performance degradation is severe (>50% slower response times)
4. Security vulnerabilities are introduced
5. The remediation process has been running for >2x the estimated time

## Rollback Authority

The following roles have authority to initiate a rollback:

1. Database Administrator
2. CTO/CIO
3. Project Manager (in consultation with technical leads)

## Rollback Procedures

### Phase-by-Phase Rollback

If issues are detected during a specific phase, you can roll back that phase using the built-in rollback functionality:

1. For Phase 1 (Security Lockdown):
   ```
   python implement_security.py --rollback
   ```

2. For Phase 2 (Schema Standardization):
   ```
   python standardize_schema.py --rollback
   ```

3. For Phase 3 (Relationship Remediation):
   ```
   python fix_relationships.py --rollback
   ```

4. For Phase 4 (Data Migration):
   No direct rollback script; restore from backup if needed

5. For Phase 5 (Performance Optimization):
   ```
   python verify_performance_indexes.py --remove
   ```

### Full Database Restore

If multiple phases have been applied or the phase-specific rollback fails:

#### Development/Staging Environment

1. Stop all applications connecting to the database
2. Drop the remediated database
3. Restore the database from the pre-remediation backup
4. Verify the restore was successful
5. Restart applications

#### Production Environment

1. Activate the maintenance page for the application
2. Stop all applications connecting to the database
3. Restore the database from the most recent pre-remediation backup
4. Verify the restore was successful using the verification queries below
5. Restart applications in the following order:
   - Backend services
   - API servers
   - Web servers
6. Deactivate the maintenance page
7. Monitor the system closely for 1 hour

## Verification Queries

After restoring the database, run these queries to verify the rollback was successful:

```sql
-- Check if tables are in their original state (singular names)
SELECT tablename FROM pg_tables WHERE schemaname = 'public';

-- Check if RLS is disabled on tables
SELECT relname, relrowsecurity 
FROM pg_class c 
JOIN pg_namespace n ON n.oid = c.relnamespace 
WHERE n.nspname = 'public' AND c.relkind = 'r';

-- Check if junction tables don't exist
SELECT EXISTS (
    SELECT FROM information_schema.tables 
    WHERE table_schema = 'public' 
    AND table_name = 'molecule_experiments'
);

-- Check if performance indexes don't exist
SELECT indexname, tablename 
FROM pg_indexes 
WHERE schemaname = 'public' 
AND indexname IN ('idx_predictions_mixture_property', 'idx_molecule_name_trgm', 'idx_mixture_name_trgm');
```

## Data Recovery Procedures

If data was modified during the remediation process:

1. Identify the affected tables
2. Extract the data from the backup
3. Compare with the current data
4. Merge the data as needed
5. Verify data integrity

## Communication Plan

### Internal Communication

1. Notify all team members of the rollback decision
2. Provide regular updates during the rollback process
3. Conduct a post-rollback meeting to discuss issues and next steps

### External Communication (if applicable)

1. Notify users of the system outage
2. Provide an estimated time for system restoration
3. Communicate when the system is back online
4. Explain the reason for the rollback (at an appropriate level of detail)

## Timeline

| Environment | Estimated Rollback Duration |
|-------------|----------------------------|
| Development | 1-2 hours |
| Staging | 2-3 hours |
| Production | 3-4 hours |

## Post-Rollback Analysis

After completing the rollback:

1. Document the issues that led to the rollback
2. Analyze the root causes
3. Update the remediation plan to address the issues
4. Adjust testing procedures to catch similar issues
5. Schedule a new remediation attempt with the improved plan

## Rollback Test

Before deploying the remediation to production, conduct a rollback test in the staging environment:

1. Apply the remediation
2. Intentionally roll back using this plan
3. Verify the rollback was successful
4. Document any issues encountered

## Contact Information

For assistance with the rollback process, contact:

- Database Administrator: [Name] - [Contact Info]
- Application Support: [Name] - [Contact Info]
- Project Manager: [Name] - [Contact Info]
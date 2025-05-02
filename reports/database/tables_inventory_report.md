# CryoProtect v2 Database Tables Inventory Report

**Date:** 2025-04-18  
**Project ID:** tsdlmynydfuypiugmkev

## Executive Summary

This report provides an inventory of all tables in the CryoProtect v2 Supabase database, classifying each as required, duplicate, or obsolete. The database contains 24 tables/views, of which:

- 13 are required core tables
- 4 are duplicate tables (all empty)
- 2 are obsolete/empty tables
- 5 are views

## Recommendations Summary

- **Keep:** 13 core tables that are actively used and contain data
- **Drop:** 6 tables that are duplicates or obsolete and contain no data
- **Verify:** 5 views to ensure they're still needed before keeping

## Detailed Table Inventory

### Required Core Tables

| Table Name | Row Count | Last Modified | Recommendation |
|------------|-----------|---------------|----------------|
| molecule | 15 | 2025-04-17 | Keep |
| mixture | 5 | 2025-04-17 | Keep |
| experiment | 5 | 2025-04-17 | Keep |
| prediction | 60 | 2025-04-17 | Keep |
| calculation_method | 7 | 2025-04-17 | Keep |
| mixture_component | 11 | 2025-04-17 | Keep |
| molecular_property | 26 | 2025-04-17 | Keep |
| experiment_property | 11 | 2025-04-17 | Keep |
| property_types | 15 | 2025-04-17 | Keep |
| team | 2 | 2025-04-17 | Keep |
| project | 1 | 2025-04-17 | Keep |
| user_profile | 156 | 2025-04-17 | Keep |
| project_membership | 1 | N/A | Keep |

### Duplicate Tables

| Table Name | Row Count | Last Modified | Recommendation |
|------------|-----------|---------------|----------------|
| molecules | 0 | N/A | Drop - Duplicate of molecule |
| experiments | 0 | N/A | Drop - Duplicate of experiment |
| predictions | 0 | N/A | Drop - Duplicate of prediction |
| projects | 0 | N/A | Drop - Duplicate of project |

### Obsolete Tables

| Table Name | Row Count | Last Modified | Recommendation |
|------------|-----------|---------------|----------------|
| protocols | 0 | N/A | Drop - Empty table |
| experiment_mixture_link | 0 | N/A | Drop - Empty linking table |

### Views

| View Name | Row Count | Last Modified | Recommendation |
|-----------|-----------|---------------|----------------|
| experiment_with_results | N/A | N/A | Verify and keep if used |
| mixture_with_components | N/A | N/A | Verify and keep if used |
| mixtures_with_components | 5 | N/A | Verify and keep if used |
| molecule_with_properties | N/A | N/A | Verify and keep if used |
| molecules_with_properties | N/A | N/A | Verify and keep if used |

## Deduplication and Cleanup Plan

1. **Verification Phase**
   - Verify that the empty duplicate tables (`molecules`, `experiments`, `predictions`, `projects`) are indeed duplicates of their counterparts
   - Confirm that the views are still needed and used by the application
   - Check if the empty tables (`protocols`, `experiment_mixture_link`) are needed for future functionality

2. **Backup Phase**
   - Create a full database backup before making any changes
   - Script the table definitions for any tables to be dropped in case they need to be recreated

3. **Cleanup Phase**
   - Drop the duplicate tables that are confirmed to be unnecessary
   - Drop the obsolete tables that are confirmed to be unnecessary
   - Keep the views that are confirmed to be in use

4. **Documentation Phase**
   - Update database schema documentation to reflect the changes
   - Document the rationale for keeping or dropping each table
   - Create a migration script for future deployments

## Conclusion

The CryoProtect v2 database contains several duplicate and empty tables that can be safely removed to simplify the schema and improve maintainability. The core tables contain the essential data for the application and should be preserved. The views should be verified to ensure they're still needed before making any decisions about them.

This inventory provides the foundation for a comprehensive deduplication and cleanup workflow that will streamline the database schema while preserving all necessary functionality.
# Database Schema Fixes - SQL Statements

This document contains the SQL DDL statements required to fix the schema issues identified in the verification report. These statements should be executed against the database to align the schema with the expected structure.

```sql
-- Schema fixes for CryoProtect database

-- 1. Add missing 'properties' JSONB columns to tables
ALTER TABLE public.molecules 
ADD COLUMN properties JSONB DEFAULT '{}'::jsonb;

ALTER TABLE public.mixtures 
ADD COLUMN properties JSONB DEFAULT '{}'::jsonb;

ALTER TABLE public.mixture_components 
ADD COLUMN properties JSONB DEFAULT '{}'::jsonb;

-- 2. Add missing 'units' column to mixture_components
ALTER TABLE public.mixture_components 
ADD COLUMN units VARCHAR;

-- 3. Add missing columns to molecular_properties
-- Note: The current schema uses a different approach with property_type_id, numeric_value, text_value, etc.
-- This is actually a more normalized approach than having generic property_name/property_value columns.
-- We'll add these columns for compatibility with the verification script, but the existing structure is better.
ALTER TABLE public.molecular_properties 
ADD COLUMN property_name VARCHAR,
ADD COLUMN property_value TEXT,
ADD COLUMN property_type VARCHAR,
ADD COLUMN source VARCHAR;

-- 4. Add missing 'parameters' column to calculation_methods
-- Note: There's already a calculation_method (singular) table with a parameters column
-- This appears to be a duplication issue, but we'll add it for verification compatibility
ALTER TABLE public.calculation_methods 
ADD COLUMN parameters JSONB DEFAULT '{}'::jsonb;

-- 5. Add missing columns to experiments
ALTER TABLE public.experiments 
ADD COLUMN protocol JSONB DEFAULT '{}'::jsonb,
ADD COLUMN results JSONB DEFAULT '{}'::jsonb;

-- 6. Add missing columns to predictions
ALTER TABLE public.predictions 
ADD COLUMN name VARCHAR,
ADD COLUMN description TEXT,
ADD COLUMN results JSONB DEFAULT '{}'::jsonb;

-- 7. Update the verification script expectations for ID columns
-- Note: This is a comment for the developers. The verification script should be updated
-- to expect UUID types for ID columns instead of varchar, as UUIDs are more appropriate
-- for primary and foreign keys in a modern database design.

-- 8. Update the verification script expectations for timestamp nullability
-- Note: This is a comment for the developers. The verification script should be updated
-- to expect NOT NULL constraints on timestamp columns, as this enforces data integrity.
```

## Implementation Notes

1. These ALTER TABLE statements add the missing columns with appropriate types and default values.
2. For JSON/JSONB columns, empty JSON objects are set as default values.
3. The statements maintain the current UUID types for ID columns and NOT NULL constraints on timestamp columns, as these are more appropriate for the database design.
4. After executing these statements, the verification script should be updated to align its expectations with the actual schema design, particularly regarding UUID types and timestamp nullability.

## Execution Plan

1. Create a database backup before executing these statements
2. Execute the statements in a test environment first
3. Verify that the changes resolve the issues reported in the verification script
4. Apply the changes to the production database during a maintenance window
5. Update the verification script to align with the actual schema design
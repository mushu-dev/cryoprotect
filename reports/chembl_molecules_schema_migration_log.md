# ChEMBL Molecules Schema Migration Log

**Date:** 2025-04-26  
**Task:** task-imp-chembl-schema-002  
**Author:** Apex Implementer  
**Migration File:** [migrations/chembl_molecules_schema_remediation.sql](../migrations/chembl_molecules_schema_remediation.sql)

## 1. Migration Summary

This migration addresses the schema drift identified in the [ChEMBL Molecules Schema Inspection Report](./chembl_molecules_schema_inspection.md) by aligning the `public.molecules` table with the canonical schema while preserving existing data and ensuring compatibility with code that references `pubchem_cid`.

### Key Changes

1. **Added `pubchem_cid` Column**
   - Added as INTEGER with UNIQUE constraint
   - Allows NULL initially to support gradual data population
   - Primary column for PubChem Compound ID storage

2. **Added `cid` as Generated Column**
   - Added as a generated column based on `pubchem_cid`
   - Ensures compatibility with canonical schema
   - Allows code to reference either column name

3. **Added `pubchem_link` Generated Column**
   - Generated from `pubchem_cid` to provide direct links to PubChem
   - Format: `https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}`

4. **Added `molecular_formula` Column**
   - Added as TEXT column
   - Copies data from existing `formula` column if present
   - Maintains canonical schema compatibility

5. **Created Indexes**
   - Created index on `cid` column: `idx_molecules_cid`
   - Created index on `pubchem_cid` column: `idx_molecules_pubchem_cid`
   - Improves query performance for lookups by PubChem ID

6. **Added Table Comment**
   - Added documentation explaining the dual column approach

## 2. Implementation Approach

### 2.1 Idempotent Design

The migration is designed to be idempotent and can be safely run multiple times:
- Each operation is wrapped in a conditional block that checks if the column/index already exists
- Uses PostgreSQL's `DO` blocks for conditional execution
- Provides informative notices about operations performed or skipped

### 2.2 Data Preservation

The migration preserves all existing data:
- No columns are dropped
- No data is modified except for copying `formula` to `molecular_formula` where applicable
- All operations are additive

### 2.3 Generated Columns

The migration uses PostgreSQL's GENERATED ALWAYS AS STORED feature:
- `cid` is generated from `pubchem_cid`
- `pubchem_link` is generated from `pubchem_cid`
- This approach ensures the columns are always in sync

## 3. Rationale for Design Decisions

### 3.1 `pubchem_cid` as Primary Column

The inspection report revealed that code references `pubchem_cid` rather than `cid`. We chose to:
- Add `pubchem_cid` as the primary storage column
- Generate `cid` from `pubchem_cid` rather than vice versa

This approach minimizes code changes while maintaining schema compatibility.

### 3.2 Allowing NULL Initially

The `pubchem_cid` column allows NULL values initially, which:
- Supports gradual data population
- Prevents issues with existing records
- Allows for a phased approach to data migration

A future migration can add a NOT NULL constraint once all records have been populated.

### 3.3 Dual Indexes

Indexes on both `cid` and `pubchem_cid` ensure:
- Optimal query performance regardless of which column is referenced
- Support for the UNIQUE constraint on `pubchem_cid`
- Compatibility with existing queries

## 4. Verification Steps

After applying this migration, the following verification steps should be performed:

1. **Schema Verification**
   ```sql
   SELECT 
       column_name, 
       data_type, 
       is_nullable,
       column_default
   FROM 
       information_schema.columns 
   WHERE 
       table_schema = 'public' 
       AND table_name = 'molecules'
   ORDER BY 
       ordinal_position;
   ```

2. **Index Verification**
   ```sql
   SELECT 
       indexname, 
       indexdef 
   FROM 
       pg_indexes 
   WHERE 
       schemaname = 'public' 
       AND tablename = 'molecules';
   ```

3. **Generated Column Verification**
   ```sql
   -- Insert a test record with pubchem_cid
   INSERT INTO public.molecules (name, pubchem_cid) 
   VALUES ('Test Molecule', 12345) 
   RETURNING id, pubchem_cid, cid, pubchem_link;
   
   -- Clean up test record
   DELETE FROM public.molecules WHERE name = 'Test Molecule';
   ```

## 5. Next Steps

1. **Data Population**
   - Populate `pubchem_cid` for existing records
   - Consider a separate migration to make `pubchem_cid` NOT NULL once populated

2. **Code Audit**
   - Audit code for references to `cid` and `pubchem_cid`
   - Document the relationship between the columns

3. **Documentation Update**
   - Update schema documentation to reflect the dual column approach
   - Add code comments explaining the generated column relationship
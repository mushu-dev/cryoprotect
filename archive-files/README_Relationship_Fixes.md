# CryoProtect v2 - Database Relationship Fixes

This document provides information about the `fix_relationships.py` script, which addresses relationship design issues in the CryoProtect Supabase database.

## Overview

The script fixes several relationship design issues:

1. **Fan Trap Replacement**: Replaces fan traps with proper junction tables
   - Creates a `molecule_proteins` junction table for many-to-many relationships
   - Creates a `molecule_experiments` junction table to directly associate molecules with experiments

2. **3NF Normalization**: Applies Third Normal Form principles
   - Ensures all tables have proper primary keys
   - Removes transitive dependencies
   - Ensures all non-key attributes are dependent on the primary key
   - Converts string-based property types to foreign key references

3. **Foreign Key Constraints**: Adds missing constraints with appropriate indexes
   - Ensures all relationships between tables are properly defined with REFERENCES constraints
   - Adds indexes on all foreign key columns for performance

4. **Data Migration**: Safely migrates data to the new structure
   - Creates new junction tables
   - Migrates data from existing tables to the new structure
   - Adds all necessary constraints and indexes

5. **Verification & Rollback**: Includes safety mechanisms
   - Verifies the new relationships after changes
   - Provides rollback capability in case of failure

## Prerequisites

- Python 3.6+
- Supabase Python client (`pip install supabase`)
- Access to the CryoProtect Supabase project (ID: tsdlmynydfuypiugmkev)
- Environment variables set in a `.env` file:
  ```
  SUPABASE_URL=https://tsdlmynydfuypiugmkev.supabase.co
  SUPABASE_KEY=your-supabase-key
  ```

## Usage

The script can be run with various command-line options:

```bash
# Show what would be done without making changes
python fix_relationships.py --dry-run

# Only verify the relationships without making changes
python fix_relationships.py --verify-only

# Apply all fixes
python fix_relationships.py

# Rollback changes using the most recent backup
python fix_relationships.py --rollback
```

## Detailed Fixes

### 1. Junction Tables

The script creates two new junction tables:

- **molecule_proteins**: Links molecules to proteins with binding information
  ```sql
  CREATE TABLE IF NOT EXISTS public.molecule_proteins (
      id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
      molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
      protein_id UUID NOT NULL REFERENCES public.proteins(id) ON DELETE CASCADE,
      binding_affinity NUMERIC,
      interaction_type VARCHAR(100),
      created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
      updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
      created_by UUID REFERENCES auth.users(id),
      UNIQUE(molecule_id, protein_id)
  );
  ```

- **molecule_experiments**: Links molecules directly to experiments
  ```sql
  CREATE TABLE IF NOT EXISTS public.molecule_experiments (
      id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
      molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
      experiment_id UUID NOT NULL REFERENCES public.experiments(id) ON DELETE CASCADE,
      role VARCHAR(100),
      concentration NUMERIC,
      concentration_unit VARCHAR(50),
      created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
      updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
      created_by UUID REFERENCES auth.users(id),
      UNIQUE(molecule_id, experiment_id)
  );
  ```

### 2. Predictions Table Fix

Modifies the predictions table to properly handle both molecule and mixture relationships:

```sql
ALTER TABLE public.predictions 
DROP CONSTRAINT IF EXISTS predictions_mixture_id_fkey,
ALTER COLUMN mixture_id DROP NOT NULL,
ADD COLUMN IF NOT EXISTS molecule_id UUID REFERENCES public.molecules(id) ON DELETE CASCADE,
DROP CONSTRAINT IF EXISTS check_prediction_target,
ADD CONSTRAINT check_prediction_target CHECK (
    (molecule_id IS NOT NULL AND mixture_id IS NULL) OR 
    (molecule_id IS NULL AND mixture_id IS NOT NULL)
);
```

### 3. Molecular Properties Fix

Converts the string-based property_type to a foreign key reference:

```sql
-- First, migrate existing property types to the property_types table
INSERT INTO public.property_types (name, data_type)
SELECT DISTINCT property_type, 'numeric'
FROM public.molecular_properties
WHERE property_type NOT IN (SELECT name FROM public.property_types)
ON CONFLICT (name) DO NOTHING;

-- Add property_type_id column
ALTER TABLE public.molecular_properties
ADD COLUMN IF NOT EXISTS property_type_id UUID;

-- Update property_type_id based on property_type
UPDATE public.molecular_properties mp
SET property_type_id = pt.id
FROM public.property_types pt
WHERE mp.property_type = pt.name
AND mp.property_type_id IS NULL;
```

## Verification

The script includes a verification step that checks:

1. The existence of the new junction tables
2. The structure of the modified tables
3. The presence of required constraints
4. The integrity of the relationships

You can run verification separately:

```bash
python fix_relationships.py --verify-only
```

## Rollback Mechanism

Before making any changes, the script creates a backup of the current database state. If something goes wrong, you can rollback:

```bash
python fix_relationships.py --rollback
```

This will:
1. Drop the newly created tables
2. Restore the original tables to their previous state

## Testing

A test script is provided to verify the fix_relationships.py script:

```bash
# Run basic tests (dry-run and verify-only)
python test_fix_relationships.py

# Run a full test including making actual changes (use with caution)
python test_fix_relationships.py --full-test
```

## After Running the Script

After successfully running the script, your database will have:

1. Proper junction tables for many-to-many relationships
2. Tables that follow 3NF normalization principles
3. Appropriate foreign key constraints with indexes
4. Improved data integrity and query performance

## Troubleshooting

If you encounter issues:

1. Check the log file (`fix_relationships.log`) for detailed error messages
2. Run with `--verify-only` to check the current state
3. Use `--rollback` to revert changes if needed
4. Ensure your Supabase credentials are correct in the `.env` file

## Contact

For assistance with this script, contact the CryoProtect development team.
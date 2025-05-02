-- ChEMBL Molecules Schema Remediation
-- Date: 2025-04-26
-- Task: task-imp-chembl-schema-002
-- Author: Apex Implementer

-- This migration is idempotent and can be safely run multiple times
-- It aligns the public.molecules table with the canonical schema while preserving existing data

-- Begin transaction
BEGIN;

-- 1. Add pubchem_cid column if it doesn't exist
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_schema = 'public' 
        AND table_name = 'molecules' 
        AND column_name = 'pubchem_cid'
    ) THEN
        ALTER TABLE public.molecules ADD COLUMN pubchem_cid INTEGER;
        -- Add unique constraint but allow NULL initially to support gradual data population
        ALTER TABLE public.molecules ADD CONSTRAINT molecules_pubchem_cid_key UNIQUE (pubchem_cid);
        RAISE NOTICE 'Added pubchem_cid column to public.molecules';
    ELSE
        RAISE NOTICE 'pubchem_cid column already exists in public.molecules';
    END IF;
END
$$;

-- 2. Add cid as a generated column based on pubchem_cid if it doesn't exist
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_schema = 'public' 
        AND table_name = 'molecules' 
        AND column_name = 'cid'
    ) THEN
        -- Add cid as a generated column based on pubchem_cid
        ALTER TABLE public.molecules ADD COLUMN cid INTEGER GENERATED ALWAYS AS (pubchem_cid) STORED;
        RAISE NOTICE 'Added cid column as a generated column based on pubchem_cid';
    ELSE
        RAISE NOTICE 'cid column already exists in public.molecules';
    END IF;
END
$$;

-- 3. Add pubchem_link as a generated column if it doesn't exist
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_schema = 'public' 
        AND table_name = 'molecules' 
        AND column_name = 'pubchem_link'
    ) THEN
        -- Add pubchem_link as a generated column based on pubchem_cid
        ALTER TABLE public.molecules ADD COLUMN pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || pubchem_cid) STORED;
        RAISE NOTICE 'Added pubchem_link column as a generated column based on pubchem_cid';
    ELSE
        RAISE NOTICE 'pubchem_link column already exists in public.molecules';
    END IF;
END
$$;

-- 4. Add molecular_formula column if it doesn't exist (noting there's already a formula column)
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_schema = 'public' 
        AND table_name = 'molecules' 
        AND column_name = 'molecular_formula'
    ) THEN
        -- Add molecular_formula column
        ALTER TABLE public.molecules ADD COLUMN molecular_formula TEXT;
        
        -- If formula column exists, copy values to molecular_formula for consistency
        IF EXISTS (
            SELECT 1 FROM information_schema.columns 
            WHERE table_schema = 'public' 
            AND table_name = 'molecules' 
            AND column_name = 'formula'
        ) THEN
            UPDATE public.molecules SET molecular_formula = formula WHERE molecular_formula IS NULL;
            RAISE NOTICE 'Copied existing formula values to molecular_formula for consistency';
        END IF;
        
        RAISE NOTICE 'Added molecular_formula column to public.molecules';
    ELSE
        RAISE NOTICE 'molecular_formula column already exists in public.molecules';
    END IF;
END
$$;

-- 5. Create indexes on cid and pubchem_cid if they don't exist
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_indexes 
        WHERE schemaname = 'public' 
        AND tablename = 'molecules' 
        AND indexname = 'idx_molecules_cid'
    ) THEN
        CREATE INDEX idx_molecules_cid ON public.molecules(cid);
        RAISE NOTICE 'Created index on cid column';
    ELSE
        RAISE NOTICE 'Index on cid column already exists';
    END IF;
    
    IF NOT EXISTS (
        SELECT 1 FROM pg_indexes 
        WHERE schemaname = 'public' 
        AND tablename = 'molecules' 
        AND indexname = 'idx_molecules_pubchem_cid'
    ) THEN
        CREATE INDEX idx_molecules_pubchem_cid ON public.molecules(pubchem_cid);
        RAISE NOTICE 'Created index on pubchem_cid column';
    ELSE
        RAISE NOTICE 'Index on pubchem_cid column already exists';
    END IF;
END
$$;

-- 6. Add a comment to the table explaining the dual column approach
COMMENT ON TABLE public.molecules IS 'Molecule information with both cid (canonical) and pubchem_cid columns for compatibility. The cid column is generated from pubchem_cid.';

-- Commit transaction
COMMIT;

-- Post-migration verification query (commented out, for manual execution)
/*
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
*/
-- Migration: unified_importer.sql
-- Description: Schema changes to support the unified chemical data importer
--
-- This migration adds tables and fields needed for the unified chemical data importer,
-- including data source tracking, molecule synonyms, import job history, and
-- performance optimizations.

-- This is a wrapper for the main 021_unified_importer.sql migration
-- that should already be in the project's migrations directory.
-- It's provided for convenience when using the unified_importer package directly.

-- Insert code to check if the migration has already been applied
DO $$
BEGIN

-- Check if migration has already been applied
IF NOT EXISTS (
    SELECT 1 FROM pg_tables 
    WHERE tablename = 'data_sources'
) THEN
    RAISE NOTICE 'The unified importer migration needs to be applied. Use the main 021_unified_importer.sql script.';
ELSE
    RAISE NOTICE 'The unified importer migration has already been applied.';
END IF;

END $$;
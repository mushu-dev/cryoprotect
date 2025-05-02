-- CryoProtect v2 Performance Optimization Indexes
-- Migration: 010_performance_indexes.sql
-- This migration adds indexes recommended by performance testing

-- Enable pg_trgm extension for text search
CREATE EXTENSION IF NOT EXISTS pg_trgm;

-- Add composite index for predictions table
-- This improves performance for queries that filter on both mixture_id and property_type_id
-- Particularly beneficial for the "compare_predictions_experiments" query
CREATE INDEX IF NOT EXISTS idx_predictions_mixture_property 
ON public.predictions(mixture_id, property_type_id);

-- Add text search index for molecule names
-- This improves performance for ILIKE/text search operations on molecule names
CREATE INDEX IF NOT EXISTS idx_molecule_name_trgm 
ON public.molecules USING gin (name gin_trgm_ops);

-- Add text search index for mixture names
-- This improves performance for ILIKE/text search operations on mixture names
CREATE INDEX IF NOT EXISTS idx_mixture_name_trgm 
ON public.mixtures USING gin (name gin_trgm_ops);

-- Add the same indexes to the scientific schema tables if they exist
DO $$
BEGIN
    -- Check if the scientific schema tables exist
    IF EXISTS (SELECT FROM pg_tables WHERE tablename = 'molecule') THEN
        -- Add text search index for molecule names in scientific schema
        IF NOT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_molecule_name_trgm') THEN
            CREATE INDEX idx_molecule_name_trgm ON public.molecule USING gin (name gin_trgm_ops);
        END IF;
    END IF;

    IF EXISTS (SELECT FROM pg_tables WHERE tablename = 'mixture') THEN
        -- Add text search index for mixture names in scientific schema
        IF NOT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_mixture_name_trgm') THEN
            CREATE INDEX idx_mixture_name_trgm ON public.mixture USING gin (name gin_trgm_ops);
        END IF;
    END IF;

    IF EXISTS (SELECT FROM pg_tables WHERE tablename = 'prediction') THEN
        -- Add composite index for prediction table in scientific schema
        IF NOT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_prediction_molecule_property') THEN
            CREATE INDEX idx_prediction_molecule_property ON public.prediction(molecule_id, property_type);
        END IF;
    END IF;
END
$$;

-- Add comment explaining the purpose of these indexes
COMMENT ON INDEX public.idx_predictions_mixture_property IS 'Improves performance for queries that filter on both mixture_id and property_type_id';
COMMENT ON INDEX public.idx_molecule_name_trgm IS 'Improves performance for text search operations on molecule names';
COMMENT ON INDEX public.idx_mixture_name_trgm IS 'Improves performance for text search operations on mixture names';
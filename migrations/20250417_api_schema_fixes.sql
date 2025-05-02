-- CryoProtect API Schema Fixes Migration
-- Date: 2025-04-17

-- Enable uuid-ossp extension if not already enabled
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Fix predictions table
ALTER TABLE public.predictions 
  ALTER COLUMN id TYPE UUID USING uuid_generate_v4(),
  ALTER COLUMN mixture_id TYPE UUID,
  ALTER COLUMN property_type_id TYPE UUID,
  ADD COLUMN IF NOT EXISTS calculation_method_id UUID,
  ADD COLUMN IF NOT EXISTS confidence NUMERIC,
  ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW(),
  ADD COLUMN IF NOT EXISTS created_by UUID;

-- Fix experiments table
ALTER TABLE public.experiments 
  ALTER COLUMN id TYPE UUID USING uuid_generate_v4(),
  ALTER COLUMN mixture_id TYPE UUID,
  ALTER COLUMN property_type_id TYPE UUID,
  ADD COLUMN IF NOT EXISTS experimental_conditions TEXT,
  ADD COLUMN IF NOT EXISTS date_performed DATE,
  ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW(),
  ADD COLUMN IF NOT EXISTS created_by UUID;

-- Add indexes
CREATE INDEX IF NOT EXISTS idx_pred_mixture ON public.predictions(mixture_id);
CREATE INDEX IF NOT EXISTS idx_exp_mixture ON public.experiments(mixture_id);
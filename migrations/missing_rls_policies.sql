-- CryoProtect Analyzer - Missing RLS Policies for Views and Tables
-- This migration adds RLS policies for missing tables and views:
-- 1. experiment_with_results
-- 2. migrations
-- 3. mixture_with_components
-- 4. molecule_with_properties

-- =========================
-- 1. molecule_with_properties view
-- =========================
-- Enable RLS on the view
ALTER VIEW public.molecule_with_properties SECURITY INVOKER;

-- Create RLS policy for molecule_with_properties (view)
-- Note: PostgreSQL RLS for views works by enforcing the policies on the underlying tables
-- The view will honor the RLS policies of the underlying molecule and molecular_property tables

-- =========================
-- 2. mixture_with_components view
-- =========================
-- Enable RLS on the view
ALTER VIEW public.mixture_with_components SECURITY INVOKER;

-- Create RLS policy for mixture_with_components (view)
-- Note: PostgreSQL RLS for views works by enforcing the policies on the underlying tables
-- The view will honor the RLS policies of the underlying mixture and mixture_component tables

-- =========================
-- 3. migrations table
-- =========================
-- Check if the migrations table exists
DO $$
BEGIN
    IF EXISTS (SELECT FROM pg_tables WHERE schemaname = 'public' AND tablename = 'migrations') THEN
        -- Enable RLS on migrations table
        ALTER TABLE public.migrations ENABLE ROW LEVEL SECURITY;

        -- Create RLS policies for migrations
        CREATE POLICY "Allow service role access to migrations" 
          ON public.migrations 
          USING (auth.role() = 'service_role');

        -- Allow authenticated users with admin role to view migrations
        CREATE POLICY "Allow admin users to view migrations" 
          ON public.migrations 
          FOR SELECT
          USING (auth.role() = 'authenticated' AND EXISTS (
              SELECT 1 FROM user_profile
              WHERE user_profile.auth_user_id = auth.uid()
              AND user_profile.role = 'admin'
          ));
    END IF;
END
$$;

-- =========================
-- 4. experiment_with_results view
-- =========================
-- Check if the experiment_with_results view exists, create if it doesn't
CREATE OR REPLACE VIEW public.experiment_with_results AS
SELECT 
    e.id,
    e.name,
    e.description,
    e.mixture_id,
    e.project_id,
    e.preparation_protocol,
    e.temperature,
    e.temperature_unit,
    e.pressure,
    e.pressure_unit,
    e.created_by,
    e.data_source,
    e.version,
    e.modification_history,
    e.created_at,
    e.updated_at,
    jsonb_agg(
        jsonb_build_object(
            'property_id', ep.id,
            'property_type', ep.property_type,
            'value', ep.value,
            'unit', ep.unit,
            'method_id', ep.method_id,
            'method_name', cm.name
        )
    ) AS results
FROM 
    public.experiment e
LEFT JOIN 
    public.experiment_property ep ON e.id = ep.experiment_id
LEFT JOIN 
    public.calculation_method cm ON ep.method_id = cm.id
GROUP BY 
    e.id;

-- Enable RLS on the view
ALTER VIEW public.experiment_with_results SECURITY INVOKER;

-- =========================
-- 5. Service Role Policies for all views
-- =========================
-- Add service role policies for views' underlying tables if they don't already exist

-- Allow service role to access molecule_with_properties
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'molecule' 
        AND policyname = 'Allow service role full access to molecule'
    ) THEN
        CREATE POLICY "Allow service role full access to molecule" 
          ON public.molecule 
          USING (auth.role() = 'service_role');
    END IF;

    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'molecular_property' 
        AND policyname = 'Allow service role full access to molecular_property'
    ) THEN
        CREATE POLICY "Allow service role full access to molecular_property" 
          ON public.molecular_property 
          USING (auth.role() = 'service_role');
    END IF;
END
$$;

-- Allow service role to access mixture_with_components
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'mixture' 
        AND policyname = 'Allow service role full access to mixture'
    ) THEN
        CREATE POLICY "Allow service role full access to mixture" 
          ON public.mixture 
          USING (auth.role() = 'service_role');
    END IF;

    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'mixture_component' 
        AND policyname = 'Allow service role full access to mixture_component'
    ) THEN
        CREATE POLICY "Allow service role full access to mixture_component" 
          ON public.mixture_component 
          USING (auth.role() = 'service_role');
    END IF;
END
$$;

-- Allow service role to access experiment_with_results
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'experiment' 
        AND policyname = 'Allow service role full access to experiment'
    ) THEN
        CREATE POLICY "Allow service role full access to experiment" 
          ON public.experiment 
          USING (auth.role() = 'service_role');
    END IF;

    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'experiment_property' 
        AND policyname = 'Allow service role full access to experiment_property'
    ) THEN
        CREATE POLICY "Allow service role full access to experiment_property" 
          ON public.experiment_property 
          USING (auth.role() = 'service_role');
    END IF;
END
$$;

-- =========================
-- 6. Comment on RLS Implementation
-- =========================
COMMENT ON VIEW public.molecule_with_properties IS 'View for molecules with their properties, honors RLS policies of underlying tables';
COMMENT ON VIEW public.mixture_with_components IS 'View for mixtures with their components, honors RLS policies of underlying tables';
COMMENT ON VIEW public.experiment_with_results IS 'View for experiments with their results, honors RLS policies of underlying tables';
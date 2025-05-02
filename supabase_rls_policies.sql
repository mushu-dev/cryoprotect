-- CryoProtect Analyzer - Missing RLS Policies for Views and Tables (Supabase SQL Editor Version)
-- This script can be directly executed in the Supabase SQL Editor to apply all missing RLS policies
-- and performance optimizations in a single transaction.

-- Start transaction to ensure atomic operation
BEGIN;

-- =========================
-- 0. Performance indexes for RLS
-- =========================
-- Adding indexes specifically for columns used in RLS policies to improve performance

-- Project membership indexes (critically important for RLS policy performance)
CREATE INDEX IF NOT EXISTS idx_user_profile_project_id ON user_profile(project_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_user_id ON user_profile(user_id);
CREATE INDEX IF NOT EXISTS idx_experiment_project_id ON experiment(project_id);
CREATE INDEX IF NOT EXISTS idx_molecule_project_id ON molecule(project_id);
CREATE INDEX IF NOT EXISTS idx_mixture_project_id ON mixture(project_id);

-- Foreign key indexes for relationship traversal in policies
CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON molecular_property(molecule_id);
CREATE INDEX IF NOT EXISTS idx_experiment_property_experiment_id ON experiment_property(experiment_id);
CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON mixture_component(mixture_id);

-- =========================
-- 1. molecule_with_properties view
-- =========================
-- Enable RLS on the view with SECURITY INVOKER to ensure it respects underlying table policies
ALTER VIEW public.molecule_with_properties SECURITY INVOKER;

-- Create an optimized policy for the underlying molecule table if not exists
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
END
$$;

-- =========================
-- 2. mixture_with_components view
-- =========================
-- Enable RLS on the view
ALTER VIEW public.mixture_with_components SECURITY INVOKER;

-- Optimized policies for underlying tables
DO $$
BEGIN
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

-- =========================
-- 3. migrations table
-- =========================
-- Check if the migrations table exists and add RLS
DO $$
BEGIN
    IF EXISTS (SELECT FROM pg_tables WHERE schemaname = 'public' AND tablename = 'migrations') THEN
        -- Enable RLS on migrations table
        ALTER TABLE public.migrations ENABLE ROW LEVEL SECURITY;

        -- Create service role policy with optimized condition
        IF NOT EXISTS (
            SELECT 1 FROM pg_policies 
            WHERE tablename = 'migrations' 
            AND policyname = 'Allow service role access to migrations'
        ) THEN
            CREATE POLICY "Allow service role access to migrations" 
              ON public.migrations 
              USING (auth.role() = 'service_role');
        END IF;

        -- Allow admin users to view migrations (optimized query)
        IF NOT EXISTS (
            SELECT 1 FROM pg_policies 
            WHERE tablename = 'migrations' 
            AND policyname = 'Allow admin users to view migrations'
        ) THEN
            CREATE POLICY "Allow admin users to view migrations" 
              ON public.migrations 
              FOR SELECT
              USING (
                  EXISTS (
                      SELECT 1 FROM user_profile
                      WHERE user_profile.auth_user_id = auth.uid()
                      AND user_profile.role = 'admin'
                  )
              );
        END IF;
    END IF;
END
$$;

-- =========================
-- 4. experiment_with_results view
-- =========================
-- Check if the experiment_with_results view exists, create if it doesn't
-- Using an optimized view definition that minimizes joins
DO $$
BEGIN
    -- Drop the view if it exists to re-create it
    DROP VIEW IF EXISTS public.experiment_with_results;
    
    -- Create the optimized view
    CREATE VIEW public.experiment_with_results AS
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
        COALESCE(
            (SELECT jsonb_agg(
                jsonb_build_object(
                    'property_id', ep.id,
                    'property_type', ep.property_type,
                    'value', ep.value,
                    'unit', ep.unit,
                    'method_id', ep.method_id,
                    'method_name', cm.name
                )
            )
            FROM public.experiment_property ep
            LEFT JOIN public.calculation_method cm ON ep.method_id = cm.id
            WHERE ep.experiment_id = e.id
            GROUP BY ep.experiment_id),
            '[]'::jsonb
        ) AS results
    FROM 
        public.experiment e;
        
    -- Set security to INVOKER
    ALTER VIEW public.experiment_with_results SECURITY INVOKER;
END
$$;

-- Create an optimized policy for the experiment table if it doesn't exist
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
END
$$;

-- =========================
-- 5. Enhanced Audit Trail for Scientific Data
-- =========================
-- Create an audit trail table for tracking access to sensitive scientific data
CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    table_name TEXT NOT NULL,
    record_id UUID NOT NULL,
    operation TEXT NOT NULL,
    user_id UUID NOT NULL,
    timestamp TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    old_data JSONB,
    new_data JSONB,
    ip_address TEXT,
    application_context TEXT
);

CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_table_record 
ON public.scientific_data_audit(table_name, record_id);

CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_user_id 
ON public.scientific_data_audit(user_id);

CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_timestamp 
ON public.scientific_data_audit(timestamp);

-- Create audit trigger function for scientific data
CREATE OR REPLACE FUNCTION public.audit_scientific_data()
RETURNS TRIGGER AS $$
DECLARE
    audit_data JSONB;
    record_jsonb JSONB;
BEGIN
    -- Skip if this is a service role operation during initial data load
    IF auth.role() = 'service_role' AND 
       coalesce(current_setting('app.bypass_audit', true), 'false') = 'true' THEN
        RETURN NEW;
    END IF;

    -- Prepare audit data
    IF TG_OP = 'DELETE' THEN
        record_jsonb = to_jsonb(OLD);
        audit_data = jsonb_build_object(
            'table_name', TG_TABLE_NAME,
            'record_id', OLD.id,
            'operation', TG_OP,
            'user_id', coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'),
            'old_data', record_jsonb,
            'new_data', null,
            'ip_address', coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
            'application_context', coalesce(current_setting('app.context', true), 'unknown')
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, OLD.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             record_jsonb, null, 
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             coalesce(current_setting('app.context', true), 'unknown'));
        RETURN OLD;
    ELSIF TG_OP = 'UPDATE' THEN
        record_jsonb = to_jsonb(NEW);
        audit_data = jsonb_build_object(
            'table_name', TG_TABLE_NAME,
            'record_id', NEW.id,
            'operation', TG_OP,
            'user_id', coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'),
            'old_data', to_jsonb(OLD),
            'new_data', record_jsonb,
            'ip_address', coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
            'application_context', coalesce(current_setting('app.context', true), 'unknown')
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, NEW.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             to_jsonb(OLD), record_jsonb,
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             coalesce(current_setting('app.context', true), 'unknown'));
        RETURN NEW;
    ELSIF TG_OP = 'INSERT' THEN
        record_jsonb = to_jsonb(NEW);
        audit_data = jsonb_build_object(
            'table_name', TG_TABLE_NAME,
            'record_id', NEW.id,
            'operation', TG_OP,
            'user_id', coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'),
            'old_data', null,
            'new_data', record_jsonb,
            'ip_address', coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
            'application_context', coalesce(current_setting('app.context', true), 'unknown')
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, NEW.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             null, record_jsonb,
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             coalesce(current_setting('app.context', true), 'unknown'));
        RETURN NEW;
    END IF;
    RETURN NULL;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Apply audit triggers to scientific data tables (conditionally)
DO $$
BEGIN
    -- Add trigger to molecule table if it doesn't exist
    IF NOT EXISTS (
        SELECT 1 FROM pg_trigger 
        WHERE tgname = 'audit_molecule_trigger' 
        AND tgrelid = 'public.molecule'::regclass
    ) THEN
        CREATE TRIGGER audit_molecule_trigger
        AFTER INSERT OR UPDATE OR DELETE ON public.molecule
        FOR EACH ROW EXECUTE FUNCTION public.audit_scientific_data();
    END IF;

    -- Add trigger to mixture table if it doesn't exist
    IF NOT EXISTS (
        SELECT 1 FROM pg_trigger 
        WHERE tgname = 'audit_mixture_trigger' 
        AND tgrelid = 'public.mixture'::regclass
    ) THEN
        CREATE TRIGGER audit_mixture_trigger
        AFTER INSERT OR UPDATE OR DELETE ON public.mixture
        FOR EACH ROW EXECUTE FUNCTION public.audit_scientific_data();
    END IF;

    -- Add trigger to experiment table if it doesn't exist
    IF NOT EXISTS (
        SELECT 1 FROM pg_trigger 
        WHERE tgname = 'audit_experiment_trigger' 
        AND tgrelid = 'public.experiment'::regclass
    ) THEN
        CREATE TRIGGER audit_experiment_trigger
        AFTER INSERT OR UPDATE OR DELETE ON public.experiment
        FOR EACH ROW EXECUTE FUNCTION public.audit_scientific_data();
    END IF;
END
$$;

-- =========================
-- 6. RLS for the scientific_data_audit table
-- =========================
ALTER TABLE public.scientific_data_audit ENABLE ROW LEVEL SECURITY;

-- Admins can view all audit records
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'scientific_data_audit' 
        AND policyname = 'Admins can view all audit records'
    ) THEN
        CREATE POLICY "Admins can view all audit records"
          ON public.scientific_data_audit
          FOR SELECT
          USING (
            EXISTS (
              SELECT 1 FROM user_profile
              WHERE user_profile.auth_user_id = auth.uid()
              AND user_profile.role = 'admin'
            )
          );
    END IF;
END
$$;

-- Users can view audit records for their own actions
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'scientific_data_audit' 
        AND policyname = 'Users can view their own audit records'
    ) THEN
        CREATE POLICY "Users can view their own audit records"
          ON public.scientific_data_audit
          FOR SELECT
          USING (user_id = auth.uid());
    END IF;
END
$$;

-- Service role can access all audit records
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_policies 
        WHERE tablename = 'scientific_data_audit' 
        AND policyname = 'Service role can access all audit records'
    ) THEN
        CREATE POLICY "Service role can access all audit records"
          ON public.scientific_data_audit
          USING (auth.role() = 'service_role');
    END IF;
END
$$;

-- =========================
-- 7. Comment on RLS Implementation
-- =========================
COMMENT ON VIEW public.molecule_with_properties IS 'View for molecules with their properties, honors RLS policies of underlying tables';
COMMENT ON VIEW public.mixture_with_components IS 'View for mixtures with their components, honors RLS policies of underlying tables';
COMMENT ON VIEW public.experiment_with_results IS 'View for experiments with their results, honors RLS policies of underlying tables';
COMMENT ON TABLE public.scientific_data_audit IS 'Audit trail for scientific data operations with RLS security';

-- =========================
-- 8. Verification Query
-- =========================
-- Create a temporary function to report on RLS status
CREATE OR REPLACE FUNCTION verify_rls_implementation()
RETURNS TEXT AS $$
DECLARE
    report TEXT := '';
    table_count INT := 0;
    policy_count INT := 0;
    table_record RECORD;
    policy_record RECORD;
    view_record RECORD;
BEGIN
    report := report || '======= CryoProtect RLS Verification Report =======\n\n';
    report := report || 'Generated: ' || NOW() || '\n\n';
    
    -- Check tables with RLS
    report := report || '1. Tables with RLS:\n\n';
    FOR table_record IN 
        SELECT c.relname AS table_name, 
               c.relrowsecurity AS rls_enabled,
               (SELECT COUNT(*) FROM pg_policy WHERE polrelid = c.oid) AS policy_count
        FROM pg_class c
        JOIN pg_namespace n ON c.relnamespace = n.oid
        WHERE n.nspname = 'public'
          AND c.relkind = 'r'
          AND c.relname IN ('molecule', 'mixture', 'experiment', 'molecular_property', 
                           'mixture_component', 'experiment_property', 'migrations',
                           'scientific_data_audit')
        ORDER BY c.relname
    LOOP
        report := report || '   - ' || table_record.table_name || ': ';
        IF table_record.rls_enabled THEN
            report := report || 'RLS ENABLED with ' || table_record.policy_count || ' policies\n';
            table_count := table_count + 1;
            policy_count := policy_count + table_record.policy_count;
        ELSE
            report := report || 'RLS DISABLED\n';
        END IF;
    END LOOP;
    
    -- Check views with SECURITY INVOKER
    report := report || '\n2. Views with SECURITY INVOKER:\n\n';
    FOR view_record IN 
        SELECT viewname, 
               pg_catalog.pg_get_viewdef(c.oid) AS view_def
        FROM pg_catalog.pg_views
        JOIN pg_catalog.pg_class c ON c.relname = viewname
        WHERE schemaname = 'public'
          AND viewname IN ('molecule_with_properties', 'mixture_with_components', 'experiment_with_results')
    LOOP
        report := report || '   - ' || view_record.viewname || ': ';
        IF view_record.view_def LIKE '%SECURITY INVOKER%' THEN
            report := report || 'SECURITY INVOKER ✓\n';
        ELSE
            report := report || 'Missing SECURITY INVOKER ✗\n';
        END IF;
    END LOOP;

    -- Summary
    report := report || '\n3. Implementation Summary:\n\n';
    report := report || '   - ' || table_count || ' tables with RLS enabled\n';
    report := report || '   - ' || policy_count || ' total RLS policies defined\n';
    
    -- Check for performance indexes
    SELECT COUNT(*) INTO table_count
    FROM pg_indexes
    WHERE indexname IN ('idx_user_profile_project_id', 'idx_experiment_project_id', 
                       'idx_molecule_project_id', 'idx_mixture_project_id');
                       
    report := report || '   - ' || table_count || ' RLS-specific performance indexes\n';
    
    -- Check for audit triggers
    SELECT COUNT(*) INTO table_count
    FROM pg_trigger
    WHERE tgname IN ('audit_molecule_trigger', 'audit_mixture_trigger', 'audit_experiment_trigger');
    
    report := report || '   - ' || table_count || ' scientific data audit triggers\n';
    
    RETURN report;
END;
$$ LANGUAGE plpgsql;

-- Function to enable bypassing audit for bulk data loading
CREATE OR REPLACE FUNCTION set_audit_bypass(bypass BOOLEAN)
RETURNS VOID AS $$
BEGIN
    IF bypass THEN
        PERFORM set_config('app.bypass_audit', 'true', false);
    ELSE
        PERFORM set_config('app.bypass_audit', 'false', false);
    END IF;
END;
$$ LANGUAGE plpgsql;

-- =========================
-- 9. Commit Transaction
-- =========================
COMMIT;

-- =========================
-- 10. Run Verification
-- =========================
-- Run the verification and show results
SELECT verify_rls_implementation();

-- Clean up (drop verification function as it's only needed during installation)
DROP FUNCTION verify_rls_implementation();
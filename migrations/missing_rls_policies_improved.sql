-- CryoProtect Analyzer - Missing RLS Policies for Views and Tables (Improved)
-- This migration adds RLS policies for missing tables and views with performance optimizations
-- and scientific data security best practices.

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

-- Enable RLS on the view
ALTER VIEW public.experiment_with_results SECURITY INVOKER;

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
       current_setting('app.bypass_audit', true) = 'true' THEN
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
            'application_context', current_setting('app.context', true)
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, OLD.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             record_jsonb, null, 
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             current_setting('app.context', true));
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
            'application_context', current_setting('app.context', true)
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, NEW.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             to_jsonb(OLD), record_jsonb,
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             current_setting('app.context', true));
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
            'application_context', current_setting('app.context', true)
        );
        INSERT INTO public.scientific_data_audit 
            (table_name, record_id, operation, user_id, old_data, new_data, ip_address, application_context)
        VALUES 
            (TG_TABLE_NAME, NEW.id, TG_OP, coalesce(auth.uid(), '00000000-0000-0000-0000-000000000000'), 
             null, record_jsonb,
             coalesce(current_setting('request.headers', true)::jsonb->>'x-forwarded-for', 'unknown'),
             current_setting('app.context', true));
        RETURN NEW;
    END IF;
    RETURN NULL;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Apply audit triggers to scientific data tables (conditionally)
DO $$
BEGIN
    -- Add trigger to molecule table
    IF NOT EXISTS (
        SELECT 1 FROM pg_trigger 
        WHERE tgname = 'audit_molecule_trigger' 
        AND tgrelid = 'public.molecule'::regclass
    ) THEN
        CREATE TRIGGER audit_molecule_trigger
        AFTER INSERT OR UPDATE OR DELETE ON public.molecule
        FOR EACH ROW EXECUTE FUNCTION public.audit_scientific_data();
    END IF;

    -- Add trigger to mixture table
    IF NOT EXISTS (
        SELECT 1 FROM pg_trigger 
        WHERE tgname = 'audit_mixture_trigger' 
        AND tgrelid = 'public.mixture'::regclass
    ) THEN
        CREATE TRIGGER audit_mixture_trigger
        AFTER INSERT OR UPDATE OR DELETE ON public.mixture
        FOR EACH ROW EXECUTE FUNCTION public.audit_scientific_data();
    END IF;

    -- Add trigger to experiment table
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

-- Users can view audit records for their own actions
CREATE POLICY "Users can view their own audit records"
  ON public.scientific_data_audit
  FOR SELECT
  USING (user_id = auth.uid());

-- Service role can access all audit records
CREATE POLICY "Service role can access all audit records"
  ON public.scientific_data_audit
  USING (auth.role() = 'service_role');

-- =========================
-- 7. Comment on RLS Implementation
-- =========================
COMMENT ON VIEW public.molecule_with_properties IS 'View for molecules with their properties, honors RLS policies of underlying tables';
COMMENT ON VIEW public.mixture_with_components IS 'View for mixtures with their components, honors RLS policies of underlying tables';
COMMENT ON VIEW public.experiment_with_results IS 'View for experiments with their results, honors RLS policies of underlying tables';
COMMENT ON TABLE public.scientific_data_audit IS 'Audit trail for scientific data operations with RLS security';
-- Migration: Simplified RLS Security Definer Functions
-- Purpose: Improve performance of RLS policies with simplified security definer functions
-- This is a simplified version that should work with Supabase

BEGIN;

-- Step 1: Create Security Definer Functions

-- Function to check if current user is a project member
CREATE OR REPLACE FUNCTION public.is_project_member(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_member.project_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.is_project_member IS
'Efficiently checks if the current user is a member of the specified project.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to get all projects for the current user
CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT project_id FROM user_profile
  WHERE user_profile.user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.user_projects IS
'Returns a set of project IDs that the current user is a member of.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to check if current user is a team member
CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.team_id = is_team_member.team_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.is_team_member IS
'Efficiently checks if the current user is a member of the specified team.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to check if current user is a project owner
CREATE OR REPLACE FUNCTION public.is_project_owner(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_owner.project_id
      AND user_profile.user_id = auth.uid()
      AND user_profile.role = 'owner'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.is_project_owner IS
'Efficiently checks if the current user is the owner of the specified project.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to check if user has access to a molecule via project membership
CREATE OR REPLACE FUNCTION public.molecule_in_user_project(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecule
    JOIN user_profile ON user_profile.project_id = molecule.project_id
    WHERE molecule.id = molecule_in_user_project.molecule_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.molecule_in_user_project IS
'Efficiently checks if the current user has access to the specified molecule via project membership.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to check if user has access to a mixture via project membership
CREATE OR REPLACE FUNCTION public.mixture_in_user_project(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixture
    JOIN user_profile ON user_profile.project_id = mixture.project_id
    WHERE mixture.id = mixture_in_user_project.mixture_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.mixture_in_user_project IS
'Efficiently checks if the current user has access to the specified mixture via project membership.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Function to check if user has access to an experiment via project membership
CREATE OR REPLACE FUNCTION public.experiment_in_user_project(experiment_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM experiment
    JOIN user_profile ON user_profile.project_id = experiment.project_id
    WHERE experiment.id = experiment_in_user_project.experiment_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION public.experiment_in_user_project IS
'Efficiently checks if the current user has access to the specified experiment via project membership.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Step 2: Create indexes
-- This will help RLS policy performance

-- User profile indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_user_id ON user_profile(user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_project_id ON user_profile(project_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_role ON user_profile(role);

-- Project keys indexes
CREATE INDEX IF NOT EXISTS idx_molecule_project_id ON molecule(project_id);
CREATE INDEX IF NOT EXISTS idx_mixture_project_id ON mixture(project_id);
CREATE INDEX IF NOT EXISTS idx_experiment_project_id ON experiment(project_id);

-- Relationship indexes
CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON mixture_component(mixture_id);
CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON molecular_property(molecule_id);
CREATE INDEX IF NOT EXISTS idx_prediction_molecule_id ON prediction(molecule_id);
CREATE INDEX IF NOT EXISTS idx_experiment_property_experiment_id ON experiment_property(experiment_id);

-- Composite indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_user_project ON user_profile(user_id, project_id);

-- Step 3: Update RLS Policies to use the Security Definer Functions

-- 1. Molecule Table Policies
DROP POLICY IF EXISTS "Select molecules for project members" ON molecule;
CREATE POLICY "Select molecules for project members"
  ON molecule
  FOR SELECT
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Insert molecules for project members" ON molecule;
CREATE POLICY "Insert molecules for project members"
  ON molecule
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

DROP POLICY IF EXISTS "Update molecules for project members" ON molecule;
CREATE POLICY "Update molecules for project members"
  ON molecule
  FOR UPDATE
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Delete molecules for project members" ON molecule;
CREATE POLICY "Delete molecules for project members"
  ON molecule
  FOR DELETE
  USING (is_project_member(project_id));

-- 2. Mixture Table Policies
DROP POLICY IF EXISTS "Select mixtures for project members" ON mixture;
CREATE POLICY "Select mixtures for project members"
  ON mixture
  FOR SELECT
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Insert mixtures for project members" ON mixture;
CREATE POLICY "Insert mixtures for project members"
  ON mixture
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

DROP POLICY IF EXISTS "Update mixtures for project members" ON mixture;
CREATE POLICY "Update mixtures for project members"
  ON mixture
  FOR UPDATE
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Delete mixtures for project members" ON mixture;
CREATE POLICY "Delete mixtures for project members"
  ON mixture
  FOR DELETE
  USING (is_project_member(project_id));

-- 3. Mixture Component Table Policies
DROP POLICY IF EXISTS "Select mixture_components for project members" ON mixture_component;
CREATE POLICY "Select mixture_components for project members"
  ON mixture_component
  FOR SELECT
  USING (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Insert mixture_components for project members" ON mixture_component;
CREATE POLICY "Insert mixture_components for project members"
  ON mixture_component
  FOR INSERT
  WITH CHECK (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Update mixture_components for project members" ON mixture_component;
CREATE POLICY "Update mixture_components for project members"
  ON mixture_component
  FOR UPDATE
  USING (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Delete mixture_components for project members" ON mixture_component;
CREATE POLICY "Delete mixture_components for project members"
  ON mixture_component
  FOR DELETE
  USING (mixture_in_user_project(mixture_id));

-- 4. Experiment Table Policies
DROP POLICY IF EXISTS "Select experiments for project members" ON experiment;
CREATE POLICY "Select experiments for project members"
  ON experiment
  FOR SELECT
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Insert experiments for project members" ON experiment;
CREATE POLICY "Insert experiments for project members"
  ON experiment
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

DROP POLICY IF EXISTS "Update experiments for project members" ON experiment;
CREATE POLICY "Update experiments for project members"
  ON experiment
  FOR UPDATE
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Delete experiments for project members" ON experiment;
CREATE POLICY "Delete experiments for project members"
  ON experiment
  FOR DELETE
  USING (is_project_member(project_id));

-- 5. Molecular Property Table Policies
DROP POLICY IF EXISTS "Select molecular_properties for project members" ON molecular_property;
CREATE POLICY "Select molecular_properties for project members"
  ON molecular_property
  FOR SELECT
  USING (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Insert molecular_properties for project members" ON molecular_property;
CREATE POLICY "Insert molecular_properties for project members"
  ON molecular_property
  FOR INSERT
  WITH CHECK (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Update molecular_properties for project members" ON molecular_property;
CREATE POLICY "Update molecular_properties for project members"
  ON molecular_property
  FOR UPDATE
  USING (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Delete molecular_properties for project members" ON molecular_property;
CREATE POLICY "Delete molecular_properties for project members"
  ON molecular_property
  FOR DELETE
  USING (molecule_in_user_project(molecule_id));

-- 6. Prediction Table Policies
DROP POLICY IF EXISTS "Select predictions for project members" ON prediction;
CREATE POLICY "Select predictions for project members"
  ON prediction
  FOR SELECT
  USING (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Insert predictions for project members" ON prediction;
CREATE POLICY "Insert predictions for project members"
  ON prediction
  FOR INSERT
  WITH CHECK (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Update predictions for project members" ON prediction;
CREATE POLICY "Update predictions for project members"
  ON prediction
  FOR UPDATE
  USING (molecule_in_user_project(molecule_id));

DROP POLICY IF EXISTS "Delete predictions for project members" ON prediction;
CREATE POLICY "Delete predictions for project members"
  ON prediction
  FOR DELETE
  USING (molecule_in_user_project(molecule_id));

-- 7. Experiment Property Table Policies
DROP POLICY IF EXISTS "Select experiment_properties for project members" ON experiment_property;
CREATE POLICY "Select experiment_properties for project members"
  ON experiment_property
  FOR SELECT
  USING (experiment_in_user_project(experiment_id));

DROP POLICY IF EXISTS "Insert experiment_properties for project members" ON experiment_property;
CREATE POLICY "Insert experiment_properties for project members"
  ON experiment_property
  FOR INSERT
  WITH CHECK (experiment_in_user_project(experiment_id));

DROP POLICY IF EXISTS "Update experiment_properties for project members" ON experiment_property;
CREATE POLICY "Update experiment_properties for project members"
  ON experiment_property
  FOR UPDATE
  USING (experiment_in_user_project(experiment_id));

DROP POLICY IF EXISTS "Delete experiment_properties for project members" ON experiment_property;
CREATE POLICY "Delete experiment_properties for project members"
  ON experiment_property
  FOR DELETE
  USING (experiment_in_user_project(experiment_id));

-- 9. Project Table Policies
DROP POLICY IF EXISTS "Select projects for project members" ON project;
CREATE POLICY "Select projects for project members"
  ON project
  FOR SELECT
  USING (is_project_member(id));

-- Original policy retained for authenticated users
-- DROP POLICY IF EXISTS "Insert projects for authenticated users" ON project;

DROP POLICY IF EXISTS "Update projects for project members" ON project;
CREATE POLICY "Update projects for project members"
  ON project
  FOR UPDATE
  USING (is_project_member(id));

DROP POLICY IF EXISTS "Delete projects for project owners" ON project;
CREATE POLICY "Delete projects for project owners"
  ON project
  FOR DELETE
  USING (is_project_owner(id));

-- Step 4: Add Service Role Policies
CREATE OR REPLACE FUNCTION create_service_role_policies()
RETURNS void AS $$
DECLARE
    table_name text;
BEGIN
    FOR table_name IN 
        SELECT tablename FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename NOT LIKE 'pg_%'
    LOOP
        EXECUTE format('
            DROP POLICY IF EXISTS "service_role_all_access" ON %I;
            CREATE POLICY "service_role_all_access" ON %I
            FOR ALL
            TO service_role
            USING (true)
            WITH CHECK (true);
        ', table_name, table_name);
    END LOOP;
END
$$ LANGUAGE plpgsql;

SELECT create_service_role_policies();

-- Log completion
DO $$
BEGIN
    RAISE NOTICE 'Simplified RLS optimization migration completed successfully';
END $$;

COMMIT;
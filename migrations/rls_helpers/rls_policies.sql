-- Row Level Security Policies using security definer functions
-- These policies have improved performance due to the use of security definer functions

-- 1. Molecules table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_molecules" ON molecules;
CREATE POLICY "users_can_view_molecules" ON molecules
FOR SELECT
TO authenticated
USING (
  has_molecule_access(id)
);

-- Insert policy
DROP POLICY IF EXISTS "users_can_insert_molecules" ON molecules;
CREATE POLICY "users_can_insert_molecules" ON molecules
FOR INSERT
TO authenticated
WITH CHECK (
  auth.uid() = created_by
);

-- Update policy
DROP POLICY IF EXISTS "users_can_update_own_molecules" ON molecules;
CREATE POLICY "users_can_update_own_molecules" ON molecules
FOR UPDATE
TO authenticated
USING (
  auth.uid() = created_by
)
WITH CHECK (
  auth.uid() = created_by
);

-- Delete policy
DROP POLICY IF EXISTS "users_can_delete_own_molecules" ON molecules;
CREATE POLICY "users_can_delete_own_molecules" ON molecules
FOR DELETE
TO authenticated
USING (
  auth.uid() = created_by
);

-- 2. Molecular Properties table policies
-- Select policy with clearance level check
DROP POLICY IF EXISTS "users_can_view_properties" ON molecular_properties;
CREATE POLICY "users_can_view_properties" ON molecular_properties
FOR SELECT
TO authenticated
USING (
  has_molecule_access(molecule_id) AND
  (sensitivity_level IS NULL OR 
   sensitivity_level <= 'medium' OR 
   user_has_clearance(sensitivity_level))
);

-- Insert policy
DROP POLICY IF EXISTS "users_can_insert_properties" ON molecular_properties;
CREATE POLICY "users_can_insert_properties" ON molecular_properties
FOR INSERT
TO authenticated
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- Update policy
DROP POLICY IF EXISTS "users_can_update_properties" ON molecular_properties;
CREATE POLICY "users_can_update_properties" ON molecular_properties
FOR UPDATE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
)
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- Delete policy
DROP POLICY IF EXISTS "users_can_delete_properties" ON molecular_properties;
CREATE POLICY "users_can_delete_properties" ON molecular_properties
FOR DELETE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- 3. Mixtures table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_mixtures" ON mixtures;
CREATE POLICY "users_can_view_mixtures" ON mixtures
FOR SELECT
TO authenticated
USING (
  has_mixture_access(id)
);

-- Insert policy
DROP POLICY IF EXISTS "users_can_insert_mixtures" ON mixtures;
CREATE POLICY "users_can_insert_mixtures" ON mixtures
FOR INSERT
TO authenticated
WITH CHECK (
  auth.uid() = created_by
);

-- Update policy
DROP POLICY IF EXISTS "users_can_update_own_mixtures" ON mixtures;
CREATE POLICY "users_can_update_own_mixtures" ON mixtures
FOR UPDATE
TO authenticated
USING (
  auth.uid() = created_by
)
WITH CHECK (
  auth.uid() = created_by
);

-- Delete policy
DROP POLICY IF EXISTS "users_can_delete_own_mixtures" ON mixtures;
CREATE POLICY "users_can_delete_own_mixtures" ON mixtures
FOR DELETE
TO authenticated
USING (
  auth.uid() = created_by
);

-- 4. Mixture Components table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_mixture_components" ON mixture_components;
CREATE POLICY "users_can_view_mixture_components" ON mixture_components
FOR SELECT
TO authenticated
USING (
  has_mixture_access(mixture_id)
);

-- Insert, Update, Delete policies
DROP POLICY IF EXISTS "users_can_modify_mixture_components" ON mixture_components;
CREATE POLICY "users_can_modify_mixture_components" ON mixture_components
FOR ALL
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM mixtures
    WHERE id = mixture_components.mixture_id AND created_by = auth.uid()
  )
)
WITH CHECK (
  EXISTS (
    SELECT 1 FROM mixtures
    WHERE id = mixture_components.mixture_id AND created_by = auth.uid()
  )
);

-- 5. Experiments table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_experiments" ON experiments;
CREATE POLICY "users_can_view_experiments" ON experiments
FOR SELECT
TO authenticated
USING (
  has_experiment_access(id)
);

-- Insert policy
DROP POLICY IF EXISTS "users_can_insert_experiments" ON experiments;
CREATE POLICY "users_can_insert_experiments" ON experiments
FOR INSERT
TO authenticated
WITH CHECK (
  auth.uid() = created_by
);

-- Update policy
DROP POLICY IF EXISTS "users_can_update_own_experiments" ON experiments;
CREATE POLICY "users_can_update_own_experiments" ON experiments
FOR UPDATE
TO authenticated
USING (
  auth.uid() = created_by
)
WITH CHECK (
  auth.uid() = created_by
);

-- Delete policy
DROP POLICY IF EXISTS "users_can_delete_own_experiments" ON experiments;
CREATE POLICY "users_can_delete_own_experiments" ON experiments
FOR DELETE
TO authenticated
USING (
  auth.uid() = created_by
);

-- 6. Team Projects table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_team_projects" ON team_projects;
CREATE POLICY "users_can_view_team_projects" ON team_projects
FOR SELECT
TO authenticated
USING (
  is_team_member(team_id)
);

-- Insert, Update, Delete policies
DROP POLICY IF EXISTS "team_admins_can_modify_projects" ON team_projects;
CREATE POLICY "team_admins_can_modify_projects" ON team_projects
FOR ALL
TO authenticated
USING (
  is_team_admin(team_id)
)
WITH CHECK (
  is_team_admin(team_id)
);

-- 7. User Profile table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_all_profiles" ON user_profile;
CREATE POLICY "users_can_view_all_profiles" ON user_profile
FOR SELECT
TO authenticated
USING (true);

-- Users can only update their own profile
DROP POLICY IF EXISTS "users_can_update_own_profile" ON user_profile;
CREATE POLICY "users_can_update_own_profile" ON user_profile
FOR UPDATE
TO authenticated
USING (
  auth_user_id = auth.uid()
)
WITH CHECK (
  auth_user_id = auth.uid()
);

-- 8. Project members table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_project_members" ON project_members;
CREATE POLICY "users_can_view_project_members" ON project_members
FOR SELECT
TO authenticated
USING (
  has_project_access(project_id)
);

-- Only project admins can modify project members
DROP POLICY IF EXISTS "admins_can_modify_project_members" ON project_members;
CREATE POLICY "admins_can_modify_project_members" ON project_members
FOR ALL
TO authenticated
USING (
  is_project_admin(project_id)
)
WITH CHECK (
  is_project_admin(project_id)
);

-- 9. Team members table policies
-- Select policy
DROP POLICY IF EXISTS "users_can_view_team_members" ON team_members;
CREATE POLICY "users_can_view_team_members" ON team_members
FOR SELECT
TO authenticated
USING (
  is_team_member(team_id)
);

-- Only team admins can modify team members
DROP POLICY IF EXISTS "admins_can_modify_team_members" ON team_members;
CREATE POLICY "admins_can_modify_team_members" ON team_members
FOR ALL
TO authenticated
USING (
  is_team_admin(team_id)
)
WITH CHECK (
  is_team_admin(team_id)
);

-- 10. Create service role policies for all tables using DO block
DO $$
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
$$;
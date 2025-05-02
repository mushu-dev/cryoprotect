-- CryoProtect Analyzer RLS Policies Migration
-- Enables Row Level Security and defines policies for collaborative scientific data protection.
-- Policies restrict access to record owners, project members, or team members as appropriate.
-- Compatible with Supabase Auth JWT claims and user_profile table.

-- =========================
-- 1. molecule table
-- =========================
ALTER TABLE molecule ENABLE ROW LEVEL SECURITY;

-- Allow SELECT/UPDATE/DELETE only for users who are project members
CREATE POLICY "Select molecules for project members"
  ON molecule
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = molecule.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select molecules for project members" ON molecule IS
  'Allows project members to view molecules in their projects.';

CREATE POLICY "Insert molecules for project members"
  ON molecule
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = molecule.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert molecules for project members" ON molecule IS
  'Allows project members to add molecules to their projects.';

CREATE POLICY "Update molecules for project members"
  ON molecule
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = molecule.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update molecules for project members" ON molecule IS
  'Allows project members to update molecules in their projects.';

CREATE POLICY "Delete molecules for project members"
  ON molecule
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = molecule.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete molecules for project members" ON molecule IS
  'Allows project members to delete molecules in their projects.';

-- =========================
-- 2. mixture table
-- =========================
ALTER TABLE mixture ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select mixtures for project members"
  ON mixture
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = mixture.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select mixtures for project members" ON mixture IS
  'Allows project members to view mixtures in their projects.';

CREATE POLICY "Insert mixtures for project members"
  ON mixture
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = mixture.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert mixtures for project members" ON mixture IS
  'Allows project members to add mixtures to their projects.';

CREATE POLICY "Update mixtures for project members"
  ON mixture
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = mixture.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update mixtures for project members" ON mixture IS
  'Allows project members to update mixtures in their projects.';

CREATE POLICY "Delete mixtures for project members"
  ON mixture
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = mixture.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete mixtures for project members" ON mixture IS
  'Allows project members to delete mixtures in their projects.';

-- =========================
-- 3. mixture_component table
-- =========================
ALTER TABLE mixture_component ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select mixture_components for project members"
  ON mixture_component
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM mixture
      JOIN project ON project.id = mixture.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE mixture.id = mixture_component.mixture_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select mixture_components for project members" ON mixture_component IS
  'Allows project members to view mixture components in their projects.';

CREATE POLICY "Insert mixture_components for project members"
  ON mixture_component
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM mixture
      JOIN project ON project.id = mixture.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE mixture.id = mixture_component.mixture_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert mixture_components for project members" ON mixture_component IS
  'Allows project members to add mixture components to their projects.';

CREATE POLICY "Update mixture_components for project members"
  ON mixture_component
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM mixture
      JOIN project ON project.id = mixture.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE mixture.id = mixture_component.mixture_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update mixture_components for project members" ON mixture_component IS
  'Allows project members to update mixture components in their projects.';

CREATE POLICY "Delete mixture_components for project members"
  ON mixture_component
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM mixture
      JOIN project ON project.id = mixture.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE mixture.id = mixture_component.mixture_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete mixture_components for project members" ON mixture_component IS
  'Allows project members to delete mixture components in their projects.';

-- =========================
-- 4. experiment table
-- =========================
ALTER TABLE experiment ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select experiments for project members"
  ON experiment
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = experiment.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select experiments for project members" ON experiment IS
  'Allows project members to view experiments in their projects.';

CREATE POLICY "Insert experiments for project members"
  ON experiment
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = experiment.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert experiments for project members" ON experiment IS
  'Allows project members to add experiments to their projects.';

CREATE POLICY "Update experiments for project members"
  ON experiment
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = experiment.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update experiments for project members" ON experiment IS
  'Allows project members to update experiments in their projects.';

CREATE POLICY "Delete experiments for project members"
  ON experiment
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = experiment.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete experiments for project members" ON experiment IS
  'Allows project members to delete experiments in their projects.';

-- =========================
-- 5. molecular_property table
-- =========================
ALTER TABLE molecular_property ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select molecular_properties for project members"
  ON molecular_property
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = molecular_property.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select molecular_properties for project members" ON molecular_property IS
  'Allows project members to view molecular properties in their projects.';

CREATE POLICY "Insert molecular_properties for project members"
  ON molecular_property
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = molecular_property.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert molecular_properties for project members" ON molecular_property IS
  'Allows project members to add molecular properties to their projects.';

CREATE POLICY "Update molecular_properties for project members"
  ON molecular_property
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = molecular_property.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update molecular_properties for project members" ON molecular_property IS
  'Allows project members to update molecular properties in their projects.';

CREATE POLICY "Delete molecular_properties for project members"
  ON molecular_property
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = molecular_property.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete molecular_properties for project members" ON molecular_property IS
  'Allows project members to delete molecular properties in their projects.';

-- =========================
-- 6. prediction table
-- =========================
ALTER TABLE prediction ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select predictions for project members"
  ON prediction
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = prediction.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select predictions for project members" ON prediction IS
  'Allows project members to view predictions in their projects.';

CREATE POLICY "Insert predictions for project members"
  ON prediction
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = prediction.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert predictions for project members" ON prediction IS
  'Allows project members to add predictions to their projects.';

CREATE POLICY "Update predictions for project members"
  ON prediction
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = prediction.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update predictions for project members" ON prediction IS
  'Allows project members to update predictions in their projects.';

CREATE POLICY "Delete predictions for project members"
  ON prediction
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM molecule
      JOIN project ON project.id = molecule.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE molecule.id = prediction.molecule_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete predictions for project members" ON prediction IS
  'Allows project members to delete predictions in their projects.';

-- =========================
-- 7. experiment_property table
-- =========================
ALTER TABLE experiment_property ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select experiment_properties for project members"
  ON experiment_property
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM experiment
      JOIN project ON project.id = experiment.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE experiment.id = experiment_property.experiment_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select experiment_properties for project members" ON experiment_property IS
  'Allows project members to view experiment properties in their projects.';

CREATE POLICY "Insert experiment_properties for project members"
  ON experiment_property
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM experiment
      JOIN project ON project.id = experiment.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE experiment.id = experiment_property.experiment_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert experiment_properties for project members" ON experiment_property IS
  'Allows project members to add experiment properties to their projects.';

CREATE POLICY "Update experiment_properties for project members"
  ON experiment_property
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM experiment
      JOIN project ON project.id = experiment.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE experiment.id = experiment_property.experiment_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update experiment_properties for project members" ON experiment_property IS
  'Allows project members to update experiment properties in their projects.';

CREATE POLICY "Delete experiment_properties for project members"
  ON experiment_property
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM experiment
      JOIN project ON project.id = experiment.project_id
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE experiment.id = experiment_property.experiment_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete experiment_properties for project members" ON experiment_property IS
  'Allows project members to delete experiment properties in their projects.';

-- =========================
-- 8. calculation_method table
-- =========================
ALTER TABLE calculation_method ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select calculation_methods for all authenticated users"
  ON calculation_method
  FOR SELECT
  USING (auth.role() = 'authenticated');
COMMENT ON POLICY "Select calculation_methods for all authenticated users" ON calculation_method IS
  'Allows all authenticated users to view calculation methods.';

-- Optionally restrict INSERT/UPDATE/DELETE to admins or project members if needed.

-- =========================
-- 9. project table
-- =========================
ALTER TABLE project ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select projects for project members"
  ON project
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.project_id = project.id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select projects for project members" ON project IS
  'Allows users to view projects they are members of.';

CREATE POLICY "Insert projects for authenticated users"
  ON project
  FOR INSERT
  WITH CHECK (auth.role() = 'authenticated');
COMMENT ON POLICY "Insert projects for authenticated users" ON project IS
  'Allows any authenticated user to create a new project.';

CREATE POLICY "Update projects for project members"
  ON project
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.project_id = project.id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update projects for project members" ON project IS
  'Allows project members to update their projects.';

CREATE POLICY "Delete projects for project owners"
  ON project
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.project_id = project.id
        AND user_profile.user_id = auth.uid()
        AND user_profile.role = 'owner'
    )
  );
COMMENT ON POLICY "Delete projects for project owners" ON project IS
  'Allows only project owners to delete their projects.';

-- =========================
-- 10. team table
-- =========================
ALTER TABLE team ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select teams for team members"
  ON team
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.team_id = team.id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select teams for team members" ON team IS
  'Allows users to view teams they are members of.';

CREATE POLICY "Insert teams for authenticated users"
  ON team
  FOR INSERT
  WITH CHECK (auth.role() = 'authenticated');
COMMENT ON POLICY "Insert teams for authenticated users" ON team IS
  'Allows any authenticated user to create a new team.';

CREATE POLICY "Update teams for team members"
  ON team
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.team_id = team.id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update teams for team members" ON team IS
  'Allows team members to update their teams.';

CREATE POLICY "Delete teams for team owners"
  ON team
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM user_profile
      WHERE user_profile.team_id = team.id
        AND user_profile.user_id = auth.uid()
        AND user_profile.role = 'owner'
    )
  );
COMMENT ON POLICY "Delete teams for team owners" ON team IS
  'Allows only team owners to delete their teams.';

-- =========================
-- 11. user_profile table
-- =========================
ALTER TABLE user_profile ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select own user_profile"
  ON user_profile
  FOR SELECT
  USING (user_id = auth.uid());
COMMENT ON POLICY "Select own user_profile" ON user_profile IS
  'Allows users to view their own profile.';

CREATE POLICY "Insert own user_profile"
  ON user_profile
  FOR INSERT
  WITH CHECK (user_id = auth.uid());
COMMENT ON POLICY "Insert own user_profile" ON user_profile IS
  'Allows users to create their own profile.';

CREATE POLICY "Update own user_profile"
  ON user_profile
  FOR UPDATE
  USING (user_id = auth.uid());
COMMENT ON POLICY "Update own user_profile" ON user_profile IS
  'Allows users to update their own profile.';

CREATE POLICY "Delete own user_profile"
  ON user_profile
  FOR DELETE
  USING (user_id = auth.uid());
COMMENT ON POLICY "Delete own user_profile" ON user_profile IS
  'Allows users to delete their own profile.';

-- =========================
-- 12. Example: graph-based relationship table (e.g., project_relationship)
-- =========================
-- Replace 'project_relationship' with actual table name if different
ALTER TABLE project_relationship ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Select project_relationships for project members"
  ON project_relationship
  FOR SELECT
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = project_relationship.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Select project_relationships for project members" ON project_relationship IS
  'Allows project members to view relationships in their projects.';

CREATE POLICY "Insert project_relationships for project members"
  ON project_relationship
  FOR INSERT
  WITH CHECK (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = project_relationship.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Insert project_relationships for project members" ON project_relationship IS
  'Allows project members to add relationships in their projects.';

CREATE POLICY "Update project_relationships for project members"
  ON project_relationship
  FOR UPDATE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = project_relationship.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Update project_relationships for project members" ON project_relationship IS
  'Allows project members to update relationships in their projects.';

CREATE POLICY "Delete project_relationships for project members"
  ON project_relationship
  FOR DELETE
  USING (
    EXISTS (
      SELECT 1 FROM project
      JOIN user_profile ON user_profile.project_id = project.id
      WHERE project.id = project_relationship.project_id
        AND user_profile.user_id = auth.uid()
    )
  );
COMMENT ON POLICY "Delete project_relationships for project members" ON project_relationship IS
  'Allows project members to delete relationships in their projects.';

-- =========================
-- Add similar policies for other graph-based relationship tables as needed.
-- =========================

-- End of RLS policies migration
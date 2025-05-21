-- Create indexes for better RLS policy performance
CREATE INDEX IF NOT EXISTS idx_team_members_user_id ON team_members (user_id);
CREATE INDEX IF NOT EXISTS idx_team_members_team_id ON team_members (team_id);
CREATE INDEX IF NOT EXISTS idx_team_members_role ON team_members (role);
CREATE INDEX IF NOT EXISTS idx_projects_team_id ON projects (team_id);
CREATE INDEX IF NOT EXISTS idx_experiments_project_id ON experiments (project_id);
CREATE INDEX IF NOT EXISTS idx_experiments_molecule_id ON experiments (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
CREATE INDEX IF NOT EXISTS idx_mixtures_project_id ON mixtures (project_id);
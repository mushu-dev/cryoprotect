-- Create indexes for better RLS policy performance
CREATE INDEX IF NOT EXISTS idx_user_teams_user_id ON user_teams (user_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_team_id ON user_teams (team_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_role ON user_teams (role);
CREATE INDEX IF NOT EXISTS idx_projects_team_id ON projects (team_id);
CREATE INDEX IF NOT EXISTS idx_experiments_project_id ON experiments (project_id);
CREATE INDEX IF NOT EXISTS idx_experiments_molecule_id ON experiments (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
CREATE INDEX IF NOT EXISTS idx_mixtures_project_id ON mixtures (project_id);
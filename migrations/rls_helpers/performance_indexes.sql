-- Performance indexes for RLS policy optimization
-- These indexes improve the performance of RLS policy evaluation and security definer functions

-- 1. User Profile Indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id ON user_profile(auth_user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_clearance_level ON user_profile(clearance_level);
CREATE INDEX IF NOT EXISTS idx_user_profile_role ON user_profile(role);

-- 2. Team Projects Indexes
CREATE INDEX IF NOT EXISTS idx_team_projects_team_id ON team_projects(team_id);
CREATE INDEX IF NOT EXISTS idx_team_projects_project_id ON team_projects(project_id);

-- 3. Project Molecule Indexes
CREATE INDEX IF NOT EXISTS idx_project_molecules_project_id ON project_molecules(project_id);
CREATE INDEX IF NOT EXISTS idx_project_molecules_molecule_id ON project_molecules(molecule_id);

-- 4. Project Mixture Indexes
CREATE INDEX IF NOT EXISTS idx_project_mixtures_project_id ON project_mixtures(project_id);
CREATE INDEX IF NOT EXISTS idx_project_mixtures_mixture_id ON project_mixtures(mixture_id);

-- 5. Project Experiment Indexes
CREATE INDEX IF NOT EXISTS idx_project_experiments_project_id ON project_experiments(project_id);
CREATE INDEX IF NOT EXISTS idx_project_experiments_experiment_id ON project_experiments(experiment_id);

-- 6. Molecule Ownership and Visibility Indexes
CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON molecules(created_by);
CREATE INDEX IF NOT EXISTS idx_molecules_is_public ON molecules(is_public);

-- 7. Mixture Ownership and Visibility Indexes
CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON mixtures(created_by);
CREATE INDEX IF NOT EXISTS idx_mixtures_is_public ON mixtures(is_public);

-- 8. Experiment Ownership and Visibility Indexes
CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON experiments(created_by);
CREATE INDEX IF NOT EXISTS idx_experiments_is_public ON experiments(is_public);

-- 9. Molecular Properties Indexes
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties(molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_name ON molecular_properties(property_name);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type ON molecular_properties(property_type);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitivity_level ON molecular_properties(sensitivity_level);

-- 10. Team Members Indexes
CREATE INDEX IF NOT EXISTS idx_team_members_team_id ON team_members(team_id);
CREATE INDEX IF NOT EXISTS idx_team_members_user_id ON team_members(user_id);
CREATE INDEX IF NOT EXISTS idx_team_members_role ON team_members(role);

-- 11. Project Members Indexes
CREATE INDEX IF NOT EXISTS idx_project_members_project_id ON project_members(project_id);
CREATE INDEX IF NOT EXISTS idx_project_members_user_id ON project_members(user_id);
CREATE INDEX IF NOT EXISTS idx_project_members_role ON project_members(role);

-- 12. Partial Indexes for Common Access Patterns

-- Partial index for public molecules (very selective)
CREATE INDEX IF NOT EXISTS idx_molecules_public ON molecules(id) WHERE is_public = true;

-- Partial index for public mixtures
CREATE INDEX IF NOT EXISTS idx_mixtures_public ON mixtures(id) WHERE is_public = true;

-- Partial index for public experiments
CREATE INDEX IF NOT EXISTS idx_experiments_public ON experiments(id) WHERE is_public = true;

-- Partial index for sensitive content
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitive 
ON molecular_properties(molecule_id, sensitivity_level)
WHERE sensitivity_level > 'low';

-- 13. Composite Indexes for Common Query Patterns

-- Composite index for team membership queries
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_team 
ON user_profile(auth_user_id, team_id);

-- Composite index for project-mixture relationships
CREATE INDEX IF NOT EXISTS idx_project_mixtures_composite
ON project_mixtures(mixture_id, project_id);

-- Composite index for project-molecule relationships
CREATE INDEX IF NOT EXISTS idx_project_molecules_composite
ON project_molecules(molecule_id, project_id);

-- Composite index for project-experiment relationships
CREATE INDEX IF NOT EXISTS idx_project_experiments_composite
ON project_experiments(experiment_id, project_id);

-- Composite index for property searches
CREATE INDEX IF NOT EXISTS idx_molecular_properties_search
ON molecular_properties(molecule_id, property_name, property_type);

-- Composite index for property value ranges
CREATE INDEX IF NOT EXISTS idx_molecular_properties_value_numeric
ON molecular_properties(property_name, (property_value::numeric))
WHERE property_value ~ '^-?[0-9]+(\.[0-9]+)?$';

-- 14. GIN indexes for text search capabilities
CREATE INDEX IF NOT EXISTS idx_molecules_name_gin
ON molecules USING GIN (to_tsvector('english', name));

CREATE INDEX IF NOT EXISTS idx_molecular_properties_text_gin
ON molecular_properties USING GIN (to_tsvector('english', property_value))
WHERE property_type = 'text';
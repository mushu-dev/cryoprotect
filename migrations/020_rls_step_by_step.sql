-- Migration: 020_rls_step_by_step.sql
-- Description: Implements improvements to RLS policies based on verification results
-- Created: 2025-05-11

-- Create individual functions with simplified syntax

-- Function: is_team_member
CREATE OR REPLACE FUNCTION auth.is_team_member(team_id UUID) RETURNS BOOLEAN 
AS 'SELECT EXISTS (SELECT 1 FROM user_teams WHERE team_id = $1 AND user_id = auth.uid());'
LANGUAGE SQL SECURITY DEFINER;

-- Function: is_team_admin
CREATE OR REPLACE FUNCTION auth.is_team_admin(team_id UUID) RETURNS BOOLEAN 
AS 'SELECT EXISTS (SELECT 1 FROM user_teams WHERE team_id = $1 AND user_id = auth.uid() AND role = ''admin'');'
LANGUAGE SQL SECURITY DEFINER;

-- Function: get_user_teams
CREATE OR REPLACE FUNCTION auth.get_user_teams() RETURNS SETOF UUID 
AS 'SELECT team_id FROM user_teams WHERE user_id = auth.uid();'
LANGUAGE SQL SECURITY DEFINER;

-- Function: can_access_project
CREATE OR REPLACE FUNCTION auth.can_access_project(project_id UUID) RETURNS BOOLEAN 
AS 'SELECT auth.is_team_member((SELECT team_id FROM projects WHERE id = $1));'
LANGUAGE SQL SECURITY DEFINER;

-- Function: can_manage_project
CREATE OR REPLACE FUNCTION auth.can_manage_project(project_id UUID) RETURNS BOOLEAN 
AS 'SELECT auth.is_team_admin((SELECT team_id FROM projects WHERE id = $1));'
LANGUAGE SQL SECURITY DEFINER;

-- Function: get_accessible_projects
CREATE OR REPLACE FUNCTION auth.get_accessible_projects() RETURNS SETOF UUID 
AS 'SELECT p.id FROM projects p JOIN teams t ON p.team_id = t.id JOIN user_teams ut ON t.id = ut.team_id WHERE ut.user_id = auth.uid();'
LANGUAGE SQL SECURITY DEFINER;

-- Function: can_access_molecule
CREATE OR REPLACE FUNCTION auth.can_access_molecule(molecule_id UUID) RETURNS BOOLEAN 
AS 'SELECT EXISTS (SELECT 1 FROM experiments e JOIN projects p ON e.project_id = p.id JOIN teams t ON p.team_id = t.id JOIN user_teams ut ON t.id = ut.team_id WHERE e.molecule_id = $1 AND ut.user_id = auth.uid());'
LANGUAGE SQL SECURITY DEFINER;

-- Now apply all the RLS policies

-- Projects table policies
DROP POLICY IF EXISTS projects_select_policy ON projects;
CREATE POLICY projects_select_policy ON projects
    FOR SELECT USING (auth.can_access_project(id));

DROP POLICY IF EXISTS projects_insert_policy ON projects;
CREATE POLICY projects_insert_policy ON projects
    FOR INSERT WITH CHECK (auth.is_team_member(team_id));

DROP POLICY IF EXISTS projects_update_policy ON projects;
CREATE POLICY projects_update_policy ON projects
    FOR UPDATE USING (auth.can_manage_project(id)) WITH CHECK (auth.can_manage_project(id));

DROP POLICY IF EXISTS projects_delete_policy ON projects;
CREATE POLICY projects_delete_policy ON projects
    FOR DELETE USING (auth.can_manage_project(id));

-- Experiments table policies
DROP POLICY IF EXISTS experiments_select_policy ON experiments;
CREATE POLICY experiments_select_policy ON experiments
    FOR SELECT USING (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_insert_policy ON experiments;
CREATE POLICY experiments_insert_policy ON experiments
    FOR INSERT WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_update_policy ON experiments;
CREATE POLICY experiments_update_policy ON experiments
    FOR UPDATE USING (auth.can_access_project(project_id)) WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_delete_policy ON experiments;
CREATE POLICY experiments_delete_policy ON experiments
    FOR DELETE USING (auth.can_access_project(project_id));

-- Molecules table policies
DROP POLICY IF EXISTS molecules_select_policy ON molecules;
CREATE POLICY molecules_select_policy ON molecules
    FOR SELECT USING (EXISTS (SELECT 1 FROM experiments e WHERE e.molecule_id = id AND auth.can_access_project(e.project_id)));

DROP POLICY IF EXISTS molecules_insert_policy ON molecules;
CREATE POLICY molecules_insert_policy ON molecules
    FOR INSERT WITH CHECK (true);  -- Allow insert, association happens in experiments

DROP POLICY IF EXISTS molecules_update_policy ON molecules;
CREATE POLICY molecules_update_policy ON molecules
    FOR UPDATE USING (auth.can_access_molecule(id)) WITH CHECK (auth.can_access_molecule(id));

DROP POLICY IF EXISTS molecules_delete_policy ON molecules;
CREATE POLICY molecules_delete_policy ON molecules
    FOR DELETE USING (auth.can_access_molecule(id));

-- Molecular properties table policies
DROP POLICY IF EXISTS molecular_properties_select_policy ON molecular_properties;
CREATE POLICY molecular_properties_select_policy ON molecular_properties
    FOR SELECT USING (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_insert_policy ON molecular_properties;
CREATE POLICY molecular_properties_insert_policy ON molecular_properties
    FOR INSERT WITH CHECK (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_update_policy ON molecular_properties;
CREATE POLICY molecular_properties_update_policy ON molecular_properties
    FOR UPDATE USING (auth.can_access_molecule(molecule_id)) WITH CHECK (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_delete_policy ON molecular_properties;
CREATE POLICY molecular_properties_delete_policy ON molecular_properties
    FOR DELETE USING (auth.can_access_molecule(molecule_id));

-- Mixtures table policies
DROP POLICY IF EXISTS mixtures_select_policy ON mixtures;
CREATE POLICY mixtures_select_policy ON mixtures
    FOR SELECT USING (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_insert_policy ON mixtures;
CREATE POLICY mixtures_insert_policy ON mixtures
    FOR INSERT WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_update_policy ON mixtures;
CREATE POLICY mixtures_update_policy ON mixtures
    FOR UPDATE USING (auth.can_access_project(project_id)) WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_delete_policy ON mixtures;
CREATE POLICY mixtures_delete_policy ON mixtures
    FOR DELETE USING (auth.can_access_project(project_id));

-- User teams table policies
DROP POLICY IF EXISTS user_teams_select_policy ON user_teams;
CREATE POLICY user_teams_select_policy ON user_teams
    FOR SELECT USING (user_id = auth.uid() OR team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_insert_policy ON user_teams;
CREATE POLICY user_teams_insert_policy ON user_teams
    FOR INSERT WITH CHECK (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_update_policy ON user_teams;
CREATE POLICY user_teams_update_policy ON user_teams
    FOR UPDATE USING (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')) 
    WITH CHECK (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_delete_policy ON user_teams;
CREATE POLICY user_teams_delete_policy ON user_teams
    FOR DELETE USING (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

-- Update indexes for better RLS policy performance
CREATE INDEX IF NOT EXISTS idx_user_teams_user_id ON user_teams (user_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_team_id ON user_teams (team_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_role ON user_teams (role);
CREATE INDEX IF NOT EXISTS idx_projects_team_id ON projects (team_id);
CREATE INDEX IF NOT EXISTS idx_experiments_project_id ON experiments (project_id);
CREATE INDEX IF NOT EXISTS idx_experiments_molecule_id ON experiments (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
CREATE INDEX IF NOT EXISTS idx_mixtures_project_id ON mixtures (project_id);

-- Set up service role access
DROP POLICY IF EXISTS projects_service_role_policy ON projects;
CREATE POLICY projects_service_role_policy ON projects FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS experiments_service_role_policy ON experiments;
CREATE POLICY experiments_service_role_policy ON experiments FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS molecules_service_role_policy ON molecules;
CREATE POLICY molecules_service_role_policy ON molecules FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS molecular_properties_service_role_policy ON molecular_properties;
CREATE POLICY molecular_properties_service_role_policy ON molecular_properties FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS mixtures_service_role_policy ON mixtures;
CREATE POLICY mixtures_service_role_policy ON mixtures FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS teams_service_role_policy ON teams;
CREATE POLICY teams_service_role_policy ON teams FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS user_teams_service_role_policy ON user_teams;
CREATE POLICY user_teams_service_role_policy ON user_teams FOR ALL USING (auth.is_service_role());

-- Ensure RLS is enabled on all tables
ALTER TABLE projects ENABLE ROW LEVEL SECURITY;
ALTER TABLE experiments ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecular_properties ENABLE ROW LEVEL SECURITY;
ALTER TABLE mixtures ENABLE ROW LEVEL SECURITY;
ALTER TABLE teams ENABLE ROW LEVEL SECURITY;
ALTER TABLE user_teams ENABLE ROW LEVEL SECURITY;
ALTER TABLE user_profile ENABLE ROW LEVEL SECURITY;
ALTER TABLE experiment_properties ENABLE ROW LEVEL SECURITY;

-- Create policies for user_profile and experiment_properties
DROP POLICY IF EXISTS user_profile_select_policy ON user_profile;
CREATE POLICY user_profile_select_policy ON user_profile
    FOR SELECT USING (auth_user_id = auth.uid() OR auth.is_service_role());

DROP POLICY IF EXISTS user_profile_update_policy ON user_profile;
CREATE POLICY user_profile_update_policy ON user_profile
    FOR UPDATE USING (auth_user_id = auth.uid()) WITH CHECK (auth_user_id = auth.uid());

DROP POLICY IF EXISTS experiment_properties_select_policy ON experiment_properties;
CREATE POLICY experiment_properties_select_policy ON experiment_properties
    FOR SELECT USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_insert_policy ON experiment_properties;
CREATE POLICY experiment_properties_insert_policy ON experiment_properties
    FOR INSERT WITH CHECK (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_update_policy ON experiment_properties;
CREATE POLICY experiment_properties_update_policy ON experiment_properties
    FOR UPDATE USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role())
    WITH CHECK (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_delete_policy ON experiment_properties;
CREATE POLICY experiment_properties_delete_policy ON experiment_properties
    FOR DELETE USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());
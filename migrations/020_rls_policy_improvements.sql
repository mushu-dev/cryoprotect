-- Migration: 020_rls_policy_improvements.sql
-- Description: Implements improvements to RLS policies based on verification results
-- Created: 2025-05-11

-- Step 1: Create missing security definer functions for team-based access

-- Function to check if a user is a member of a team
CREATE OR REPLACE FUNCTION auth.is_team_member(team_id UUID)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 FROM user_teams
        WHERE team_id = $1
        AND user_id = auth.uid()
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to check if a user is an admin of a team
CREATE OR REPLACE FUNCTION auth.is_team_admin(team_id UUID)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 FROM user_teams
        WHERE team_id = $1
        AND user_id = auth.uid()
        AND role = 'admin'
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to get teams that a user belongs to
CREATE OR REPLACE FUNCTION auth.get_user_teams()
RETURNS SETOF UUID AS $$
BEGIN
    RETURN QUERY
    SELECT team_id FROM user_teams
    WHERE user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to check if a project belongs to a team the user is a member of
CREATE OR REPLACE FUNCTION auth.can_access_project(project_id UUID)
RETURNS BOOLEAN AS $$
DECLARE
    project_team_id UUID;
BEGIN
    -- Get the team_id for the project
    SELECT team_id INTO project_team_id
    FROM projects
    WHERE id = project_id;
    
    -- Check if user is a member of this team
    RETURN auth.is_team_member(project_team_id);
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to check if user can manage a project
CREATE OR REPLACE FUNCTION auth.can_manage_project(project_id UUID)
RETURNS BOOLEAN AS $$
DECLARE
    project_team_id UUID;
BEGIN
    -- Get the team_id for the project
    SELECT team_id INTO project_team_id
    FROM projects
    WHERE id = project_id;
    
    -- Check if user is an admin of this team
    RETURN auth.is_team_admin(project_team_id);
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to get accessible projects for a user
CREATE OR REPLACE FUNCTION auth.get_accessible_projects()
RETURNS SETOF UUID AS $$
BEGIN
    RETURN QUERY
    SELECT p.id
    FROM projects p
    JOIN teams t ON p.team_id = t.id
    JOIN user_teams ut ON t.id = ut.team_id
    WHERE ut.user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to check if a molecule is accessible to the user through projects
CREATE OR REPLACE FUNCTION auth.can_access_molecule(molecule_id UUID)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM experiments e
        JOIN projects p ON e.project_id = p.id
        JOIN teams t ON p.team_id = t.id
        JOIN user_teams ut ON t.id = ut.team_id
        WHERE e.molecule_id = molecule_id
        AND ut.user_id = auth.uid()
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 2: Update RLS policies on core tables to use the security definer functions

-- Drop existing RLS policies on projects table
DROP POLICY IF EXISTS projects_select_policy ON projects;
DROP POLICY IF EXISTS projects_insert_policy ON projects;
DROP POLICY IF EXISTS projects_update_policy ON projects;
DROP POLICY IF EXISTS projects_delete_policy ON projects;

-- Create improved RLS policies on projects table
CREATE POLICY projects_select_policy ON projects
    FOR SELECT USING (
        auth.can_access_project(id)
    );

CREATE POLICY projects_insert_policy ON projects
    FOR INSERT WITH CHECK (
        auth.is_team_member(team_id)
    );

CREATE POLICY projects_update_policy ON projects
    FOR UPDATE USING (
        auth.can_manage_project(id)
    ) WITH CHECK (
        auth.can_manage_project(id)
    );

CREATE POLICY projects_delete_policy ON projects
    FOR DELETE USING (
        auth.can_manage_project(id)
    );

-- Drop existing RLS policies on experiments table
DROP POLICY IF EXISTS experiments_select_policy ON experiments;
DROP POLICY IF EXISTS experiments_insert_policy ON experiments;
DROP POLICY IF EXISTS experiments_update_policy ON experiments;
DROP POLICY IF EXISTS experiments_delete_policy ON experiments;

-- Create improved RLS policies on experiments table
CREATE POLICY experiments_select_policy ON experiments
    FOR SELECT USING (
        auth.can_access_project(project_id)
    );

CREATE POLICY experiments_insert_policy ON experiments
    FOR INSERT WITH CHECK (
        auth.can_access_project(project_id)
    );

CREATE POLICY experiments_update_policy ON experiments
    FOR UPDATE USING (
        auth.can_access_project(project_id)
    ) WITH CHECK (
        auth.can_access_project(project_id)
    );

CREATE POLICY experiments_delete_policy ON experiments
    FOR DELETE USING (
        auth.can_access_project(project_id)
    );

-- Drop existing RLS policies on molecules table
DROP POLICY IF EXISTS molecules_select_policy ON molecules;
DROP POLICY IF EXISTS molecules_insert_policy ON molecules;
DROP POLICY IF EXISTS molecules_update_policy ON molecules;
DROP POLICY IF EXISTS molecules_delete_policy ON molecules;

-- Create improved RLS policies on molecules table based on experiment association
CREATE POLICY molecules_select_policy ON molecules
    FOR SELECT USING (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.molecule_id = id
            AND auth.can_access_project(e.project_id)
        )
    );

CREATE POLICY molecules_insert_policy ON molecules
    FOR INSERT WITH CHECK (true);  -- Allow insert, association happens in experiments

CREATE POLICY molecules_update_policy ON molecules
    FOR UPDATE USING (
        auth.can_access_molecule(id)
    ) WITH CHECK (
        auth.can_access_molecule(id)
    );

CREATE POLICY molecules_delete_policy ON molecules
    FOR DELETE USING (
        auth.can_access_molecule(id)
    );

-- Drop existing RLS policies on molecular_properties table
DROP POLICY IF EXISTS molecular_properties_select_policy ON molecular_properties;
DROP POLICY IF EXISTS molecular_properties_insert_policy ON molecular_properties;
DROP POLICY IF EXISTS molecular_properties_update_policy ON molecular_properties;
DROP POLICY IF EXISTS molecular_properties_delete_policy ON molecular_properties;

-- Create improved RLS policies on molecular_properties table
CREATE POLICY molecular_properties_select_policy ON molecular_properties
    FOR SELECT USING (
        auth.can_access_molecule(molecule_id)
    );

CREATE POLICY molecular_properties_insert_policy ON molecular_properties
    FOR INSERT WITH CHECK (
        auth.can_access_molecule(molecule_id)
    );

CREATE POLICY molecular_properties_update_policy ON molecular_properties
    FOR UPDATE USING (
        auth.can_access_molecule(molecule_id)
    ) WITH CHECK (
        auth.can_access_molecule(molecule_id)
    );

CREATE POLICY molecular_properties_delete_policy ON molecular_properties
    FOR DELETE USING (
        auth.can_access_molecule(molecule_id)
    );

-- Drop existing RLS policies on mixtures table
DROP POLICY IF EXISTS mixtures_select_policy ON mixtures;
DROP POLICY IF EXISTS mixtures_insert_policy ON mixtures;
DROP POLICY IF EXISTS mixtures_update_policy ON mixtures;
DROP POLICY IF EXISTS mixtures_delete_policy ON mixtures;

-- Create improved RLS policies on mixtures table
CREATE POLICY mixtures_select_policy ON mixtures
    FOR SELECT USING (
        auth.can_access_project(project_id)
    );

CREATE POLICY mixtures_insert_policy ON mixtures
    FOR INSERT WITH CHECK (
        auth.can_access_project(project_id)
    );

CREATE POLICY mixtures_update_policy ON mixtures
    FOR UPDATE USING (
        auth.can_access_project(project_id)
    ) WITH CHECK (
        auth.can_access_project(project_id)
    );

CREATE POLICY mixtures_delete_policy ON mixtures
    FOR DELETE USING (
        auth.can_access_project(project_id)
    );

-- Create policies for user_teams table
DROP POLICY IF EXISTS user_teams_select_policy ON user_teams;
DROP POLICY IF EXISTS user_teams_insert_policy ON user_teams;
DROP POLICY IF EXISTS user_teams_update_policy ON user_teams;
DROP POLICY IF EXISTS user_teams_delete_policy ON user_teams;

CREATE POLICY user_teams_select_policy ON user_teams
    FOR SELECT USING (
        user_id = auth.uid() OR
        team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')
    );

CREATE POLICY user_teams_insert_policy ON user_teams
    FOR INSERT WITH CHECK (
        team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')
    );

CREATE POLICY user_teams_update_policy ON user_teams
    FOR UPDATE USING (
        team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')
    ) WITH CHECK (
        team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')
    );

CREATE POLICY user_teams_delete_policy ON user_teams
    FOR DELETE USING (
        team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')
    );

-- Step 3: Update indices for better RLS policy performance

-- Add index for faster user team membership lookup
CREATE INDEX IF NOT EXISTS idx_user_teams_user_id ON user_teams (user_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_team_id ON user_teams (team_id);
CREATE INDEX IF NOT EXISTS idx_user_teams_role ON user_teams (role);

-- Add index for projects team lookup
CREATE INDEX IF NOT EXISTS idx_projects_team_id ON projects (team_id);

-- Add index for experiments project lookup
CREATE INDEX IF NOT EXISTS idx_experiments_project_id ON experiments (project_id);
CREATE INDEX IF NOT EXISTS idx_experiments_molecule_id ON experiments (molecule_id);

-- Add index for faster molecule access checks
CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);

-- Ensure mixture project relationship is indexed
CREATE INDEX IF NOT EXISTS idx_mixtures_project_id ON mixtures (project_id);

-- Step 4: Set up service role access correctly

-- Create policies for service role access on critical tables
-- These use the special 'service_role' check from auth.enable_service_role()

-- Projects service role policy
DROP POLICY IF EXISTS projects_service_role_policy ON projects;
CREATE POLICY projects_service_role_policy ON projects
    FOR ALL USING (auth.is_service_role());

-- Experiments service role policy
DROP POLICY IF EXISTS experiments_service_role_policy ON experiments;
CREATE POLICY experiments_service_role_policy ON experiments
    FOR ALL USING (auth.is_service_role());

-- Molecules service role policy
DROP POLICY IF EXISTS molecules_service_role_policy ON molecules;
CREATE POLICY molecules_service_role_policy ON molecules
    FOR ALL USING (auth.is_service_role());

-- Molecular properties service role policy
DROP POLICY IF EXISTS molecular_properties_service_role_policy ON molecular_properties;
CREATE POLICY molecular_properties_service_role_policy ON molecular_properties
    FOR ALL USING (auth.is_service_role());

-- Mixtures service role policy
DROP POLICY IF EXISTS mixtures_service_role_policy ON mixtures;
CREATE POLICY mixtures_service_role_policy ON mixtures
    FOR ALL USING (auth.is_service_role());

-- Teams service role policy
DROP POLICY IF EXISTS teams_service_role_policy ON teams;
CREATE POLICY teams_service_role_policy ON teams
    FOR ALL USING (auth.is_service_role());

-- User teams service role policy
DROP POLICY IF EXISTS user_teams_service_role_policy ON user_teams;
CREATE POLICY user_teams_service_role_policy ON user_teams
    FOR ALL USING (auth.is_service_role());

-- Step 5: Verify all tables have RLS enabled

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

-- Create policies for user_profile and experiment_properties if not already existing
DROP POLICY IF EXISTS user_profile_select_policy ON user_profile;
CREATE POLICY user_profile_select_policy ON user_profile
    FOR SELECT USING (auth_user_id = auth.uid() OR auth.is_service_role());

DROP POLICY IF EXISTS user_profile_update_policy ON user_profile;
CREATE POLICY user_profile_update_policy ON user_profile
    FOR UPDATE USING (auth_user_id = auth.uid()) WITH CHECK (auth_user_id = auth.uid());

DROP POLICY IF EXISTS experiment_properties_select_policy ON experiment_properties;
CREATE POLICY experiment_properties_select_policy ON experiment_properties
    FOR SELECT USING (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.id = experiment_id
            AND auth.can_access_project(e.project_id)
        ) OR auth.is_service_role()
    );

DROP POLICY IF EXISTS experiment_properties_insert_policy ON experiment_properties;
CREATE POLICY experiment_properties_insert_policy ON experiment_properties
    FOR INSERT WITH CHECK (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.id = experiment_id
            AND auth.can_access_project(e.project_id)
        ) OR auth.is_service_role()
    );

DROP POLICY IF EXISTS experiment_properties_update_policy ON experiment_properties;
CREATE POLICY experiment_properties_update_policy ON experiment_properties
    FOR UPDATE USING (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.id = experiment_id
            AND auth.can_access_project(e.project_id)
        ) OR auth.is_service_role()
    ) WITH CHECK (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.id = experiment_id
            AND auth.can_access_project(e.project_id)
        ) OR auth.is_service_role()
    );

DROP POLICY IF EXISTS experiment_properties_delete_policy ON experiment_properties;
CREATE POLICY experiment_properties_delete_policy ON experiment_properties
    FOR DELETE USING (
        EXISTS (
            SELECT 1 FROM experiments e
            WHERE e.id = experiment_id
            AND auth.can_access_project(e.project_id)
        ) OR auth.is_service_role()
    );
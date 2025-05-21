-- Migration: 030_team_based_access_control.sql
-- Purpose: Implement proper team-based access control with RLS
-- This migration builds on previous RLS enhancements

-- Step 1: Create team-related helper functions

-- Function to check if a user is a member of a team
CREATE OR REPLACE FUNCTION auth.is_team_member(team_id UUID) RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM public.team_members tm
        WHERE tm.team_id = team_id
        AND tm.user_id = auth.uid()
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.is_team_member IS 'Check if the current user is a member of the specified team';

-- Function to check if a user is a team admin
CREATE OR REPLACE FUNCTION auth.is_team_admin(team_id UUID) RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM public.team_members tm
        WHERE tm.team_id = team_id
        AND tm.user_id = auth.uid()
        AND tm.role = 'admin'
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.is_team_admin IS 'Check if the current user is an admin of the specified team';

-- Function to get all teams a user belongs to
CREATE OR REPLACE FUNCTION auth.get_user_teams() RETURNS TABLE(team_id UUID, role TEXT) AS $$
BEGIN
    RETURN QUERY
    SELECT tm.team_id, tm.role
    FROM public.team_members tm
    WHERE tm.user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.get_user_teams IS 'Get all teams the current user belongs to';

-- Function to check if a molecule belongs to a team the user is a member of
CREATE OR REPLACE FUNCTION auth.can_access_team_molecule(molecule_id UUID) RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM public.molecules m
        JOIN public.team_members tm ON m.team_id = tm.team_id
        WHERE m.id = molecule_id
        AND tm.user_id = auth.uid()
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.can_access_team_molecule IS 'Check if the current user can access a molecule through team membership';

-- Function to check if a mixture belongs to a team the user is a member of
CREATE OR REPLACE FUNCTION auth.can_access_team_mixture(mixture_id UUID) RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM public.mixtures m
        JOIN public.team_members tm ON m.team_id = tm.team_id
        WHERE m.id = mixture_id
        AND tm.user_id = auth.uid()
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.can_access_team_mixture IS 'Check if the current user can access a mixture through team membership';

-- Function to check if a user can modify a team resource
CREATE OR REPLACE FUNCTION auth.can_modify_team_resource(team_id UUID) RETURNS BOOLEAN AS $$
DECLARE
    user_role TEXT;
BEGIN
    -- Get the user's role in the team
    SELECT tm.role INTO user_role
    FROM public.team_members tm
    WHERE tm.team_id = team_id
    AND tm.user_id = auth.uid();
    
    -- Allow if user is a team admin or has 'write' role
    RETURN user_role IN ('admin', 'write');
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.can_modify_team_resource IS 'Check if the current user has permission to modify a team resource';

-- Step 2: Update the existing access control functions to properly incorporate team access

-- Update the molecule access function to use the new team functions
CREATE OR REPLACE FUNCTION auth.has_access_to_molecule(molecule_id UUID) RETURNS BOOLEAN AS $$
DECLARE
    is_public BOOLEAN;
    owner_id UUID;
    team_id UUID;
BEGIN
    -- Get molecule privacy setting, owner and team
    SELECT m.is_public, m.created_by, m.team_id INTO is_public, owner_id, team_id
    FROM public.molecules m
    WHERE m.id = molecule_id;
    
    -- Allow if molecule is public
    IF is_public THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is the owner
    IF owner_id = auth.uid() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is an admin
    IF auth.is_admin() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role_cached() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a member of the team that owns the molecule
    IF team_id IS NOT NULL AND auth.is_team_member(team_id) THEN
        RETURN TRUE;
    END IF;
    
    -- By default, deny access
    RETURN FALSE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Update the mixture access function to use the new team functions
CREATE OR REPLACE FUNCTION auth.has_access_to_mixture(mixture_id UUID) RETURNS BOOLEAN AS $$
DECLARE
    is_public BOOLEAN;
    owner_id UUID;
    project_id UUID;
    team_id UUID;
BEGIN
    -- Get mixture privacy setting, owner, project and team
    SELECT m.is_public, m.created_by, m.project_id, m.team_id INTO is_public, owner_id, project_id, team_id
    FROM public.mixtures m
    WHERE m.id = mixture_id;
    
    -- Allow if mixture is public
    IF is_public THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is the owner
    IF owner_id = auth.uid() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is an admin
    IF auth.is_admin() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role_cached() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a member of the team that owns the mixture
    IF team_id IS NOT NULL AND auth.is_team_member(team_id) THEN
        RETURN TRUE;
    END IF;
    
    -- Check if the user has access to the project
    IF project_id IS NOT NULL AND auth.can_access_project(project_id) THEN
        RETURN TRUE;
    END IF;
    
    -- By default, deny access
    RETURN FALSE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 3: Create RLS policies specifically for team access

-- Update molecules table policies
DROP POLICY IF EXISTS molecules_team_access_policy ON molecules;
CREATE POLICY molecules_team_access_policy ON molecules
    FOR SELECT
    USING (
        team_id IS NOT NULL AND auth.is_team_member(team_id)
    );

DROP POLICY IF EXISTS molecules_team_modify_policy ON molecules;
CREATE POLICY molecules_team_modify_policy ON molecules
    FOR ALL
    USING (
        team_id IS NOT NULL AND auth.can_modify_team_resource(team_id)
    );

-- Update mixtures table policies
DROP POLICY IF EXISTS mixtures_team_access_policy ON mixtures;
CREATE POLICY mixtures_team_access_policy ON mixtures
    FOR SELECT
    USING (
        team_id IS NOT NULL AND auth.is_team_member(team_id)
    );

DROP POLICY IF EXISTS mixtures_team_modify_policy ON mixtures;
CREATE POLICY mixtures_team_modify_policy ON mixtures
    FOR ALL
    USING (
        team_id IS NOT NULL AND auth.can_modify_team_resource(team_id)
    );

-- Step 4: Create policies for team tables

-- Team table
ALTER TABLE IF EXISTS teams ENABLE ROW LEVEL SECURITY;

DROP POLICY IF EXISTS teams_access_policy ON teams;
CREATE POLICY teams_access_policy ON teams
    FOR SELECT
    USING (
        id IN (SELECT team_id FROM auth.get_user_teams())
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

DROP POLICY IF EXISTS teams_modify_policy ON teams;
CREATE POLICY teams_modify_policy ON teams
    FOR ALL
    USING (
        auth.is_team_admin(id)
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

-- Team members table
ALTER TABLE IF EXISTS team_members ENABLE ROW LEVEL SECURITY;

DROP POLICY IF EXISTS team_members_access_policy ON team_members;
CREATE POLICY team_members_access_policy ON team_members
    FOR SELECT
    USING (
        auth.is_team_member(team_id)
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

DROP POLICY IF EXISTS team_members_modify_policy ON team_members;
CREATE POLICY team_members_modify_policy ON team_members
    FOR ALL
    USING (
        auth.is_team_admin(team_id)
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

-- Step 5: Create team invitations table and RLS policies

-- Create team invitations table if it doesn't exist
CREATE TABLE IF NOT EXISTS team_invitations (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    team_id UUID NOT NULL REFERENCES teams(id) ON DELETE CASCADE,
    email TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'read',
    invited_by UUID NOT NULL REFERENCES auth.users(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    expires_at TIMESTAMPTZ NOT NULL DEFAULT NOW() + INTERVAL '7 days',
    token TEXT NOT NULL UNIQUE DEFAULT encode(gen_random_bytes(32), 'hex'),
    CONSTRAINT valid_role CHECK (role IN ('admin', 'write', 'read'))
);

-- Enable RLS on team invitations
ALTER TABLE team_invitations ENABLE ROW LEVEL SECURITY;

-- Team invitations policies
DROP POLICY IF EXISTS team_invitations_access_policy ON team_invitations;
CREATE POLICY team_invitations_access_policy ON team_invitations
    FOR SELECT
    USING (
        auth.is_team_member(team_id)
        OR email = (SELECT email FROM auth.users WHERE id = auth.uid())
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

DROP POLICY IF EXISTS team_invitations_modify_policy ON team_invitations;
CREATE POLICY team_invitations_modify_policy ON team_invitations
    FOR ALL
    USING (
        auth.is_team_admin(team_id)
        OR auth.is_admin()
        OR auth.is_service_role_cached()
    );

-- Step 6: Create team resources view

CREATE OR REPLACE VIEW auth.team_resources AS
SELECT
    'molecule'::TEXT AS resource_type,
    m.id AS resource_id,
    m.name AS resource_name,
    m.team_id,
    t.name AS team_name,
    m.created_by,
    u.email AS owner_email,
    m.created_at,
    m.is_public
FROM
    public.molecules m
JOIN
    public.teams t ON m.team_id = t.id
JOIN
    auth.users u ON m.created_by = u.id
WHERE
    m.team_id IS NOT NULL

UNION ALL

SELECT
    'mixture'::TEXT AS resource_type,
    m.id AS resource_id,
    m.name AS resource_name,
    m.team_id,
    t.name AS team_name,
    m.created_by,
    u.email AS owner_email,
    m.created_at,
    m.is_public
FROM
    public.mixtures m
JOIN
    public.teams t ON m.team_id = t.id
JOIN
    auth.users u ON m.created_by = u.id
WHERE
    m.team_id IS NOT NULL;

COMMENT ON VIEW auth.team_resources IS 'View of all resources (molecules, mixtures) that belong to teams';

-- Step 7: Create functions for team management

-- Function to create a team and add the current user as admin
CREATE OR REPLACE FUNCTION auth.create_team(team_name TEXT, team_description TEXT DEFAULT NULL) 
RETURNS UUID AS $$
DECLARE
    new_team_id UUID;
BEGIN
    -- Insert the new team
    INSERT INTO public.teams (name, description, created_by)
    VALUES (team_name, team_description, auth.uid())
    RETURNING id INTO new_team_id;
    
    -- Add the current user as an admin
    INSERT INTO public.team_members (team_id, user_id, role)
    VALUES (new_team_id, auth.uid(), 'admin');
    
    RETURN new_team_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.create_team IS 'Create a new team and add the current user as admin';

-- Function to invite a user to a team
CREATE OR REPLACE FUNCTION auth.invite_to_team(team_id UUID, email TEXT, role TEXT DEFAULT 'read')
RETURNS UUID AS $$
DECLARE
    invitation_id UUID;
BEGIN
    -- Check if the current user is a team admin
    IF NOT auth.is_team_admin(team_id) AND NOT auth.is_admin() THEN
        RAISE EXCEPTION 'Permission denied: Only team admins can invite users';
    END IF;
    
    -- Insert the invitation
    INSERT INTO public.team_invitations (team_id, email, role, invited_by)
    VALUES (team_id, email, role, auth.uid())
    RETURNING id INTO invitation_id;
    
    RETURN invitation_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.invite_to_team IS 'Invite a user to a team by email';

-- Function to accept a team invitation
CREATE OR REPLACE FUNCTION auth.accept_team_invitation(invitation_token TEXT)
RETURNS BOOLEAN AS $$
DECLARE
    invite_record RECORD;
BEGIN
    -- Get the invitation
    SELECT * INTO invite_record
    FROM public.team_invitations
    WHERE token = invitation_token
    AND expires_at > NOW();
    
    -- Check if invitation exists and is valid
    IF invite_record IS NULL THEN
        RAISE EXCEPTION 'Invalid or expired invitation';
    END IF;
    
    -- Check if the invitation is for the current user
    IF invite_record.email <> (SELECT email FROM auth.users WHERE id = auth.uid()) THEN
        RAISE EXCEPTION 'This invitation is not for you';
    END IF;
    
    -- Add the user to the team
    INSERT INTO public.team_members (team_id, user_id, role)
    VALUES (invite_record.team_id, auth.uid(), invite_record.role)
    ON CONFLICT (team_id, user_id) 
    DO UPDATE SET role = invite_record.role;
    
    -- Delete the invitation
    DELETE FROM public.team_invitations
    WHERE id = invite_record.id;
    
    RETURN TRUE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.accept_team_invitation IS 'Accept a team invitation using the invitation token';

-- Step 8: Create triggers to manage team-based resources

-- Function to check if user can assign a resource to a team
CREATE OR REPLACE FUNCTION auth.check_team_resource_permission()
RETURNS TRIGGER AS $$
BEGIN
    -- If team_id is being set or changed
    IF TG_OP = 'INSERT' AND NEW.team_id IS NOT NULL THEN
        -- Check if user is a member of the team with write permissions
        IF NOT auth.can_modify_team_resource(NEW.team_id) AND NOT auth.is_admin() THEN
            RAISE EXCEPTION 'You cannot assign resources to this team';
        END IF;
    ELSIF TG_OP = 'UPDATE' AND NEW.team_id IS DISTINCT FROM OLD.team_id AND NEW.team_id IS NOT NULL THEN
        -- Check if user is a member of the new team with write permissions
        IF NOT auth.can_modify_team_resource(NEW.team_id) AND NOT auth.is_admin() THEN
            RAISE EXCEPTION 'You cannot assign resources to this team';
        END IF;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Apply the trigger to molecules table
DROP TRIGGER IF EXISTS check_team_molecule_permission ON molecules;
CREATE TRIGGER check_team_molecule_permission
    BEFORE INSERT OR UPDATE ON molecules
    FOR EACH ROW
    EXECUTE FUNCTION auth.check_team_resource_permission();

-- Apply the trigger to mixtures table
DROP TRIGGER IF EXISTS check_team_mixture_permission ON mixtures;
CREATE TRIGGER check_team_mixture_permission
    BEFORE INSERT OR UPDATE ON mixtures
    FOR EACH ROW
    EXECUTE FUNCTION auth.check_team_resource_permission();

-- Step 9: Create a function to transfer resource ownership
CREATE OR REPLACE FUNCTION auth.transfer_resource_to_team(
    resource_type TEXT, 
    resource_id UUID, 
    team_id UUID
) RETURNS BOOLEAN AS $$
DECLARE
    can_transfer BOOLEAN;
BEGIN
    -- Check if the current user can transfer this resource
    IF resource_type = 'molecule' THEN
        SELECT created_by = auth.uid() INTO can_transfer
        FROM public.molecules
        WHERE id = resource_id;
    ELSIF resource_type = 'mixture' THEN
        SELECT created_by = auth.uid() INTO can_transfer
        FROM public.mixtures
        WHERE id = resource_id;
    ELSE
        RAISE EXCEPTION 'Invalid resource type: %', resource_type;
    END IF;
    
    -- Check if the user is a team admin and can receive resources
    IF NOT can_transfer AND NOT auth.is_admin() THEN
        RAISE EXCEPTION 'You do not have permission to transfer this resource';
    END IF;
    
    IF NOT auth.can_modify_team_resource(team_id) AND NOT auth.is_admin() THEN
        RAISE EXCEPTION 'You cannot transfer resources to this team';
    END IF;
    
    -- Perform the transfer
    IF resource_type = 'molecule' THEN
        UPDATE public.molecules
        SET team_id = team_id
        WHERE id = resource_id;
    ELSIF resource_type = 'mixture' THEN
        UPDATE public.mixtures
        SET team_id = team_id
        WHERE id = resource_id;
    END IF;
    
    RETURN TRUE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.transfer_resource_to_team IS 'Transfer ownership of a resource to a team';

-- Step 10: Create indexes to improve performance
CREATE INDEX IF NOT EXISTS idx_team_members_user_id ON team_members(user_id);
CREATE INDEX IF NOT EXISTS idx_team_members_team_id ON team_members(team_id);
CREATE INDEX IF NOT EXISTS idx_team_invitations_email ON team_invitations(email);
CREATE INDEX IF NOT EXISTS idx_team_invitations_token ON team_invitations(token);
CREATE INDEX IF NOT EXISTS idx_molecules_team_id ON molecules(team_id);
CREATE INDEX IF NOT EXISTS idx_mixtures_team_id ON mixtures(team_id);
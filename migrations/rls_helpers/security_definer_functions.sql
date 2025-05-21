-- RLS security definer functions for CryoProtect
-- These functions encapsulate common access patterns for RLS policies
-- They run with the privileges of the function creator (DB owner)
-- This improves performance of RLS policy evaluation significantly

-- 1. User Team Membership Function
CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() AND user_profile.team_id = $1
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 2. User Teams Function
CREATE OR REPLACE FUNCTION public.user_teams()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT team_id FROM user_profile
  WHERE auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 3. User Projects Function
CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT tp.project_id 
  FROM team_projects tp
  JOIN user_profile up ON tp.team_id = up.team_id
  WHERE up.auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 4. Molecule Access Function
CREATE OR REPLACE FUNCTION public.has_molecule_access(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecules m
    WHERE m.id = molecule_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = molecule_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 5. Mixture Access Function
CREATE OR REPLACE FUNCTION public.has_mixture_access(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixtures m
    WHERE m.id = mixture_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_mixtures pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.mixture_id = mixture_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 6. Clearance Level Function
CREATE OR REPLACE FUNCTION public.user_has_clearance(required_level text)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() 
    AND (clearance_level >= required_level OR required_level IS NULL)
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 7. Experiment Access Function
CREATE OR REPLACE FUNCTION public.has_experiment_access(experiment_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM experiments e
    WHERE e.id = experiment_id AND (
      e.is_public = true OR
      e.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_experiments pe
        JOIN team_projects tp ON pe.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pe.experiment_id = experiment_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 8. Project Access Function
CREATE OR REPLACE FUNCTION public.has_project_access(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM team_projects tp
    JOIN user_profile up ON tp.team_id = up.team_id
    WHERE tp.project_id = project_id AND up.auth_user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 9. Property Access Function
CREATE OR REPLACE FUNCTION public.has_property_access(property_id uuid)
RETURNS boolean AS $$
DECLARE
  molecule_id uuid;
BEGIN
  -- First, get the molecule ID for this property
  SELECT mp.molecule_id INTO molecule_id
  FROM molecular_properties mp
  WHERE mp.id = property_id;
  
  -- If molecule ID was found, check access to that molecule
  IF molecule_id IS NOT NULL THEN
    RETURN public.has_molecule_access(molecule_id);
  END IF;
  
  -- No molecule ID means no access
  RETURN FALSE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 10. Is Admin Function
CREATE OR REPLACE FUNCTION public.is_admin()
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() AND role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 11. Is Project Admin Function
CREATE OR REPLACE FUNCTION public.is_project_admin(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM project_members pm
    JOIN user_profile up ON pm.user_id = up.id
    WHERE pm.project_id = project_id 
    AND up.auth_user_id = auth.uid() 
    AND pm.role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 12. Is Team Admin Function
CREATE OR REPLACE FUNCTION public.is_team_admin(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM team_members tm
    JOIN user_profile up ON tm.user_id = up.id
    WHERE tm.team_id = team_id 
    AND up.auth_user_id = auth.uid() 
    AND tm.role = 'admin'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create a function to get the current user's profile ID
CREATE OR REPLACE FUNCTION public.current_user_profile_id()
RETURNS uuid AS $$
DECLARE
  profile_id uuid;
BEGIN
  SELECT id INTO profile_id FROM user_profile WHERE auth_user_id = auth.uid();
  RETURN profile_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create a function to check if a resource belongs to the current user
CREATE OR REPLACE FUNCTION public.is_owner(resource_table text, resource_id uuid)
RETURNS boolean AS $$
DECLARE
  owner_id text;
BEGIN
  EXECUTE format('SELECT created_by FROM %I WHERE id = $1', resource_table)
  INTO owner_id
  USING resource_id;
  
  RETURN owner_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
-- Migration: 028_enhance_rls_policies.sql
-- Purpose: Enhance RLS policies for better access control
-- This migration consolidates and improves the Row Level Security policies

-- Step 1: Create helper functions for RLS policy checks
CREATE OR REPLACE FUNCTION auth.has_access_to_molecule(molecule_id UUID) RETURNS BOOLEAN AS $$
DECLARE
    is_public BOOLEAN;
    owner_id UUID;
BEGIN
    -- Get molecule privacy setting and owner
    SELECT m.is_public, m.created_by INTO is_public, owner_id
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
    IF EXISTS (
        SELECT 1 FROM auth.users
        WHERE id = auth.uid() AND role = 'admin'
    ) THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role() THEN
        RETURN TRUE;
    END IF;
    
    -- Check if the user is part of a team that can access this molecule
    IF EXISTS (
        SELECT 1
        FROM public.team_members tm
        JOIN public.molecules m ON tm.team_id = m.team_id
        WHERE tm.user_id = auth.uid() AND m.id = molecule_id
    ) THEN
        RETURN TRUE;
    END IF;
    
    -- By default, deny access
    RETURN FALSE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION auth.has_access_to_mixture(mixture_id UUID) RETURNS BOOLEAN AS $$
DECLARE
    is_public BOOLEAN;
    owner_id UUID;
    project_id UUID;
BEGIN
    -- Get mixture privacy setting, owner and project
    SELECT m.is_public, m.created_by, m.project_id INTO is_public, owner_id, project_id
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
    IF EXISTS (
        SELECT 1 FROM auth.users
        WHERE id = auth.uid() AND role = 'admin'
    ) THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role() THEN
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

-- Step 2: Create a function to get user's clearance level
CREATE OR REPLACE FUNCTION auth.get_user_clearance_level() RETURNS TEXT AS $$
DECLARE
    user_clearance TEXT;
BEGIN
    -- Get the user's clearance level from user_profile
    SELECT up.clearance_level INTO user_clearance
    FROM public.user_profile up
    WHERE up.auth_user_id = auth.uid();
    
    -- If not found, assume 'low' clearance
    IF user_clearance IS NULL THEN
        RETURN 'low';
    END IF;
    
    -- Return the clearance level
    RETURN user_clearance;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Function to convert clearance level to numeric value for comparison
CREATE OR REPLACE FUNCTION auth.clearance_level_value(level TEXT) RETURNS INTEGER AS $$
BEGIN
    RETURN CASE
        WHEN level = 'admin' THEN 1000
        WHEN level = 'high' THEN 100
        WHEN level = 'medium' THEN 50
        WHEN level = 'low' THEN 10
        WHEN level IS NULL THEN 0
        ELSE 0
    END;
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Step 3: Create a function to check if user has required clearance
CREATE OR REPLACE FUNCTION auth.has_clearance(required_level TEXT) RETURNS BOOLEAN AS $$
DECLARE
    user_level TEXT;
    user_value INTEGER;
    required_value INTEGER;
BEGIN
    -- Get the user's clearance level
    user_level := auth.get_user_clearance_level();
    
    -- Convert clearance levels to numeric values
    user_value := auth.clearance_level_value(user_level);
    required_value := auth.clearance_level_value(required_level);
    
    -- Compare values
    RETURN user_value >= required_value;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 4: Drop existing duplicate policies
DO $$
DECLARE
    policy_record RECORD;
BEGIN
    -- Find and drop duplicate service_role policies
    FOR policy_record IN 
        SELECT policyname, tablename
        FROM pg_policies
        WHERE schemaname = 'public'
          AND policyname LIKE '%service_role%'
          AND tablename IN (
              SELECT tablename 
              FROM pg_policies 
              WHERE schemaname = 'public' 
                AND policyname LIKE '%service_role%'
              GROUP BY tablename 
              HAVING COUNT(*) > 1
          )
    LOOP
        EXECUTE format('DROP POLICY IF EXISTS %I ON %I', 
                      policy_record.policyname, 
                      policy_record.tablename);
    END LOOP;
END $$;

-- Step 5: Create consolidated RLS policies for core tables

-- Molecules table
DROP POLICY IF EXISTS molecules_access_policy ON molecules;
CREATE POLICY molecules_access_policy ON molecules
    FOR SELECT
    USING (
        is_public
        OR created_by = auth.uid()
        OR auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
        OR EXISTS (
            SELECT 1
            FROM team_members tm
            WHERE tm.user_id = auth.uid() AND tm.team_id = molecules.team_id
        )
    );

DROP POLICY IF EXISTS molecules_modify_policy ON molecules;
CREATE POLICY molecules_modify_policy ON molecules
    FOR ALL
    USING (
        created_by = auth.uid()
        OR auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
    );

-- Mixtures table
DROP POLICY IF EXISTS mixtures_access_policy ON mixtures;
CREATE POLICY mixtures_access_policy ON mixtures
    FOR SELECT
    USING (
        is_public
        OR created_by = auth.uid()
        OR auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
        OR (project_id IS NOT NULL AND auth.can_access_project(project_id))
    );

DROP POLICY IF EXISTS mixtures_modify_policy ON mixtures;
CREATE POLICY mixtures_modify_policy ON mixtures
    FOR ALL
    USING (
        created_by = auth.uid()
        OR auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
        OR (project_id IS NOT NULL AND auth.can_manage_project(project_id))
    );

-- Molecular properties table
DROP POLICY IF EXISTS molecular_properties_access_policy ON molecular_properties;
CREATE POLICY molecular_properties_access_policy ON molecular_properties
    FOR SELECT
    USING (
        auth.has_access_to_molecule(molecule_id)
        AND (
            sensitivity_level IS NULL
            OR auth.clearance_level_value(sensitivity_level) <= auth.clearance_level_value(auth.get_user_clearance_level())
            OR auth.is_service_role()
        )
    );

DROP POLICY IF EXISTS molecular_properties_modify_policy ON molecular_properties;
CREATE POLICY molecular_properties_modify_policy ON molecular_properties
    FOR ALL
    USING (
        auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM molecules m
            WHERE m.id = molecular_properties.molecule_id
            AND m.created_by = auth.uid()
        )
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
    );

-- Mixture components table
DROP POLICY IF EXISTS mixture_components_access_policy ON mixture_components;
CREATE POLICY mixture_components_access_policy ON mixture_components
    FOR SELECT
    USING (
        auth.has_access_to_mixture(mixture_id)
        AND auth.has_access_to_molecule(molecule_id)
    );

DROP POLICY IF EXISTS mixture_components_modify_policy ON mixture_components;
CREATE POLICY mixture_components_modify_policy ON mixture_components
    FOR ALL
    USING (
        auth.is_service_role()
        OR EXISTS (
            SELECT 1 FROM mixtures m
            WHERE m.id = mixture_components.mixture_id
            AND m.created_by = auth.uid()
        )
        OR EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = auth.uid() AND role = 'admin'
        )
    );

-- Step 6: Create an RLS test function to validate policies
CREATE OR REPLACE FUNCTION auth.test_rls_access(
    entity_type TEXT,
    entity_id UUID,
    test_user_id UUID DEFAULT NULL
) RETURNS TABLE (
    access_type TEXT,
    has_access BOOLEAN,
    reason TEXT
) AS $$
DECLARE
    current_user_id UUID;
    was_role_changed BOOLEAN := FALSE;
BEGIN
    -- Save current user ID
    current_user_id := auth.uid();
    
    -- Set the user if provided
    IF test_user_id IS NOT NULL AND test_user_id != current_user_id THEN
        -- This requires the calling user to have admin privileges
        IF NOT EXISTS (
            SELECT 1 FROM auth.users
            WHERE id = current_user_id AND role = 'admin'
        ) AND NOT auth.is_service_role() THEN
            RETURN QUERY
            SELECT 'error'::TEXT, FALSE, 'Permission denied: Only admins can test as other users'::TEXT;
            RETURN;
        END IF;
        
        -- Temporarily set role to the test user
        PERFORM set_config('request.jwt.claims', json_build_object('sub', test_user_id)::TEXT, TRUE);
        was_role_changed := TRUE;
    END IF;
    
    -- Test access based on entity type
    CASE entity_type
        WHEN 'molecule' THEN
            -- Test SELECT access
            RETURN QUERY
            SELECT 
                'select'::TEXT,
                auth.has_access_to_molecule(entity_id),
                CASE 
                    WHEN auth.has_access_to_molecule(entity_id) THEN 'Access granted'
                    ELSE 'Access denied'
                END;
                
            -- Test UPDATE access
            RETURN QUERY
            SELECT 
                'update'::TEXT,
                CASE 
                    WHEN EXISTS (
                        SELECT 1 FROM molecules m
                        WHERE m.id = entity_id AND m.created_by = auth.uid()
                    ) OR auth.is_service_role() OR EXISTS (
                        SELECT 1 FROM auth.users
                        WHERE id = auth.uid() AND role = 'admin'
                    ) THEN TRUE
                    ELSE FALSE
                END,
                CASE 
                    WHEN EXISTS (
                        SELECT 1 FROM molecules m
                        WHERE m.id = entity_id AND m.created_by = auth.uid()
                    ) THEN 'User is owner'
                    WHEN auth.is_service_role() THEN 'User is service role'
                    WHEN EXISTS (
                        SELECT 1 FROM auth.users
                        WHERE id = auth.uid() AND role = 'admin'
                    ) THEN 'User is admin'
                    ELSE 'User does not have permission to update'
                END;
                
        WHEN 'mixture' THEN
            -- Test SELECT access
            RETURN QUERY
            SELECT 
                'select'::TEXT,
                auth.has_access_to_mixture(entity_id),
                CASE 
                    WHEN auth.has_access_to_mixture(entity_id) THEN 'Access granted'
                    ELSE 'Access denied'
                END;
                
            -- Test UPDATE access
            RETURN QUERY
            SELECT 
                'update'::TEXT,
                CASE 
                    WHEN EXISTS (
                        SELECT 1 FROM mixtures m
                        WHERE m.id = entity_id AND m.created_by = auth.uid()
                    ) OR auth.is_service_role() OR EXISTS (
                        SELECT 1 FROM auth.users
                        WHERE id = auth.uid() AND role = 'admin'
                    ) THEN TRUE
                    ELSE FALSE
                END,
                CASE 
                    WHEN EXISTS (
                        SELECT 1 FROM mixtures m
                        WHERE m.id = entity_id AND m.created_by = auth.uid()
                    ) THEN 'User is owner'
                    WHEN auth.is_service_role() THEN 'User is service role'
                    WHEN EXISTS (
                        SELECT 1 FROM auth.users
                        WHERE id = auth.uid() AND role = 'admin'
                    ) THEN 'User is admin'
                    ELSE 'User does not have permission to update'
                END;
                
        WHEN 'property' THEN
            -- Test SELECT access (assumes entity_id is a molecular_properties.id)
            RETURN QUERY
            WITH prop AS (
                SELECT molecule_id, sensitivity_level FROM molecular_properties WHERE id = entity_id
            )
            SELECT 
                'select'::TEXT,
                CASE 
                    WHEN EXISTS (
                        SELECT 1 FROM prop p
                        WHERE auth.has_access_to_molecule(p.molecule_id)
                        AND (
                            p.sensitivity_level IS NULL
                            OR auth.clearance_level_value(p.sensitivity_level) <= auth.clearance_level_value(auth.get_user_clearance_level())
                            OR auth.is_service_role()
                        )
                    ) THEN TRUE
                    ELSE FALSE
                END,
                CASE 
                    WHEN NOT EXISTS (SELECT 1 FROM prop) THEN 'Property not found'
                    WHEN EXISTS (
                        SELECT 1 FROM prop p
                        WHERE auth.has_access_to_molecule(p.molecule_id)
                        AND (
                            p.sensitivity_level IS NULL
                            OR auth.clearance_level_value(p.sensitivity_level) <= auth.clearance_level_value(auth.get_user_clearance_level())
                            OR auth.is_service_role()
                        )
                    ) THEN 'Access granted'
                    ELSE 'Access denied due to insufficient clearance or molecule access'
                END;
                
        ELSE
            -- Unknown entity type
            RETURN QUERY
            SELECT 'error'::TEXT, FALSE, 'Unknown entity type'::TEXT;
    END CASE;
    
    -- Restore original user if changed
    IF was_role_changed THEN
        PERFORM set_config('request.jwt.claims', json_build_object('sub', current_user_id)::TEXT, TRUE);
    END IF;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 7: Add comments to the functions for better documentation
COMMENT ON FUNCTION auth.has_access_to_molecule IS 'Checks if the current user has access to the specified molecule';
COMMENT ON FUNCTION auth.has_access_to_mixture IS 'Checks if the current user has access to the specified mixture';
COMMENT ON FUNCTION auth.get_user_clearance_level IS 'Gets the clearance level of the current user';
COMMENT ON FUNCTION auth.clearance_level_value IS 'Converts a clearance level string to a numeric value for comparison';
COMMENT ON FUNCTION auth.has_clearance IS 'Checks if the current user has the required clearance level';
COMMENT ON FUNCTION auth.test_rls_access IS 'Tests RLS access for a given entity and user';
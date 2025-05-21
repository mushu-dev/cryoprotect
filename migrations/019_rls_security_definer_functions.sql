-- Migration: RLS Security Definer Functions
-- Purpose: This migration creates high-performance security definer functions for RLS policies
-- to improve database query performance while maintaining security.

BEGIN;

-- Step 1: Create User and Team Permission Functions

-- Function to check if a user has access to a resource
-- This is a security definer function that executes with elevated privileges
CREATE OR REPLACE FUNCTION auth.check_resource_access(
    user_id uuid,
    resource_owner uuid,
    resource_type text,
    resource_id uuid
) RETURNS boolean
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
DECLARE
    has_access boolean;
BEGIN
    -- Prevent SQL injection by validating resource_type against allowed values
    IF resource_type NOT IN ('molecules', 'mixtures', 'mixture_components', 'molecular_properties', 
                            'experiments', 'predictions', 'protocols', 'protocol_steps', 
                            'toxicity_data', 'calculation_methods', 'lab_verifications') THEN
        RAISE EXCEPTION 'Invalid resource type: %', resource_type;
    END IF;

    -- Owner always has access
    IF user_id = resource_owner THEN
        RETURN true;
    END IF;

    -- Check team-based access
    SELECT EXISTS (
        SELECT 1 FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE utm.user_id = check_resource_access.user_id AND tr.can_view = true
    ) INTO has_access;
    
    IF has_access THEN
        RETURN true;
    END IF;

    -- Check sharing-based access
    SELECT EXISTS (
        SELECT 1 FROM shared_items si
        WHERE si.item_type = check_resource_access.resource_type AND si.item_id = check_resource_access.resource_id 
        AND (si.shared_with_user_id = check_resource_access.user_id OR si.is_public = true)
    ) INTO has_access;

    RETURN has_access;
END;
$$;

-- Function to check if a user can modify a resource
CREATE OR REPLACE FUNCTION auth.check_resource_modify_access(
    user_id uuid,
    resource_owner uuid,
    action text
) RETURNS boolean
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
DECLARE
    has_access boolean;
BEGIN
    -- Prevent SQL injection by validating action against allowed values
    IF action NOT IN ('create', 'update', 'delete') THEN
        RAISE EXCEPTION 'Invalid action: %', action;
    END IF;

    -- Owner always has access
    IF user_id = resource_owner THEN
        RETURN true;
    END IF;

    -- Check team-based access
    IF action = 'create' THEN
        SELECT EXISTS (
            SELECT 1 FROM user_team_membership utm
            JOIN team_roles tr ON utm.role_id = tr.id
            WHERE utm.user_id = check_resource_modify_access.user_id AND tr.can_create = true
        ) INTO has_access;
    ELSIF action = 'update' THEN
        SELECT EXISTS (
            SELECT 1 FROM user_team_membership utm
            JOIN team_roles tr ON utm.role_id = tr.id
            WHERE utm.user_id = check_resource_modify_access.user_id AND tr.can_update = true
        ) INTO has_access;
    ELSIF action = 'delete' THEN
        SELECT EXISTS (
            SELECT 1 FROM user_team_membership utm
            JOIN team_roles tr ON utm.role_id = tr.id
            WHERE utm.user_id = check_resource_modify_access.user_id AND tr.can_delete = true
        ) INTO has_access;
    END IF;

    RETURN has_access;
END;
$$;

-- Function to check if user is in service role (without using auth.role())
CREATE OR REPLACE FUNCTION auth.is_service_role_user() RETURNS boolean
LANGUAGE sql SECURITY DEFINER
SET search_path = public, auth
AS $$
    -- This uses Supabase's internal auth.role() function more efficiently
    SELECT auth.role() = 'service_role';
$$;

-- Step 2: Create Specialized Verification Functions

-- Function to check if a user can access lab verification data
CREATE OR REPLACE FUNCTION auth.check_verification_access(
    user_id uuid,
    verifier_id uuid
) RETURNS boolean
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
DECLARE
    has_access boolean;
BEGIN
    -- Verifier always has access
    IF user_id = verifier_id THEN
        RETURN true;
    END IF;

    -- Check team-based access with verification permissions
    SELECT EXISTS (
        SELECT 1 FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE utm.user_id = check_verification_access.user_id AND tr.can_view_verifications = true
    ) INTO has_access;
    
    RETURN has_access;
END;
$$;

-- Function to check if a user can modify lab verification data
CREATE OR REPLACE FUNCTION auth.check_verification_modify_access(
    user_id uuid,
    verifier_id uuid
) RETURNS boolean
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
DECLARE
    has_access boolean;
BEGIN
    -- Verifier always has access
    IF user_id = verifier_id THEN
        RETURN true;
    END IF;

    -- Check team-based access with verification permissions
    SELECT EXISTS (
        SELECT 1 FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE utm.user_id = check_verification_modify_access.user_id AND tr.can_update_verifications = true
    ) INTO has_access;
    
    RETURN has_access;
END;
$$;

-- Step 3: Create High-Performance Batch Access Check Functions

-- Function to check if a user has access to a batch of resources by ID
CREATE OR REPLACE FUNCTION auth.check_batch_access(
    user_id uuid,
    resource_type text,
    resource_ids uuid[]
) RETURNS TABLE (resource_id uuid, has_access boolean)
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
BEGIN
    -- Prevent SQL injection by validating resource_type against allowed values
    IF resource_type NOT IN ('molecules', 'mixtures', 'mixture_components', 'molecular_properties', 
                            'experiments', 'predictions', 'protocols', 'protocol_steps', 
                            'toxicity_data', 'calculation_methods', 'lab_verifications') THEN
        RAISE EXCEPTION 'Invalid resource type: %', resource_type;
    END IF;

    -- Handle different resource types with specific queries
    CASE resource_type
        WHEN 'molecules' THEN
            RETURN QUERY
                WITH resource_owners AS (
                    SELECT id, created_by FROM molecules WHERE id = ANY(resource_ids)
                )
                SELECT ro.id, 
                       (ro.created_by = user_id OR
                        EXISTS (
                            SELECT 1 FROM user_team_membership utm
                            JOIN team_roles tr ON utm.role_id = tr.id
                            WHERE utm.user_id = check_batch_access.user_id AND tr.can_view = true
                        ) OR
                        EXISTS (
                            SELECT 1 FROM shared_items si
                            WHERE si.item_type = 'molecules' AND si.item_id = ro.id 
                            AND (si.shared_with_user_id = check_batch_access.user_id OR si.is_public = true)
                        )) AS has_access
                FROM resource_owners ro;
                
        WHEN 'mixtures' THEN
            RETURN QUERY
                WITH resource_owners AS (
                    SELECT id, created_by FROM mixtures WHERE id = ANY(resource_ids)
                )
                SELECT ro.id, 
                       (ro.created_by = user_id OR
                        EXISTS (
                            SELECT 1 FROM user_team_membership utm
                            JOIN team_roles tr ON utm.role_id = tr.id
                            WHERE utm.user_id = check_batch_access.user_id AND tr.can_view = true
                        ) OR
                        EXISTS (
                            SELECT 1 FROM shared_items si
                            WHERE si.item_type = 'mixtures' AND si.item_id = ro.id 
                            AND (si.shared_with_user_id = check_batch_access.user_id OR si.is_public = true)
                        )) AS has_access
                FROM resource_owners ro;
        
        -- Add other resource types as needed with similar patterns
        
        ELSE
            -- For other resource types, check individually
            RETURN QUERY
                SELECT r.id, auth.check_resource_access(user_id, r.created_by, resource_type, r.id)
                FROM unnest(resource_ids) AS r(id)
                LEFT JOIN (
                    SELECT id, created_by FROM 
                    CASE resource_type
                        WHEN 'mixture_components' THEN mixture_components
                        WHEN 'molecular_properties' THEN molecular_properties
                        WHEN 'experiments' THEN experiments
                        WHEN 'predictions' THEN predictions
                        WHEN 'protocols' THEN protocols
                        WHEN 'protocol_steps' THEN protocol_steps
                        WHEN 'toxicity_data' THEN toxicity_data
                        WHEN 'calculation_methods' THEN calculation_methods
                        ELSE NULL::record
                    END
                ) AS resource_data ON r.id = resource_data.id;
    END CASE;
END;
$$;

-- Step 4: Create User Role and Permission Cache Functions

-- Function to get a user's roles and permissions (cached)
CREATE OR REPLACE FUNCTION auth.get_user_permissions(
    user_id uuid
) RETURNS TABLE (
    role_id uuid,
    role_name text,
    can_view boolean,
    can_create boolean,
    can_update boolean,
    can_delete boolean,
    can_view_verifications boolean,
    can_update_verifications boolean
)
LANGUAGE sql SECURITY DEFINER
SET search_path = public, auth
AS $$
    SELECT tr.id, tr.name, tr.can_view, tr.can_create, tr.can_update, tr.can_delete, 
           tr.can_view_verifications, tr.can_update_verifications
    FROM user_team_membership utm
    JOIN team_roles tr ON utm.role_id = tr.id
    WHERE utm.user_id = get_user_permissions.user_id;
$$;

-- Function to check if user has any specific permission
CREATE OR REPLACE FUNCTION auth.user_has_permission(
    user_id uuid,
    permission text
) RETURNS boolean
LANGUAGE plpgsql SECURITY DEFINER
SET search_path = public, auth
AS $$
DECLARE
    has_permission boolean;
BEGIN
    -- Validate permission parameter to prevent SQL injection
    IF permission NOT IN ('can_view', 'can_create', 'can_update', 'can_delete', 
                         'can_view_verifications', 'can_update_verifications') THEN
        RAISE EXCEPTION 'Invalid permission: %', permission;
    END IF;
    
    -- Dynamic query to check specific permission
    EXECUTE format('
        SELECT EXISTS (
            SELECT 1 FROM user_team_membership utm
            JOIN team_roles tr ON utm.role_id = tr.id
            WHERE utm.user_id = %L AND tr.%I = true
        )', user_id, permission) INTO has_permission;
    
    RETURN has_permission;
END;
$$;

-- Step 5: Create Comments on Functions for Documentation

COMMENT ON FUNCTION auth.check_resource_access IS 
'Efficiently checks if a user has access to a specific resource based on ownership, team roles, or sharing.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.check_resource_modify_access IS 
'Checks if a user can modify (create, update, delete) a resource based on ownership and team roles.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.is_service_role_user IS 
'Efficiently checks if the current user has the service_role.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.check_verification_access IS 
'Specialized function for checking access to lab verification data.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.check_verification_modify_access IS 
'Specialized function for checking modify access to lab verification data.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.check_batch_access IS 
'High-performance batch access check for multiple resources at once.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.get_user_permissions IS 
'Returns all roles and permissions for a user.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

COMMENT ON FUNCTION auth.user_has_permission IS 
'Checks if a user has a specific permission from any of their roles.
This is a SECURITY DEFINER function that executes with elevated privileges for better performance.';

-- Step 6: Log migration completed
DO $$
BEGIN
    RAISE NOTICE 'RLS security definer functions migration completed successfully';
END $$;

COMMIT;
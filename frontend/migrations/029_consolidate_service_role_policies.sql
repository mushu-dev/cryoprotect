-- Migration: 029_consolidate_service_role_policies.sql
-- Purpose: Consolidate duplicative RLS policies, with a focus on service_role policies
-- This migration builds on the previous RLS enhancements

-- Step 1: Create a function to efficiently check if a user is a service role
-- This improves performance by caching the result within a transaction
CREATE OR REPLACE FUNCTION auth.is_service_role_cached() RETURNS BOOLEAN AS $$
DECLARE
    is_service BOOLEAN;
    setting_name TEXT := 'service_role_check_' || auth.uid()::TEXT;
BEGIN
    -- Check if we've already calculated this in the current transaction
    is_service := current_setting(setting_name, TRUE)::BOOLEAN;
    
    -- If not found, calculate and store it
    IF is_service IS NULL THEN
        is_service := auth.is_service_role();
        PERFORM set_config(setting_name, is_service::TEXT, TRUE);
    END IF;
    
    RETURN is_service;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.is_service_role_cached IS 'Performance-optimized version of is_service_role() that caches the result within a transaction';

-- Step 2: Create a function to efficiently get the user's role
-- This improves performance by caching the result within a transaction
CREATE OR REPLACE FUNCTION auth.get_user_role_cached() RETURNS TEXT AS $$
DECLARE
    user_role TEXT;
    setting_name TEXT := 'user_role_check_' || auth.uid()::TEXT;
BEGIN
    -- Check if we've already calculated this in the current transaction
    BEGIN
        user_role := current_setting(setting_name, TRUE);
    EXCEPTION WHEN OTHERS THEN
        user_role := NULL;
    END;
    
    -- If not found, calculate and store it
    IF user_role IS NULL THEN
        SELECT role INTO user_role FROM auth.users WHERE id = auth.uid();
        PERFORM set_config(setting_name, COALESCE(user_role, 'none'), TRUE);
    END IF;
    
    RETURN user_role;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.get_user_role_cached IS 'Performance-optimized version of getting user role that caches the result within a transaction';

-- Step 3: Create a function to check if a user is an admin
CREATE OR REPLACE FUNCTION auth.is_admin() RETURNS BOOLEAN AS $$
BEGIN
    RETURN auth.get_user_role_cached() = 'admin';
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

COMMENT ON FUNCTION auth.is_admin IS 'Check if the current user has the admin role';

-- Step 4: Update the access check functions to use the cached versions
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
    IF auth.is_admin() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role_cached() THEN
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
    IF auth.is_admin() THEN
        RETURN TRUE;
    END IF;
    
    -- Allow if user is a service role
    IF auth.is_service_role_cached() THEN
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

-- Step 5: Function to get all tables with RLS enabled
CREATE OR REPLACE FUNCTION auth.get_tables_with_rls() RETURNS TABLE(
    schema_name TEXT,
    table_name TEXT,
    has_rls BOOLEAN
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        nsp.nspname::TEXT AS schema_name,
        cls.relname::TEXT AS table_name,
        cls.relrowsecurity AS has_rls
    FROM
        pg_class cls
    JOIN
        pg_namespace nsp ON cls.relnamespace = nsp.oid
    WHERE
        cls.relkind = 'r'  -- Only regular tables
        AND nsp.nspname = 'public'  -- Only public schema
    ORDER BY
        nsp.nspname, cls.relname;
END;
$$ LANGUAGE plpgsql;

COMMENT ON FUNCTION auth.get_tables_with_rls IS 'Get all tables in the public schema with their RLS status';

-- Step 6: Drop duplicative policies from all tables
DO $$
DECLARE
    policy_record RECORD;
    policy_count INT;
    table_record RECORD;
BEGIN
    -- First, find all tables with RLS enabled
    FOR table_record IN 
        SELECT table_name FROM auth.get_tables_with_rls() WHERE has_rls = TRUE
    LOOP
        -- Check for service_role policies
        SELECT COUNT(*) INTO policy_count
        FROM pg_policies
        WHERE schemaname = 'public'
          AND tablename = table_record.table_name
          AND policyname LIKE '%service_role%';
        
        -- If there are multiple service_role policies, drop them
        IF policy_count > 1 THEN
            RAISE NOTICE 'Table % has % service_role policies - consolidating...', 
                table_record.table_name, policy_count;
            
            -- Drop each service_role policy
            FOR policy_record IN 
                SELECT policyname
                FROM pg_policies
                WHERE schemaname = 'public'
                  AND tablename = table_record.table_name
                  AND policyname LIKE '%service_role%'
            LOOP
                EXECUTE format('DROP POLICY IF EXISTS %I ON %I', 
                              policy_record.policyname, 
                              table_record.table_name);
                              
                RAISE NOTICE 'Dropped policy % on table %', 
                    policy_record.policyname, table_record.table_name;
            END LOOP;
            
            -- Create a single consolidated service_role policy
            EXECUTE format(
                'CREATE POLICY %I ON %I FOR ALL USING (auth.is_service_role_cached())',
                'service_role_access_policy',
                table_record.table_name
            );
            
            RAISE NOTICE 'Created consolidated service_role policy on table %', 
                table_record.table_name;
        END IF;
        
        -- Also check for duplicate admin policies
        SELECT COUNT(*) INTO policy_count
        FROM pg_policies
        WHERE schemaname = 'public'
          AND tablename = table_record.table_name
          AND qual LIKE '%role = ''admin''%';
        
        -- If there are multiple admin-related policies, consolidate them
        IF policy_count > 1 THEN
            RAISE NOTICE 'Table % has % admin-related policies - consolidating...', 
                table_record.table_name, policy_count;
            
            -- Drop each admin-related policy
            FOR policy_record IN 
                SELECT policyname
                FROM pg_policies
                WHERE schemaname = 'public'
                  AND tablename = table_record.table_name
                  AND qual LIKE '%role = ''admin''%'
            LOOP
                EXECUTE format('DROP POLICY IF EXISTS %I ON %I', 
                              policy_record.policyname, 
                              table_record.table_name);
                              
                RAISE NOTICE 'Dropped policy % on table %', 
                    policy_record.policyname, table_record.table_name;
            END LOOP;
            
            -- Create a single consolidated admin policy
            EXECUTE format(
                'CREATE POLICY %I ON %I FOR ALL USING (auth.is_admin())',
                'admin_access_policy',
                table_record.table_name
            );
            
            RAISE NOTICE 'Created consolidated admin policy on table %', 
                table_record.table_name;
        END IF;
    END LOOP;
END $$;

-- Step 7: Identify tables without service_role or admin policies and add them
DO $$
DECLARE
    table_record RECORD;
    has_service_role_policy BOOLEAN;
    has_admin_policy BOOLEAN;
BEGIN
    -- Find all tables with RLS enabled
    FOR table_record IN 
        SELECT table_name FROM auth.get_tables_with_rls() WHERE has_rls = TRUE
    LOOP
        -- Check if the table has a service_role policy
        SELECT COUNT(*) > 0 INTO has_service_role_policy
        FROM pg_policies
        WHERE schemaname = 'public'
          AND tablename = table_record.table_name
          AND (
              policyname LIKE '%service_role%' 
              OR qual LIKE '%service_role%'
          );
        
        -- Check if the table has an admin policy
        SELECT COUNT(*) > 0 INTO has_admin_policy
        FROM pg_policies
        WHERE schemaname = 'public'
          AND tablename = table_record.table_name
          AND (
              policyname LIKE '%admin%' 
              OR qual LIKE '%role = ''admin''%'
          );
        
        -- If no service_role policy, add one
        IF NOT has_service_role_policy THEN
            EXECUTE format(
                'CREATE POLICY %I ON %I FOR ALL USING (auth.is_service_role_cached())',
                'service_role_access_policy',
                table_record.table_name
            );
            
            RAISE NOTICE 'Added missing service_role policy to table %', 
                table_record.table_name;
        END IF;
        
        -- If no admin policy, add one
        IF NOT has_admin_policy THEN
            EXECUTE format(
                'CREATE POLICY %I ON %I FOR ALL USING (auth.is_admin())',
                'admin_access_policy',
                table_record.table_name
            );
            
            RAISE NOTICE 'Added missing admin policy to table %', 
                table_record.table_name;
        END IF;
    END LOOP;
END $$;

-- Step 8: Update existing RLS policies to use cached functions
DO $$
DECLARE
    policy_record RECORD;
    updated_qual TEXT;
BEGIN
    -- Find policies with service_role checks
    FOR policy_record IN 
        SELECT 
            schemaname, 
            tablename, 
            policyname, 
            qual
        FROM 
            pg_policies
        WHERE 
            schemaname = 'public'
            AND qual LIKE '%auth.is_service_role()%'
    LOOP
        -- Replace is_service_role() with is_service_role_cached()
        updated_qual := REPLACE(
            policy_record.qual, 
            'auth.is_service_role()', 
            'auth.is_service_role_cached()'
        );
        
        -- Update the policy
        EXECUTE format(
            'ALTER POLICY %I ON %I USING (%s)',
            policy_record.policyname,
            policy_record.tablename,
            updated_qual
        );
        
        RAISE NOTICE 'Updated policy % on table % to use cached service role function', 
            policy_record.policyname, policy_record.tablename;
    END LOOP;
    
    -- Find policies with admin role checks
    FOR policy_record IN 
        SELECT 
            schemaname, 
            tablename, 
            policyname, 
            qual
        FROM 
            pg_policies
        WHERE 
            schemaname = 'public'
            AND qual LIKE '%role = ''admin''%'
            AND qual NOT LIKE '%auth.is_admin()%'
    LOOP
        -- Replace admin check with is_admin() function
        updated_qual := REGEXP_REPLACE(
            policy_record.qual,
            'EXISTS \\(.*?role = ''admin''.*?\\)',
            'auth.is_admin()',
            'g'
        );
        
        -- Only update if the replacement was successful
        IF updated_qual <> policy_record.qual THEN
            -- Update the policy
            EXECUTE format(
                'ALTER POLICY %I ON %I USING (%s)',
                policy_record.policyname,
                policy_record.tablename,
                updated_qual
            );
            
            RAISE NOTICE 'Updated policy % on table % to use admin function', 
                policy_record.policyname, policy_record.tablename;
        END IF;
    END LOOP;
END $$;

-- Step 9: Create a monitoring view for RLS policies
CREATE OR REPLACE VIEW auth.rls_policy_overview AS
SELECT
    p.schemaname,
    p.tablename,
    p.policyname,
    p.cmd,
    p.roles,
    CASE
        WHEN p.qual LIKE '%auth.is_service_role_cached()%' THEN 'Service Role (cached)'
        WHEN p.qual LIKE '%auth.is_service_role()%' THEN 'Service Role (uncached)'
        WHEN p.qual LIKE '%auth.is_admin()%' THEN 'Admin'
        WHEN p.qual LIKE '%role = ''admin''%' THEN 'Admin (direct check)'
        WHEN p.policyname LIKE '%service_role%' THEN 'Service Role (by name)'
        WHEN p.policyname LIKE '%admin%' THEN 'Admin (by name)'
        WHEN p.qual LIKE '%created_by = auth.uid()%' THEN 'Owner'
        WHEN p.qual LIKE '%is_public%' THEN 'Public'
        ELSE 'Other'
    END AS policy_type,
    LENGTH(p.qual) AS qual_length,
    p.qual
FROM
    pg_policies p
WHERE
    p.schemaname = 'public'
ORDER BY
    p.tablename,
    p.policyname;

COMMENT ON VIEW auth.rls_policy_overview IS 'Overview of all RLS policies for monitoring and administration';

-- Step 10: Create a function to analyze RLS policy coverage
CREATE OR REPLACE FUNCTION auth.analyze_rls_coverage() RETURNS TABLE (
    tablename TEXT,
    has_rls BOOLEAN,
    total_policies INTEGER,
    has_service_role_policy BOOLEAN,
    has_admin_policy BOOLEAN,
    has_owner_policy BOOLEAN,
    has_public_policy BOOLEAN,
    coverage_score INTEGER,
    recommendations TEXT
) AS $$
BEGIN
    RETURN QUERY
    WITH table_stats AS (
        SELECT
            t.table_name,
            t.has_rls,
            COUNT(p.policyname) AS total_policies,
            BOOL_OR(p.qual LIKE '%auth.is_service_role_cached()%' OR 
                    p.qual LIKE '%auth.is_service_role()%' OR
                    p.policyname LIKE '%service_role%') AS has_service_role_policy,
            BOOL_OR(p.qual LIKE '%auth.is_admin()%' OR 
                    p.qual LIKE '%role = ''admin''%' OR
                    p.policyname LIKE '%admin%') AS has_admin_policy,
            BOOL_OR(p.qual LIKE '%created_by = auth.uid()%') AS has_owner_policy,
            BOOL_OR(p.qual LIKE '%is_public%') AS has_public_policy
        FROM
            auth.get_tables_with_rls() t
        LEFT JOIN
            pg_policies p ON p.schemaname = 'public' AND p.tablename = t.table_name
        GROUP BY
            t.table_name, t.has_rls
    )
    SELECT
        ts.table_name,
        ts.has_rls,
        ts.total_policies,
        ts.has_service_role_policy,
        ts.has_admin_policy,
        ts.has_owner_policy,
        ts.has_public_policy,
        -- Calculate coverage score (0-100)
        CASE WHEN ts.has_rls THEN
            (CASE WHEN ts.has_service_role_policy THEN 25 ELSE 0 END) +
            (CASE WHEN ts.has_admin_policy THEN 25 ELSE 0 END) +
            (CASE WHEN ts.has_owner_policy THEN 25 ELSE 0 END) +
            (CASE WHEN ts.has_public_policy THEN 25 ELSE 0 END)
        ELSE 0 END AS coverage_score,
        -- Generate recommendations
        CASE
            WHEN NOT ts.has_rls THEN 'Enable RLS for this table'
            WHEN ts.total_policies = 0 THEN 'Add RLS policies to this table'
            WHEN NOT ts.has_service_role_policy THEN 'Add service_role policy'
            WHEN NOT ts.has_admin_policy THEN 'Add admin policy'
            WHEN NOT ts.has_owner_policy THEN 'Add owner policy'
            WHEN NOT ts.has_public_policy THEN 'Consider adding public access policy'
            WHEN ts.total_policies > 4 THEN 'Consider consolidating policies'
            ELSE 'Policies look good'
        END AS recommendations
    FROM
        table_stats ts
    ORDER BY
        -- Sort by RLS status, then by coverage score (low to high)
        ts.has_rls DESC,
        coverage_score ASC,
        ts.table_name;
END;
$$ LANGUAGE plpgsql;

COMMENT ON FUNCTION auth.analyze_rls_coverage IS 'Analyze RLS policy coverage and provide recommendations';
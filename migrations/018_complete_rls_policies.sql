-- Migration: Complete RLS Policies
-- Purpose: This migration ensures that all tables have appropriate Row Level Security (RLS) policies
-- to control data access while providing service role bypass for admin functions.

BEGIN;

-- Step 1: Ensure RLS is enabled on all relevant tables

-- Function to enable RLS on a table if it exists
CREATE OR REPLACE FUNCTION enable_rls_if_exists(tbl text) RETURNS void AS $$
BEGIN
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = tbl
    ) THEN
        EXECUTE format('ALTER TABLE public.%I ENABLE ROW LEVEL SECURITY', tbl);
        RAISE NOTICE 'Enabled RLS on table: %', tbl;
    ELSE
        RAISE NOTICE 'Table % does not exist, skipping', tbl;
    END IF;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error enabling RLS on %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Enable RLS on core tables
SELECT enable_rls_if_exists('molecules');
SELECT enable_rls_if_exists('mixtures');
SELECT enable_rls_if_exists('mixture_components');
SELECT enable_rls_if_exists('molecular_properties');
SELECT enable_rls_if_exists('experiments');
SELECT enable_rls_if_exists('property_types');
SELECT enable_rls_if_exists('predictions');
SELECT enable_rls_if_exists('lab_verifications');

-- Enable RLS on additional tables that may exist
SELECT enable_rls_if_exists('teams');
SELECT enable_rls_if_exists('team_roles');
SELECT enable_rls_if_exists('user_team_membership');
SELECT enable_rls_if_exists('shared_items');
SELECT enable_rls_if_exists('user_profiles');
SELECT enable_rls_if_exists('protocols');
SELECT enable_rls_if_exists('protocol_steps');
SELECT enable_rls_if_exists('toxicity_data');
SELECT enable_rls_if_exists('calculation_methods');

-- Step 2: Create or replace standard RLS policy functions

-- Function to create the standard select policy if it doesn't exist
CREATE OR REPLACE FUNCTION create_standard_select_policy(tbl text, owner_column text DEFAULT 'created_by') RETURNS void AS $$
DECLARE
    policy_name text := tbl || '_select_policy';
BEGIN
    -- First drop policy if it exists
    EXECUTE format('DROP POLICY IF EXISTS %I ON public.%I', policy_name, tbl);
    
    -- Create the policy
    IF EXISTS (
        SELECT FROM information_schema.columns 
        WHERE table_schema = 'public' AND table_name = tbl AND column_name = owner_column
    ) THEN
        -- Create policy with owner column
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR SELECT
            USING (
                auth.uid() = %I
                OR
                EXISTS (
                    SELECT 1 FROM user_team_membership utm
                    JOIN team_roles tr ON utm.role_id = tr.id
                    WHERE utm.user_id = auth.uid() AND tr.can_view = true
                )
                OR
                EXISTS (
                    SELECT 1 FROM shared_items si
                    WHERE si.item_type = %L AND si.item_id = id AND (si.shared_with_user_id = auth.uid() OR si.is_public = true)
                )
            )', policy_name, tbl, owner_column, tbl);
    ELSE
        -- Create policy without owner column (public read)
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR SELECT
            USING (true)', policy_name, tbl);
    END IF;
    
    RAISE NOTICE 'Created select policy for table: %', tbl;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error creating select policy for %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Function to create the standard insert policy if it doesn't exist
CREATE OR REPLACE FUNCTION create_standard_insert_policy(tbl text, owner_column text DEFAULT 'created_by') RETURNS void AS $$
DECLARE
    policy_name text := tbl || '_insert_policy';
BEGIN
    -- First drop policy if it exists
    EXECUTE format('DROP POLICY IF EXISTS %I ON public.%I', policy_name, tbl);
    
    -- Create the policy
    IF EXISTS (
        SELECT FROM information_schema.columns 
        WHERE table_schema = 'public' AND table_name = tbl AND column_name = owner_column
    ) THEN
        -- Create policy with owner column
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR INSERT
            WITH CHECK (
                %I = auth.uid()
                OR
                EXISTS (
                    SELECT 1 FROM user_team_membership utm
                    JOIN team_roles tr ON utm.role_id = tr.id
                    WHERE utm.user_id = auth.uid() AND tr.can_create = true
                )
            )', policy_name, tbl, owner_column);
    ELSE
        -- Create policy without owner column (authenticated insert)
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR INSERT
            WITH CHECK (auth.role() = ''authenticated'')', policy_name, tbl);
    END IF;
    
    RAISE NOTICE 'Created insert policy for table: %', tbl;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error creating insert policy for %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Function to create the standard update policy if it doesn't exist
CREATE OR REPLACE FUNCTION create_standard_update_policy(tbl text, owner_column text DEFAULT 'created_by') RETURNS void AS $$
DECLARE
    policy_name text := tbl || '_update_policy';
BEGIN
    -- First drop policy if it exists
    EXECUTE format('DROP POLICY IF EXISTS %I ON public.%I', policy_name, tbl);
    
    -- Create the policy
    IF EXISTS (
        SELECT FROM information_schema.columns 
        WHERE table_schema = 'public' AND table_name = tbl AND column_name = owner_column
    ) THEN
        -- Create policy with owner column
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR UPDATE
            USING (
                %I = auth.uid()
                OR
                EXISTS (
                    SELECT 1 FROM user_team_membership utm
                    JOIN team_roles tr ON utm.role_id = tr.id
                    WHERE utm.user_id = auth.uid() AND tr.can_update = true
                )
            )', policy_name, tbl, owner_column);
    ELSE
        -- Create policy without owner column (no update allowed)
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR UPDATE
            USING (false)', policy_name, tbl);
    END IF;
    
    RAISE NOTICE 'Created update policy for table: %', tbl;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error creating update policy for %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Function to create the standard delete policy if it doesn't exist
CREATE OR REPLACE FUNCTION create_standard_delete_policy(tbl text, owner_column text DEFAULT 'created_by') RETURNS void AS $$
DECLARE
    policy_name text := tbl || '_delete_policy';
BEGIN
    -- First drop policy if it exists
    EXECUTE format('DROP POLICY IF EXISTS %I ON public.%I', policy_name, tbl);
    
    -- Create the policy
    IF EXISTS (
        SELECT FROM information_schema.columns 
        WHERE table_schema = 'public' AND table_name = tbl AND column_name = owner_column
    ) THEN
        -- Create policy with owner column
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR DELETE
            USING (
                %I = auth.uid()
                OR
                EXISTS (
                    SELECT 1 FROM user_team_membership utm
                    JOIN team_roles tr ON utm.role_id = tr.id
                    WHERE utm.user_id = auth.uid() AND tr.can_delete = true
                )
            )', policy_name, tbl, owner_column);
    ELSE
        -- Create policy without owner column (no delete allowed)
        EXECUTE format('
            CREATE POLICY %I ON public.%I
            FOR DELETE
            USING (false)', policy_name, tbl);
    END IF;
    
    RAISE NOTICE 'Created delete policy for table: %', tbl;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error creating delete policy for %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Function to create the service role bypass policy
CREATE OR REPLACE FUNCTION create_service_role_policy(tbl text) RETURNS void AS $$
DECLARE
    policy_name text := tbl || '_service_role_policy';
BEGIN
    -- First drop policy if it exists
    EXECUTE format('DROP POLICY IF EXISTS %I ON public.%I', policy_name, tbl);
    
    -- Create the policy
    EXECUTE format('
        CREATE POLICY %I ON public.%I
        USING (auth.role() = ''service_role'')
        WITH CHECK (auth.role() = ''service_role'')',
        policy_name, tbl);
    
    RAISE NOTICE 'Created service role bypass policy for table: %', tbl;
EXCEPTION WHEN OTHERS THEN
    RAISE NOTICE 'Error creating service role policy for %: %', tbl, SQLERRM;
END;
$$ LANGUAGE plpgsql;

-- Step 3: Apply policies to all tables

-- Apply policies to core tables
SELECT create_standard_select_policy('molecules');
SELECT create_standard_insert_policy('molecules');
SELECT create_standard_update_policy('molecules');
SELECT create_standard_delete_policy('molecules');
SELECT create_service_role_policy('molecules');

SELECT create_standard_select_policy('mixtures');
SELECT create_standard_insert_policy('mixtures');
SELECT create_standard_update_policy('mixtures');
SELECT create_standard_delete_policy('mixtures');
SELECT create_service_role_policy('mixtures');

SELECT create_standard_select_policy('mixture_components');
SELECT create_standard_insert_policy('mixture_components');
SELECT create_standard_update_policy('mixture_components');
SELECT create_standard_delete_policy('mixture_components');
SELECT create_service_role_policy('mixture_components');

SELECT create_standard_select_policy('molecular_properties');
SELECT create_standard_insert_policy('molecular_properties');
SELECT create_standard_update_policy('molecular_properties');
SELECT create_standard_delete_policy('molecular_properties');
SELECT create_service_role_policy('molecular_properties');

SELECT create_standard_select_policy('experiments');
SELECT create_standard_insert_policy('experiments');
SELECT create_standard_update_policy('experiments');
SELECT create_standard_delete_policy('experiments');
SELECT create_service_role_policy('experiments');

SELECT create_standard_select_policy('property_types');
SELECT create_standard_insert_policy('property_types');
SELECT create_standard_update_policy('property_types');
SELECT create_standard_delete_policy('property_types');
SELECT create_service_role_policy('property_types');

SELECT create_standard_select_policy('predictions');
SELECT create_standard_insert_policy('predictions');
SELECT create_standard_update_policy('predictions');
SELECT create_standard_delete_policy('predictions');
SELECT create_service_role_policy('predictions');

-- Lab verifications have special policies - reapply them
DO $$
BEGIN
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'lab_verifications'
    ) THEN
        -- Drop existing policies
        DROP POLICY IF EXISTS lab_verifications_select ON public.lab_verifications;
        DROP POLICY IF EXISTS lab_verifications_insert ON public.lab_verifications;
        DROP POLICY IF EXISTS lab_verifications_update ON public.lab_verifications;
        DROP POLICY IF EXISTS lab_verifications_service_role_policy ON public.lab_verifications;
        
        -- Create specialized policies for lab_verifications
        CREATE POLICY lab_verifications_select ON public.lab_verifications
            FOR SELECT USING (auth.uid() = verifier OR auth.uid() IN (
                SELECT user_id FROM user_team_membership utm
                JOIN team_roles tr ON utm.role_id = tr.id
                WHERE tr.can_view_verifications = true
            ));

        CREATE POLICY lab_verifications_insert ON public.lab_verifications
            FOR INSERT WITH CHECK (auth.uid() = verifier);

        CREATE POLICY lab_verifications_update ON public.lab_verifications
            FOR UPDATE USING (auth.uid() = verifier OR auth.uid() IN (
                SELECT user_id FROM user_team_membership utm
                JOIN team_roles tr ON utm.role_id = tr.id
                WHERE tr.can_update_verifications = true
            ));
            
        CREATE POLICY lab_verifications_service_role_policy ON public.lab_verifications
            USING (auth.role() = 'service_role')
            WITH CHECK (auth.role() = 'service_role');
        
        RAISE NOTICE 'Created specialized policies for lab_verifications table';
    END IF;
END $$;

-- Apply policies to additional tables that may exist
DO $$
BEGIN
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'teams'
    ) THEN
        PERFORM create_standard_select_policy('teams');
        PERFORM create_standard_insert_policy('teams');
        PERFORM create_standard_update_policy('teams');
        PERFORM create_standard_delete_policy('teams');
        PERFORM create_service_role_policy('teams');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'team_roles'
    ) THEN
        PERFORM create_standard_select_policy('team_roles');
        PERFORM create_standard_insert_policy('team_roles');
        PERFORM create_standard_update_policy('team_roles');
        PERFORM create_standard_delete_policy('team_roles');
        PERFORM create_service_role_policy('team_roles');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'user_team_membership'
    ) THEN
        PERFORM create_standard_select_policy('user_team_membership');
        PERFORM create_standard_insert_policy('user_team_membership');
        PERFORM create_standard_update_policy('user_team_membership');
        PERFORM create_standard_delete_policy('user_team_membership');
        PERFORM create_service_role_policy('user_team_membership');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'shared_items'
    ) THEN
        PERFORM create_standard_select_policy('shared_items');
        PERFORM create_standard_insert_policy('shared_items');
        PERFORM create_standard_update_policy('shared_items');
        PERFORM create_standard_delete_policy('shared_items');
        PERFORM create_service_role_policy('shared_items');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'user_profiles'
    ) THEN
        PERFORM create_standard_select_policy('user_profiles');
        PERFORM create_standard_insert_policy('user_profiles');
        PERFORM create_standard_update_policy('user_profiles');
        PERFORM create_standard_delete_policy('user_profiles');
        PERFORM create_service_role_policy('user_profiles');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'protocols'
    ) THEN
        PERFORM create_standard_select_policy('protocols');
        PERFORM create_standard_insert_policy('protocols');
        PERFORM create_standard_update_policy('protocols');
        PERFORM create_standard_delete_policy('protocols');
        PERFORM create_service_role_policy('protocols');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'protocol_steps'
    ) THEN
        PERFORM create_standard_select_policy('protocol_steps');
        PERFORM create_standard_insert_policy('protocol_steps');
        PERFORM create_standard_update_policy('protocol_steps');
        PERFORM create_standard_delete_policy('protocol_steps');
        PERFORM create_service_role_policy('protocol_steps');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'toxicity_data'
    ) THEN
        PERFORM create_standard_select_policy('toxicity_data');
        PERFORM create_standard_insert_policy('toxicity_data');
        PERFORM create_standard_update_policy('toxicity_data');
        PERFORM create_standard_delete_policy('toxicity_data');
        PERFORM create_service_role_policy('toxicity_data');
    END IF;
    
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' AND table_name = 'calculation_methods'
    ) THEN
        PERFORM create_standard_select_policy('calculation_methods');
        PERFORM create_standard_insert_policy('calculation_methods');
        PERFORM create_standard_update_policy('calculation_methods');
        PERFORM create_standard_delete_policy('calculation_methods');
        PERFORM create_service_role_policy('calculation_methods');
    END IF;
END $$;

-- Step 4: Clean up temporary functions
DROP FUNCTION IF EXISTS enable_rls_if_exists;
DROP FUNCTION IF EXISTS create_standard_select_policy;
DROP FUNCTION IF EXISTS create_standard_insert_policy;
DROP FUNCTION IF EXISTS create_standard_update_policy;
DROP FUNCTION IF EXISTS create_standard_delete_policy;
DROP FUNCTION IF EXISTS create_service_role_policy;

-- Log migration completed
DO $$
BEGIN
    RAISE NOTICE 'RLS policy completion migration completed successfully';
END $$;

COMMIT;
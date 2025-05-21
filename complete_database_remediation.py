#!/usr/bin/env python3
"""
CryoProtect v2 - Complete Database Remediation

This script implements the comprehensive database remediation plan for CryoProtect v2,
addressing all critical issues in the following order:

1. SECURITY: Enable Row Level Security (RLS)
2. STRUCTURE: Standardize Schema & Fix Relationships
3. PERFORMANCE: Add Missing Indexes
4. ROLES: Create Application-Specific Roles
5. DATA: Consolidate Duplicate Tables

Usage:
    python complete_database_remediation.py [--dry-run] [--phase PHASE]
"""

import os
import sys
import argparse
import logging
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"complete_remediation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def get_supabase_client():
    """Get a Supabase client with service role key."""
    try:
        from supabase import create_client, Client
        
        SUPABASE_URL = os.getenv("SUPABASE_URL")
        SUPABASE_KEY = os.getenv("SUPABASE_KEY")
        
        if not SUPABASE_URL or not SUPABASE_KEY:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        return create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        sys.exit(1)

def execute_sql(supabase, sql, description, dry_run=False):
    """Execute SQL using the Supabase client."""
    if dry_run:
        logger.info(f"DRY RUN: Would execute SQL for: {description}")
        logger.info(f"SQL: {sql[:200]}...")  # Log first 200 chars of SQL
        return True, None
    
    try:
        logger.info(f"Executing SQL: {description}")
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing SQL: {response.error}")
            return False, response.error
        
        logger.info(f"SQL executed successfully: {description}")
        return True, response.data
    except Exception as e:
        logger.error(f"Error executing SQL ({description}): {str(e)}")
        return False, str(e)

def implement_security(supabase, dry_run=False):
    """
    Phase 1: Implement security by enabling RLS on all tables and creating policies.
    This is the HIGHEST PRIORITY phase.
    """
    logger.info("Phase 1: Implementing Security (RLS)")
    
    # SQL to enable RLS on all tables
    enable_rls_sql = """
    -- Enable RLS on all tables
    DO $$
    DECLARE r RECORD;
    BEGIN
      FOR r IN SELECT tablename FROM pg_tables WHERE schemaname = 'public'
      LOOP
        EXECUTE format('ALTER TABLE public.%I ENABLE ROW LEVEL SECURITY;', r.tablename);
      END LOOP;
    END $$;

    -- Restrict anonymous access
    REVOKE ALL ON SCHEMA public FROM anon;
    REVOKE ALL ON ALL TABLES IN SCHEMA public FROM anon;
    GRANT USAGE ON SCHEMA public TO anon;
    GRANT SELECT ON public.property_types TO anon; -- Only non-sensitive tables

    -- Base RLS policies for authenticated users
    DO $$
    DECLARE r RECORD;
    BEGIN
      FOR r IN SELECT tablename FROM pg_tables WHERE schemaname = 'public' 
               AND tablename NOT IN ('property_types')
      LOOP
        EXECUTE format('
          CREATE POLICY "Auth users access own data" ON public.%I
          FOR ALL TO authenticated
          USING (created_by = (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));
        ', r.tablename);
      END LOOP;
    END $$;

    -- Add RLS performance indexes
    CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id ON public.user_profile(auth_user_id);
    CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON public.molecules(created_by);
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON public.mixtures(created_by);
    """
    
    success, result = execute_sql(supabase, enable_rls_sql, "Enable RLS on all tables", dry_run)
    
    if not success and not dry_run:
        logger.error("Failed to implement security")
        return False
    
    logger.info("Security implementation completed successfully")
    return True

def standardize_schema(supabase, dry_run=False):
    """
    Phase 2: Standardize schema by renaming tables to plural form and fixing relationships.
    """
    logger.info("Phase 2: Standardizing Schema")
    
    # SQL to standardize schema
    standardize_schema_sql = """
    -- Standardize to plural names
    ALTER TABLE IF EXISTS public.molecule RENAME TO molecules;
    ALTER TABLE IF EXISTS public.experiment RENAME TO experiments;
    ALTER TABLE IF EXISTS public.prediction RENAME TO predictions;
    ALTER TABLE IF EXISTS public.project RENAME TO projects;

    -- Update foreign keys after renaming
    ALTER TABLE public.mixture_component
    DROP CONSTRAINT IF EXISTS mixture_component_molecule_id_fkey,
    ADD CONSTRAINT mixture_component_molecule_id_fkey 
    FOREIGN KEY (molecule_id) REFERENCES public.molecules(id);

    -- Fix fan trap with junction table
    CREATE TABLE IF NOT EXISTS public.experiment_mixtures (
      id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
      experiment_id UUID NOT NULL REFERENCES public.experiments(id),
      mixture_id UUID NOT NULL REFERENCES public.mixtures(id),
      relationship_type TEXT,
      created_at TIMESTAMPTZ NOT NULL DEFAULT now(),
      updated_at TIMESTAMPTZ NOT NULL DEFAULT now(),
      created_by UUID REFERENCES public.user_profile(id),
      UNIQUE(experiment_id, mixture_id)
    );

    -- Secure new junction table
    ALTER TABLE public.experiment_mixtures ENABLE ROW LEVEL SECURITY;
    CREATE POLICY "Users access their experiment mixtures" ON public.experiment_mixtures
    FOR ALL TO authenticated USING (created_by = (SELECT id FROM user_profile WHERE auth_user_id = auth.uid()));
    CREATE INDEX idx_experiment_mixtures_experiment_id ON public.experiment_mixtures(experiment_id);
    CREATE INDEX idx_experiment_mixtures_mixture_id ON public.experiment_mixtures(mixture_id);
    """
    
    success, result = execute_sql(supabase, standardize_schema_sql, "Standardize schema", dry_run)
    
    if not success and not dry_run:
        logger.error("Failed to standardize schema")
        return False
    
    logger.info("Schema standardization completed successfully")
    return True

def add_performance_indexes(supabase, dry_run=False):
    """
    Phase 3: Add missing indexes for performance optimization.
    """
    logger.info("Phase 3: Adding Performance Indexes")
    
    # SQL to add missing indexes
    add_indexes_sql = """
    -- Index all foreign keys
    DO $$
    DECLARE
        r RECORD;
        index_name TEXT;
        index_sql TEXT;
    BEGIN
        FOR r IN
            SELECT
                tc.table_name, kcu.column_name
            FROM
                information_schema.table_constraints tc
                JOIN information_schema.key_column_usage kcu ON tc.constraint_name = kcu.constraint_name
                JOIN information_schema.constraint_column_usage ccu ON ccu.constraint_name = tc.constraint_name
            WHERE
                tc.constraint_type = 'FOREIGN KEY'
                AND tc.table_schema = 'public'
        LOOP
            index_name := 'idx_' || r.table_name || '_' || r.column_name;
            
            -- Check if index already exists
            IF NOT EXISTS (
                SELECT 1 FROM pg_indexes 
                WHERE tablename = r.table_name 
                AND indexname = index_name
            ) THEN
                index_sql := 'CREATE INDEX ' || index_name || ' ON public.' || r.table_name || '(' || r.column_name || ');';
                EXECUTE index_sql;
                RAISE NOTICE 'Created index: %', index_name;
            ELSE
                RAISE NOTICE 'Index already exists: %', index_name;
            END IF;
        END LOOP;
    END $$;
    """
    
    success, result = execute_sql(supabase, add_indexes_sql, "Add performance indexes", dry_run)
    
    if not success and not dry_run:
        logger.error("Failed to add performance indexes")
        return False
    
    logger.info("Performance indexes added successfully")
    return True

def create_application_roles(supabase, dry_run=False):
    """
    Phase 4: Create application-specific roles with appropriate permissions.
    """
    logger.info("Phase 4: Creating Application-Specific Roles")
    
    # SQL to create application roles
    create_roles_sql = """
    -- Create application roles
    CREATE ROLE app_readonly;
    GRANT USAGE ON SCHEMA public TO app_readonly;
    GRANT SELECT ON ALL TABLES IN SCHEMA public TO app_readonly;

    CREATE ROLE app_readwrite;
    GRANT USAGE ON SCHEMA public TO app_readwrite;
    GRANT SELECT, INSERT, UPDATE ON ALL TABLES IN SCHEMA public TO app_readwrite;

    -- Secure admin operations with SECURITY DEFINER functions
    CREATE OR REPLACE FUNCTION admin_create_user(
      display_name TEXT, email TEXT, password TEXT
    ) RETURNS UUID LANGUAGE plpgsql SECURITY DEFINER AS $$
    DECLARE user_id UUID;
    BEGIN
      INSERT INTO auth.users (email, password) VALUES (email, password) RETURNING id INTO user_id;
      INSERT INTO public.user_profile (auth_user_id, display_name, email) VALUES (user_id, display_name, email);
      RETURN user_id;
    END; $$;
    """
    
    success, result = execute_sql(supabase, create_roles_sql, "Create application roles", dry_run)
    
    if not success and not dry_run:
        logger.error("Failed to create application roles")
        return False
    
    logger.info("Application roles created successfully")
    return True

def consolidate_duplicate_tables(supabase, dry_run=False):
    """
    Phase 5: Consolidate duplicate tables.
    """
    logger.info("Phase 5: Consolidating Duplicate Tables")
    
    # SQL to consolidate duplicate tables
    consolidate_tables_sql = """
    -- Migrate data from prediction to predictions
    INSERT INTO predictions (id, mixture_id, property_type_id, calculation_method_id, numeric_value, created_by)
    SELECT id, mixture_id, property_type_id, method_id, value, created_by
    FROM prediction
    ON CONFLICT (id) DO UPDATE SET updated_at = NOW();
    """
    
    success, result = execute_sql(supabase, consolidate_tables_sql, "Consolidate duplicate tables", dry_run)
    
    if not success and not dry_run:
        logger.error("Failed to consolidate duplicate tables")
        return False
    
    logger.info("Duplicate tables consolidated successfully")
    return True

def verify_remediation(supabase, dry_run=False):
    """Verify that all remediation steps were applied correctly."""
    logger.info("Verifying remediation...")
    
    if dry_run:
        logger.info("DRY RUN: Would verify remediation")
        return True
    
    # Verify RLS enabled
    verify_rls_sql = """
    SELECT relname, relrowsecurity FROM pg_class c 
    JOIN pg_namespace n ON n.oid = c.relnamespace 
    WHERE n.nspname = 'public' AND c.relkind = 'r';
    """
    
    success, rls_results = execute_sql(supabase, verify_rls_sql, "Verify RLS enabled")
    
    if success:
        tables_without_rls = [r['relname'] for r in rls_results if not r['relrowsecurity']]
        if tables_without_rls:
            logger.warning(f"Tables without RLS: {tables_without_rls}")
            return False
    else:
        return False
    
    # Verify anon permissions removed
    verify_anon_sql = """
    SELECT table_name, has_table_privilege('anon', 'public.' || table_name, 'SELECT') as anon_select,
           has_table_privilege('anon', 'public.' || table_name, 'INSERT') as anon_insert
    FROM information_schema.tables WHERE table_schema = 'public';
    """
    
    success, anon_results = execute_sql(supabase, verify_anon_sql, "Verify anon permissions removed")
    
    if success:
        tables_with_anon_access = [r['table_name'] for r in anon_results if r['anon_select'] or r['anon_insert']]
        if len(tables_with_anon_access) > 1:  # Allow property_types
            logger.warning(f"Tables with anon access: {tables_with_anon_access}")
            return False
    else:
        return False
    
    logger.info("All remediation steps verified successfully")
    return True

def main():
    """Main function to run the complete database remediation."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 Complete Database Remediation")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    parser.add_argument("--phase", type=int, choices=range(1, 6), help="Run a specific phase (1-5)")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Complete Database Remediation")
    
    # Connect to Supabase
    supabase = get_supabase_client()
    
    # Run a specific phase or all phases
    success = True
    
    if args.phase:
        logger.info(f"Running Phase {args.phase} only")
        
        if args.phase == 1:
            success = implement_security(supabase, args.dry_run)
        elif args.phase == 2:
            success = standardize_schema(supabase, args.dry_run)
        elif args.phase == 3:
            success = add_performance_indexes(supabase, args.dry_run)
        elif args.phase == 4:
            success = create_application_roles(supabase, args.dry_run)
        elif args.phase == 5:
            success = consolidate_duplicate_tables(supabase, args.dry_run)
    else:
        logger.info("Running all phases in sequence")
        
        # Phase 1: Security (HIGHEST PRIORITY)
        if not implement_security(supabase, args.dry_run):
            logger.error("Phase 1 (Security) failed")
            success = False
        
        # Phase 2: Schema Standardization
        if success and not standardize_schema(supabase, args.dry_run):
            logger.error("Phase 2 (Schema Standardization) failed")
            success = False
        
        # Phase 3: Performance Indexes
        if success and not add_performance_indexes(supabase, args.dry_run):
            logger.error("Phase 3 (Performance Indexes) failed")
            success = False
        
        # Phase 4: Application Roles
        if success and not create_application_roles(supabase, args.dry_run):
            logger.error("Phase 4 (Application Roles) failed")
            success = False
        
        # Phase 5: Consolidate Duplicate Tables
        if success and not consolidate_duplicate_tables(supabase, args.dry_run):
            logger.error("Phase 5 (Consolidate Duplicate Tables) failed")
            success = False
    
    # Verify the remediation
    if success and not args.dry_run:
        if not verify_remediation(supabase):
            logger.warning("Verification found issues after remediation")
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("CryoProtect v2 Database Remediation Summary")
    logger.info("=" * 80)
    
    if success:
        logger.info("\nStatus: SUCCESS")
        logger.info("All specified phases completed successfully.")
    else:
        logger.info("\nStatus: FAILED")
        logger.info("One or more phases failed. Check the log for details.")
    
    logger.info("\nFor detailed information, check the log file.")
    logger.info("=" * 80 + "\n")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())

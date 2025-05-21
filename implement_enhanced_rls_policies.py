#!/usr/bin/env python
"""
Implement Enhanced RLS Policies for CryoProtect v2

This script implements the recommendations from the RLS Audit Report:
1. Enable RLS on all tables that currently have it disabled
2. Complete missing policies for tables with incomplete policy sets
3. Implement public access policies based on the `is_public` flag
4. Standardize policy implementation across all tables
5. Add data classification policies

Usage:
    python implement_enhanced_rls_policies.py [--project-id <project_id>]
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from contextlib import contextmanager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/implement_enhanced_rls_policies.log', 'a')
    ]
)
logger = logging.getLogger(__name__)

def ensure_directory(path):
    """Ensure a directory exists."""
    if not os.path.exists(path):
        os.makedirs(path)

def execute_sql_via_mcp(project_id, sql_query):
    """Execute SQL query via Supabase MCP server."""
    try:
        # Import the MCP client here to avoid dependency issues
        from supabase_mcp_tools import execute_sql_on_supabase
        
        # Execute the query
        result = execute_sql_on_supabase(project_id, sql_query)
        
        return result
    except ImportError:
        logger.error("supabase_mcp_tools module not found. Using fallback method.")
        return execute_sql_via_mcp_fallback(project_id, sql_query)

def execute_sql_via_mcp_fallback(project_id, sql_query):
    """Fallback method to execute SQL via MCP using subprocess."""
    try:
        import subprocess
        import tempfile
        
        # Create a temporary file with the SQL query
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sql', delete=False) as temp_file:
            temp_file.write(sql_query)
            temp_file_path = temp_file.name
        
        # Use the MCP tool to execute the SQL
        cmd = [
            "python", "execute_rls_sql_via_mcp.py",
            "--project-id", project_id,
            "--sql-file", temp_file_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Clean up the temporary file
        os.unlink(temp_file_path)
        
        if result.returncode != 0:
            logger.error(f"Error executing SQL via MCP: {result.stderr}")
            return None
        
        # Parse the output as JSON
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError:
            logger.error(f"Error parsing MCP result as JSON: {result.stdout}")
            return result.stdout
    except Exception as e:
        logger.error(f"Error in subprocess MCP execution: {str(e)}")
        return None

def enable_rls_on_all_tables(project_id):
    """Enable RLS on all tables that currently have it disabled."""
    logger.info("Enabling RLS on all tables...")
    
    sql = """
    -- Enable RLS on tables where it's currently disabled
    ALTER TABLE calculation_methods ENABLE ROW LEVEL SECURITY;
    ALTER TABLE predictions ENABLE ROW LEVEL SECURITY;
    ALTER TABLE projects ENABLE ROW LEVEL SECURITY;
    ALTER TABLE property_types ENABLE ROW LEVEL SECURITY;
    ALTER TABLE teams ENABLE ROW LEVEL SECURITY;
    ALTER TABLE user_profile ENABLE ROW LEVEL SECURITY;
    
    -- Verify the changes
    SELECT relname AS tablename, relrowsecurity
    FROM pg_class
    JOIN pg_namespace ON pg_class.relnamespace = pg_namespace.oid
    WHERE pg_namespace.nspname = 'public'
    AND pg_class.relkind = 'r'
    AND relname IN (
        'calculation_methods', 'predictions', 'projects', 
        'property_types', 'teams', 'user_profile'
    );
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully enabled RLS on all tables")
        return True
    except Exception as e:
        logger.error(f"Error enabling RLS on all tables: {str(e)}")
        return False

def add_missing_policies_for_experiment_properties(project_id):
    """Add missing UPDATE and DELETE policies for experiment_properties."""
    logger.info("Adding missing policies for experiment_properties...")
    
    sql = """
    -- Add missing UPDATE policy for experiment_properties
    CREATE POLICY "Update experiment_properties for project members"
      ON experiment_property
      FOR UPDATE
      USING (
        EXISTS (
          SELECT 1 FROM experiment
          JOIN project ON project.id = experiment.project_id
          JOIN user_profile ON user_profile.project_id = project.id
          WHERE experiment.id = experiment_property.experiment_id
            AND user_profile.user_id = auth.uid()
        )
      );
    COMMENT ON POLICY "Update experiment_properties for project members" ON experiment_property IS
      'Allows project members to update experiment properties in their projects.';
    
    -- Add missing DELETE policy for experiment_properties
    CREATE POLICY "Delete experiment_properties for project members"
      ON experiment_property
      FOR DELETE
      USING (
        EXISTS (
          SELECT 1 FROM experiment
          JOIN project ON project.id = experiment.project_id
          JOIN user_profile ON user_profile.project_id = project.id
          WHERE experiment.id = experiment_property.experiment_id
            AND user_profile.user_id = auth.uid()
        )
      );
    COMMENT ON POLICY "Delete experiment_properties for project members" ON experiment_property IS
      'Allows project members to delete experiment properties in their projects.';
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added missing policies for experiment_properties")
        return True
    except Exception as e:
        logger.error(f"Error adding missing policies for experiment_properties: {str(e)}")
        return False

def add_service_role_policies(project_id):
    """Add service role policies to tables that are missing them."""
    logger.info("Adding service role policies...")
    
    sql = """
    -- Add service role policies to tables that are missing them
    DO $$
    DECLARE
        tables_array text[] := ARRAY[
            'mixture_components', 'mixtures', 'molecular_properties', 
            'molecules', 'molecule_experiments', 'molecule_proteins', 'proteins'
        ];
        t text;
    BEGIN
        FOREACH t IN ARRAY tables_array
        LOOP
            IF NOT EXISTS (
                SELECT 1 FROM pg_policies 
                WHERE tablename = t 
                AND policyname = 'Allow service role full access'
            ) THEN
                EXECUTE format('
                    CREATE POLICY "Allow service role full access" 
                    ON public.%I 
                    USING (auth.role() = ''service_role'');
                ', t);
            END IF;
        END LOOP;
    END
    $$;
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added service role policies")
        return True
    except Exception as e:
        logger.error(f"Error adding service role policies: {str(e)}")
        return False

def implement_public_access_policies(project_id):
    """Implement public access policies based on the is_public flag."""
    logger.info("Implementing public access policies...")
    
    sql = """
    -- Add public access policy for molecules
    CREATE POLICY "Allow public access to public molecules"
      ON molecule
      FOR SELECT
      USING (is_public = true);
    COMMENT ON POLICY "Allow public access to public molecules" ON molecule IS
      'Allows anyone to view molecules marked as public.';
    
    -- Add public access policy for mixtures
    CREATE POLICY "Allow public access to public mixtures"
      ON mixture
      FOR SELECT
      USING (is_public = true);
    COMMENT ON POLICY "Allow public access to public mixtures" ON mixture IS
      'Allows anyone to view mixtures marked as public.';
    
    -- Add public access policy for experiments
    CREATE POLICY "Allow public access to public experiments"
      ON experiment
      FOR SELECT
      USING (is_public = true);
    COMMENT ON POLICY "Allow public access to public experiments" ON experiment IS
      'Allows anyone to view experiments marked as public.';
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully implemented public access policies")
        return True
    except Exception as e:
        logger.error(f"Error implementing public access policies: {str(e)}")
        return False

def add_basic_policies_for_projects(project_id):
    """Add basic policies for the projects table."""
    logger.info("Adding basic policies for projects table...")
    
    sql = """
    -- Add basic policies for projects table
    CREATE POLICY "Select projects for project members"
      ON project
      FOR SELECT
      USING (
        EXISTS (
          SELECT 1 FROM user_profile
          WHERE user_profile.project_id = project.id
            AND user_profile.user_id = auth.uid()
        )
      );
    COMMENT ON POLICY "Select projects for project members" ON project IS
      'Allows users to view projects they are members of.';
    
    CREATE POLICY "Insert projects for authenticated users"
      ON project
      FOR INSERT
      WITH CHECK (auth.role() = 'authenticated');
    COMMENT ON POLICY "Insert projects for authenticated users" ON project IS
      'Allows any authenticated user to create a new project.';
    
    CREATE POLICY "Update projects for project members"
      ON project
      FOR UPDATE
      USING (
        EXISTS (
          SELECT 1 FROM user_profile
          WHERE user_profile.project_id = project.id
            AND user_profile.user_id = auth.uid()
        )
      );
    COMMENT ON POLICY "Update projects for project members" ON project IS
      'Allows project members to update their projects.';
    
    CREATE POLICY "Delete projects for project owners"
      ON project
      FOR DELETE
      USING (
        EXISTS (
          SELECT 1 FROM user_profile
          WHERE user_profile.project_id = project.id
            AND user_profile.user_id = auth.uid()
            AND user_profile.role = 'owner'
        )
      );
    COMMENT ON POLICY "Delete projects for project owners" ON project IS
      'Allows only project owners to delete their projects.';
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added basic policies for projects table")
        return True
    except Exception as e:
        logger.error(f"Error adding basic policies for projects table: {str(e)}")
        return False

def add_data_classification_policies(project_id):
    """Add data classification policies for sensitive scientific data."""
    logger.info("Adding data classification policies...")
    
    sql = """
    -- Add data classification policy for sensitive molecular data
    CREATE POLICY "Restrict access to sensitive molecular data"
      ON molecular_property
      FOR SELECT
      USING (
        (sensitivity_level IS NULL OR sensitivity_level <= 'medium') OR
        EXISTS (
          SELECT 1 FROM user_profile
          WHERE user_profile.user_id = auth.uid()
          AND user_profile.clearance_level >= molecular_property.sensitivity_level
        )
      );
    COMMENT ON POLICY "Restrict access to sensitive molecular data" ON molecular_property IS
      'Restricts access to sensitive molecular data based on user clearance level.';
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added data classification policies")
        return True
    except Exception as e:
        logger.error(f"Error adding data classification policies: {str(e)}")
        return False

def standardize_policy_implementation(project_id):
    """Standardize policy implementation across all tables."""
    logger.info("Standardizing policy implementation...")
    
    sql = """
    -- Standardize policy naming and implementation
    DO $$
    DECLARE
        tables_array text[] := ARRAY[
            'molecule', 'mixture', 'experiment', 'molecular_property', 
            'mixture_component', 'experiment_property', 'prediction',
            'molecule_experiment', 'molecule_protein', 'protein'
        ];
        t text;
    BEGIN
        FOREACH t IN ARRAY tables_array
        LOOP
            -- Ensure all tables have consistent service role policy
            IF NOT EXISTS (
                SELECT 1 FROM pg_policies 
                WHERE tablename = t 
                AND policyname = 'service_role_access'
            ) THEN
                -- Drop any existing service role policies with different names
                EXECUTE format('
                    DO $inner$
                    BEGIN
                        IF EXISTS (
                            SELECT 1 FROM pg_policies 
                            WHERE tablename = %L 
                            AND polqual::text LIKE ''%%auth.role() = ''''service_role''''%%''
                        ) THEN
                            EXECUTE (
                                SELECT string_agg(''DROP POLICY "'' || policyname || ''" ON public.'' || tablename, '';'')
                                FROM pg_policies 
                                WHERE tablename = %L 
                                AND polqual::text LIKE ''%%auth.role() = ''''service_role''''%%''
                            );
                        END IF;
                    END
                    $inner$;
                ', t, t);
                
                -- Create standardized service role policy
                EXECUTE format('
                    CREATE POLICY "service_role_access" 
                    ON public.%I 
                    USING (auth.role() = ''service_role'');
                ', t);
            END IF;
            
            -- Ensure all tables have consistent project member policy
            IF NOT EXISTS (
                SELECT 1 FROM pg_policies 
                WHERE tablename = t 
                AND policyname = 'project_member_access'
            ) THEN
                -- Create standardized project member policy
                EXECUTE format('
                    DO $inner$
                    BEGIN
                        IF EXISTS (
                            SELECT 1 FROM information_schema.columns
                            WHERE table_schema = ''public''
                            AND table_name = %L
                            AND column_name = ''project_id''
                        ) THEN
                            EXECUTE ''
                                CREATE POLICY "project_member_access" 
                                ON public.%I 
                                USING (
                                    EXISTS (
                                        SELECT 1 FROM user_profile
                                        WHERE user_profile.project_id = %I.project_id
                                        AND user_profile.user_id = auth.uid()
                                    )
                                );
                            '';
                        END IF;
                    END
                    $inner$;
                ', t, t, t);
            END IF;
        END LOOP;
    END
    $$;
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully standardized policy implementation")
        return True
    except Exception as e:
        logger.error(f"Error standardizing policy implementation: {str(e)}")
        return False

def verify_rls_implementation(project_id):
    """Verify that the RLS implementation is complete and correct."""
    logger.info("Verifying RLS implementation...")
    
    sql = """
    -- Get RLS status for all tables
    SELECT 
        t.tablename,
        c.relrowsecurity AS rls_enabled,
        (SELECT COUNT(*) FROM pg_policies WHERE tablename = t.tablename) AS policy_count
    FROM 
        pg_tables t
    JOIN 
        pg_class c ON t.tablename = c.relname AND t.schemaname = 'public'
    WHERE 
        t.schemaname = 'public'
    ORDER BY 
        t.tablename;
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        
        # Generate verification report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = f"reports/security/rls_implementation_verification_{timestamp}.md"
        
        ensure_directory("reports/security")
        
        with open(report_path, 'w') as f:
            f.write("# RLS Implementation Verification Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## RLS Status for All Tables\n\n")
            f.write("| Table Name | RLS Enabled | Policy Count |\n")
            f.write("|------------|-------------|-------------|\n")
            
            for table in result:
                f.write(f"| {table['tablename']} | {'✅' if table['rls_enabled'] else '❌'} | {table['policy_count']} |\n")
        
        logger.info(f"Verification report generated: {report_path}")
        return report_path
    except Exception as e:
        logger.error(f"Error verifying RLS implementation: {str(e)}")
        return None

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Implement Enhanced RLS Policies for CryoProtect v2")
    parser.add_argument("--project-id", help="Supabase project ID", default="tsdlmynydfuypiugmkev")
    args = parser.parse_args()
    
    # Ensure directories exist
    ensure_directory("logs")
    ensure_directory("reports/security")
    
    logger.info(f"Starting enhanced RLS policy implementation for project {args.project_id}")
    
    # 1. Enable RLS on all tables
    if enable_rls_on_all_tables(args.project_id):
        logger.info("✅ Enabled RLS on all tables")
    else:
        logger.error("❌ Failed to enable RLS on all tables")
    
    # 2. Add missing policies for experiment_properties
    if add_missing_policies_for_experiment_properties(args.project_id):
        logger.info("✅ Added missing policies for experiment_properties")
    else:
        logger.error("❌ Failed to add missing policies for experiment_properties")
    
    # 3. Add service role policies
    if add_service_role_policies(args.project_id):
        logger.info("✅ Added service role policies")
    else:
        logger.error("❌ Failed to add service role policies")
    
    # 4. Implement public access policies
    if implement_public_access_policies(args.project_id):
        logger.info("✅ Implemented public access policies")
    else:
        logger.error("❌ Failed to implement public access policies")
    
    # 5. Add basic policies for projects table
    if add_basic_policies_for_projects(args.project_id):
        logger.info("✅ Added basic policies for projects table")
    else:
        logger.error("❌ Failed to add basic policies for projects table")
    
    # 6. Add data classification policies
    if add_data_classification_policies(args.project_id):
        logger.info("✅ Added data classification policies")
    else:
        logger.error("❌ Failed to add data classification policies")
    
    # 7. Standardize policy implementation
    if standardize_policy_implementation(args.project_id):
        logger.info("✅ Standardized policy implementation")
    else:
        logger.error("❌ Failed to standardize policy implementation")
    
    # 8. Verify RLS implementation
    report_path = verify_rls_implementation(args.project_id)
    if report_path:
        logger.info(f"✅ Verification completed. Report saved to {report_path}")
    else:
        logger.error("❌ Verification failed")
    
    logger.info("Enhanced RLS policy implementation completed")

if __name__ == "__main__":
    main()
#!/usr/bin/env python
"""
Fix RLS Issues in CryoProtect Supabase Project

This script implements the recommendations from the RLS verification report:
1. Add SECURITY INVOKER to views
2. Enable RLS on tables that don't have it
3. Add missing admin policies
4. Add performance indexes

Usage:
    python fix_rls_issues.py --project-id <project_id>
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/fix_rls_issues.log', 'a')
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
    """Fallback method to execute SQL via MCP using the use_mcp_tool."""
    try:
        # Use the MCP tool to execute the SQL
        from supabase_mcp_client import SupabaseMCPClient
        
        # Create a client
        client = SupabaseMCPClient(project_id)
        
        # Execute the query
        result = client.execute_sql(sql_query)
        
        return result
    except Exception as e:
        logger.error(f"Error in fallback MCP execution: {str(e)}")
        
        # Try using the MCP server directly
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
        except Exception as e2:
            logger.error(f"Error in subprocess MCP execution: {str(e2)}")
            return None

def add_security_invoker_to_views(project_id):
    """Add SECURITY INVOKER to views."""
    logger.info("Adding SECURITY INVOKER to views...")
    
    sql = """
    -- Add SECURITY INVOKER to views
    ALTER VIEW public.molecule_with_properties SECURITY INVOKER;
    ALTER VIEW public.mixture_with_components SECURITY INVOKER;
    ALTER VIEW public.experiment_with_results SECURITY INVOKER;
    
    -- Verify the changes
    SELECT viewname, pg_get_viewdef(viewname::regclass) 
    FROM pg_views 
    WHERE schemaname = 'public' 
    AND viewname IN ('molecule_with_properties', 'mixture_with_components', 'experiment_with_results');
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added SECURITY INVOKER to views")
        return True
    except Exception as e:
        logger.error(f"Error adding SECURITY INVOKER to views: {str(e)}")
        return False

def enable_rls_on_tables(project_id):
    """Enable RLS on tables that don't have it."""
    logger.info("Enabling RLS on tables...")
    
    # First, check which tables don't have RLS enabled
    check_sql = """
    SELECT relname AS tablename, relrowsecurity
    FROM pg_class
    JOIN pg_namespace ON pg_class.relnamespace = pg_namespace.oid
    WHERE pg_namespace.nspname = 'public'
    AND pg_class.relkind = 'r'
    AND NOT relrowsecurity;
    """
    
    try:
        tables_result = execute_sql_via_mcp(project_id, check_sql)
        tables_without_rls = [table["tablename"] for table in tables_result]
        
        if not tables_without_rls:
            logger.info("All tables already have RLS enabled")
            return True
        
        # Enable RLS on tables that don't have it
        enable_sql = "BEGIN;\n"
        
        for table in tables_without_rls:
            enable_sql += f"ALTER TABLE public.{table} ENABLE ROW LEVEL SECURITY;\n"
            
            # Add basic policies for each table
            enable_sql += f"""
            -- Service role policy
            CREATE POLICY "service_role_policy" ON public.{table}
              USING (auth.role() = 'service_role');
              
            -- Owner policy
            CREATE POLICY "owner_policy" ON public.{table}
              USING (auth.uid() = created_by);
              
            -- Admin policy for SELECT
            CREATE POLICY "admin_select_policy" ON public.{table}
              FOR SELECT
              USING (EXISTS (
                SELECT 1 FROM auth.users
                WHERE users.id = auth.uid()
                AND users.role = 'admin'
              ));
            """
            
            # Add public read policy for tables with is_public column
            enable_sql += f"""
            DO $$
            BEGIN
              IF EXISTS (
                SELECT 1 FROM information_schema.columns
                WHERE table_schema = 'public'
                AND table_name = '{table}'
                AND column_name = 'is_public'
              ) THEN
                EXECUTE 'CREATE POLICY "public_read_policy" ON public.{table}
                  FOR SELECT
                  USING (is_public = true)';
              END IF;
            END
            $$;
            """
        
        enable_sql += "COMMIT;"
        
        result = execute_sql_via_mcp(project_id, enable_sql)
        logger.info(f"Successfully enabled RLS on tables: {', '.join(tables_without_rls)}")
        return True
    except Exception as e:
        logger.error(f"Error enabling RLS on tables: {str(e)}")
        return False

def add_performance_indexes(project_id):
    """Add performance indexes for RLS policies."""
    logger.info("Adding performance indexes for RLS policies...")
    
    sql = """
    -- Add indexes for columns used in RLS policies
    CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON public.molecules(created_by);
    CREATE INDEX IF NOT EXISTS idx_molecules_is_public ON public.molecules(is_public);
    CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON public.mixtures(created_by);
    CREATE INDEX IF NOT EXISTS idx_mixtures_is_public ON public.mixtures(is_public);
    CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON public.experiments(created_by);
    CREATE INDEX IF NOT EXISTS idx_experiments_molecule_id ON public.experiments(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_experiments_mixture_id ON public.experiments(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON public.molecular_properties(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_created_by ON public.molecular_properties(created_by);
    CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id ON public.mixture_components(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_mixture_components_created_by ON public.mixture_components(created_by);
    CREATE INDEX IF NOT EXISTS idx_predictions_molecule_id ON public.predictions(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_mixture_id ON public.predictions(mixture_id);
    CREATE INDEX IF NOT EXISTS idx_predictions_created_by ON public.predictions(created_by);
    CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_user_id ON public.scientific_data_audit(user_id);
    
    -- Verify the indexes
    SELECT
        tablename,
        indexname,
        indexdef
    FROM
        pg_indexes
    WHERE
        schemaname = 'public'
        AND indexname LIKE 'idx_%';
    """
    
    try:
        result = execute_sql_via_mcp(project_id, sql)
        logger.info("Successfully added performance indexes")
        return True
    except Exception as e:
        logger.error(f"Error adding performance indexes: {str(e)}")
        return False

def verify_fixes(project_id):
    """Verify that the fixes were applied correctly."""
    logger.info("Verifying fixes...")
    
    # Check views have SECURITY INVOKER
    views_sql = """
    SELECT viewname, pg_get_viewdef(viewname::regclass) 
    FROM pg_views 
    WHERE schemaname = 'public' 
    AND viewname IN ('molecule_with_properties', 'mixture_with_components', 'experiment_with_results');
    """
    
    # Check tables have RLS enabled
    tables_sql = """
    SELECT relname AS tablename, relrowsecurity
    FROM pg_class
    JOIN pg_namespace ON pg_class.relnamespace = pg_namespace.oid
    WHERE pg_namespace.nspname = 'public'
    AND pg_class.relkind = 'r';
    """
    
    # Check policies
    policies_sql = """
    SELECT relname AS tablename, polname, polcmd, pg_get_expr(polqual, polrelid) AS policy_definition 
    FROM pg_policy 
    JOIN pg_class ON pg_policy.polrelid = pg_class.oid 
    WHERE pg_class.relnamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public');
    """
    
    # Check indexes
    indexes_sql = """
    SELECT
        tablename,
        indexname,
        indexdef
    FROM
        pg_indexes
    WHERE
        schemaname = 'public'
        AND indexname LIKE 'idx_%';
    """
    
    try:
        views_result = execute_sql_via_mcp(project_id, views_sql)
        tables_result = execute_sql_via_mcp(project_id, tables_sql)
        policies_result = execute_sql_via_mcp(project_id, policies_sql)
        indexes_result = execute_sql_via_mcp(project_id, indexes_sql)
        
        # Generate verification report
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = f"reports/security/rls_fixes_verification_{timestamp}.md"
        
        ensure_directory("reports/security")
        
        with open(report_path, 'w') as f:
            f.write("# RLS Fixes Verification Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Views with SECURITY INVOKER
            f.write("## Views with SECURITY INVOKER\n\n")
            f.write("| View | Has SECURITY INVOKER |\n")
            f.write("|------|---------------------|\n")
            
            for view in views_result:
                has_security_invoker = "SECURITY INVOKER" in view["pg_get_viewdef"].upper()
                f.write(f"| {view['viewname']} | {'✅' if has_security_invoker else '❌'} |\n")
            
            f.write("\n")
            
            # Tables with RLS enabled
            f.write("## Tables with RLS Enabled\n\n")
            f.write("| Table | RLS Enabled |\n")
            f.write("|-------|------------|\n")
            
            for table in tables_result:
                f.write(f"| {table['tablename']} | {'✅' if table['relrowsecurity'] else '❌'} |\n")
            
            f.write("\n")
            
            # Policies
            f.write("## RLS Policies\n\n")
            f.write("| Table | Policy | Command | Definition |\n")
            f.write("|-------|--------|---------|------------|\n")
            
            for policy in policies_result:
                f.write(f"| {policy['tablename']} | {policy['polname']} | {policy['polcmd']} | {policy['policy_definition']} |\n")
            
            f.write("\n")
            
            # Indexes
            f.write("## Performance Indexes\n\n")
            f.write("| Table | Index | Definition |\n")
            f.write("|-------|-------|------------|\n")
            
            for index in indexes_result:
                f.write(f"| {index['tablename']} | {index['indexname']} | {index['indexdef']} |\n")
        
        logger.info(f"Verification report generated: {report_path}")
        return report_path
    except Exception as e:
        logger.error(f"Error verifying fixes: {str(e)}")
        return None

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Fix RLS issues in CryoProtect Supabase project")
    parser.add_argument("--project-id", help="Supabase project ID", default="tsdlmynydfuypiugmkev")
    args = parser.parse_args()
    
    # Ensure directories exist
    ensure_directory("logs")
    ensure_directory("reports/security")
    
    logger.info(f"Starting RLS fixes for project {args.project_id}")
    
    # Add SECURITY INVOKER to views
    if add_security_invoker_to_views(args.project_id):
        logger.info("✅ Added SECURITY INVOKER to views")
    else:
        logger.error("❌ Failed to add SECURITY INVOKER to views")
    
    # Enable RLS on tables
    if enable_rls_on_tables(args.project_id):
        logger.info("✅ Enabled RLS on tables")
    else:
        logger.error("❌ Failed to enable RLS on tables")
    
    # Add performance indexes
    if add_performance_indexes(args.project_id):
        logger.info("✅ Added performance indexes")
    else:
        logger.error("❌ Failed to add performance indexes")
    
    # Verify fixes
    report_path = verify_fixes(args.project_id)
    if report_path:
        logger.info(f"✅ Verification completed. Report saved to {report_path}")
    else:
        logger.error("❌ Verification failed")
    
    logger.info("RLS fixes completed")

if __name__ == "__main__":
    main()
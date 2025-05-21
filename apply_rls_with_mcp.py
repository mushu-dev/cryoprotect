#!/usr/bin/env python3
"""
Apply RLS Optimizations - Using Supabase MCP

This script applies the optimized RLS policies migration to improve
performance of Row Level Security policies in the database using
the Supabase MCP.
"""

import os
import sys
import argparse
import time
import json
from pathlib import Path
from typing import Optional, Dict, Any

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_migration_file(file_path: str) -> str:
    """
    Read the migration SQL file.
    
    Args:
        file_path: Path to the migration file
        
    Returns:
        SQL content as string
    """
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading migration file {file_path}: {str(e)}")
        raise

def execute_sql_via_mcp(project_id: str, sql: str) -> bool:
    """
    Execute SQL via Supabase MCP.
    
    Args:
        project_id: Supabase project ID
        sql: SQL to execute
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Import the MCP execute_sql function
        try:
            from supabase_mcp_tools import execute_sql_on_supabase
            
            logger.info("Executing SQL via MCP...")
            result = execute_sql_on_supabase(project_id, sql)
            
            # Check if result indicates success or error
            if isinstance(result, dict) and result.get('error'):
                logger.error(f"MCP execution error: {result.get('error')}")
                return False
            
            logger.info("SQL executed successfully via MCP")
            return True
        except ImportError:
            logger.warning("supabase_mcp_tools not found, trying alternate method...")
            
            # Fallback to mcp__supabase__execute_sql
            try:
                import mcp__supabase__execute_sql
                
                logger.info("Executing SQL via direct MCP function...")
                result = mcp__supabase__execute_sql.execute({
                    "project_id": project_id,
                    "query": sql
                })
                
                # Check if result indicates success or error
                if isinstance(result, dict) and result.get('error'):
                    logger.error(f"MCP execution error: {result.get('error')}")
                    return False
                
                logger.info("SQL executed successfully via direct MCP function")
                return True
            except ImportError:
                logger.error("Could not import MCP modules, falling back to command line execution...")
                
                # Last resort: Use command line execution
                return execute_sql_via_command_line(project_id, sql)
    except Exception as e:
        logger.error(f"Error executing SQL via MCP: {str(e)}")
        return False

def execute_sql_via_command_line(project_id: str, sql: str) -> bool:
    """
    Execute SQL via command line MCP tool.
    
    Args:
        project_id: Supabase project ID
        sql: SQL to execute
        
    Returns:
        True if successful, False otherwise
    """
    try:
        import subprocess
        import tempfile
        
        # Create temporary file for SQL
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sql', delete=False) as temp_file:
            temp_file.write(sql)
            temp_file_path = temp_file.name
        
        try:
            # Try execute_rls_sql_via_mcp.py first
            if os.path.exists("execute_rls_sql_via_mcp.py"):
                logger.info("Executing SQL via execute_rls_sql_via_mcp.py...")
                cmd = [
                    "python", "execute_rls_sql_via_mcp.py",
                    "--project-id", project_id,
                    "--sql-file", temp_file_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    logger.error(f"Command failed: {result.stderr}")
                    return False
                
                logger.info("SQL executed successfully via execute_rls_sql_via_mcp.py")
                return True
            else:
                # Try mcp_supabase_execute_sql.py next
                if os.path.exists("mcp_supabase_execute_sql.py"):
                    logger.info("Executing SQL via mcp_supabase_execute_sql.py...")
                    cmd = [
                        "python", "mcp_supabase_execute_sql.py",
                        "--project-id", project_id,
                        "--sql-file", temp_file_path
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode != 0:
                        logger.error(f"Command failed: {result.stderr}")
                        return False
                    
                    logger.info("SQL executed successfully via mcp_supabase_execute_sql.py")
                    return True
                else:
                    # Try generic approach
                    logger.info("Executing SQL via mcp__supabase__execute_sql...")
                    cmd = [
                        "python", "-c",
                        f"import json; import mcp__supabase__execute_sql; "
                        f"print(json.dumps(mcp__supabase__execute_sql.execute({{'project_id': '{project_id}', 'query': open('{temp_file_path}').read()}})))"
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode != 0:
                        logger.error(f"Command failed: {result.stderr}")
                        return False
                    
                    logger.info("SQL executed successfully via mcp__supabase__execute_sql")
                    return True
        finally:
            # Clean up temporary file
            os.unlink(temp_file_path)
    except Exception as e:
        logger.error(f"Error executing SQL via command line: {str(e)}")
        return False

def generate_report(success: bool, migration_file: str, project_id: str) -> str:
    """
    Generate a simple report file.
    
    Args:
        success: Whether the migration was successful
        migration_file: Path to the migration file
        project_id: Supabase project ID
        
    Returns:
        Path to the report file
    """
    # Ensure directory exists
    os.makedirs("reports", exist_ok=True)
    
    # Generate timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/rls_optimization_report_{timestamp}.md"
    
    with open(report_path, 'w') as f:
        f.write("# RLS Optimization Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Applied migration: {migration_file}\n\n")
        f.write(f"Supabase Project ID: {project_id}\n\n")
        f.write(f"Status: {'✅ Success' if success else '❌ Failed'}\n\n")
        
        if success:
            f.write("## Next Steps\n\n")
            f.write("1. Verify the RLS policies are working as expected\n")
            f.write("2. Monitor query performance to ensure the optimizations are effective\n")
            f.write("3. Update application code to use any new optimized views or functions\n")
        else:
            f.write("## Troubleshooting\n\n")
            f.write("1. Check for syntax errors in the migration file\n")
            f.write("2. Ensure the Supabase project ID is correct\n")
            f.write("3. Check if the database user has sufficient privileges\n")
            f.write("4. Examine any error messages in the logs\n")
            
    logger.info(f"Report generated: {report_path}")
    return report_path

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Apply Optimized RLS Policies Migration via Supabase MCP")
    # Migration file parameters
    parser.add_argument("--security-definer-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/019_rls_security_definer_functions.sql",
                    help="Path to the security definer functions migration file")
    parser.add_argument("--optimized-policies-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/020_optimized_rls_policies.sql",
                    help="Path to the optimized RLS policies migration file")
    
    # Mode parameters
    parser.add_argument("--security-definer-only", action="store_true",
                    help="Only apply the security definer functions")
    parser.add_argument("--optimized-policies-only", action="store_true",
                    help="Only apply the optimized RLS policies")
    
    # Supabase project ID
    parser.add_argument("--project-id", required=True,
                    help="Supabase project ID")
    
    args = parser.parse_args()
    
    # Validate files exist
    if not args.optimized_policies_only and not os.path.exists(args.security_definer_file):
        logger.error(f"Security definer functions migration file not found: {args.security_definer_file}")
        return 1
    
    if not args.security_definer_only and not os.path.exists(args.optimized_policies_file):
        logger.error(f"Optimized RLS policies migration file not found: {args.optimized_policies_file}")
        return 1
    
    # Apply security definer functions migration
    if not args.optimized_policies_only:
        logger.info(f"Applying security definer functions migration: {args.security_definer_file}")
        sql = read_migration_file(args.security_definer_file)
        success = execute_sql_via_mcp(args.project_id, sql)
        if not success:
            logger.error("Failed to apply security definer functions migration")
            generate_report(False, args.security_definer_file, args.project_id)
            return 1
        logger.info("Security definer functions migration applied successfully")
    
    # Apply optimized RLS policies migration
    if not args.security_definer_only:
        logger.info(f"Applying optimized RLS policies migration: {args.optimized_policies_file}")
        sql = read_migration_file(args.optimized_policies_file)
        success = execute_sql_via_mcp(args.project_id, sql)
        if not success:
            logger.error("Failed to apply optimized RLS policies migration")
            generate_report(False, args.optimized_policies_file, args.project_id)
            return 1
        logger.info("Optimized RLS policies migration applied successfully")
    
    # Generate a simple success report
    report_path = generate_report(
        True, 
        args.optimized_policies_file if not args.security_definer_only else args.security_definer_file,
        args.project_id
    )
    
    logger.info(f"RLS optimization complete. Report saved to: {report_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
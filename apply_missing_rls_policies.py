#!/usr/bin/env python
"""
CryoProtect RLS Policy Applier

This script applies comprehensive Row Level Security (RLS) policies to all tables:
1. Core tables (molecules, mixtures, mixture_components, experiments, etc.)
2. Support tables (teams, team_roles, user_team_membership, etc.)
3. Special tables with custom policies (lab_verifications)
4. Views with SECURITY INVOKER

The script is idempotent and safe to run multiple times.
"""

import os
import sys
import time
import json
import logging
import argparse
from datetime import datetime
from dotenv import load_dotenv
import psycopg2
from psycopg2 import sql
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/apply_rls_policies.log', 'a')
    ]
)
logger = logging.getLogger('apply_rls_policies')

def ensure_logs_dir():
    """Ensure logs directory exists."""
    os.makedirs('logs', exist_ok=True)

def ensure_reports_dir():
    """Ensure reports directory exists."""
    os.makedirs('reports/auth', exist_ok=True)

def get_db_connection():
    """Create database connection from environment variables."""
    load_dotenv()  # Load environment variables from .env file
    
    # Get database connection parameters from environment variables
    db_host = os.getenv('SUPABASE_DB_HOST', 'localhost')
    db_port = os.getenv('SUPABASE_DB_PORT', '5432')
    db_name = os.getenv('SUPABASE_DB_NAME', 'postgres')
    db_user = os.getenv('SUPABASE_DB_USER', 'postgres')
    db_pass = os.getenv('SUPABASE_DB_PASSWORD', '')
    
    # Create and return the connection
    conn = psycopg2.connect(
        host=db_host,
        port=db_port,
        dbname=db_name,
        user=db_user,
        password=db_pass
    )
    conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    return conn

def execute_sql_file(conn, file_path):
    """Execute SQL commands from a file."""
    logger.info(f"Executing SQL file: {file_path}")
    
    # Read SQL file
    with open(file_path, 'r') as file:
        sql_commands = file.read()
    
    cursor = conn.cursor()
    try:
        cursor.execute(sql_commands)
        logger.info("SQL execution completed successfully")
    except Exception as e:
        logger.error(f"Error executing SQL: {e}")
        raise
    finally:
        cursor.close()

def get_all_public_tables(conn):
    """Get a list of all tables in the public schema."""
    cursor = conn.cursor()
    try:
        cursor.execute(
            "SELECT tablename FROM pg_tables WHERE schemaname = 'public'"
        )
        tables = [row[0] for row in cursor.fetchall()]
        return tables
    finally:
        cursor.close()

def get_all_public_views(conn):
    """Get a list of all views in the public schema."""
    cursor = conn.cursor()
    try:
        cursor.execute(
            "SELECT viewname FROM pg_views WHERE schemaname = 'public'"
        )
        views = [row[0] for row in cursor.fetchall()]
        return views
    finally:
        cursor.close()

def verify_rls_policies(conn):
    """Verify that RLS policies have been applied correctly."""
    logger.info("Verifying RLS policies...")
    
    # Get all tables and views
    tables = get_all_public_tables(conn)
    views = get_all_public_views(conn)
    
    results = {
        "timestamp": datetime.now().isoformat(),
        "tables": {},
        "views": {}
    }
    
    cursor = conn.cursor()
    
    # Check tables
    for table in tables:
        cursor.execute(
            "SELECT COUNT(*) FROM pg_policies WHERE tablename = %s", 
            (table,)
        )
        policy_count = cursor.fetchone()[0]
        
        cursor.execute(
            "SELECT count(*) FROM pg_tables WHERE tablename = %s AND rowsecurity = true",
            (table,)
        )
        rls_enabled = cursor.fetchone()[0] > 0
        
        results["tables"][table] = {
            "rls_enabled": rls_enabled,
            "policy_count": policy_count,
            "policies": []
        }
        
        if policy_count > 0:
            # Get policy details
            cursor.execute(
                """
                SELECT policyname, permissive, cmd, 
                       pg_catalog.pg_get_expr(qual, tar.oid) AS using_expr,
                       pg_catalog.pg_get_expr(with_check, tar.oid) AS check_expr
                FROM pg_catalog.pg_policies pol
                JOIN pg_catalog.pg_class tar ON pol.tablename = tar.relname
                WHERE pol.tablename = %s
                """, 
                (table,)
            )
            policies = []
            for row in cursor.fetchall():
                policies.append({
                    "name": row[0],
                    "permissive": row[1],
                    "command": row[2],
                    "using_expr": row[3],
                    "check_expr": row[4]
                })
            results["tables"][table]["policies"] = policies
            
            logger.info(f"✓ Table {table} has {policy_count} RLS policies and RLS is {'enabled' if rls_enabled else 'disabled'}")
        else:
            logger.warning(f"✗ Table {table} has no RLS policies")
    
    # Check views
    for view in views:
        cursor.execute(
            """
            SELECT pg_catalog.pg_get_viewdef(c.oid) 
            FROM pg_catalog.pg_class c 
            JOIN pg_catalog.pg_namespace n ON c.relnamespace = n.oid 
            WHERE c.relname = %s 
              AND n.nspname = 'public'
            """, 
            (view,)
        )
        view_def = cursor.fetchone()[0] if cursor.rowcount > 0 else ""
        security_invoker = 'SECURITY INVOKER' in view_def.upper()
        
        results["views"][view] = {
            "security_invoker": security_invoker,
            "definition": view_def
        }
        
        if security_invoker:
            logger.info(f"✓ View {view} has SECURITY INVOKER set")
        else:
            logger.warning(f"✗ View {view} does not have SECURITY INVOKER set")
    
    cursor.close()
    
    # Save verification results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/auth/rls_verification_summary_{timestamp}.json"
    
    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Generate human-readable summary
    generate_verification_summary(results, timestamp)
    
    logger.info(f"RLS verification results saved to {report_path}")
    
    return results

def generate_verification_summary(results, timestamp):
    """Generate a human-readable verification summary."""
    report_path = f"reports/auth/RLS_Verification_Summary_{timestamp}.md"
    
    with open(report_path, 'w') as f:
        f.write("# CryoProtect RLS Verification Summary\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Summary statistics
        table_count = len(results["tables"])
        tables_with_rls = sum(1 for table in results["tables"].values() if table["rls_enabled"])
        tables_with_policies = sum(1 for table in results["tables"].values() if table["policy_count"] > 0)
        
        view_count = len(results["views"])
        views_with_security_invoker = sum(1 for view in results["views"].values() if view["security_invoker"])
        
        f.write("## Summary\n\n")
        f.write(f"- **Tables**: {table_count} total, {tables_with_rls} with RLS enabled, {tables_with_policies} with policies\n")
        f.write(f"- **Views**: {view_count} total, {views_with_security_invoker} with SECURITY INVOKER\n\n")
        
        # Table details
        f.write("## Tables\n\n")
        f.write("| Table | RLS Enabled | Policy Count | Status |\n")
        f.write("|-------|-------------|--------------|--------|\n")
        
        for table_name, table in sorted(results["tables"].items()):
            status = "✅" if table["rls_enabled"] and table["policy_count"] > 0 else "❌"
            f.write(f"| {table_name} | {'Yes' if table['rls_enabled'] else 'No'} | {table['policy_count']} | {status} |\n")
        
        # View details
        f.write("\n## Views\n\n")
        f.write("| View | SECURITY INVOKER | Status |\n")
        f.write("|------|-----------------|--------|\n")
        
        for view_name, view in sorted(results["views"].items()):
            status = "✅" if view["security_invoker"] else "❌"
            f.write(f"| {view_name} | {'Yes' if view['security_invoker'] else 'No'} | {status} |\n")
        
        # Issues and recommendations
        f.write("\n## Issues and Recommendations\n\n")
        
        # Tables without RLS
        tables_without_rls = [name for name, table in results["tables"].items() if not table["rls_enabled"]]
        if tables_without_rls:
            f.write("### Tables without RLS Enabled\n\n")
            for table in tables_without_rls:
                f.write(f"- `{table}`: Enable RLS with `ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;`\n")
        
        # Tables without policies
        tables_without_policies = [name for name, table in results["tables"].items() if table["policy_count"] == 0]
        if tables_without_policies:
            f.write("\n### Tables without RLS Policies\n\n")
            for table in tables_without_policies:
                f.write(f"- `{table}`: Add appropriate RLS policies\n")
        
        # Views without SECURITY INVOKER
        views_without_security_invoker = [name for name, view in results["views"].items() if not view["security_invoker"]]
        if views_without_security_invoker:
            f.write("\n### Views without SECURITY INVOKER\n\n")
            for view in views_without_security_invoker:
                f.write(f"- `{view}`: Recreate with SECURITY INVOKER option\n")
        
        # If no issues
        if not tables_without_rls and not tables_without_policies and not views_without_security_invoker:
            f.write("No issues found. All tables have RLS enabled with appropriate policies, and all views use SECURITY INVOKER.\n")
    
    logger.info(f"Verification summary saved to {report_path}")
    return report_path

def main():
    """Main function to apply RLS policies."""
    parser = argparse.ArgumentParser(description="Apply RLS policies to CryoProtect database")
    parser.add_argument("--skip-verification", action="store_true", help="Skip RLS policy verification")
    parser.add_argument("--sql-file", default="migrations/018_complete_rls_policies.sql", help="SQL file containing RLS policies")
    args = parser.parse_args()
    
    ensure_logs_dir()
    ensure_reports_dir()
    
    logger.info("Starting comprehensive RLS policy application")
    
    try:
        # Get database connection
        conn = get_db_connection()
        logger.info("Successfully connected to the database")
        
        # Execute SQL file with RLS policies
        execute_sql_file(conn, args.sql_file)
        
        # Verify RLS policies
        if not args.skip_verification:
            verify_rls_policies(conn)
        
        conn.close()
        logger.info("RLS policy application completed successfully")
        
    except Exception as e:
        logger.error(f"Error applying RLS policies: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
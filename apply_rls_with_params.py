#!/usr/bin/env python3
"""
Apply RLS Optimizations - With Parameters

This script applies the optimized RLS policies migration to improve
performance of Row Level Security policies in the database.
This version takes database connection parameters directly.
"""

import os
import sys
import argparse
import time
import psycopg2
from psycopg2.extras import RealDictCursor
from pathlib import Path
from typing import Optional, List, Dict, Any

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

def execute_sql(conn_params: Dict[str, str], sql: str) -> bool:
    """
    Execute SQL directly using psycopg2.

    Args:
        conn_params: Database connection parameters
        sql: SQL to execute

    Returns:
        True if successful, False otherwise
    """
    try:
        # Log connection details (without password)
        safe_params = conn_params.copy()
        if 'password' in safe_params:
            safe_params['password'] = '********'
        logger.info(f"Connecting to database with params: {safe_params}")

        # Connect to the database
        conn = psycopg2.connect(**conn_params)
        conn.autocommit = False

        # Create cursor
        cursor = conn.cursor()

        try:
            # Execute the SQL
            logger.info("Executing SQL...")
            cursor.execute(sql)

            # Commit the transaction
            conn.commit()
            logger.info("SQL executed successfully")
            return True
        except Exception as e:
            # Rollback on error
            conn.rollback()
            logger.error(f"Error executing SQL: {str(e)}")
            return False
        finally:
            # Close cursor
            cursor.close()
    except Exception as e:
        logger.error(f"Error connecting to database: {str(e)}")
        return False
    finally:
        # Close connection
        if 'conn' in locals() and conn:
            conn.close()

def generate_report(success: bool, migration_file: str, conn_params: Dict[str, str]) -> str:
    """
    Generate a simple report file.
    
    Args:
        success: Whether the migration was successful
        migration_file: Path to the migration file
        conn_params: Database connection parameters (sensitive info redacted)
        
    Returns:
        Path to the report file
    """
    # Ensure directory exists
    os.makedirs("reports", exist_ok=True)
    
    # Generate timestamp
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    report_path = f"reports/rls_optimization_report_{timestamp}.md"
    
    # Create a sanitized copy of connection params
    safe_params = conn_params.copy()
    if 'password' in safe_params:
        safe_params['password'] = '********'
    
    with open(report_path, 'w') as f:
        f.write("# RLS Optimization Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"Applied migration: {migration_file}\n\n")
        f.write(f"Status: {'✅ Success' if success else '❌ Failed'}\n\n")
        
        f.write("## Connection Details\n\n")
        f.write("```\n")
        for key, value in safe_params.items():
            f.write(f"{key}: {value}\n")
        f.write("```\n\n")
            
    logger.info(f"Report generated: {report_path}")
    return report_path

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Apply Optimized RLS Policies Migration with Parameters")
    # Migration file parameters
    parser.add_argument("--security-definer-file", type=str,
                    default="/home/mushu/Projects/CryoProtect/migrations/simplified_rls_functions.sql",
                    help="Path to the security definer functions migration file")
    parser.add_argument("--optimized-policies-file", type=str,
                    default="/home/mushu/Projects/CryoProtect/migrations/simplified_rls_functions.sql",
                    help="Path to the optimized RLS policies migration file")
    
    # Mode parameters
    parser.add_argument("--security-definer-only", action="store_true",
                    help="Only apply the security definer functions")
    parser.add_argument("--optimized-policies-only", action="store_true",
                    help="Only apply the optimized RLS policies")
    
    # Database connection parameters
    parser.add_argument("--host", default="localhost",
                    help="Database host")
    parser.add_argument("--port", default="5432",
                    help="Database port")
    parser.add_argument("--dbname", default="postgres",
                    help="Database name")
    parser.add_argument("--user", default="postgres",
                    help="Database user")
    parser.add_argument("--password",
                    help="Database password")
    parser.add_argument("--sslmode", default="require",
                    help="SSL mode (disable, allow, prefer, require, verify-ca, verify-full)")
    
    args = parser.parse_args()
    
    # Validate files exist
    if not args.optimized_policies_only and not os.path.exists(args.security_definer_file):
        logger.error(f"Security definer functions migration file not found: {args.security_definer_file}")
        return 1
    
    if not args.security_definer_only and not os.path.exists(args.optimized_policies_file):
        logger.error(f"Optimized RLS policies migration file not found: {args.optimized_policies_file}")
        return 1
    
    # Prepare connection parameters
    conn_params = {
        'host': args.host,
        'port': args.port,
        'dbname': args.dbname,
        'user': args.user,
        'sslmode': args.sslmode  # SSL mode from command line
    }

    # Add password if provided
    if args.password:
        conn_params['password'] = args.password
    
    # Apply security definer functions migration
    if not args.optimized_policies_only:
        logger.info(f"Applying security definer functions migration: {args.security_definer_file}")
        sql = read_migration_file(args.security_definer_file)
        success = execute_sql(conn_params, sql)
        if not success:
            logger.error("Failed to apply security definer functions migration")
            generate_report(False, args.security_definer_file, conn_params)
            return 1
        logger.info("Security definer functions migration applied successfully")
    
    # Apply optimized RLS policies migration
    if not args.security_definer_only:
        logger.info(f"Applying optimized RLS policies migration: {args.optimized_policies_file}")
        sql = read_migration_file(args.optimized_policies_file)
        success = execute_sql(conn_params, sql)
        if not success:
            logger.error("Failed to apply optimized RLS policies migration")
            generate_report(False, args.optimized_policies_file, conn_params)
            return 1
        logger.info("Optimized RLS policies migration applied successfully")
    
    # Generate a simple success report
    report_path = generate_report(
        True, 
        args.optimized_policies_file if not args.security_definer_only else args.security_definer_file,
        conn_params
    )
    
    logger.info(f"RLS optimization complete. Report saved to: {report_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
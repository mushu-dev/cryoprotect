#!/usr/bin/env python3
"""
Apply RLS Optimizations - Simple Version

This script applies the optimized RLS policies migration to improve
performance of Row Level Security policies in the database.
This is a simplified version that doesn't rely on connection pools.
"""

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path
from typing import Optional, List, Dict, Any

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import from the new configuration system
try:
    from database.connection_config import (
        validate_config,
        get_connection_config,
        test_adapter_configuration
    )
    HAS_CENTRAL_CONFIG = True
except ImportError:
    HAS_CENTRAL_CONFIG = False

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

def apply_migration(migration_file: str) -> bool:
    """
    Apply a migration file using psql directly.
    
    Args:
        migration_file: Path to migration file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # First try to use the new configuration system if available
        if HAS_CENTRAL_CONFIG:
            try:
                # Validate the configuration
                validate_config()
                
                # Get configuration for local adapter (for psql connection)
                config = get_connection_config('local')
                
                # Extract connection details
                dbname = config.get('database', 'postgres')
                user = config.get('user', 'postgres')
                host = config.get('host', 'localhost')
                port = str(config.get('port', '5432'))
                password = config.get('password', '')
                
                logger.info("Using database configuration from central configuration system")
            except Exception as e:
                logger.warning(f"Error getting database configuration from central config: {str(e)}")
                logger.info("Falling back to environment variables")
                # Fall back to environment variables
                dbname = os.environ.get('PGDATABASE', 'postgres')
                user = os.environ.get('PGUSER', 'postgres')
                host = os.environ.get('PGHOST', 'localhost')
                port = os.environ.get('PGPORT', '5432')
                password = os.environ.get('PGPASSWORD', '')
        else:
            # Get database connection details from environment variables
            dbname = os.environ.get('PGDATABASE', 'postgres')
            user = os.environ.get('PGUSER', 'postgres')
            host = os.environ.get('PGHOST', 'localhost')
            port = os.environ.get('PGPORT', '5432')
            password = os.environ.get('PGPASSWORD', '')
        
        # Set PGPASSWORD environment variable for psql
        env = os.environ.copy()
        if password:
            env['PGPASSWORD'] = password
            
        # Use psql to apply migration
        cmd = [
            "psql",
            "-h", host,
            "-p", port,
            "-U", user,
            "-d", dbname,
            "-f", migration_file
        ]
        
        logger.info(f"Applying migration: {migration_file}")
        process = subprocess.run(
            cmd, 
            env=env, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True
        )
        
        if process.returncode == 0:
            logger.info(f"Migration applied successfully: {migration_file}")
            return True
        else:
            logger.error(f"Migration failed: {process.stderr}")
            return False
            
    except Exception as e:
        logger.error(f"Error applying migration: {str(e)}")
        return False

def generate_report(success: bool, migration_file: str, output: str) -> str:
    """
    Generate a simple report file.
    
    Args:
        success: Whether the migration was successful
        migration_file: Path to the migration file
        output: Command output to include in the report
        
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
        f.write(f"Status: {'✅ Success' if success else '❌ Failed'}\n\n")
        
        if output:
            f.write("## Command Output\n\n")
            f.write("```\n")
            f.write(output)
            f.write("\n```\n\n")
            
    logger.info(f"Report generated: {report_path}")
    return report_path

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Apply Optimized RLS Policies Migration - Simple Version")
    parser.add_argument("--security-definer-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/019_rls_security_definer_functions.sql",
                    help="Path to the security definer functions migration file")
    parser.add_argument("--optimized-policies-file", type=str, 
                    default="/home/mushu/Projects/CryoProtect/migrations/020_optimized_rls_policies.sql",
                    help="Path to the optimized RLS policies migration file")
    parser.add_argument("--security-definer-only", action="store_true",
                    help="Only apply the security definer functions")
    parser.add_argument("--optimized-policies-only", action="store_true",
                    help="Only apply the optimized RLS policies")
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
        success = apply_migration(args.security_definer_file)
        if not success:
            logger.error("Failed to apply security definer functions migration")
            return 1
        logger.info("Security definer functions migration applied successfully")
    
    # Apply optimized RLS policies migration
    if not args.security_definer_only:
        logger.info(f"Applying optimized RLS policies migration: {args.optimized_policies_file}")
        success = apply_migration(args.optimized_policies_file)
        if not success:
            logger.error("Failed to apply optimized RLS policies migration")
            return 1
        logger.info("Optimized RLS policies migration applied successfully")
    
    # Generate a simple success report
    report_path = generate_report(
        True, 
        args.optimized_policies_file if not args.security_definer_only else args.security_definer_file,
        "Migrations applied successfully"
    )
    
    logger.info(f"RLS optimization complete. Report saved to: {report_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
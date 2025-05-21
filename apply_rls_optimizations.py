#!/usr/bin/env python3
"""
Apply RLS Optimizations Migration

This script applies the RLS security definer functions migration to improve
performance of Row Level Security policies in the database.
"""

import os
import sys
import logging
import argparse
import time
from pathlib import Path
from typing import Optional, List, Dict, Any

# Import optimized connection pool
from optimized_connection_pool import (
    execute_query_with_retry,
    transaction_context,
    ConnectionManager
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
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

def execute_migration(sql: str) -> bool:
    """
    Execute the migration SQL.
    
    Args:
        sql: SQL migration to execute
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with transaction_context() as conn:
            with conn.cursor() as cursor:
                cursor.execute(sql)
            
        logger.info("Migration executed successfully")
        return True
    except Exception as e:
        logger.error(f"Error executing migration: {str(e)}")
        return False

def verify_migration() -> Dict[str, Any]:
    """
    Verify the migration was applied correctly by checking for the new functions.
    
    Returns:
        Dictionary with verification results
    """
    results = {
        'success': False,
        'functions_found': 0,
        'functions_expected': 8,
        'function_list': []
    }
    
    # List of functions to check for
    functions_to_check = [
        'auth.check_resource_access',
        'auth.check_resource_modify_access',
        'auth.is_service_role_user',
        'auth.check_verification_access',
        'auth.check_verification_modify_access',
        'auth.check_batch_access',
        'auth.get_user_permissions',
        'auth.user_has_permission'
    ]
    
    try:
        # Query to check for the functions
        sql = """
        SELECT 
            n.nspname as schema,
            p.proname as name,
            pg_get_functiondef(p.oid) as definition
        FROM 
            pg_proc p
            JOIN pg_namespace n ON p.pronamespace = n.oid
        WHERE 
            n.nspname = 'auth' AND 
            p.proname IN (
                'check_resource_access',
                'check_resource_modify_access',
                'is_service_role_user',
                'check_verification_access',
                'check_verification_modify_access',
                'check_batch_access',
                'get_user_permissions',
                'user_has_permission'
            );
        """
        
        result = execute_query_with_retry(sql)
        
        if result:
            # Process results
            for row in result:
                function_name = f"{row['schema']}.{row['name']}"
                results['function_list'].append(function_name)
            
            results['functions_found'] = len(results['function_list'])
            results['success'] = results['functions_found'] == results['functions_expected']
            
            logger.info(f"Found {results['functions_found']} of {results['functions_expected']} expected functions")
            
            # Log specific results
            for function in functions_to_check:
                if function in results['function_list']:
                    logger.info(f"✓ Function {function} was created successfully")
                else:
                    logger.warning(f"✗ Function {function} was not found!")
        else:
            logger.error("Failed to query database for functions")
    
    except Exception as e:
        logger.error(f"Error verifying migration: {str(e)}")
    
    return results

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Apply RLS Optimizations Migration")
    parser.add_argument("--migration-file", type=str, 
                        default="/home/mushu/Projects/CryoProtect/migrations/019_rls_security_definer_functions.sql",
                        help="Path to the migration file")
    parser.add_argument("--verify-only", action="store_true", 
                        help="Only verify the migration (don't execute it)")
    args = parser.parse_args()
    
    migration_file = args.migration_file
    
    logger.info(f"Using migration file: {migration_file}")
    
    if not os.path.exists(migration_file):
        logger.error(f"Migration file not found: {migration_file}")
        return 1
    
    # Verify only mode
    if args.verify_only:
        logger.info("Verify-only mode: Checking if migration has been applied")
        results = verify_migration()
        
        if results['success']:
            logger.info("Migration verification successful: All functions found")
            return 0
        else:
            logger.warning(f"Migration verification failed: Found {results['functions_found']} "
                         f"of {results['functions_expected']} expected functions")
            return 1
    
    # Execute migration mode
    try:
        # Read migration file
        sql = read_migration_file(migration_file)
        
        # Execute migration
        logger.info("Executing migration...")
        success = execute_migration(sql)
        
        if not success:
            logger.error("Migration failed")
            return 1
        
        # Verify migration
        logger.info("Verifying migration...")
        results = verify_migration()
        
        if results['success']:
            logger.info("Migration verification successful: All functions created")
            return 0
        else:
            logger.warning(f"Migration verification partial: Found {results['functions_found']} "
                         f"of {results['functions_expected']} expected functions")
            return 1
        
    except Exception as e:
        logger.error(f"Error in migration process: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
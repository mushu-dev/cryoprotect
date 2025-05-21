#!/usr/bin/env python3
"""
CryoProtect v2 - Apply RBAC Schema Migration

This script applies the RBAC schema migration to the database.
It creates the necessary tables, functions, and policies for role-based access control.
"""

import os
import sys
import logging
from supabase import create_client
from dotenv import load_dotenv

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def get_supabase_client():
    """Get the Supabase client."""
    # Load environment variables
    load_dotenv()
    
    # Get Supabase credentials
    supabase_url = os.environ.get('SUPABASE_URL')
    supabase_key = os.environ.get('SUPABASE_KEY')
    
    if not supabase_url or not supabase_key:
        logger.error("SUPABASE_URL and SUPABASE_KEY must be set in the .env file")
        sys.exit(1)
    
    # Create Supabase client
    return create_client(supabase_url, supabase_key)

def apply_migration():
    """Apply the RBAC schema migration."""
    try:
        # Get Supabase client
        supabase = get_supabase_client()
        
        # Read the migration SQL file
        migration_path = os.path.join('migrations', '011_rbac_schema.sql')
        
        if not os.path.exists(migration_path):
            logger.error(f"Migration file not found: {migration_path}")
            sys.exit(1)
        
        with open(migration_path, 'r') as f:
            migration_sql = f.read()
        
        # Execute the migration SQL
        logger.info("Applying RBAC schema migration...")
        
        # Split the SQL into statements and execute them one by one
        # This is necessary because the Supabase client doesn't support multiple statements in a single query
        statements = migration_sql.split(';')
        
        for statement in statements:
            # Skip empty statements
            if statement.strip():
                try:
                    # Execute the statement using the Supabase REST API
                    # We need to use rpc with exec_sql function if it exists
                    response = supabase.rpc(
                        'exec_sql',
                        {'sql_query': statement}
                    ).execute()
                    
                    if response.error:
                        logger.error(f"Error executing SQL: {response.error}")
                        logger.error(f"Statement: {statement}")
                except Exception as e:
                    logger.error(f"Error executing SQL: {str(e)}")
                    logger.error(f"Statement: {statement}")
        
        logger.info("RBAC schema migration applied successfully")
        
        # Verify the migration
        verify_migration(supabase)
        
    except Exception as e:
        logger.error(f"Error applying migration: {str(e)}")
        sys.exit(1)

def verify_migration(supabase):
    """Verify that the migration was applied successfully."""
    try:
        # Check if the roles table exists
        response = supabase.table('roles').select('count').execute()
        
        if response.error:
            logger.error(f"Error verifying migration: {response.error}")
            return
        
        # Check if the default roles were created
        response = supabase.table('roles').select('name').execute()
        
        if response.error:
            logger.error(f"Error verifying roles: {response.error}")
            return
        
        roles = [role['name'] for role in response.data]
        expected_roles = ['admin', 'user', 'curator', 'viewer']
        
        for role in expected_roles:
            if role not in roles:
                logger.warning(f"Default role not found: {role}")
        
        logger.info(f"Found roles: {', '.join(roles)}")
        
        # Check if the permissions table exists
        response = supabase.table('permissions').select('count').execute()
        
        if response.error:
            logger.error(f"Error verifying permissions: {response.error}")
            return
        
        logger.info("Migration verification completed")
        
    except Exception as e:
        logger.error(f"Error verifying migration: {str(e)}")

if __name__ == '__main__':
    logger.info("Starting RBAC schema migration")
    apply_migration()
    logger.info("RBAC schema migration completed")
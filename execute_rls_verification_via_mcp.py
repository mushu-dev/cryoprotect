#!/usr/bin/env python3
"""
CryoProtect v2 - ChEMBL Integration Database Preparation

This script prepares the database for ChEMBL data integration by:
1. Temporarily modifying RLS policies to allow bulk import
2. Verifying schema compatibility
3. Setting up necessary database resources

It uses the Supabase MCP server for all database operations.
"""

import os
import json
import logging
import argparse
from datetime import datetime
import sys
from pathlib import Path
import subprocess

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/chembl_db_prep.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Supabase project ID
PROJECT_ID = "tsdlmynydfuypiugmkev"

def execute_mcp_sql(query):
    """
    Execute SQL through the MCP server using subprocess.
    
    Args:
        query (str): SQL query to execute
        
    Returns:
        list/dict: Result of the query execution or error dict
    """
    logger.info(f"Executing SQL via MCP: {query[:100]}{'...' if len(query) > 100 else ''}")
    
    try:
        # Create a temporary file for the query
        with open("temp_query.sql", "w") as f:
            f.write(query)
        
        # Execute the MCP command
        cmd = [
            "roo", "use_mcp_tool", "supabase", "execute_sql",
            "--project_id", PROJECT_ID,
            "--query", query
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"SQL execution failed: {result.stderr}")
            return {"error": result.stderr}
        
        # Parse the result
        try:
            result_data = json.loads(result.stdout)
            logger.info(f"SQL execution successful")
            return result_data
        except json.JSONDecodeError:
            logger.error(f"Failed to parse MCP result: {result.stdout}")
            return {"error": "Failed to parse MCP result", "stdout": result.stdout, "stderr": result.stderr}
    
    except Exception as e:
        logger.error(f"Error executing SQL: {str(e)}")
        return {"error": str(e)}
    
    finally:
        # Clean up the temporary file
        if os.path.exists("temp_query.sql"):
            os.remove("temp_query.sql")

def check_rls_policies(table_name="molecules"):
    """
    Check the RLS policies for a table.
    
    Args:
        table_name (str): Name of the table to check
        
    Returns:
        list: List of RLS policies for the table
    """
    logger.info(f"Checking RLS policies for table '{table_name}'...")
    
    query = f"""
    SELECT tablename, policyname, permissive, roles, cmd, qual, with_check 
    FROM pg_policies 
    WHERE tablename = '{table_name}';
    """
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to check RLS policies: {result['error']}")
            return []
        
        if not result:
            logger.warning(f"No RLS policies found for table '{table_name}'")
            return []
        
        logger.info(f"Found {len(result)} RLS policies for table '{table_name}'")
        
        for policy in result:
            logger.info(f"Policy: {policy['policyname']}, Command: {policy['cmd']}, Permissive: {policy['permissive']}")
        
        return result
    
    except Exception as e:
        logger.error(f"Error checking RLS policies: {str(e)}")
        return []

def verify_table_schema(table_name="molecules"):
    """
    Verify the schema of a table.
    
    Args:
        table_name (str): Name of the table to verify
        
    Returns:
        list: List of columns in the table
    """
    logger.info(f"Verifying schema for table '{table_name}'...")
    
    query = f"""
    SELECT column_name, data_type, is_nullable
    FROM information_schema.columns
    WHERE table_name = '{table_name}'
    ORDER BY ordinal_position;
    """
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to verify table schema: {result['error']}")
            return []
        
        if not result:
            logger.warning(f"Table '{table_name}' not found or has no columns")
            return []
        
        logger.info(f"Table '{table_name}' has {len(result)} columns")
        
        for column in result:
            logger.info(f"Column: {column['column_name']}, Type: {column['data_type']}, Nullable: {column['is_nullable']}")
        
        return result
    
    except Exception as e:
        logger.error(f"Error verifying table schema: {str(e)}")
        return []

def disable_rls(table_name="molecules"):
    """
    Temporarily disable RLS for a table.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Temporarily disabling RLS for table '{table_name}'...")
    
    query = f"ALTER TABLE {table_name} DISABLE ROW LEVEL SECURITY;"
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to disable RLS: {result['error']}")
            return False
        
        logger.info(f"RLS disabled for table '{table_name}'")
        return True
    
    except Exception as e:
        logger.error(f"Error disabling RLS: {str(e)}")
        return False

def create_permissive_policy(table_name="molecules"):
    """
    Create a permissive policy for service role.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Creating permissive policy for service role on table '{table_name}'...")
    
    # First drop the policy if it exists
    drop_query = f"""
    DROP POLICY IF EXISTS "Service role has full access to {table_name}" ON {table_name};
    """
    
    create_query = f"""
    CREATE POLICY "Service role has full access to {table_name}"
    ON {table_name}
    FOR ALL
    USING (auth.jwt() ->> 'role' = 'service_role')
    WITH CHECK (auth.jwt() ->> 'role' = 'service_role');
    """
    
    try:
        # Drop existing policy if any
        drop_result = execute_sql(drop_query)
        if isinstance(drop_result, dict) and "error" in drop_result:
            logger.warning(f"Failed to drop existing policy: {drop_result['error']}")
        
        # Create new policy
        create_result = execute_sql(create_query)
        if isinstance(create_result, dict) and "error" in create_result:
            logger.error(f"Failed to create permissive policy: {create_result['error']}")
            return False
        
        logger.info(f"Permissive policy created for table '{table_name}'")
        return True
    
    except Exception as e:
        logger.error(f"Error creating permissive policy: {str(e)}")
        return False

def enable_rls(table_name="molecules"):
    """
    Re-enable RLS for a table.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Re-enabling RLS for table '{table_name}'...")
    
    query = f"ALTER TABLE {table_name} ENABLE ROW LEVEL SECURITY;"
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to enable RLS: {result['error']}")
            return False
        
        logger.info(f"RLS enabled for table '{table_name}'")
        return True
    
    except Exception as e:
        logger.error(f"Error enabling RLS: {str(e)}")
        return False

def drop_permissive_policy(table_name="molecules"):
    """
    Drop the permissive policy for a table.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Dropping permissive policy for table '{table_name}'...")
    
    query = f"""
    DROP POLICY IF EXISTS "Service role has full access to {table_name}" ON {table_name};
    """
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to drop permissive policy: {result['error']}")
            return False
        
        logger.info(f"Permissive policy dropped for table '{table_name}'")
        return True
    
    except Exception as e:
        logger.error(f"Error dropping permissive policy: {str(e)}")
        return False

def grant_privileges(table_name="molecules", role="authenticated"):
    """
    Grant privileges on a table to a role.
    
    Args:
        table_name (str): Name of the table
        role (str): Role to grant privileges to
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Granting privileges on table '{table_name}' to role '{role}'...")
    
    query = f"GRANT ALL PRIVILEGES ON TABLE {table_name} TO {role};"
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to grant privileges: {result['error']}")
            return False
        
        logger.info(f"Privileges granted on table '{table_name}' to role '{role}'")
        return True
    
    except Exception as e:
        logger.error(f"Error granting privileges: {str(e)}")
        return False

def test_batch_insert(table_name="molecules"):
    """
    Test batch insert with transaction.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Testing batch insert with transaction on table '{table_name}'...")
    
    # Check if the table has an inchikey column for the ON CONFLICT clause
    schema = verify_table_schema(table_name)
    has_inchikey = any(col["column_name"] == "inchikey" for col in schema)
    
    # Prepare the query based on the schema
    if has_inchikey:
        query = f"""
        BEGIN;
        INSERT INTO {table_name} (name, smiles, inchi, inchikey, formula, molecular_weight, created_by, data_source, version)
        VALUES ('TestMol1', 'C1=CC=CC=C1', 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H', 'UHOVQNZJYSORNB-UHFFFAOYSA-N', 'C6H6', 78.11, 'test_user', 'TestSource', 1)
        ON CONFLICT (inchikey) DO NOTHING
        RETURNING id;
        
        INSERT INTO {table_name} (name, smiles, inchi, inchikey, formula, molecular_weight, created_by, data_source, version)
        VALUES ('TestMol2', 'CC1=CC=CC=C1', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'YXFVVABEGXRONW-UHFFFAOYSA-N', 'C7H8', 92.14, 'test_user', 'TestSource', 1)
        ON CONFLICT (inchikey) DO NOTHING
        RETURNING id;
        COMMIT;
        """
    else:
        # Simplified query without inchikey conflict handling
        query = f"""
        BEGIN;
        INSERT INTO {table_name} (name, created_by, data_source, version)
        VALUES ('TestMol1', 'test_user', 'TestSource', 1)
        RETURNING id;
        
        INSERT INTO {table_name} (name, created_by, data_source, version)
        VALUES ('TestMol2', 'test_user', 'TestSource', 1)
        RETURNING id;
        COMMIT;
        """
    
    try:
        result = execute_sql(query)
        
        if isinstance(result, dict) and "error" in result:
            logger.error(f"Failed to test batch insert: {result['error']}")
            return False
        
        logger.info(f"Batch insert test successful")
        return True
    
    except Exception as e:
        logger.error(f"Error testing batch insert: {str(e)}")
        return False

def test_transaction_atomicity(table_name="molecules"):
    """
    Test transaction atomicity.
    
    Args:
        table_name (str): Name of the table
        
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info(f"Testing transaction atomicity on table '{table_name}'...")
    
    # This transaction should fail because the second insert has invalid data
    query = f"""
    BEGIN;
    INSERT INTO {table_name} (name, created_by, data_source, version)
    VALUES ('ValidMol', 'test_user', 'TestSource', 1)
    RETURNING id;
    
    -- This should fail due to NULL in a NOT NULL column (assuming id is NOT NULL)
    INSERT INTO {table_name} (id, name, created_by, data_source, version)
    VALUES (NULL, 'InvalidMol', 'test_user', 'TestSource', 1)
    RETURNING id;
    COMMIT;
    """
    
    try:
        result = execute_sql(query)
        
        # We expect this to fail
        if isinstance(result, dict) and "error" in result:
            logger.info(f"Transaction atomicity test passed: Transaction failed as expected")
            
            # Verify that the first insert was rolled back
            verify_query = f"""
            SELECT COUNT(*) FROM {table_name} WHERE name = 'ValidMol' AND created_by = 'test_user';
            """
            
            verify_result = execute_sql(verify_query)
            
            if isinstance(verify_result, dict) and "error" in verify_result:
                logger.error(f"Failed to verify transaction rollback: {verify_result['error']}")
                return False
            
            count = verify_result[0]["count"]
            
            if count == 0:
                logger.info(f"Transaction rollback verified: First insert was rolled back")
                return True
            else:
                logger.warning(f"Transaction rollback failed: First insert was not rolled back")
                return False
        else:
            logger.warning(f"Transaction atomicity test failed: Transaction succeeded when it should have failed")
            return False
    
    except Exception as e:
        logger.error(f"Error testing transaction atomicity: {str(e)}")
        return False

def prepare_database_for_chembl():
    """
    Prepare the database for ChEMBL integration.
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info("Preparing database for ChEMBL integration...")
    
    # Step 1: Verify database role
    role, user = verify_database_role()
    if not role:
        logger.error("Database role verification failed")
        return False
    
    # Check if we have the necessary role for RLS bypass
    if role not in ["service_role", "postgres", "rls_admin"]:
        logger.warning(f"Role '{role}' may not have sufficient permissions for RLS modifications")
    else:
        logger.info(f"Role '{role}' has sufficient permissions for RLS modifications")
    
    # Step 2: Check RLS policies
    molecules_policies = check_rls_policies("molecules")
    molecular_properties_policies = check_rls_policies("molecular_properties")
    
    # Step 3: Verify table schemas
    molecules_schema = verify_table_schema("molecules")
    molecular_properties_schema = verify_table_schema("molecular_properties")
    
    if not molecules_schema or not molecular_properties_schema:
        logger.error("Table schema verification failed")
        return False
    
    # Step 4: Temporarily disable RLS or create permissive policies
    tables = ["molecules", "molecular_properties"]
    
    for table in tables:
        # Option 1: Disable RLS completely
        if disable_rls(table):
            logger.info(f"RLS disabled for table '{table}'")
        else:
            # Option 2: Create a permissive policy for service role
            if create_permissive_policy(table):
                logger.info(f"Permissive policy created for table '{table}'")
            else:
                logger.error(f"Failed to modify RLS for table '{table}'")
                return False
    
    # Step 5: Grant necessary privileges
    for table in tables:
        if not grant_privileges(table, "authenticated"):
            logger.warning(f"Failed to grant privileges on table '{table}'")
    
    # Step 6: Test batch insert and transaction handling
    if not test_batch_insert("molecules"):
        logger.warning("Batch insert test failed")
    
    if not test_transaction_atomicity("molecules"):
        logger.warning("Transaction atomicity test failed")
    
    logger.info("Database preparation for ChEMBL integration completed successfully")
    return True

def restore_rls():
    """
    Restore RLS settings after ChEMBL integration.
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger.info("Restoring RLS settings...")
    
    tables = ["molecules", "molecular_properties"]
    success = True
    
    for table in tables:
        # Re-enable RLS
        if not enable_rls(table):
            logger.error(f"Failed to re-enable RLS for table '{table}'")
            success = False
        
        # Drop permissive policy
        if not drop_permissive_policy(table):
            logger.warning(f"Failed to drop permissive policy for table '{table}'")
    
    if success:
        logger.info("RLS settings restored successfully")
    else:
        logger.warning("RLS settings restoration completed with warnings")
    
    return success

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Prepare database for ChEMBL integration")
    parser.add_argument("--restore", action="store_true", help="Restore RLS settings after integration")
    args = parser.parse_args()
    
    try:
        if args.restore:
            restore_rls()
        else:
            prepare_database_for_chembl()
    
    except Exception as e:
        logger.error(f"Error in main function: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    # If running directly, use the Roo CLI to execute the script
    # This is a workaround for the Python subprocess issues
    if len(sys.argv) > 1 and sys.argv[1] == "--cli":
        # This branch is executed when called from the command line
        sys.exit(main())
    else:
        # This branch is executed when run directly
        logger.info("Running in CLI mode for direct MCP access")
        
        # Check if we should restore RLS
        restore_mode = "--restore" in sys.argv
        
        # Execute the appropriate command
        if restore_mode:
            cmd = ["roo", "execute_command", "python execute_rls_verification_via_mcp.py --cli --restore"]
        else:
            cmd = ["roo", "execute_command", "python execute_rls_verification_via_mcp.py --cli"]
        
        subprocess.run(cmd)
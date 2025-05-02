#!/usr/bin/env python3
"""
supabase_mcp_tools.py - Helper functions for interacting with Supabase via MCP

This module provides utility functions to execute SQL queries and other operations
on a Supabase project using the MCP tools.
"""

import json
import logging
import subprocess
from typing import Dict, List, Tuple, Union, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_sql_on_supabase(project_id: str, query: str) -> List[Dict[str, Any]]:
    """
    Execute a SQL query on a Supabase project using MCP tools
    
    Args:
        project_id: The Supabase project ID
        query: The SQL query to execute
        
    Returns:
        The query results as a list of dictionaries
        
    Raises:
        Exception: If the query execution fails
    """
    try:
        # Create a temporary JSON file with the MCP tool arguments
        args = {
            "project_id": project_id,
            "query": query
        }
        
        # Use the MCP tool to execute the query
        result = subprocess.run(
            ["npx", "-y", "@supabase/mcp-server-supabase@latest", "execute_sql", 
             "--access-token", "sbp_2ef753d5e351cd40412c70c7ed9852c59f18a559",
             "--args", json.dumps(args)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse the result
        try:
            result_json = json.loads(result.stdout)
            return result_json
        except json.JSONDecodeError:
            logger.error(f"Error parsing result: {result.stdout}")
            raise Exception(f"Failed to parse query result: {result.stdout}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing query: {e.stderr}")
        raise Exception(f"Failed to execute query: {e.stderr}")
    
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        raise

def get_project_details(project_id: str) -> Dict[str, Any]:
    """
    Get details about a Supabase project
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        Project details as a dictionary
    """
    try:
        result = subprocess.run(
            ["npx", "-y", "@supabase/mcp-server-supabase@latest", "get_project", 
             "--access-token", "sbp_2ef753d5e351cd40412c70c7ed9852c59f18a559",
             "--args", json.dumps({"id": project_id})],
            capture_output=True,
            text=True,
            check=True
        )
        
        return json.loads(result.stdout)
    
    except Exception as e:
        logger.error(f"Error getting project details: {str(e)}")
        raise Exception(f"Failed to get project details: {str(e)}")

def list_tables(project_id: str, schemas: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    List all tables in a Supabase project
    
    Args:
        project_id: The Supabase project ID
        schemas: Optional list of schemas to include
        
    Returns:
        List of tables with their details
    """
    try:
        args = {
            "project_id": project_id
        }
        
        if schemas:
            args["schemas"] = schemas
        
        result = subprocess.run(
            ["npx", "-y", "@supabase/mcp-server-supabase@latest", "list_tables", 
             "--access-token", "sbp_2ef753d5e351cd40412c70c7ed9852c59f18a559",
             "--args", json.dumps(args)],
            capture_output=True,
            text=True,
            check=True
        )
        
        return json.loads(result.stdout)
    
    except Exception as e:
        logger.error(f"Error listing tables: {str(e)}")
        raise Exception(f"Failed to list tables: {str(e)}")

def apply_migration(project_id: str, name: str, query: str) -> Dict[str, Any]:
    """
    Apply a migration to the database
    
    Args:
        project_id: The Supabase project ID
        name: The name of the migration
        query: The SQL query to execute
        
    Returns:
        Migration result details
    """
    try:
        args = {
            "project_id": project_id,
            "name": name,
            "query": query
        }
        
        result = subprocess.run(
            ["npx", "-y", "@supabase/mcp-server-supabase@latest", "apply_migration", 
             "--access-token", "sbp_2ef753d5e351cd40412c70c7ed9852c59f18a559",
             "--args", json.dumps(args)],
            capture_output=True,
            text=True,
            check=True
        )
        
        return json.loads(result.stdout)
    
    except Exception as e:
        logger.error(f"Error applying migration: {str(e)}")
        raise Exception(f"Failed to apply migration: {str(e)}")

def backup_table(project_id: str, table_name: str) -> str:
    """
    Create a backup of a table before modifying it
    
    Args:
        project_id: The Supabase project ID
        table_name: The name of the table to backup
        
    Returns:
        The name of the backup table
    """
    import time
    
    backup_table_name = f"{table_name}_backup_{int(time.time())}"
    
    query = f"""
    CREATE TABLE public.{backup_table_name} AS 
    SELECT * FROM public.{table_name};
    """
    
    try:
        execute_sql_on_supabase(project_id, query)
        logger.info(f"Created backup of {table_name} as {backup_table_name}")
        return backup_table_name
    except Exception as e:
        logger.error(f"Failed to create backup of {table_name}: {str(e)}")
        raise

def restore_from_backup(project_id: str, original_table: str, backup_table: str) -> bool:
    """
    Restore a table from its backup
    
    Args:
        project_id: The Supabase project ID
        original_table: The name of the original table
        backup_table: The name of the backup table
        
    Returns:
        True if successful, False otherwise
    """
    # First drop the original table
    drop_query = f"DROP TABLE IF EXISTS public.{original_table};"
    
    # Then recreate it from the backup
    restore_query = f"""
    CREATE TABLE public.{original_table} AS 
    SELECT * FROM public.{backup_table};
    """
    
    try:
        execute_sql_on_supabase(project_id, drop_query)
        execute_sql_on_supabase(project_id, restore_query)
        logger.info(f"Restored {original_table} from {backup_table}")
        return True
    except Exception as e:
        logger.error(f"Failed to restore {original_table} from {backup_table}: {str(e)}")
        return False

if __name__ == "__main__":
    # Example usage
    PROJECT_ID = "tsdlmynydfuypiugmkev"
    
    try:
        # Get project details
        project = get_project_details(PROJECT_ID)
        print(f"Project: {project['name']} (Region: {project['region']})")
        
        # List tables
        tables = list_tables(PROJECT_ID, ["public"])
        print(f"Found {len(tables)} tables in the public schema")
        
        # Example query
        result = execute_sql_on_supabase(PROJECT_ID, "SELECT COUNT(*) FROM public.molecules;")
        print(f"Query result: {result}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
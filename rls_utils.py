#!/usr/bin/env python3
"""
CryoProtect v2 - RLS Verification and Remediation Utilities

This module provides utilities for verifying and remediating Row Level Security (RLS)
settings for the molecule and molecular_property tables. It ensures that RLS is always
enabled and all required policies are present after any import or DDL operation.

Usage:
    from rls_utils import verify_rls, remediate_rls, ensure_rls_restored

    # Verify RLS settings
    issues = verify_rls()
    if issues:
        # Remediate issues
        remediate_rls(issues)

    # Or use the decorator to automatically verify and remediate
    @ensure_rls_restored
    def import_compounds_to_database(...):
        ...
"""

import os
import json
import logging
import functools
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional, Callable

# Import the MCP tool helper for SQL execution
from use_mcp_tool import execute_sql
import os
from dotenv import load_dotenv
from supabase import create_client, Client

# Load environment variables
load_dotenv()

# Initialize Supabase client
supabase_client = None

def get_supabase_client() -> Client:
    """
    Get a Supabase client using the service role key.
    
    Returns:
        Client: Supabase client with service role permissions
    """
    global supabase_client
    
    if supabase_client is not None:
        return supabase_client
        
    try:
        # Get Supabase credentials from environment
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY") or os.getenv("SUPABASE_SERVICE_KEY") or os.getenv("SUPABASE_SERVICE_ROLE_KEY")
        
        if not supabase_url or not supabase_key:
            raise Exception("Missing Supabase credentials in .env file")
        
        # Create Supabase client
        supabase_client = create_client(supabase_url, supabase_key)
        return supabase_client
            
    except Exception as e:
        logger.error(f"Error creating Supabase client: {str(e)}")
        raise

def execute_ddl(query: str) -> bool:
    """
    Execute a DDL query that doesn't return results.
    
    Args:
        query: The DDL query to execute
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        client = get_supabase_client()
        
        # Use the rpc function with a different approach for DDL
        # The exec_sql_ddl function should be created in the database
        # to handle DDL statements that don't return results
        result = client.rpc('exec_sql', {'query': query}).execute()
        
        # For DDL statements, we don't expect data to be returned
        # So we consider it successful if no error occurred
        return True
            
    except Exception as e:
        logger.error(f"Error executing DDL: {str(e)}")
        return False

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/rls_audit.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Define the tables to verify
TABLES_TO_VERIFY = ["molecule", "molecular_property"]

# Define the canonical RLS policies from migration files
# These are the policies that should be present on the tables
CANONICAL_POLICIES = {
    "molecule": {
        "Select molecules for project members": {
            "operation": "SELECT",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM project
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE project.id = molecule.project_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to view molecules in their projects."
        },
        "Insert molecules for project members": {
            "operation": "INSERT",
            "using_clause": None,
            "with_check_clause": """
                EXISTS (
                  SELECT 1 FROM project
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE project.id = molecule.project_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "comment": "Allows project members to add molecules to their projects."
        },
        "Update molecules for project members": {
            "operation": "UPDATE",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM project
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE project.id = molecule.project_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to update molecules in their projects."
        },
        "Delete molecules for project members": {
            "operation": "DELETE",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM project
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE project.id = molecule.project_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to delete molecules in their projects."
        },
        "Allow service role inserts on molecule": {
            "operation": "INSERT",
            "using_clause": None,
            "with_check_clause": "auth.role() = 'service_role'",
            "comment": "Allows service role to insert molecules regardless of project membership."
        }
    },
    "molecular_property": {
        "Select molecular_properties for project members": {
            "operation": "SELECT",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM molecule
                  JOIN project ON project.id = molecule.project_id
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE molecule.id = molecular_property.molecule_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to view molecular properties in their projects."
        },
        "Insert molecular_properties for project members": {
            "operation": "INSERT",
            "using_clause": None,
            "with_check_clause": """
                EXISTS (
                  SELECT 1 FROM molecule
                  JOIN project ON project.id = molecule.project_id
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE molecule.id = molecular_property.molecule_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "comment": "Allows project members to add molecular properties to their projects."
        },
        "Update molecular_properties for project members": {
            "operation": "UPDATE",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM molecule
                  JOIN project ON project.id = molecule.project_id
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE molecule.id = molecular_property.molecule_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to update molecular properties in their projects."
        },
        "Delete molecular_properties for project members": {
            "operation": "DELETE",
            "using_clause": """
                EXISTS (
                  SELECT 1 FROM molecule
                  JOIN project ON project.id = molecule.project_id
                  JOIN user_profile ON user_profile.project_id = project.id
                  WHERE molecule.id = molecular_property.molecule_id
                    AND user_profile.user_id = auth.uid()
                )
            """,
            "with_check_clause": None,
            "comment": "Allows project members to delete molecular properties in their projects."
        },
        "Allow service role inserts on molecular_property": {
            "operation": "INSERT",
            "using_clause": None,
            "with_check_clause": "auth.role() = 'service_role'",
            "comment": "Allows service role to insert molecular properties regardless of project membership."
        }
    }
}

def log_rls_audit(action: str, table: str, details: Dict[str, Any]) -> None:
    """
    Log an RLS audit entry to the audit log.
    
    Args:
        action: The action being performed (verify, remediate)
        table: The table being affected
        details: Details of the action
    """
    timestamp = datetime.now().isoformat()
    
    # Create the audit log entry
    entry = {
        "timestamp": timestamp,
        "action": action,
        "table": table,
        "details": details,
        "user": os.environ.get("USER", "unknown"),
        "script": os.path.basename(__file__)
    }
    
    # Log to the audit log file
    audit_log_path = Path("logs/rls_audit.jsonl")
    with open(audit_log_path, "a") as f:
        f.write(json.dumps(entry) + "\n")
    
    # Also log to the regular logger
    logger.info(f"RLS {action} - {table}: {json.dumps(details)}")

def is_rls_enabled(table: str) -> bool:
    """
    Check if RLS is enabled for a table.
    
    Args:
        table: The table to check
        
    Returns:
        bool: True if RLS is enabled, False otherwise
    """
    try:
        client = get_supabase_client()
        
        query = f"""
        SELECT relrowsecurity
        FROM pg_class
        JOIN pg_namespace ON pg_namespace.oid = pg_class.relnamespace
        WHERE pg_namespace.nspname = 'public'
        AND pg_class.relname = '{table}';
        """
        
        result = client.rpc('exec_sql', {'query': query}).execute()
        
        if hasattr(result, 'data') and result.data and len(result.data) > 0:
            return result.data[0]['relrowsecurity']
        else:
            logger.error(f"Table {table} not found")
            return False
            
    except Exception as e:
        logger.error(f"Error checking if RLS is enabled for {table}: {str(e)}")
        return False

def get_table_policies(table: str) -> Dict[str, Dict[str, Any]]:
    """
    Get all policies for a table.
    
    Args:
        table: The table to get policies for
        
    Returns:
        Dict[str, Dict[str, Any]]: Dictionary of policy names to policy details
    """
    try:
        client = get_supabase_client()
        
        query = f"""
        SELECT
            p.polname AS policy_name,
            p.polcmd AS operation,
            pg_get_expr(p.polqual, p.polrelid) AS using_clause,
            pg_get_expr(p.polwithcheck, p.polrelid) AS with_check_clause,
            obj_description(p.oid, 'pg_policy') AS comment
        FROM pg_policy p
        JOIN pg_class c ON p.polrelid = c.oid
        JOIN pg_namespace n ON c.relnamespace = n.oid
        WHERE n.nspname = 'public'
        AND c.relname = '{table}';
        """
        
        result = client.rpc('exec_sql', {'query': query}).execute()
        
        policies = {}
        if hasattr(result, 'data') and result.data:
            for policy in result.data:
                policies[policy["policy_name"]] = {
                    "operation": policy["operation"],
                    "using_clause": policy["using_clause"],
                    "with_check_clause": policy["with_check_clause"],
                    "comment": policy["comment"]
                }
        
        return policies
            
    except Exception as e:
        logger.error(f"Error getting policies for {table}: {str(e)}")
        return {}

def normalize_sql(sql: str) -> str:
    """
    Normalize SQL for comparison by removing whitespace and comments.
    
    Args:
        sql: The SQL to normalize
        
    Returns:
        str: Normalized SQL
    """
    if not sql:
        return ""
    
    # Remove comments
    lines = sql.split("\n")
    lines = [line.split("--")[0] for line in lines]
    
    # Remove whitespace
    sql = " ".join(lines)
    sql = " ".join(sql.split())
    
    return sql.lower()

def compare_policies(actual_policies: Dict[str, Dict[str, Any]], 
                    canonical_policies: Dict[str, Dict[str, Any]]) -> Dict[str, List[str]]:
    """
    Compare actual policies with canonical policies.
    
    Args:
        actual_policies: The actual policies on the table
        canonical_policies: The canonical policies that should be present
        
    Returns:
        Dict[str, List[str]]: Dictionary of issues found
    """
    issues = {
        "missing_policies": [],
        "incorrect_policies": [],
        "unauthorized_policies": []
    }
    
    # Check for missing or incorrect policies
    for policy_name, canonical_policy in canonical_policies.items():
        if policy_name not in actual_policies:
            issues["missing_policies"].append(policy_name)
            continue
        
        actual_policy = actual_policies[policy_name]
        
        # Compare operation
        if actual_policy["operation"] != canonical_policy["operation"]:
            issues["incorrect_policies"].append(f"{policy_name} (operation mismatch)")
            continue
        
        # Compare using clause
        canonical_using = normalize_sql(canonical_policy["using_clause"])
        actual_using = normalize_sql(actual_policy["using_clause"])
        if canonical_using != actual_using:
            issues["incorrect_policies"].append(f"{policy_name} (using clause mismatch)")
            continue
        
        # Compare with check clause
        canonical_with_check = normalize_sql(canonical_policy["with_check_clause"])
        actual_with_check = normalize_sql(actual_policy["with_check_clause"])
        if canonical_with_check != actual_with_check:
            issues["incorrect_policies"].append(f"{policy_name} (with check clause mismatch)")
            continue
    
    # Check for unauthorized policies
    for policy_name in actual_policies:
        if policy_name not in canonical_policies:
            issues["unauthorized_policies"].append(policy_name)
    
    return issues

def verify_rls() -> Dict[str, Dict[str, Any]]:
    """
    Verify RLS settings for all tables.
    
    Returns:
        Dict[str, Dict[str, Any]]: Dictionary of issues found for each table
    """
    all_issues = {}
    
    for table in TABLES_TO_VERIFY:
        table_issues = {}
        
        # Check if RLS is enabled
        rls_enabled = is_rls_enabled(table)
        if not rls_enabled:
            table_issues["rls_disabled"] = True
        
        # Get actual policies
        actual_policies = get_table_policies(table)
        
        # Compare with canonical policies
        policy_issues = compare_policies(actual_policies, CANONICAL_POLICIES[table])
        
        # Add policy issues to table issues
        for issue_type, issues in policy_issues.items():
            if issues:
                table_issues[issue_type] = issues
        
        # Log audit entry
        if table_issues:
            log_rls_audit("verify", table, table_issues)
            all_issues[table] = table_issues
    
    return all_issues

def enable_rls(table: str) -> bool:
    """
    Enable RLS for a table.
    
    Args:
        table: The table to enable RLS for
        
    Returns:
        bool: True if successful, False otherwise
    """
    query = f"ALTER TABLE {table} ENABLE ROW LEVEL SECURITY;"
    
    success = execute_ddl(query)
    
    if success:
        log_rls_audit("remediate", table, {"action": "enable_rls", "success": True})
    else:
        logger.error(f"Error enabling RLS for {table}")
    
    return success

def create_policy(table: str, policy_name: str, policy_details: Dict[str, Any]) -> bool:
    """
    Create a policy on a table.
    
    Args:
        table: The table to create the policy on
        policy_name: The name of the policy
        policy_details: The details of the policy
        
    Returns:
        bool: True if successful, False otherwise
    """
    operation = policy_details["operation"]
    using_clause = policy_details["using_clause"]
    with_check_clause = policy_details["with_check_clause"]
    comment = policy_details["comment"]
    
    # Build the query
    query = f"CREATE POLICY \"{policy_name}\" ON {table} FOR {operation} "
    
    if using_clause:
        query += f"USING ({using_clause.strip()})"
    
    if with_check_clause:
        query += f"WITH CHECK ({with_check_clause.strip()})"
    
    query += ";"
    
    # Add comment if available
    if comment:
        query += f"\nCOMMENT ON POLICY \"{policy_name}\" ON {table} IS '{comment}';"
    
    success = execute_ddl(query)
    
    if success:
        log_rls_audit("remediate", table, {
            "action": "create_policy",
            "policy_name": policy_name,
            "success": True
        })
    else:
        logger.error(f"Error creating policy {policy_name} on {table}")
    
    return success

def drop_policy(table: str, policy_name: str) -> bool:
    """
    Drop a policy from a table.
    
    Args:
        table: The table to drop the policy from
        policy_name: The name of the policy
        
    Returns:
        bool: True if successful, False otherwise
    """
    query = f"DROP POLICY \"{policy_name}\" ON {table};"
    
    success = execute_ddl(query)
    
    if success:
        log_rls_audit("remediate", table, {
            "action": "drop_policy",
            "policy_name": policy_name,
            "success": True
        })
    else:
        logger.error(f"Error dropping policy {policy_name} from {table}")
    
    return success

def remediate_rls(issues: Dict[str, Dict[str, Any]]) -> bool:
    """
    Remediate RLS issues.
    
    Args:
        issues: Dictionary of issues found for each table
        
    Returns:
        bool: True if all issues were remediated, False otherwise
    """
    all_success = True
    
    for table, table_issues in issues.items():
        # Enable RLS if disabled
        if table_issues.get("rls_disabled", False):
            success = enable_rls(table)
            all_success = all_success and success
        
        # Create missing policies
        for policy_name in table_issues.get("missing_policies", []):
            policy_details = CANONICAL_POLICIES[table][policy_name]
            success = create_policy(table, policy_name, policy_details)
            all_success = all_success and success
        
        # Fix incorrect policies by dropping and recreating them
        for policy_name_with_issue in table_issues.get("incorrect_policies", []):
            policy_name = policy_name_with_issue.split(" (")[0]
            
            # Drop the policy
            success = drop_policy(table, policy_name)
            all_success = all_success and success
            
            # Recreate the policy
            policy_details = CANONICAL_POLICIES[table][policy_name]
            success = create_policy(table, policy_name, policy_details)
            all_success = all_success and success
        
        # Drop unauthorized policies
        for policy_name in table_issues.get("unauthorized_policies", []):
            success = drop_policy(table, policy_name)
            all_success = all_success and success
    
    return all_success

def ensure_rls_restored(func: Callable) -> Callable:
    """
    Decorator to ensure RLS is restored before and after a function call.
    
    Args:
        func: The function to decorate
        
    Returns:
        Callable: The decorated function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Verify and remediate RLS before function call
        logger.info(f"Verifying RLS before executing {func.__name__}")
        issues_before = verify_rls()
        if issues_before:
            logger.warning(f"RLS issues found before executing {func.__name__}: {issues_before}")
            remediate_rls(issues_before)
        
        try:
            # Call the function
            result = func(*args, **kwargs)
            
            # Verify and remediate RLS after function call
            logger.info(f"Verifying RLS after executing {func.__name__}")
            issues_after = verify_rls()
            if issues_after:
                logger.warning(f"RLS issues found after executing {func.__name__}: {issues_after}")
                remediate_rls(issues_after)
            
            return result
            
        except Exception as e:
            # Verify and remediate RLS even if an exception occurs
            logger.error(f"Exception in {func.__name__}: {str(e)}")
            logger.info(f"Verifying RLS after exception in {func.__name__}")
            issues_after_exception = verify_rls()
            if issues_after_exception:
                logger.warning(f"RLS issues found after exception in {func.__name__}: {issues_after_exception}")
                remediate_rls(issues_after_exception)
            
            # Re-raise the exception
            raise
    
    return wrapper

if __name__ == "__main__":
    # Test the RLS verification and remediation
    print("Verifying RLS settings...")
    issues = verify_rls()
    
    if not issues:
        print("No RLS issues found.")
    else:
        print(f"RLS issues found: {json.dumps(issues, indent=2)}")
        
        print("\nRemediating RLS issues...")
        success = remediate_rls(issues)
        
        if success:
            print("All RLS issues remediated successfully.")
        else:
            print("Some RLS issues could not be remediated.")
        
        print("\nVerifying RLS settings again...")
        issues_after = verify_rls()
        
        if not issues_after:
            print("No RLS issues found after remediation.")
        else:
            print(f"RLS issues still present after remediation: {json.dumps(issues_after, indent=2)}")
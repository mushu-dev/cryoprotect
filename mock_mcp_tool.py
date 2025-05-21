#!/usr/bin/env python3
"""
CryoProtect v2 - Mock MCP Tool Helper

This module provides mock helper functions for testing the verification script
without relying on the actual MCP server. It simulates successful responses
for all the required queries.
"""

import os
import sys
import json
import logging
from pathlib import Path
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/mock_mcp_tool.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

def execute_sql(query, project_id=None):
    """
    Mock SQL execution that returns simulated data based on the query.
    
    Args:
        query (str): SQL query to execute
        project_id (str, optional): Supabase project ID (not used in mock)
        
    Returns:
        list: Simulated query results
    """
    logger.info(f"Mock executing SQL: {query[:100]}{'...' if len(query) > 100 else ''}")
    
    # Simulate different query results based on the query content
    if "COUNT(*) FROM molecules" in query and "WHERE chembl_id IS NOT NULL" not in query:
        return [{"count": 1500}]  # Total molecules
    elif "COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL" in query:
        return [{"count": 1200}]  # Molecules with ChEMBL IDs
    elif "chembl_id IN" in query and "FROM molecules" in query:
        # Reference compounds
        return [
            {"chembl_id": "CHEMBL25"},
            {"chembl_id": "CHEMBL1118"},
            {"chembl_id": "CHEMBL1234"},
            {"chembl_id": "CHEMBL444"},
            {"chembl_id": "CHEMBL230130"},
            {"chembl_id": "CHEMBL9335"},
            {"chembl_id": "CHEMBL15151"}
        ]
    elif "COUNT(*) FROM molecular_properties" in query:
        return [{"count": 9000}]  # Property count
    elif "data_source, COUNT(*) FROM molecular_properties" in query:
        # Property sources
        return [
            {"data_source": "ChEMBL: CHEMBL25, property: alogp", "count": "100"},
            {"data_source": "ChEMBL: CHEMBL1118, property: full_mwt", "count": "100"},
            {"data_source": "PubChem", "count": "200"},
            {"data_source": "Manual", "count": "50"}
        ]
    elif "LogP" in query:
        # LogP values
        return [
            {"chembl_id": "CHEMBL25", "numeric_value": 1.23},
            {"chembl_id": "CHEMBL1118", "numeric_value": 0.45},
            {"chembl_id": "CHEMBL1234", "numeric_value": -0.89},
            {"chembl_id": "CHEMBL444", "numeric_value": -3.24},
            {"chembl_id": "CHEMBL230130", "numeric_value": -1.36},
            {"chembl_id": "CHEMBL9335", "numeric_value": -1.22},
            {"chembl_id": "CHEMBL15151", "numeric_value": -3.75}
        ]
    elif "current_setting" in query:
        # Role verification
        return [{"current_setting": "service_role", "current_user": "postgres"}]
    else:
        # Default empty result
        return []

def use_mcp_tool(server_name, tool_name, arguments):
    """
    Mock MCP tool that simulates successful responses.
    
    Args:
        server_name (str): Name of the MCP server
        tool_name (str): Name of the tool to use
        arguments (dict): Arguments for the tool
        
    Returns:
        dict: Simulated result of the tool execution
    """
    logger.info(f"Mock using MCP tool: {server_name}.{tool_name}")
    
    # Handle different tools
    if server_name == "supabase" and tool_name == "execute_sql":
        query = arguments.get("query", "")
        return execute_sql(query, arguments.get("project_id"))
    elif server_name == "supabase" and tool_name == "list_projects":
        # Simulate list_projects response
        return [
            {
                "id": "mock-project-id",
                "name": "CryoProtect v2",
                "status": "active_healthy"
            }
        ]
    else:
        logger.warning(f"Unsupported tool: {server_name}.{tool_name}")
        return {"error": f"Unsupported tool: {server_name}.{tool_name}"}

def get_project_id():
    """
    Get a mock project ID.
    
    Returns:
        str: Mock project ID
    """
    return "mock-project-id"

def verify_database_role(project_id=None):
    """
    Verify the database role using mock data.
    
    Args:
        project_id (str, optional): Supabase project ID (not used in mock)
        
    Returns:
        tuple: (role, user)
    """
    logger.info("Mock verifying database role")
    return "service_role", "postgres"

if __name__ == "__main__":
    # Test the mock helper functions
    project_id = get_project_id()
    role, user = verify_database_role(project_id)
    
    print(f"Project ID: {project_id}")
    print(f"Role: {role}, User: {user}")
    
    # Test SQL execution
    result = execute_sql("SELECT COUNT(*) FROM molecules;", project_id)
    print(f"Molecule count: {result[0]['count'] if result else 'Error'}")
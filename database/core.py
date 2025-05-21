#!/usr/bin/env python3
"""
Database Core Module for CryoProtect on Fedora

This module provides a consolidated approach to database connectivity,
combining direct REST API access with MCP SQL execution capabilities.
It serves as the primary entry point for database operations in the
CryoProtect application running on Fedora Linux.
"""

import os
import sys
import logging
import requests
import json
import time
from typing import Any, Dict, List, Optional, Union, Tuple
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def load_config() -> Dict[str, Any]:
    """
    Load database configuration from environment variables.
    
    Returns:
        Dict with configuration values
    """
    return {
        'supabase_url': os.environ.get('SUPABASE_URL'),
        'supabase_key': os.environ.get('SUPABASE_KEY'),
        'supabase_service_key': os.environ.get('SUPABASE_SERVICE_KEY', os.environ.get('SUPABASE_KEY')),
        'project_id': os.environ.get('MCP_PROJECT_ID'),
        'debug': os.environ.get('DB_DEBUG', 'false').lower() == 'true'
    }

def validate_config(config: Dict[str, Any] = None) -> Tuple[bool, str]:
    """
    Validate database configuration.
    
    Args:
        config: Optional config dict (loaded from environment if not provided)
        
    Returns:
        Tuple of (is_valid: bool, message: str)
    """
    if config is None:
        config = load_config()
    
    missing = []
    
    if not config.get('supabase_url'):
        missing.append('SUPABASE_URL')
    
    if not config.get('supabase_key'):
        missing.append('SUPABASE_KEY')
    
    # Project ID is only required for MCP operations
    if not config.get('project_id'):
        logger.warning("MCP_PROJECT_ID is not set. MCP operations will not be available.")
    
    if missing:
        return False, f"Missing required configuration values: {', '.join(missing)}"
    
    return True, "Configuration is valid"

# REST API functions
def rest_request(
    method: str, 
    path: str, 
    data: Any = None, 
    params: Dict[str, Any] = None, 
    headers: Dict[str, str] = None,
    use_service_role: bool = False
) -> requests.Response:
    """
    Make a request to Supabase REST API.
    
    Args:
        method: HTTP method (GET, POST, PUT, PATCH, DELETE)
        path: API path (without /rest/v1/ prefix)
        data: Request data (for POST, PUT, PATCH)
        params: Query parameters
        headers: Additional headers
        use_service_role: Whether to use service role key instead of anon key
        
    Returns:
        requests.Response object
    """
    config = load_config()
    
    # Validate config
    is_valid, message = validate_config(config)
    if not is_valid:
        raise ValueError(f"Invalid configuration: {message}")
    
    # Build URL (remove leading slash if present)
    if path.startswith('/'):
        path = path[1:]
    
    url = f"{config['supabase_url']}/rest/v1/{path}"
    
    # Prepare headers
    if headers is None:
        headers = {}
    
    # Use service role key if requested
    key = config['supabase_service_key'] if use_service_role else config['supabase_key']
    
    headers.update({
        'apikey': key,
        'Authorization': f"Bearer {key}"
    })
    
    # Log request details if debug is enabled
    if config.get('debug'):
        safe_headers = headers.copy()
        if 'Authorization' in safe_headers:
            safe_headers['Authorization'] = 'Bearer ********'
        if 'apikey' in safe_headers:
            safe_headers['apikey'] = '********'
            
        logger.debug(f"Making {method} request to {url}")
        logger.debug(f"Headers: {json.dumps(safe_headers)}")
        logger.debug(f"Params: {json.dumps(params) if params else None}")
        if data and method in ['POST', 'PUT', 'PATCH']:
            logger.debug(f"Data: {json.dumps(data)}")
    
    start_time = time.time()
    try:
        response = requests.request(
            method=method,
            url=url,
            json=data,
            params=params,
            headers=headers,
            timeout=10  # 10 second timeout
        )
        
        elapsed = time.time() - start_time
        
        # Log response details if debug is enabled
        if config.get('debug'):
            logger.debug(f"Response status: {response.status_code} ({elapsed:.3f}s)")
            logger.debug(f"Response headers: {json.dumps(dict(response.headers))}")
            
            # Try to log response as JSON if possible
            try:
                resp_data = response.json()
                logger.debug(f"Response data: {json.dumps(resp_data)[:500]}...")
            except:
                # If not JSON, log first part of text
                logger.debug(f"Response text: {response.text[:500]}...")
        
        return response
    except requests.RequestException as e:
        elapsed = time.time() - start_time
        logger.error(f"REST request failed after {elapsed:.3f}s: {str(e)}")
        raise

# MCP SQL execution
def execute_sql(
    query: str, 
    params: Optional[Union[Tuple, Dict]] = None, 
    project_id: str = None
) -> List[Dict[str, Any]]:
    """
    Execute SQL query through MCP.
    
    Args:
        query: SQL query to execute
        params: Query parameters
        project_id: Supabase project ID (loaded from environment if not provided)
        
    Returns:
        List of result rows as dictionaries
    """
    config = load_config()
    proj_id = project_id or config.get('project_id')
    
    if not proj_id:
        raise ValueError("Project ID not provided and not found in environment")
    
    # Log query if debug is enabled
    if config.get('debug'):
        logger.debug(f"Executing SQL query through MCP: {query}")
        if params:
            logger.debug(f"Query params: {params}")
    
    # Handle query parameters by formatting the query if provided
    if params:
        if isinstance(params, dict):
            # Convert dict params to SQL
            for key, value in params.items():
                placeholder = f"%({key})s"
                if isinstance(value, str):
                    # Escape single quotes for SQL
                    escaped_value = value.replace("'", "''")
                    formatted_value = f"'{escaped_value}'"
                elif value is None:
                    formatted_value = "NULL"
                elif isinstance(value, bool):
                    formatted_value = "TRUE" if value else "FALSE"
                else:
                    formatted_value = str(value)
                query = query.replace(placeholder, formatted_value)
        else:
            # Convert tuple params to SQL
            for i, value in enumerate(params):
                placeholder = "%s"
                if isinstance(value, str):
                    # Escape single quotes for SQL
                    escaped_value = value.replace("'", "''")
                    formatted_value = f"'{escaped_value}'"
                elif value is None:
                    formatted_value = "NULL"
                elif isinstance(value, bool):
                    formatted_value = "TRUE" if value else "FALSE"
                else:
                    formatted_value = str(value)
                query = query.replace(placeholder, formatted_value, 1)
    
    start_time = time.time()
    try:
        # Import the appropriate MCP function
        try:
            # First try to import from the supabase_mcp_tools module
            from supabase_mcp_tools import execute_sql_through_mcp
            result = execute_sql_through_mcp(query, proj_id)
        except ImportError:
            # Fallback to use_mcp_tool if available
            try:
                from use_mcp_tool import execute_sql
                result = execute_sql(query, proj_id)
            except ImportError:
                logger.error("Failed to import MCP tools. MCP SQL execution will not work.")
                raise ValueError("MCP tools not available. Install supabase_mcp_tools or use_mcp_tool.")
        
        elapsed = time.time() - start_time
        
        # Log result if debug is enabled
        if config.get('debug'):
            logger.debug(f"SQL query completed in {elapsed:.3f}s")
            if result:
                logger.debug(f"Returned {len(result)} rows")
                if len(result) > 0:
                    logger.debug(f"First row: {json.dumps(result[0])}")
        
        return result
    except Exception as e:
        elapsed = time.time() - start_time
        logger.error(f"SQL execution failed after {elapsed:.3f}s: {str(e)}")
        raise

# Connection testing
def test_connection() -> Dict[str, Any]:
    """
    Test database connection using both REST API and MCP (if available).
    
    Returns:
        Dict with test results
    """
    results = {
        'rest_api': {
            'status': 'untested',
            'message': '',
            'elapsed': 0
        },
        'mcp': {
            'status': 'untested',
            'message': '',
            'elapsed': 0
        },
        'overall': {
            'status': 'untested',
            'message': ''
        }
    }
    
    # Test REST API
    start_time = time.time()
    try:
        response = rest_request('GET', '')
        elapsed = time.time() - start_time
        results['rest_api']['elapsed'] = elapsed
        
        if response.status_code < 300:
            results['rest_api']['status'] = 'success'
            results['rest_api']['message'] = f"Connected in {elapsed:.3f}s"
        else:
            results['rest_api']['status'] = 'error'
            results['rest_api']['message'] = f"Failed with status {response.status_code}"
    except Exception as e:
        elapsed = time.time() - start_time
        results['rest_api']['elapsed'] = elapsed
        results['rest_api']['status'] = 'error'
        results['rest_api']['message'] = str(e)
    
    # Test MCP SQL
    config = load_config()
    if config.get('project_id'):
        start_time = time.time()
        try:
            result = execute_sql("SELECT 1 as test")
            elapsed = time.time() - start_time
            results['mcp']['elapsed'] = elapsed
            
            if result and result[0].get('test') == 1:
                results['mcp']['status'] = 'success'
                results['mcp']['message'] = f"Connected in {elapsed:.3f}s"
            else:
                results['mcp']['status'] = 'error'
                results['mcp']['message'] = "Unexpected result"
        except Exception as e:
            elapsed = time.time() - start_time
            results['mcp']['elapsed'] = elapsed
            results['mcp']['status'] = 'error'
            results['mcp']['message'] = str(e)
    else:
        results['mcp']['status'] = 'skipped'
        results['mcp']['message'] = "No project ID configured"
    
    # Determine overall status
    if results['rest_api']['status'] == 'success' or results['mcp']['status'] == 'success':
        results['overall']['status'] = 'success'
        results['overall']['message'] = "At least one connection method succeeded"
    else:
        results['overall']['status'] = 'error'
        results['overall']['message'] = "All connection methods failed"
    
    return results

# Convenience functions for common operations
def get_table_data(
    table_name: str, 
    filters: Dict[str, Any] = None, 
    limit: int = 100, 
    offset: int = 0,
    order_by: str = None,
    select: str = "*"
) -> List[Dict[str, Any]]:
    """
    Get data from a table with optional filters.
    
    Args:
        table_name: Name of the table
        filters: Postgrest filter conditions
        limit: Maximum number of records to return
        offset: Number of records to skip
        order_by: Order by clause
        select: Columns to select
        
    Returns:
        List of records as dictionaries
    """
    params = {
        'limit': limit,
        'offset': offset
    }
    
    # Add select parameter if not default
    if select != "*":
        params['select'] = select
    
    # Add order parameter if provided
    if order_by:
        params['order'] = order_by
    
    # Add filter parameters
    if filters:
        params.update(filters)
    
    response = rest_request('GET', table_name, params=params)
    response.raise_for_status()
    return response.json()

def insert_data(
    table_name: str, 
    data: Union[Dict[str, Any], List[Dict[str, Any]]],
    upsert: bool = False
) -> List[Dict[str, Any]]:
    """
    Insert data into a table.
    
    Args:
        table_name: Name of the table
        data: Data to insert (single record or list of records)
        upsert: Whether to perform upsert operation
        
    Returns:
        List of inserted records
    """
    headers = {'Prefer': 'return=representation'}
    
    if upsert:
        headers['Prefer'] += ', resolution=merge-duplicates'
    
    response = rest_request(
        'POST', 
        table_name, 
        data=data, 
        headers=headers
    )
    response.raise_for_status()
    return response.json()

def update_data(
    table_name: str, 
    filters: Dict[str, Any],
    data: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Update data in a table.
    
    Args:
        table_name: Name of the table
        filters: Filters to identify records to update
        data: Data to update
        
    Returns:
        List of updated records
    """
    path = table_name
    headers = {'Prefer': 'return=representation'}
    
    response = rest_request(
        'PATCH', 
        path, 
        data=data, 
        params=filters, 
        headers=headers
    )
    response.raise_for_status()
    return response.json()

def delete_data(
    table_name: str, 
    filters: Dict[str, Any],
    return_records: bool = False
) -> Union[bool, List[Dict[str, Any]]]:
    """
    Delete data from a table.
    
    Args:
        table_name: Name of the table
        filters: Filters to identify records to delete
        return_records: Whether to return deleted records
        
    Returns:
        True if successful or list of deleted records if return_records=True
    """
    path = table_name
    headers = {}
    
    if return_records:
        headers['Prefer'] = 'return=representation'
    
    response = rest_request(
        'DELETE', 
        path, 
        params=filters, 
        headers=headers
    )
    
    response.raise_for_status()
    
    if return_records and response.text:
        return response.json()
    
    return True

def count_records(
    table_name: str, 
    filters: Dict[str, Any] = None
) -> int:
    """
    Count records in a table.
    
    Args:
        table_name: Name of the table
        filters: Filters to apply
        
    Returns:
        Number of records
    """
    headers = {'Prefer': 'count=exact'}
    params = {}
    
    if filters:
        params.update(filters)
    
    response = rest_request(
        'GET', 
        table_name, 
        params=params, 
        headers=headers
    )
    
    response.raise_for_status()
    
    # Extract count from headers
    count = response.headers.get('content-range', '').split('/')[-1]
    
    try:
        return int(count)
    except ValueError:
        # If header parsing fails, fall back to counting results
        return len(response.json())

def run_transaction(queries: List[str]) -> List[Dict[str, Any]]:
    """
    Run a series of queries as a transaction.
    
    Args:
        queries: List of SQL queries to execute in transaction
        
    Returns:
        Results of the transaction
    """
    # Convert queries to a single transaction
    transaction = [
        "BEGIN;",
        *queries,
        "COMMIT;"
    ]
    
    combined_query = "\n".join(transaction)
    
    try:
        return execute_sql(combined_query)
    except Exception as e:
        logger.error(f"Transaction failed: {str(e)}")
        # Try to execute rollback
        try:
            execute_sql("ROLLBACK;")
        except:
            pass
        raise

def get_schema_info(table_name: str = None) -> Dict[str, Any]:
    """
    Get schema information for a table or all tables.

    Args:
        table_name: Optional table name to get schema for

    Returns:
        Schema information
    """
    config = load_config()

    # If MCP is available, use SQL to get schema info
    if config.get('project_id'):
        try:
            query = """
            SELECT
                table_schema,
                table_name,
                column_name,
                data_type,
                is_nullable,
                column_default
            FROM
                information_schema.columns
            WHERE
                table_schema = 'public'
            """

            if table_name:
                query += f" AND table_name = '{table_name}'"

            query += " ORDER BY table_name, ordinal_position"

            results = execute_sql(query)

            # Organize results by table
            schema = {}

            for row in results:
                table = row['table_name']

                if table not in schema:
                    schema[table] = {
                        'schema': row['table_schema'],
                        'columns': []
                    }

                schema[table]['columns'].append({
                    'name': row['column_name'],
                    'type': row['data_type'],
                    'nullable': row['is_nullable'] == 'YES',
                    'default': row['column_default']
                })

            if table_name and table_name in schema:
                return schema[table_name]

            return schema

        except Exception as e:
            logger.warning(f"Failed to get schema using SQL, falling back to REST API: {str(e)}")

    # Fallback to REST API
    try:
        response = rest_request('GET', '')
        response.raise_for_status()

        # Parse OpenAPI schema to get tables
        api_schema = response.json()
        paths = api_schema.get('paths', {})

        # Extract tables from paths
        tables = []
        for path in paths:
            parts = path.strip('/').split('/')
            if len(parts) >= 3 and parts[0] == 'rest' and parts[1] == 'v1':
                table_name_from_path = parts[2]
                if table_name_from_path and '.' not in table_name_from_path and '{' not in table_name_from_path:
                    if table_name and table_name_from_path != table_name:
                        continue
                    tables.append(table_name_from_path)

        # Get schema for each table
        schema = {}
        for table in tables:
            try:
                # Get a sample record to infer schema
                sample = get_table_data(table, limit=1)

                schema[table] = {
                    'schema': 'public',
                    'columns': []
                }

                if sample:
                    # Extract columns from first record
                    first_record = sample[0]
                    for col_name, value in first_record.items():
                        type_str = 'unknown'
                        if value is None:
                            type_str = 'unknown'
                        elif isinstance(value, int):
                            type_str = 'integer'
                        elif isinstance(value, float):
                            type_str = 'number'
                        elif isinstance(value, bool):
                            type_str = 'boolean'
                        elif isinstance(value, str):
                            type_str = 'string'
                        elif isinstance(value, dict):
                            type_str = 'object'
                        elif isinstance(value, list):
                            type_str = 'array'

                        schema[table]['columns'].append({
                            'name': col_name,
                            'type': type_str,
                            'nullable': True,  # Assume nullable without schema info
                            'default': None
                        })
                else:
                    # No sample data available
                    schema[table]['columns'] = []
            except Exception as e:
                logger.warning(f"Error getting schema for table {table}: {str(e)}")
                # Still include the table but without columns
                schema[table] = {
                    'schema': 'public',
                    'columns': []
                }

        if table_name and table_name in schema:
            return schema[table_name]

        return schema

    except Exception as e:
        logger.error(f"Failed to get schema using REST API: {str(e)}")
        return {}

# Main functions for testing
if __name__ == "__main__":
    print("CryoProtect Database Core Module")
    print("--------------------------------")

    # Load and validate configuration
    config = load_config()
    valid, message = validate_config(config)

    print(f"Configuration: {'Valid' if valid else 'Invalid'}")
    print(f"Message: {message}")

    if valid:
        print("\nTesting connection...")
        results = test_connection()

        print(f"\nREST API: {results['rest_api']['status']}")
        print(f"Message: {results['rest_api']['message']}")

        print(f"\nMCP SQL: {results['mcp']['status']}")
        print(f"Message: {results['mcp']['message']}")

        print(f"\nOverall: {results['overall']['status']}")
        print(f"Message: {results['overall']['message']}")

        if results['overall']['status'] == 'success':
            try:
                print("\nGetting schema information...")
                schema = get_schema_info()

                if schema:
                    print(f"\nFound {len(schema)} tables:")
                    for table, info in schema.items():
                        print(f"- {table} ({len(info['columns'])} columns)")
                else:
                    print("\nNo schema information available")
            except Exception as e:
                print(f"\nError getting schema information: {str(e)}")

            # Test getting table data if REST API is working
            if results['rest_api']['status'] == 'success':
                try:
                    print("\nGetting table list via REST API...")
                    response = rest_request('GET', '')
                    openapi = response.json()

                    paths = openapi.get('paths', {})
                    tables = []

                    for path in paths:
                        parts = path.strip('/').split('/')
                        if len(parts) >= 3 and parts[0] == 'rest' and parts[1] == 'v1':
                            table_name = parts[2]
                            if table_name and '.' not in table_name and '{' not in table_name:
                                tables.append(table_name)

                    if tables:
                        print(f"\nFound {len(tables)} tables via REST API:")
                        for table in sorted(tables):
                            print(f"- {table}")

                            # Try to get record count for this table
                            try:
                                count = count_records(table)
                                print(f"  - Records: {count}")
                            except:
                                print(f"  - Records: Unknown")
                except Exception as e:
                    print(f"\nError getting table list: {str(e)}")

    print("\nDone.")
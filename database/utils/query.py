"""
Database query utilities.

This module provides functions for executing database queries and
processing query results.
"""

import logging
from typing import Dict, List, Any, Optional, Union, Tuple

from .connection import supabase_connection, get_supabase_service_client

# Configure logging
logger = logging.getLogger(__name__)

def execute_query(table_name: str, query_type: str, 
                  data: Optional[Dict[str, Any]] = None,
                  columns: Optional[List[str]] = None,
                  filters: Optional[Dict[str, Any]] = None,
                  order_by: Optional[str] = None,
                  limit: Optional[int] = None,
                  use_service_role: bool = False) -> Dict[str, Any]:
    """
    Execute a database query on the specified table.
    
    Args:
        table_name: Name of the table to query
        query_type: Type of query (select, insert, update, delete)
        data: Data to insert or update (for insert/update queries)
        columns: Columns to select (for select queries)
        filters: Filter conditions (for select/update/delete queries)
        order_by: Column to order by (for select queries)
        limit: Maximum number of results to return (for select queries)
        use_service_role: Whether to use service role permissions
        
    Returns:
        Dictionary with query results
    """
    try:
        if use_service_role:
            supabase = get_supabase_service_client()
        else:
            with supabase_connection() as supabase:
                return _execute_query_internal(
                    supabase, table_name, query_type, data, 
                    columns, filters, order_by, limit
                )
        
        return _execute_query_internal(
            supabase, table_name, query_type, data, 
            columns, filters, order_by, limit
        )
    except Exception as e:
        logger.error(f"Query execution error ({query_type} on {table_name}): {str(e)}")
        return {"error": str(e), "data": None}

def _execute_query_internal(supabase, table_name: str, query_type: str, 
                           data: Optional[Dict[str, Any]],
                           columns: Optional[List[str]],
                           filters: Optional[Dict[str, Any]],
                           order_by: Optional[str],
                           limit: Optional[int]) -> Dict[str, Any]:
    """
    Internal function to execute a database query.
    
    Args:
        supabase: Supabase client
        table_name: Name of the table to query
        query_type: Type of query (select, insert, update, delete)
        data: Data to insert or update
        columns: Columns to select
        filters: Filter conditions
        order_by: Column to order by
        limit: Maximum number of results to return
        
    Returns:
        Dictionary with query results
    """
    query = supabase.table(table_name)
    
    if query_type == 'select':
        # Select specific columns or all columns
        if columns:
            select_str = ','.join(columns)
            query = query.select(select_str)
        else:
            query = query.select('*')
        
        # Apply filters
        if filters:
            for key, value in filters.items():
                if isinstance(value, dict) and 'operator' in value:
                    query = query.filter(key, value['operator'], value['value'])
                else:
                    query = query.eq(key, value)
        
        # Apply order by
        if order_by:
            query = query.order(order_by)
        
        # Apply limit
        if limit:
            query = query.limit(limit)
        
        # Execute query
        result = query.execute()
        return {"data": result.data, "error": None}
    
    elif query_type == 'insert':
        if not data:
            return {"error": "No data provided for insert", "data": None}
        
        # Execute query
        result = query.insert(data).execute()
        return {"data": result.data, "error": None}
    
    elif query_type == 'update':
        if not data:
            return {"error": "No data provided for update", "data": None}
        
        # Apply filters
        if filters:
            for key, value in filters.items():
                if isinstance(value, dict) and 'operator' in value:
                    query = query.filter(key, value['operator'], value['value'])
                else:
                    query = query.eq(key, value)
        
        # Execute query
        result = query.update(data).execute()
        return {"data": result.data, "error": None}
    
    elif query_type == 'delete':
        # Apply filters
        if filters:
            for key, value in filters.items():
                if isinstance(value, dict) and 'operator' in value:
                    query = query.filter(key, value['operator'], value['value'])
                else:
                    query = query.eq(key, value)
        
        # Execute query
        result = query.delete().execute()
        return {"data": result.data, "error": None}
    
    else:
        return {"error": f"Invalid query type: {query_type}", "data": None}

def execute_rpc(function_name: str, params: Optional[Dict[str, Any]] = None,
               use_service_role: bool = False) -> Dict[str, Any]:
    """
    Execute a stored procedure or function using RPC.
    
    Args:
        function_name: Name of the function to execute
        params: Parameters to pass to the function
        use_service_role: Whether to use service role permissions
        
    Returns:
        Dictionary with RPC results
    """
    try:
        if use_service_role:
            supabase = get_supabase_service_client()
        else:
            with supabase_connection() as supabase:
                return _execute_rpc_internal(supabase, function_name, params)
        
        return _execute_rpc_internal(supabase, function_name, params)
    except Exception as e:
        logger.error(f"RPC execution error ({function_name}): {str(e)}")
        return {"error": str(e), "data": None}

def _execute_rpc_internal(supabase, function_name: str, 
                          params: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Internal function to execute a stored procedure or function using RPC.
    
    Args:
        supabase: Supabase client
        function_name: Name of the function to execute
        params: Parameters to pass to the function
        
    Returns:
        Dictionary with RPC results
    """
    if params:
        result = supabase.rpc(function_name, params).execute()
    else:
        result = supabase.rpc(function_name).execute()
    
    return {"data": result.data, "error": None}

def execute_raw_sql(sql: str, params: Optional[Dict[str, Any]] = None,
                   use_service_role: bool = True) -> Dict[str, Any]:
    """
    Execute a raw SQL query using the exec_sql function.
    
    Note: This requires the exec_sql function to be created in the database.
    
    Args:
        sql: SQL query to execute
        params: Parameters to pass to the query
        use_service_role: Whether to use service role permissions
        
    Returns:
        Dictionary with query results
    """
    rpc_params = {"sql_query": sql}
    if params:
        rpc_params["params"] = params
    
    return execute_rpc("exec_sql", rpc_params, use_service_role=use_service_role)

def batch_insert(table_name: str, records: List[Dict[str, Any]],
                batch_size: int = 100, use_service_role: bool = False) -> Dict[str, Any]:
    """
    Insert records in batches to avoid hitting size limits.
    
    Args:
        table_name: Name of the table to insert into
        records: List of records to insert
        batch_size: Number of records per batch
        use_service_role: Whether to use service role permissions
        
    Returns:
        Dictionary with insert results
    """
    if not records:
        return {"data": [], "error": None}
    
    all_inserted = []
    error = None
    
    for i in range(0, len(records), batch_size):
        batch = records[i:i+batch_size]
        result = execute_query(
            table_name=table_name,
            query_type='insert',
            data=batch,
            use_service_role=use_service_role
        )
        
        if result.get('error'):
            error = result['error']
            logger.error(f"Batch insert error (batch {i//batch_size}): {error}")
            break
        
        if result.get('data'):
            all_inserted.extend(result['data'])
    
    return {"data": all_inserted, "error": error}

def get_table_column_names(table_name: str) -> List[str]:
    """
    Get the column names for a specific table.
    
    Args:
        table_name: Name of the table
        
    Returns:
        List of column names
    """
    sql = f"""
    SELECT column_name 
    FROM information_schema.columns 
    WHERE table_schema = 'public' 
    AND table_name = '{table_name}'
    ORDER BY ordinal_position
    """
    
    result = execute_raw_sql(sql)
    
    if result.get('error') or not result.get('data'):
        logger.error(f"Error getting column names for table {table_name}")
        return []
    
    return [row['column_name'] for row in result['data']]

def table_exists(table_name: str) -> bool:
    """
    Check if a table exists in the database.
    
    Args:
        table_name: Name of the table to check
        
    Returns:
        True if the table exists, False otherwise
    """
    sql = f"""
    SELECT EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' 
        AND table_name = '{table_name}'
    )
    """
    
    result = execute_raw_sql(sql)
    
    if result.get('error') or not result.get('data'):
        logger.error(f"Error checking if table {table_name} exists")
        return False
    
    return result['data'][0]['exists']
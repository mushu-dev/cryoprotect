"""
CryoProtect v2 - API Compatibility Layer

This module provides helper functions and decorators to ensure backward compatibility
during the migration from singular to plural table names.

Key features:
1. Table name mapping between singular and plural forms
2. Logging of deprecated endpoint usage
3. Deprecation notice injection into API responses
4. Error handling for both table formats

Author: CryoProtect Team
Date: 2025-04-18
"""

import functools
import logging
import time
import warnings
from datetime import datetime
from flask import request, g, current_app, jsonify
import json
from typing import Dict, Any, Callable, List, Union, Optional

# Configure logging
logger = logging.getLogger("api_compatibility")
if not logger.handlers:
    handler = logging.FileHandler(f"deprecated_api_usage_{datetime.now().strftime('%Y%m%d')}.log")
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# Table name mappings
TABLE_MAPPINGS = {
    "molecule": "molecules",
    "mixture": "mixtures",
    "experiment": "experiments",
    "prediction": "predictions",
    "project": "projects"
}

# Reverse mappings for convenience
REVERSE_MAPPINGS = {v: k for k, v in TABLE_MAPPINGS.items()}

# Deprecation notice template
DEPRECATION_NOTICE = {
    "deprecation_warning": "This endpoint uses a deprecated singular table name. Please update your code to use the plural form.",
    "migration_date": "2025-04-18",
    "end_of_support_date": "2025-07-18",
    "documentation_url": "https://docs.cryoprotect.com/api/migration-to-plural"
}

def is_singular_endpoint(path: str) -> bool:
    """
    Check if the current API endpoint uses a singular table name.
    
    Args:
        path (str): The request path
        
    Returns:
        bool: True if the endpoint uses a singular table name
    """
    parts = path.strip('/').split('/')
    if len(parts) < 2:
        return False
    
    # Check if the endpoint starts with a singular table name
    # Example: /api/v1/molecule/123 -> True
    # Example: /api/v1/molecules/123 -> False
    return parts[1] in TABLE_MAPPINGS.keys()

def get_plural_table_name(singular_name: str) -> str:
    """
    Get the plural form of a table name.
    
    Args:
        singular_name (str): Singular table name
        
    Returns:
        str: Plural table name
    """
    return TABLE_MAPPINGS.get(singular_name, singular_name)

def get_singular_table_name(plural_name: str) -> str:
    """
    Get the singular form of a table name.
    
    Args:
        plural_name (str): Plural table name
        
    Returns:
        str: Singular table name
    """
    return REVERSE_MAPPINGS.get(plural_name, plural_name)

def log_deprecated_endpoint_usage(func: Callable) -> Callable:
    """
    Decorator to log usage of deprecated singular endpoints.
    
    Args:
        func (Callable): The API endpoint function
        
    Returns:
        Callable: Wrapped function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if is_singular_endpoint(request.path):
            # Log the deprecated endpoint usage
            logger.warning(
                "Deprecated endpoint used: %s, IP: %s, User-Agent: %s",
                request.path,
                request.remote_addr,
                request.headers.get('User-Agent', 'Unknown')
            )
            
            # If user ID is available, log it too
            if hasattr(g, 'user_id'):
                logger.warning("User ID: %s", g.user_id)
        
        return func(*args, **kwargs)
    
    return wrapper

def add_deprecation_notice(response_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Add a deprecation notice to the API response.
    
    Args:
        response_data (Dict[str, Any]): The original response data
        
    Returns:
        Dict[str, Any]: Response data with deprecation notice
    """
    if isinstance(response_data, dict):
        # Add deprecation notice to the response
        response_data["_deprecation"] = DEPRECATION_NOTICE
    elif isinstance(response_data, list):
        # For list responses, wrap in a container object
        response_data = {
            "data": response_data,
            "_deprecation": DEPRECATION_NOTICE
        }
    else:
        # For other types, convert to a dict
        response_data = {
            "data": response_data,
            "_deprecation": DEPRECATION_NOTICE
        }
    
    return response_data

def inject_deprecation_notice(func: Callable) -> Callable:
    """
    Decorator to inject deprecation notices into API responses for singular endpoints.
    
    Args:
        func (Callable): The API endpoint function
        
    Returns:
        Callable: Wrapped function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        
        # If this is a singular endpoint, add deprecation notice
        if is_singular_endpoint(request.path):
            if isinstance(result, tuple) and len(result) >= 2:
                # Handle (data, status_code) or (data, status_code, headers) responses
                data, status_code = result[0], result[1]
                
                # Add deprecation notice to the data
                modified_data = add_deprecation_notice(data)
                
                # Reconstruct the response tuple
                if len(result) > 2:
                    return (modified_data, status_code) + result[2:]
                else:
                    return modified_data, status_code
            else:
                # Handle direct response objects
                return add_deprecation_notice(result)
        
        return result
    
    return wrapper

def handle_table_name_compatibility(func: Callable) -> Callable:
    """
    Decorator to handle both singular and plural table names in database queries.
    
    This decorator modifies the function arguments to use the correct table name
    based on the current database schema.
    
    Args:
        func (Callable): The database query function
        
    Returns:
        Callable: Wrapped function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Check if 'table_name' is in kwargs
        if 'table_name' in kwargs:
            original_table = kwargs['table_name']
            
            # If it's a singular name, convert to plural
            if original_table in TABLE_MAPPINGS:
                kwargs['table_name'] = TABLE_MAPPINGS[original_table]
                
                # Log the table name conversion
                logger.debug(
                    "Converting table name from %s to %s",
                    original_table,
                    kwargs['table_name']
                )
        
        try:
            return func(*args, **kwargs)
        except Exception as e:
            # If the error might be related to table names, try with the alternative name
            if 'relation' in str(e) and 'does not exist' in str(e):
                if 'table_name' in kwargs:
                    current_table = kwargs['table_name']
                    
                    # If we're using plural, try singular as fallback
                    if current_table in REVERSE_MAPPINGS:
                        kwargs['table_name'] = REVERSE_MAPPINGS[current_table]
                        logger.warning(
                            "Table %s not found, trying fallback to %s",
                            current_table,
                            kwargs['table_name']
                        )
                        return func(*args, **kwargs)
            
            # Re-raise the original exception if we couldn't handle it
            raise
    
    return wrapper

def compatibility_layer(func: Callable) -> Callable:
    """
    Combined decorator that applies all compatibility features:
    - Logs deprecated endpoint usage
    - Injects deprecation notices
    - Handles table name compatibility
    
    Args:
        func (Callable): The API endpoint function
        
    Returns:
        Callable: Wrapped function with full compatibility support
    """
    @functools.wraps(func)
    @log_deprecated_endpoint_usage
    @inject_deprecation_notice
    @handle_table_name_compatibility
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    
    return wrapper

def update_query_table_names(query: str) -> str:
    """
    Update SQL query to use plural table names.
    
    Args:
        query (str): Original SQL query
        
    Returns:
        str: Updated SQL query
    """
    updated_query = query
    
    # Replace table names in FROM and JOIN clauses
    for singular, plural in TABLE_MAPPINGS.items():
        # Match "FROM singular" or "JOIN singular"
        updated_query = updated_query.replace(f" FROM {singular}", f" FROM {plural}")
        updated_query = updated_query.replace(f" JOIN {singular}", f" JOIN {plural}")
        
        # Match "FROM singular " or "JOIN singular "
        updated_query = updated_query.replace(f" FROM {singular} ", f" FROM {plural} ")
        updated_query = updated_query.replace(f" JOIN {singular} ", f" JOIN {plural} ")
        
        # Match "FROM singular\n" or "JOIN singular\n"
        updated_query = updated_query.replace(f" FROM {singular}\n", f" FROM {plural}\n")
        updated_query = updated_query.replace(f" JOIN {singular}\n", f" JOIN {plural}\n")
        
        # Match "table.column" references
        updated_query = updated_query.replace(f"{singular}.", f"{plural}.")
    
    return updated_query

def get_table_with_compatibility(table_name: str) -> str:
    """
    Get the appropriate table name with compatibility handling.
    
    This function is used when constructing Supabase queries to ensure
    the correct table name is used regardless of whether the code
    is using singular or plural names.
    
    Args:
        table_name (str): Original table name (singular or plural)
        
    Returns:
        str: The correct table name to use
    """
    # If it's a singular name, convert to plural
    if table_name in TABLE_MAPPINGS:
        return TABLE_MAPPINGS[table_name]
    
    # If it's already plural, return as is
    return table_name

def track_deprecated_usage(endpoint: str, user_id: Optional[str] = None) -> None:
    """
    Track usage of deprecated endpoints for analytics and monitoring.
    
    Args:
        endpoint (str): The API endpoint
        user_id (Optional[str]): User ID if available
    """
    try:
        # Log to the dedicated deprecated usage log
        usage_data = {
            "timestamp": datetime.now().isoformat(),
            "endpoint": endpoint,
            "user_id": user_id or "anonymous",
            "ip_address": request.remote_addr,
            "user_agent": request.headers.get('User-Agent', 'Unknown')
        }
        
        logger.info("DEPRECATED_USAGE: %s", json.dumps(usage_data))
        
        # Optionally, you could also send this data to a monitoring service
        # or store it in a database for further analysis
    except Exception as e:
        # Don't let tracking errors affect the API
        logger.error("Error tracking deprecated usage: %s", str(e))
"""
CryoProtect Analyzer API - Pagination Utilities

This module provides utilities for implementing standardized pagination across API endpoints.
It includes functions for:
- Parsing pagination parameters from requests
- Calculating pagination metadata (total pages, next/prev links)
- Formatting paginated responses with consistent metadata
- Handling edge cases and validation

All API endpoints that return potentially large lists should use these utilities
to ensure consistent pagination behavior and response formats.
"""

from flask import request, url_for
from typing import Dict, List, Any, Tuple, Optional, Union
import math
from urllib.parse import urlencode


def parse_pagination_params(
    default_limit: int = 100, 
    max_limit: int = 500
) -> Tuple[int, int]:
    """
    Parse pagination parameters (limit and offset) from the request query string.
    
    Args:
        default_limit (int): Default number of items per page if not specified
        max_limit (int): Maximum allowed limit to prevent excessive server load
        
    Returns:
        Tuple[int, int]: A tuple containing (limit, offset)
        
    Raises:
        ValueError: If pagination parameters are invalid
    """
    try:
        # Get pagination parameters with defaults
        limit = request.args.get("limit", default_limit, type=int)
        offset = request.args.get("offset", 0, type=int)
        page = request.args.get("page", None, type=int)
        per_page = request.args.get("per_page", None, type=int)
        
        # Handle both limit/offset and page/per_page styles
        if page is not None and per_page is not None:
            # Convert page/per_page to limit/offset
            limit = min(per_page, max_limit)
            offset = (page - 1) * limit if page > 0 else 0
        else:
            # Use limit/offset directly
            limit = min(limit, max_limit)
        
        # Validate parameters
        if limit <= 0:
            raise ValueError("Limit must be a positive integer")
        if offset < 0:
            raise ValueError("Offset must be a non-negative integer")
            
        return limit, offset
        
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid pagination parameters: {str(e)}")


def get_pagination_metadata(
    total_count: int,
    limit: int,
    offset: int,
    resource_endpoint: str,
    query_params: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Generate pagination metadata including total items, pages, and navigation links.
    
    Args:
        total_count (int): Total number of items available
        limit (int): Number of items per page
        offset (int): Current offset
        resource_endpoint (str): API endpoint for generating links
        query_params (Dict[str, Any], optional): Additional query parameters to include in links
        
    Returns:
        Dict[str, Any]: Pagination metadata dictionary
    """
    # Calculate current page and total pages
    current_page = math.floor(offset / limit) + 1 if limit > 0 else 1
    total_pages = math.ceil(total_count / limit) if limit > 0 else 1
    
    # Prepare base query parameters
    base_params = query_params.copy() if query_params else {}
    base_params["limit"] = limit
    
    # Generate links
    links = {}
    
    # First page
    first_params = base_params.copy()
    first_params["offset"] = 0
    links["first"] = f"{resource_endpoint}?{urlencode(first_params)}"
    
    # Previous page
    if offset > 0:
        prev_params = base_params.copy()
        prev_params["offset"] = max(0, offset - limit)
        links["prev"] = f"{resource_endpoint}?{urlencode(prev_params)}"
    
    # Next page
    if offset + limit < total_count:
        next_params = base_params.copy()
        next_params["offset"] = offset + limit
        links["next"] = f"{resource_endpoint}?{urlencode(next_params)}"
    
    # Last page
    last_params = base_params.copy()
    last_offset = (total_pages - 1) * limit if total_pages > 0 else 0
    last_params["offset"] = last_offset
    links["last"] = f"{resource_endpoint}?{urlencode(last_params)}"
    
    # Build metadata object
    metadata = {
        "total_items": total_count,
        "total_pages": total_pages,
        "current_page": current_page,
        "per_page": limit,
        "offset": offset,
        "links": links
    }
    
    return metadata


def paginate_query_response(
    query_result: List[Dict[str, Any]],
    total_count: int,
    limit: int,
    offset: int,
    resource_endpoint: str,
    query_params: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Format a paginated response with data and pagination metadata.
    
    Args:
        query_result (List[Dict[str, Any]]): The query results to paginate
        total_count (int): Total number of items available
        limit (int): Number of items per page
        offset (int): Current offset
        resource_endpoint (str): API endpoint for generating links
        query_params (Dict[str, Any], optional): Additional query parameters to include in links
        
    Returns:
        Dict[str, Any]: Formatted response with data and pagination metadata
    """
    # Get pagination metadata
    pagination = get_pagination_metadata(
        total_count, limit, offset, resource_endpoint, query_params
    )
    
    # Format response
    response = {
        "data": query_result,
        "pagination": pagination
    }
    
    return response


def get_total_count(supabase, table_name: str, filters: Optional[Dict[str, Any]] = None) -> int:
    """
    Get the total count of items in a table, optionally with filters.
    
    Args:
        supabase: Supabase client instance
        table_name (str): Name of the table or view to count
        filters (Dict[str, Any], optional): Filters to apply before counting
        
    Returns:
        int: Total count of items
    """
    try:
        # Start with a count query
        query = supabase.table(table_name).select("*", count="exact")
        
        # Apply filters if provided
        if filters:
            for key, value in filters.items():
                if isinstance(value, list):
                    query = query.in_(key, value)
                else:
                    query = query.eq(key, value)
        
        # Execute query
        response = query.execute()
        
        # Extract count from response
        count = response.count
        
        return count
    except Exception as e:
        # Log the error and return 0 as fallback
        print(f"Error getting count from {table_name}: {str(e)}")
        return 0
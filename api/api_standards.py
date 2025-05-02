"""
CryoProtect Analyzer API - Standardization Utilities

This module provides utilities for standardizing API responses, error handling,
and HTTP status codes across all API endpoints. It ensures consistent response
formats and proper status code usage throughout the application.
"""

from typing import Dict, Any, Tuple, Optional, List, Union
from datetime import datetime
import logging
from flask import jsonify, current_app, g, request
from http import HTTPStatus

# Set up logging
logger = logging.getLogger(__name__)

# Standard HTTP status codes with descriptions
HTTP_STATUS_CODES = {
    # 2xx - Success
    HTTPStatus.OK: "Request succeeded",
    HTTPStatus.CREATED: "Resource created successfully",
    HTTPStatus.ACCEPTED: "Request accepted for processing",
    HTTPStatus.NO_CONTENT: "Request succeeded with no content to return",
    
    # 4xx - Client errors
    HTTPStatus.BAD_REQUEST: "Invalid request parameters or data",
    HTTPStatus.UNAUTHORIZED: "Authentication required",
    HTTPStatus.FORBIDDEN: "Permission denied",
    HTTPStatus.NOT_FOUND: "Resource not found",
    HTTPStatus.METHOD_NOT_ALLOWED: "HTTP method not allowed for this resource",
    HTTPStatus.CONFLICT: "Request conflicts with current state of the resource",
    HTTPStatus.GONE: "Resource no longer available",
    HTTPStatus.UNPROCESSABLE_ENTITY: "Request data validation failed",
    HTTPStatus.TOO_MANY_REQUESTS: "Rate limit exceeded",
    
    # 5xx - Server errors
    HTTPStatus.INTERNAL_SERVER_ERROR: "Internal server error",
    HTTPStatus.NOT_IMPLEMENTED: "Functionality not implemented",
    HTTPStatus.BAD_GATEWAY: "Invalid response from upstream server",
    HTTPStatus.SERVICE_UNAVAILABLE: "Service temporarily unavailable",
    HTTPStatus.GATEWAY_TIMEOUT: "Upstream server timeout"
}

# Map common exceptions to appropriate HTTP status codes
EXCEPTION_STATUS_MAPPING = {
    ValueError: HTTPStatus.BAD_REQUEST,
    TypeError: HTTPStatus.BAD_REQUEST,
    KeyError: HTTPStatus.BAD_REQUEST,
    AttributeError: HTTPStatus.BAD_REQUEST,
    LookupError: HTTPStatus.NOT_FOUND,
    FileNotFoundError: HTTPStatus.NOT_FOUND,
    PermissionError: HTTPStatus.FORBIDDEN,
    NotImplementedError: HTTPStatus.NOT_IMPLEMENTED,
    TimeoutError: HTTPStatus.GATEWAY_TIMEOUT,
    ConnectionError: HTTPStatus.SERVICE_UNAVAILABLE
}

def format_timestamp(timestamp=None):
    """
    Format a timestamp in ISO format.
    
    Args:
        timestamp: Timestamp to format (default: current time)
        
    Returns:
        Formatted timestamp string
    """
    if timestamp is None:
        timestamp = datetime.now()
    return timestamp.isoformat()

def create_standard_response(
    data: Any = None,
    message: str = None,
    status_code: int = HTTPStatus.OK,
    meta: Dict[str, Any] = None,
    errors: List[Dict[str, Any]] = None,
    pagination: Dict[str, Any] = None
) -> Tuple[Dict[str, Any], int]:
    """
    Create a standardized API response.
    
    Args:
        data: The response data
        message: A human-readable message
        status_code: HTTP status code
        meta: Additional metadata
        errors: List of error details
        pagination: Pagination information
        
    Returns:
        Tuple of (response_dict, status_code)
    """
    # Initialize response
    response = {
        "status": "success" if status_code < 400 else "error",
        "timestamp": format_timestamp(),
        "code": status_code,
        "message": message or HTTP_STATUS_CODES.get(status_code, "")
    }
    
    # Add data if provided
    if data is not None:
        response["data"] = data
    
    # Add errors if provided
    if errors:
        response["errors"] = errors
        
    # Add metadata if provided
    if meta:
        response["meta"] = meta
        
    # Add pagination if provided
    if pagination:
        response["pagination"] = pagination
        
    return response, status_code

def create_error_response(
    error: Union[Exception, str],
    status_code: int = None,
    context: str = None,
    details: Dict[str, Any] = None
) -> Tuple[Dict[str, Any], int]:
    """
    Create a standardized error response.
    
    Args:
        error: Exception or error message
        status_code: HTTP status code (will be inferred from exception if not provided)
        context: Additional context about where the error occurred
        details: Additional error details
        
    Returns:
        Tuple of (response_dict, status_code)
    """
    # Determine error message and type
    if isinstance(error, Exception):
        error_message = str(error)
        error_type = error.__class__.__name__
        
        # Determine status code from exception if not provided
        if status_code is None:
            status_code = EXCEPTION_STATUS_MAPPING.get(
                error.__class__, HTTPStatus.INTERNAL_SERVER_ERROR
            )
    else:
        error_message = str(error)
        error_type = "Error"
        
        # Default to internal server error if status code not provided
        if status_code is None:
            status_code = HTTPStatus.INTERNAL_SERVER_ERROR
    
    # Create error details
    error_details = {
        "type": error_type,
        "message": error_message
    }
    
    # Add context if provided
    if context:
        error_details["context"] = context
        
    # Add additional details if provided
    if details:
        error_details.update(details)
        
    # Log the error
    log_message = f"{error_type}: {error_message}"
    if context:
        log_message = f"{context} - {log_message}"
        
    if status_code >= 500:
        logger.error(log_message, exc_info=isinstance(error, Exception))
    elif status_code >= 400:
        logger.warning(log_message)
    else:
        logger.info(log_message)
        
    # Create and return the error response
    return create_standard_response(
        message=HTTP_STATUS_CODES.get(status_code, "An error occurred"),
        status_code=status_code,
        errors=[error_details]
    )

def create_success_response(
    data: Any = None,
    message: str = None,
    status_code: int = HTTPStatus.OK,
    meta: Dict[str, Any] = None,
    pagination: Dict[str, Any] = None
) -> Tuple[Dict[str, Any], int]:
    """
    Create a standardized success response.
    
    Args:
        data: The response data
        message: A human-readable message
        status_code: HTTP status code
        meta: Additional metadata
        pagination: Pagination information
        
    Returns:
        Tuple of (response_dict, status_code)
    """
    return create_standard_response(
        data=data,
        message=message or HTTP_STATUS_CODES.get(status_code, "Request succeeded"),
        status_code=status_code,
        meta=meta,
        pagination=pagination
    )

def create_pagination_info(
    page: int,
    per_page: int,
    total_count: int,
    base_url: str = None
) -> Dict[str, Any]:
    """
    Create standardized pagination information.
    
    Args:
        page: Current page number
        per_page: Number of items per page
        total_count: Total number of items
        base_url: Base URL for pagination links
        
    Returns:
        Dictionary with pagination information
    """
    total_pages = (total_count + per_page - 1) // per_page if per_page > 0 else 0
    
    pagination = {
        "page": page,
        "per_page": per_page,
        "total_items": total_count,
        "total_pages": total_pages,
        "has_next": page < total_pages,
        "has_prev": page > 1
    }
    
    # Add links if base_url is provided
    if base_url:
        links = {}
        
        # Current page
        links["self"] = f"{base_url}?page={page}&per_page={per_page}"
        
        # First page
        links["first"] = f"{base_url}?page=1&per_page={per_page}"
        
        # Last page
        links["last"] = f"{base_url}?page={total_pages}&per_page={per_page}"
        
        # Next page
        if pagination["has_next"]:
            links["next"] = f"{base_url}?page={page + 1}&per_page={per_page}"
            
        # Previous page
        if pagination["has_prev"]:
            links["prev"] = f"{base_url}?page={page - 1}&per_page={per_page}"
            
        pagination["links"] = links
        
    return pagination

def get_pagination_params(request):
    """
    Extract pagination parameters from request.
    
    Args:
        request: Flask request object
        
    Returns:
        Dictionary with page and per_page values
    """
    try:
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 20))
        
        # Validate and constrain values
        page = max(1, page)  # Minimum page is 1
        per_page = max(1, min(100, per_page))  # Between 1 and 100
        
        return {
            'page': page,
            'per_page': per_page,
            'offset': (page - 1) * per_page,
            'limit': per_page
        }
    except (ValueError, TypeError):
        # Default values if parsing fails
        return {
            'page': 1,
            'per_page': 20,
            'offset': 0,
            'limit': 20
        }

def jsonify_standard_response(response_data, status_code):
    """
    Convert a standard response to a Flask JSON response.
    
    Args:
        response_data: Response data dictionary
        status_code: HTTP status code
        
    Returns:
        Flask JSON response with appropriate status code
    """
    return jsonify(response_data), status_code
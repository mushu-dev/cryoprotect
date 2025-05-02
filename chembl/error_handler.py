from enum import Enum
from typing import Dict, Any, Optional, Tuple
import requests
import time
import logging

logger = logging.getLogger(__name__)

class ErrorCategory(Enum):
    """Categories of errors that can occur during ChEMBL API operations"""
    API_RATE_LIMIT = "rate_limit"       # API rate limiting (429 responses)
    CONNECTION_ERROR = "connection"     # Network/connection issues
    API_SERVER_ERROR = "server_error"   # Server-side errors (5xx)
    API_CLIENT_ERROR = "client_error"   # Client-side errors (4xx except 429)
    DATA_VALIDATION = "validation"      # Data validation errors
    TRANSFORMATION = "transformation"   # Data transformation errors
    PARSING_ERROR = "parsing"           # JSON/XML parsing errors
    DATABASE_ERROR = "database"         # Database operation errors
    UNKNOWN = "unknown"                 # Uncategorized errors

def classify_error(exception: Exception) -> Tuple[ErrorCategory, str]:
    """
    Classify an exception into an error category and provide a description.
    
    Args:
        exception: The exception to classify
        
    Returns:
        Tuple containing (ErrorCategory, description)
    """
    # Handle requests/HTTP exceptions
    if isinstance(exception, requests.exceptions.RequestException):
        if isinstance(exception, requests.exceptions.Timeout):
            return ErrorCategory.CONNECTION_ERROR, "Request timed out"
        
        if isinstance(exception, requests.exceptions.ConnectionError):
            return ErrorCategory.CONNECTION_ERROR, "Connection error"
            
        if isinstance(exception, requests.exceptions.HTTPError):
            response = exception.response
            if response is not None:
                status_code = response.status_code
                
                if status_code == 429:
                    return ErrorCategory.API_RATE_LIMIT, f"Rate limit exceeded: {status_code}"
                    
                if 400 <= status_code < 500:
                    return ErrorCategory.API_CLIENT_ERROR, f"Client error: {status_code}"
                    
                if 500 <= status_code < 600:
                    return ErrorCategory.API_SERVER_ERROR, f"Server error: {status_code}"
    
    # Handle JSON parsing errors
    if isinstance(exception, (ValueError, TypeError)) and "JSON" in str(exception):
        return ErrorCategory.PARSING_ERROR, f"JSON parsing error: {str(exception)}"
    
    # Handle database errors (psycopg2)
    try:
        import psycopg2
        if isinstance(exception, psycopg2.Error):
            return ErrorCategory.DATABASE_ERROR, f"Database error: {str(exception)}"
    except ImportError:
        # psycopg2 not available, skip this check
        pass
    
    # Handle validation errors (often custom exceptions)
    if "validation" in str(exception).lower() or "invalid" in str(exception).lower():
        return ErrorCategory.DATA_VALIDATION, f"Validation error: {str(exception)}"
    
    # Handle transformation errors
    if "transform" in str(exception).lower() or "conversion" in str(exception).lower():
        return ErrorCategory.TRANSFORMATION, f"Transformation error: {str(exception)}"
    
    # Default case
    return ErrorCategory.UNKNOWN, f"Unclassified error: {str(exception)}"


def get_recovery_strategy(category: ErrorCategory, error_details: str = "", attempt: int = 0) -> str:
    """
    Determine appropriate recovery action for an error category.
    
    Args:
        category: The error category
        error_details: Details about the error (optional)
        attempt: Current attempt number (0-based, optional)
        
    Returns:
        Recovery action: 'RETRY', 'SKIP', 'ABORT', or 'LOG_ONLY'
    """
    max_retries = 5
    
    # Check if max retries exceeded for retryable errors
    if attempt >= max_retries:
        return "ABORT"
    
    # Determine strategy based on category
    if category == ErrorCategory.API_RATE_LIMIT:
        return "RETRY"
        
    elif category == ErrorCategory.CONNECTION_ERROR:
        return "RETRY"
        
    elif category == ErrorCategory.API_SERVER_ERROR:
        return "RETRY"
        
    elif category == ErrorCategory.API_CLIENT_ERROR:
        # Client errors are typically not recoverable
        return "SKIP"
        
    elif category == ErrorCategory.PARSING_ERROR:
        # Parsing errors might be temporary or due to malformed responses
        return "RETRY"
        
    elif category == ErrorCategory.DATA_VALIDATION:
        # Data validation errors indicate issues with the input data
        return "SKIP"
        
    elif category == ErrorCategory.TRANSFORMATION:
        # Transformation errors might be due to unexpected data formats
        return "SKIP"
        
    elif category == ErrorCategory.DATABASE_ERROR:
        # Database errors might be temporary
        return "RETRY"
        
    # Default strategy for unknown errors
    return "LOG_ONLY"
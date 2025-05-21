#!/usr/bin/env python3
"""
CryoProtect v2 - Supabase Adapter

This module provides a SupabaseAdapter class that normalizes Supabase SQL query responses
into a standardized format, making it easier to handle different response structures
consistently throughout the application.
"""

import logging
from typing import Any, Dict, List, Optional, Union

# Import logging configuration
from logging_config import setup_logging

# Set up logging
setup_logging()
logger = logging.getLogger(__name__)


class SupabaseAdapter:
    """
    Adapter class for normalizing Supabase responses into a standardized format.
    
    This class provides methods to convert various Supabase response formats into
    a consistent structure with 'success', 'data', and 'error' keys, making it
    easier to handle responses throughout the application.
    """
    
    @staticmethod
    def normalize_response(raw_response: Any) -> Dict[str, Any]:
        """
        Normalize a raw Supabase response into a standardized format.
        
        Args:
            raw_response: The raw response from a Supabase query, which could be
                          in various formats depending on the query and response type.
        
        Returns:
            A dictionary with the following keys:
                - 'success': bool - True if the query succeeded, False otherwise
                - 'data': list or None - The result rows, or None if error
                - 'error': str or None - Error message if any, else None
                
        Examples:
            >>> adapter = SupabaseAdapter()
            >>> result = adapter.normalize_response({'data': [{'id': 1}], 'error': None})
            >>> result
            {'success': True, 'data': [{'id': 1}], 'error': None}
            
            >>> result = adapter.normalize_response({'error': {'message': 'Not found'}})
            >>> result
            {'success': False, 'data': None, 'error': 'Not found'}
        """
        logger.debug(f"Normalizing Supabase response: {raw_response}")
        
        # Initialize the normalized response
        normalized = {
            'success': False,
            'data': None,
            'error': None
        }
        
        # Handle None response
        if raw_response is None:
            normalized['error'] = "Empty response received"
            logger.warning("Empty response received from Supabase")
            return normalized
            
        # Handle response object with data/error attributes (common Supabase pattern)
        if hasattr(raw_response, 'data') and hasattr(raw_response, 'error'):
            data = getattr(raw_response, 'data')
            error = getattr(raw_response, 'error')
            
            if error:
                # Extract error message if it's in a nested structure
                error_msg = SupabaseAdapter._extract_error_message(error)
                normalized['error'] = error_msg
                logger.error(f"Supabase query error: {error_msg}")
            else:
                normalized['success'] = True
                normalized['data'] = data
            
            return normalized
        
        # Handle dictionary response
        if isinstance(raw_response, dict):
            # Case: {'data': [...], 'error': None} or {'data': [...], 'error': {...}}
            if 'data' in raw_response:
                data = raw_response.get('data')
                error = raw_response.get('error')
                
                if error:
                    error_msg = SupabaseAdapter._extract_error_message(error)
                    normalized['error'] = error_msg
                    logger.error(f"Supabase query error: {error_msg}")
                else:
                    normalized['success'] = True
                    normalized['data'] = data
                
                return normalized
            
            # Case: {'result': [...], 'status': 'ok'}
            if 'result' in raw_response:
                status = raw_response.get('status')
                result = raw_response.get('result')
                
                if status == 'ok' or status is None:
                    normalized['success'] = True
                    normalized['data'] = result
                else:
                    normalized['error'] = f"Query failed with status: {status}"
                    logger.error(f"Supabase query failed with status: {status}")
                
                return normalized
            
            # Case: {'error': {...}} or {'message': 'error', ...}
            if 'error' in raw_response or 'message' in raw_response:
                error = raw_response.get('error')
                message = raw_response.get('message')
                
                if error:
                    error_msg = SupabaseAdapter._extract_error_message(error)
                    normalized['error'] = error_msg
                elif message:
                    normalized['error'] = message
                
                logger.error(f"Supabase query error: {normalized['error']}")
                return normalized
            
            # If it's a dictionary but doesn't match known patterns, assume it's data
            logger.warning(f"Unrecognized response format, treating as data: {raw_response}")
            normalized['success'] = True
            normalized['data'] = [raw_response]  # Wrap in list for consistency
            return normalized
        
        # Handle list response (assume it's data)
        if isinstance(raw_response, list):
            normalized['success'] = True
            normalized['data'] = raw_response
            return normalized
        
        # Handle any other type of response
        logger.warning(f"Unexpected response type: {type(raw_response)}")
        normalized['error'] = f"Unexpected response type: {type(raw_response)}"
        return normalized
    
    @staticmethod
    def _extract_error_message(error: Any) -> str:
        """
        Extract a readable error message from various error formats.
        
        Args:
            error: The error object or message from a Supabase response
            
        Returns:
            A string containing the error message
        """
        if error is None:
            return "Unknown error"
            
        if isinstance(error, str):
            return error
            
        if isinstance(error, dict):
            # Try common error message fields
            for field in ['message', 'msg', 'detail', 'description', 'error']:
                if field in error and error[field]:
                    return str(error[field])
            
            # If no specific message field found, return the whole dict as string
            return str(error)
        
        # For any other type, convert to string
        return str(error)
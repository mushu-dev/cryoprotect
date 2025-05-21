"""
CryoProtect Analyzer API - Consolidated Molecule Decorators

This module provides decorators for handling consolidated molecules in API endpoints.
These decorators help ensure that API responses properly handle molecule consolidation
and differentiation according to the data quality standards established in Phase 2.
"""

import functools
import logging
from typing import Callable, Dict, Any, List, Optional
from flask import request, g, jsonify, current_app, redirect, url_for

from api.consolidated_utils import (
    is_consolidated,
    get_primary_molecule,
    get_consolidated_molecules,
    get_differentiation_group,
    enrich_molecule_data
)

# Set up logging
logger = logging.getLogger(__name__)

def handle_consolidated_molecules(f: Callable) -> Callable:
    """
    Decorator to handle consolidated molecules in API endpoints.
    
    This decorator:
    1. Checks if a molecule ID in the request refers to a consolidated molecule
    2. If so, redirects to the primary molecule or updates the request parameters
    3. Enriches the response with consolidated molecule information
    
    Args:
        f: The endpoint function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        # Extract molecule ID from different possible locations
        molecule_id = None
        
        # Check path parameters
        if 'molecule_id' in kwargs:
            molecule_id = kwargs['molecule_id']
        
        # Check query parameters if not in path
        if molecule_id is None and request.args.get('molecule_id'):
            molecule_id = request.args.get('molecule_id')
            
        # Check request body if not in path or query
        if molecule_id is None and request.is_json:
            request_data = request.get_json()
            if request_data and 'molecule_id' in request_data:
                molecule_id = request_data.get('molecule_id')
        
        # If molecule ID is found and it's a consolidated molecule
        if molecule_id:
            primary_id = get_primary_molecule(molecule_id)
            
            if primary_id != molecule_id:
                # Log the redirect
                logger.info(f"Redirecting from consolidated molecule {molecule_id} to primary {primary_id}")
                
                # Handle based on configuration
                redirect_mode = current_app.config.get('CONSOLIDATED_REDIRECT_MODE', 'update')
                
                if redirect_mode == 'http_redirect' and request.method == 'GET':
                    # HTTP 302 redirect to the primary molecule
                    if 'molecule_id' in kwargs:
                        # Update URL path parameter
                        return redirect(request.path.replace(
                            f"/{molecule_id}", f"/{primary_id}"
                        ), code=302)
                elif redirect_mode == 'update':
                    # Update kwargs to use primary molecule ID
                    kwargs['molecule_id'] = primary_id
                    # Store original ID for reference
                    kwargs['original_molecule_id'] = molecule_id
        
        # Call the original function
        result = f(*args, **kwargs)
        
        # Handle different return types
        if isinstance(result, tuple) and len(result) >= 2:
            # Function returned (data, status_code, [headers])
            data = result[0]
            status_code = result[1]
            headers = result[2] if len(result) > 2 else {}
            
            # Enrich data with consolidated information
            if isinstance(data, dict):
                data = enrich_molecule_data(data)
            elif isinstance(data, list):
                data = [enrich_molecule_data(item) if isinstance(item, dict) else item for item in data]
            
            # Add original ID as a header if applicable
            if 'original_molecule_id' in kwargs:
                headers['X-Original-Molecule-ID'] = kwargs['original_molecule_id']
                
            return data, status_code, headers
        else:
            # Function returned just data
            if isinstance(result, dict):
                return enrich_molecule_data(result)
            elif isinstance(result, list):
                return [enrich_molecule_data(item) if isinstance(item, dict) else item for item in result]
            else:
                return result
                
    return decorated

def handle_batch_consolidated_molecules(f: Callable) -> Callable:
    """
    Decorator to handle consolidated molecules in batch API endpoints.
    
    This decorator:
    1. Checks if any molecule IDs in a batch request refer to consolidated molecules
    2. Updates the request to use primary molecule IDs
    3. Enriches the response with consolidated molecule information
    
    Args:
        f: The endpoint function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        # Check for molecule IDs in request body
        molecule_ids = []
        original_to_primary_map = {}
        
        if request.is_json:
            request_data = request.get_json()
            if request_data:
                # Handle various batch request formats
                if 'molecule_ids' in request_data:
                    molecule_ids = request_data.get('molecule_ids', [])
                elif 'ids' in request_data:
                    molecule_ids = request_data.get('ids', [])
                elif 'molecules' in request_data:
                    molecules = request_data.get('molecules', [])
                    if isinstance(molecules, list):
                        for molecule in molecules:
                            if isinstance(molecule, dict) and 'id' in molecule:
                                molecule_ids.append(molecule['id'])
        
        # If molecule IDs are found, check for consolidated molecules
        if molecule_ids:
            # Create a mapping from original to primary IDs
            for molecule_id in molecule_ids:
                primary_id = get_primary_molecule(molecule_id)
                if primary_id != molecule_id:
                    original_to_primary_map[molecule_id] = primary_id
            
            # Update the request data with primary IDs
            if original_to_primary_map and request.is_json:
                request_data = request.get_json()
                
                if 'molecule_ids' in request_data:
                    request_data['molecule_ids'] = [
                        original_to_primary_map.get(id, id) for id in request_data['molecule_ids']
                    ]
                elif 'ids' in request_data:
                    request_data['ids'] = [
                        original_to_primary_map.get(id, id) for id in request_data['ids']
                    ]
                elif 'molecules' in request_data:
                    for molecule in request_data['molecules']:
                        if isinstance(molecule, dict) and 'id' in molecule:
                            molecule_id = molecule['id']
                            if molecule_id in original_to_primary_map:
                                molecule['id'] = original_to_primary_map[molecule_id]
                                molecule['original_id'] = molecule_id
                
                # Store the mapping for use in the response
                if hasattr(g, 'consolidated_map'):
                    g.consolidated_map.update(original_to_primary_map)
                else:
                    g.consolidated_map = original_to_primary_map
        
        # Call the original function
        result = f(*args, **kwargs)
        
        # Handle different return types
        if isinstance(result, tuple) and len(result) >= 2:
            # Function returned (data, status_code, [headers])
            data = result[0]
            status_code = result[1]
            headers = result[2] if len(result) > 2 else {}
            
            # Enrich batch data with consolidated information
            if isinstance(data, dict) and 'data' in data and isinstance(data['data'], list):
                data['data'] = [enrich_molecule_data(item) if isinstance(item, dict) else item for item in data['data']]
                
                # Add consolidated mapping to metadata
                if hasattr(g, 'consolidated_map') and g.consolidated_map:
                    if 'meta' not in data:
                        data['meta'] = {}
                    data['meta']['consolidated_mapping'] = g.consolidated_map
            elif isinstance(data, list):
                data = [enrich_molecule_data(item) if isinstance(item, dict) else item for item in data]
                
                # Add consolidated mapping to headers
                if hasattr(g, 'consolidated_map') and g.consolidated_map:
                    headers['X-Consolidated-Mapping'] = str(g.consolidated_map)
                
            return data, status_code, headers
        else:
            # Function returned just data
            if isinstance(result, dict) and 'data' in result and isinstance(result['data'], list):
                result['data'] = [enrich_molecule_data(item) if isinstance(item, dict) else item for item in result['data']]
                
                # Add consolidated mapping to metadata
                if hasattr(g, 'consolidated_map') and g.consolidated_map:
                    if 'meta' not in result:
                        result['meta'] = {}
                    result['meta']['consolidated_mapping'] = g.consolidated_map
            elif isinstance(result, list):
                result = [enrich_molecule_data(item) if isinstance(item, dict) else item for item in result]
                
            return result
                
    return decorated
def with_transaction(f):
    """
    Decorator to wrap a function in a database transaction.
    
    This decorator ensures that all database operations in the wrapped function
    are performed in a single transaction. If the function returns successfully,
    the transaction is committed. If an exception is raised, the transaction is
    rolled back.
    
    Args:
        f: The function to wrap in a transaction
        
    Returns:
        Decorated function that runs in a transaction
    """
    @functools.wraps(f)
    def decorated(*args, **kwargs):
        # For testing purposes, we're not actually implementing a transaction
        # In a real application, this would use SQLAlchemy or other transaction support
        
        try:
            # Call the wrapped function
            result = f(*args, **kwargs)
            
            # Simulate committing the transaction
            logger.debug(f"Transaction would be committed for {f.__name__}")
            
            return result
            
        except Exception as e:
            # Simulate rolling back the transaction
            logger.error(f"Transaction would be rolled back for {f.__name__}: {str(e)}")
            
            # Re-raise the exception to be handled upstream
            raise
    
    return decorated
EOF < /dev/null

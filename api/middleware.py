"""
CryoProtect Analyzer API - Middleware

This module provides middleware functionality for the API, including functions
for resolving molecule IDs to their primary molecules and retrieving molecule
data with consolidated information.
"""

from flask import g, current_app, request
import logging
from typing import Tuple, List, Dict, Optional, Any, Union
import uuid

from database.adapter import get_connection

# Set up logging
logger = logging.getLogger(__name__)

def resolve_molecule_id(molecule_id: str) -> Tuple[str, bool]:
    """
    Resolve a molecule ID to its primary molecule ID if it's consolidated.
    
    Args:
        molecule_id: The molecule ID to resolve
        
    Returns:
        Tuple of (resolved_id, is_consolidated)
        - resolved_id: The primary molecule ID if consolidated, otherwise the original ID
        - is_consolidated: True if the molecule is consolidated, False otherwise
    """
    if not molecule_id:
        return None, False
        
    try:
        # Validate UUID format
        uuid.UUID(molecule_id)
    except ValueError:
        logger.warning(f"Invalid molecule ID format: {molecule_id}")
        return molecule_id, False
    
    try:
        conn = get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                id, 
                consolidated_to,
                consolidated_to IS NOT NULL as is_consolidated
            FROM molecules
            WHERE id = %s
        """
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        if not result:
            return molecule_id, False
            
        if result[2]:  # Is consolidated
            logger.debug(f"Molecule {molecule_id} is consolidated to {result[1]}")
            return result[1], True
        
        return result[0], False
    except Exception as e:
        logger.error(f"Error resolving molecule ID {molecule_id}: {str(e)}", exc_info=True)
        return molecule_id, False
    finally:
        if 'cursor' in locals() and cursor:
            cursor.close()

def resolve_molecule_ids(molecule_ids: List[str]) -> Dict[str, Tuple[str, bool]]:
    """
    Resolve multiple molecule IDs to their primary molecule IDs.
    
    Args:
        molecule_ids: List of molecule IDs to resolve
        
    Returns:
        Dictionary mapping original IDs to tuples of (resolved_id, is_consolidated)
    """
    results = {}
    
    # Default to original IDs if resolution fails
    for molecule_id in molecule_ids:
        results[molecule_id] = (molecule_id, False)
    
    if not molecule_ids:
        return results
    
    # Filter valid UUIDs
    valid_ids = []
    for molecule_id in molecule_ids:
        try:
            uuid.UUID(molecule_id)
            valid_ids.append(molecule_id)
        except ValueError:
            logger.warning(f"Invalid molecule ID format: {molecule_id}")
    
    if not valid_ids:
        return results
    
    try:
        conn = get_connection()
        cursor = conn.cursor()
        
        placeholders = ','.join(['%s'] * len(valid_ids))
        query = f"""
            SELECT 
                id, 
                consolidated_to,
                consolidated_to IS NOT NULL as is_consolidated
            FROM molecules
            WHERE id IN ({placeholders})
        """
        
        cursor.execute(query, valid_ids)
        for row in cursor.fetchall():
            molecule_id = row[0]
            consolidated_to = row[1]
            is_consolidated = row[2]
            
            if is_consolidated:
                results[molecule_id] = (consolidated_to, True)
            else:
                results[molecule_id] = (molecule_id, False)
        
        return results
    except Exception as e:
        logger.error(f"Error resolving molecule IDs: {str(e)}", exc_info=True)
        return results
    finally:
        if 'cursor' in locals() and cursor:
            cursor.close()

def get_molecule_with_consolidated_info(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule with consolidated molecule information.
    
    Args:
        molecule_id: The molecule ID to get
        
    Returns:
        Dictionary containing molecule data with consolidated information,
        or None if the molecule doesn't exist
    """
    resolved_id, is_consolidated = resolve_molecule_id(molecule_id)
    if not resolved_id:
        return None
    
    try:
        conn = get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                id, 
                name,
                molecular_formula,
                smiles,
                inchi_key,
                consolidated_to,
                (SELECT jsonb_agg(id) 
                 FROM molecules 
                 WHERE consolidated_to = m.id) as consolidated_molecules,
                created_at,
                updated_at
            FROM molecules m
            WHERE id = %s
        """
        
        cursor.execute(query, (resolved_id,))
        result = cursor.fetchone()
        
        if not result:
            return None
            
        # Convert to dictionary
        columns = [desc[0] for desc in cursor.description]
        molecule_data = dict(zip(columns, result))
        
        # Get differentiation information
        diff_query = """
            SELECT 
                pt.name as property_type,
                mp.property_value
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE mp.molecule_id = %s AND pt.name IN ('differentiationGroup', 'differentiationDescription')
        """
        
        cursor.execute(diff_query, (resolved_id,))
        diff_results = cursor.fetchall()
        
        for row in diff_results:
            property_type = row[0]
            property_value = row[1]
            
            if property_type == 'differentiationGroup':
                molecule_data['differentiation_group'] = property_value
            elif property_type == 'differentiationDescription':
                molecule_data['differentiation_description'] = property_value
        
        # Add consolidated info
        molecule_data['is_consolidated'] = is_consolidated
        
        # If this is a consolidated molecule (secondary), include the original ID
        if is_consolidated:
            molecule_data['original_molecule_id'] = molecule_id
        
        # Ensure consolidated_molecules is not null
        if molecule_data['consolidated_molecules'] is None:
            molecule_data['consolidated_molecules'] = []
        
        return molecule_data
    except Exception as e:
        logger.error(f"Error getting molecule with consolidated info {molecule_id}: {str(e)}", exc_info=True)
        return None
    finally:
        if 'cursor' in locals() and cursor:
            cursor.close()

def get_differentiation_group_members(group_id: str) -> List[Dict[str, Any]]:
    """
    Get all molecules in a differentiation group.
    
    Args:
        group_id: The differentiation group ID
        
    Returns:
        List of molecules in the differentiation group
    """
    if not group_id:
        return []
    
    try:
        conn = get_connection()
        cursor = conn.cursor()
        
        query = """
            SELECT 
                m.id,
                m.name,
                m.molecular_formula,
                m.smiles,
                m.inchi_key,
                mp2.property_value as differentiation_description
            FROM molecular_properties mp
            JOIN molecules m ON mp.molecule_id = m.id
            JOIN property_types pt ON mp.property_type_id = pt.id
            LEFT JOIN molecular_properties mp2 ON mp2.molecule_id = m.id
            LEFT JOIN property_types pt2 ON mp2.property_type_id = pt2.id AND pt2.name = 'differentiationDescription'
            WHERE pt.name = 'differentiationGroup'
            AND mp.property_value = %s
            ORDER BY m.name
        """
        
        cursor.execute(query, (group_id,))
        results = cursor.fetchall()
        
        members = []
        for row in results:
            # Convert to dictionary
            columns = [desc[0] for desc in cursor.description]
            member_data = dict(zip(columns, row))
            members.append(member_data)
        
        return members
    except Exception as e:
        logger.error(f"Error getting differentiation group members {group_id}: {str(e)}", exc_info=True)
        return []
    finally:
        if 'cursor' in locals() and cursor:
            cursor.close()

def handle_consolidated_molecules(view_func):
    """
    Decorator for Flask view functions to handle consolidated molecules.
    
    This decorator checks if a molecule ID in the request parameters or URL
    refers to a consolidated molecule, and if so, replaces it with the primary
    molecule ID. It also adds information about the consolidation to the response.
    
    Args:
        view_func: The Flask view function to decorate
        
    Returns:
        Decorated function
    """
    def wrapper(*args, **kwargs):
        # Check if molecule_id is in URL parameters
        molecule_id = kwargs.get('molecule_id')
        
        # If not in URL, check request parameters
        if not molecule_id and request.args.get('molecule_id'):
            molecule_id = request.args.get('molecule_id')
            
        # If not in parameters, check request body
        if not molecule_id and request.is_json:
            try:
                request_data = request.get_json()
                if request_data and 'molecule_id' in request_data:
                    molecule_id = request_data.get('molecule_id')
            except Exception:
                pass
        
        # If molecule ID is found, resolve it
        original_id = None
        if molecule_id:
            resolved_id, is_consolidated = resolve_molecule_id(molecule_id)
            
            if is_consolidated:
                # Save the original ID
                original_id = molecule_id
                
                # Replace with resolved ID
                if 'molecule_id' in kwargs:
                    kwargs['molecule_id'] = resolved_id
                elif request.args.get('molecule_id'):
                    request.args = request.args.copy()
                    request.args['molecule_id'] = resolved_id
                    request.args['original_molecule_id'] = original_id
        
        # Store the original ID for later use
        if original_id:
            g.original_molecule_id = original_id
        
        # Call the original function
        response = view_func(*args, **kwargs)
        
        # If the response is a tuple with data and status code
        if isinstance(response, tuple) and len(response) >= 2:
            data, status_code = response[0], response[1]
            
            # If response is JSON and original_id exists
            if isinstance(data, dict) and original_id:
                data.setdefault('meta', {})
                data['meta']['consolidated_from'] = original_id
            
            return response
        
        return response
    
    return wrapper

def molecule_batch_middleware(view_func):
    """
    Middleware for batch operations with molecules.
    
    This middleware resolves multiple molecule IDs in batch operations
    and ensures that consolidated molecules are correctly handled.
    
    Args:
        view_func: The Flask view function to decorate
        
    Returns:
        Decorated function
    """
    def wrapper(*args, **kwargs):
        if not request.is_json:
            return view_func(*args, **kwargs)
            
        try:
            request_data = request.get_json()
            
            # Look for molecule IDs in various formats
            molecule_ids = []
            
            if 'molecule_ids' in request_data:
                molecule_ids = request_data['molecule_ids']
            elif 'ids' in request_data:
                molecule_ids = request_data['ids']
            elif 'molecules' in request_data and isinstance(request_data['molecules'], list):
                for molecule in request_data['molecules']:
                    if isinstance(molecule, dict) and 'id' in molecule:
                        molecule_ids.append(molecule['id'])
            
            if not molecule_ids:
                return view_func(*args, **kwargs)
            
            # Resolve molecule IDs
            resolution_map = resolve_molecule_ids(molecule_ids)
            
            # Check if any molecules were resolved to a different ID
            if any(original_id != resolved[0] for original_id, resolved in resolution_map.items()):
                # Create a copy of the request data
                request._cached_json = request_data.copy()
                
                # Update molecule IDs in the request data
                if 'molecule_ids' in request_data:
                    request._cached_json['molecule_ids'] = [resolution_map[id][0] for id in request_data['molecule_ids']]
                elif 'ids' in request_data:
                    request._cached_json['ids'] = [resolution_map[id][0] for id in request_data['ids']]
                elif 'molecules' in request_data and isinstance(request_data['molecules'], list):
                    for i, molecule in enumerate(request_data['molecules']):
                        if isinstance(molecule, dict) and 'id' in molecule:
                            original_id = molecule['id']
                            resolved_id, is_consolidated = resolution_map.get(original_id, (original_id, False))
                            request._cached_json['molecules'][i]['id'] = resolved_id
                            if is_consolidated:
                                request._cached_json['molecules'][i]['original_id'] = original_id
                
                # Store the resolution map for use in the response
                g.molecule_resolution_map = {
                    original_id: resolved_id
                    for original_id, (resolved_id, is_consolidated) in resolution_map.items()
                    if original_id != resolved_id
                }
        except Exception as e:
            logger.error(f"Error in molecule_batch_middleware: {str(e)}", exc_info=True)
        
        # Call the original function
        response = view_func(*args, **kwargs)
        
        # Add resolution map to the response
        if hasattr(g, 'molecule_resolution_map') and g.molecule_resolution_map:
            if isinstance(response, tuple) and len(response) >= 2:
                data, status_code = response[0], response[1]
                
                if isinstance(data, dict):
                    data.setdefault('meta', {})
                    data['meta']['consolidated_molecules'] = g.molecule_resolution_map
                
                return response
        
        return response
    
    return wrapper
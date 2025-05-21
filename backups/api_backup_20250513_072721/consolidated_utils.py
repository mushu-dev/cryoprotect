"""
CryoProtect Analyzer API - Consolidated Molecule Utilities

This module provides utilities for handling consolidated molecules in the API.
It helps maintain data consistency when working with molecules that have been
consolidated or differentiated during the data quality enhancement process.
"""

import logging
from typing import Dict, List, Optional, Union, Any, Tuple
from flask import current_app, g
import uuid

# Set up logging
logger = logging.getLogger(__name__)

def get_db_connection():
    """
    Get a database connection.
    
    Returns:
        Database connection from the current context
    """
    if hasattr(g, 'db_conn'):
        return g.db_conn
    
    # Import here to avoid circular imports
    from database.adapter import get_connection
    g.db_conn = get_connection()
    return g.db_conn

def is_consolidated(molecule_id: str) -> bool:
    """
    Check if a molecule is consolidated (secondary molecule).
    
    Args:
        molecule_id: The molecule ID to check
        
    Returns:
        True if the molecule is consolidated (points to a primary molecule),
        False otherwise
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to check if molecule is consolidated
        query = """
        SELECT consolidated_to IS NOT NULL
        FROM molecules
        WHERE id = %s
        """
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        return result[0] if result is not None else False
        
    except Exception as e:
        logger.error(f"Error checking if molecule is consolidated: {str(e)}")
        return False

def get_primary_molecule(molecule_id: str) -> str:
    """
    Get the primary molecule ID for a given molecule ID.
    
    If the molecule is a secondary (consolidated) molecule, return the ID
    of its primary molecule. Otherwise, return the original ID.
    
    Args:
        molecule_id: The molecule ID to check
        
    Returns:
        Primary molecule ID if consolidated, otherwise the original ID
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get primary molecule
        query = """
        SELECT COALESCE(consolidated_to, id) as primary_id
        FROM molecules
        WHERE id = %s
        """
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        return result[0] if result is not None else molecule_id
        
    except Exception as e:
        logger.error(f"Error getting primary molecule: {str(e)}")
        return molecule_id

def get_consolidated_molecules(primary_id: str) -> List[str]:
    """
    Get the list of consolidated (secondary) molecules for a primary molecule.
    
    Args:
        primary_id: The primary molecule ID
        
    Returns:
        List of consolidated molecule IDs
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get consolidated molecules
        query = """
        SELECT id
        FROM molecules
        WHERE consolidated_to = %s
        """
        
        cursor.execute(query, (primary_id,))
        results = cursor.fetchall()
        
        return [row[0] for row in results] if results else []
        
    except Exception as e:
        logger.error(f"Error getting consolidated molecules: {str(e)}")
        return []

def get_differentiation_group(molecule_id: str) -> Optional[str]:
    """
    Get the differentiation group for a molecule.
    
    Args:
        molecule_id: The molecule ID to check
        
    Returns:
        Differentiation group ID if available, otherwise None
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get differentiation group
        query = """
        SELECT property_value
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s AND pt.name = 'differentiationGroup'
        """
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        return result[0] if result is not None else None
        
    except Exception as e:
        logger.error(f"Error getting differentiation group: {str(e)}")
        return None

def get_differentiation_group_members(group_id: str) -> List[str]:
    """
    Get all molecules in a differentiation group.
    
    Args:
        group_id: The differentiation group ID
        
    Returns:
        List of molecule IDs in the group
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get molecules in differentiation group
        query = """
        SELECT molecule_id
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE pt.name = 'differentiationGroup' AND mp.property_value = %s
        """
        
        cursor.execute(query, (group_id,))
        results = cursor.fetchall()
        
        return [row[0] for row in results] if results else []
        
    except Exception as e:
        logger.error(f"Error getting differentiation group members: {str(e)}")
        return []

def enrich_molecule_data(molecule_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enrich molecule data with consolidated and differentiation information.
    
    Args:
        molecule_data: Dictionary containing molecule data
        
    Returns:
        Enriched molecule data dictionary
    """
    if not molecule_data or 'id' not in molecule_data:
        return molecule_data
    
    molecule_id = molecule_data['id']
    
    # Add consolidated information
    molecule_data['is_consolidated'] = is_consolidated(molecule_id)
    molecule_data['consolidated_to'] = get_primary_molecule(molecule_id)
    
    # Only add consolidated_molecules list for primary molecules
    if not molecule_data['is_consolidated']:
        molecule_data['consolidated_molecules'] = get_consolidated_molecules(molecule_id)
    else:
        molecule_data['consolidated_molecules'] = []
    
    # Add differentiation information
    diff_group = get_differentiation_group(molecule_id)
    if diff_group:
        molecule_data['differentiation_group'] = diff_group
        molecule_data['differentiation_description'] = get_differentiation_description(molecule_id)
    
    return molecule_data

def get_differentiation_description(molecule_id: str) -> Optional[str]:
    """
    Get the differentiation description for a molecule.
    
    Args:
        molecule_id: The molecule ID to check
        
    Returns:
        Differentiation description if available, otherwise None
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get differentiation description
        query = """
        SELECT property_value
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s AND pt.name = 'differentiationDescription'
        """
        
        cursor.execute(query, (molecule_id,))
        result = cursor.fetchone()
        
        return result[0] if result is not None else None
        
    except Exception as e:
        logger.error(f"Error getting differentiation description: {str(e)}")
        return None

def get_all_differentiation_groups() -> List[Dict[str, Any]]:
    """
    Get all differentiation groups and their members.
    
    Returns:
        List of dictionaries containing group information
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get all differentiation groups
        query = """
        SELECT DISTINCT property_value as group_id
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE pt.name = 'differentiationGroup'
        """
        
        cursor.execute(query)
        groups = cursor.fetchall()
        
        result = []
        for row in groups:
            group_id = row[0]
            members = get_differentiation_group_members(group_id)
            
            # Get group name (use the first molecule's name as basis)
            if members:
                name_query = """
                SELECT name
                FROM molecules
                WHERE id = %s
                """
                cursor.execute(name_query, (members[0],))
                name_result = cursor.fetchone()
                group_name = f"Differentiated: {name_result[0]}" if name_result else f"Group {group_id}"
            else:
                group_name = f"Group {group_id}"
            
            result.append({
                'id': group_id,
                'name': group_name,
                'member_count': len(members),
                'members': members
            })
        
        return result
        
    except Exception as e:
        logger.error(f"Error getting all differentiation groups: {str(e)}")
        return []
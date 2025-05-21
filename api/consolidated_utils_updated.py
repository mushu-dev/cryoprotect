"""
CryoProtect Analyzer API - Consolidated Molecule Utilities

This module provides utilities for handling consolidated molecules in the API.
It helps maintain data consistency when working with molecules that have been
consolidated during the data quality enhancement process.
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
        SELECT molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL
        FROM consolidated_molecules
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
        SELECT 
            CASE 
                WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL 
                THEN primary_molecule_id
                ELSE id
            END as primary_id
        FROM consolidated_molecules
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
        FROM consolidated_molecules
        WHERE primary_molecule_id = %s AND molecule_status = 'duplicate'
        """
        
        cursor.execute(query, (primary_id,))
        results = cursor.fetchall()
        
        return [row[0] for row in results] if results else []
        
    except Exception as e:
        logger.error(f"Error getting consolidated molecules: {str(e)}")
        return []

def enrich_molecule_data(molecule_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enrich molecule data with consolidated information.
    
    Args:
        molecule_data: Dictionary containing molecule data
        
    Returns:
        Enriched molecule data dictionary
    """
    if not molecule_data or 'id' not in molecule_data:
        return molecule_data
    
    molecule_id = molecule_data['id']
    
    # Add consolidated information
    is_duplicate = is_consolidated(molecule_id)
    molecule_data['is_consolidated'] = is_duplicate
    
    if is_duplicate:
        # This is a duplicate molecule
        primary_id = get_primary_molecule(molecule_id)
        molecule_data['primary_molecule_id'] = primary_id
        
        # Get primary molecule name
        try:
            conn = get_db_connection()
            cursor = conn.cursor()
            
            # Query to get primary molecule name
            query = """
            SELECT name
            FROM consolidated_molecules
            WHERE id = %s
            """
            
            cursor.execute(query, (primary_id,))
            result = cursor.fetchone()
            
            if result:
                molecule_data['primary_molecule_name'] = result[0]
        except Exception as e:
            logger.error(f"Error getting primary molecule name: {str(e)}")
    
    elif molecule_data.get('molecule_status') == 'primary':
        # This is a primary molecule - get its duplicates
        duplicates = get_consolidated_molecules(molecule_id)
        molecule_data['duplicate_molecules'] = duplicates
        molecule_data['duplicate_count'] = len(duplicates)
        
        # Get duplicate names
        if duplicates:
            try:
                conn = get_db_connection()
                cursor = conn.cursor()
                
                placeholders = ', '.join(['%s'] * len(duplicates))
                query = f"""
                SELECT id, name
                FROM consolidated_molecules
                WHERE id IN ({placeholders})
                """
                
                cursor.execute(query, duplicates)
                results = cursor.fetchall()
                
                duplicate_names = {row[0]: row[1] for row in results}
                molecule_data['duplicate_names'] = duplicate_names
            except Exception as e:
                logger.error(f"Error getting duplicate molecule names: {str(e)}")
    
    return molecule_data

def get_molecule_audit_history(molecule_id: str) -> List[Dict[str, Any]]:
    """
    Get the audit history for a molecule.
    
    Args:
        molecule_id: The molecule ID to check
        
    Returns:
        List of audit events
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Query to get audit history
        query = """
        SELECT 
            id, operation, table_name, 
            old_value, new_value, timestamp
        FROM scientific_data_audit
        WHERE record_id = %s
        ORDER BY timestamp DESC
        """
        
        cursor.execute(query, (molecule_id,))
        results = cursor.fetchall()
        
        if not results:
            return []
        
        history = []
        for row in results:
            audit_id, operation, table_name, old_value, new_value, timestamp = row
            
            history.append({
                'id': audit_id,
                'operation': operation,
                'table_name': table_name,
                'old_value': old_value,
                'new_value': new_value,
                'timestamp': timestamp.isoformat() if timestamp else None
            })
        
        return history
        
    except Exception as e:
        logger.error(f"Error getting audit history: {str(e)}")
        return []
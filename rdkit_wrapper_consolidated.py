#!/usr/bin/env python3
"""
CryoProtect Analyzer - RDKit Wrapper for Consolidated Molecules

This module extends the RDKit wrapper to support consolidated molecules.
It ensures that property calculations are only performed on primary molecules,
and provides functions for migrating properties between consolidated molecules.
"""

import os
import sys
import logging
from typing import Dict, List, Optional, Union, Any, Tuple
import json

# Add the parent directory to the path to import the rdkit_wrapper
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from rdkit_wrapper import (
    calculate_properties,
    generate_visualization,
    search_substructure,
    search_similar_molecules,
    calculate_molecular_descriptors,
    check_rdkit_available
)

# Set up logging
logger = logging.getLogger(__name__)

def get_db_connection():
    """
    Get a database connection.
    
    Returns:
        Database connection from the current context
    """
    # Import here to avoid circular imports
    from database.adapter import get_connection
    return get_connection()

def get_primary_molecule_id(molecule_id: str, conn=None) -> str:
    """
    Get the primary molecule ID for a given molecule.
    
    If the molecule is a secondary (consolidated) molecule, return the ID
    of its primary molecule. Otherwise, return the original ID.
    
    Args:
        molecule_id: The molecule ID to check
        conn: Optional database connection
        
    Returns:
        Primary molecule ID if consolidated, otherwise the original ID
    """
    close_conn = False
    try:
        if conn is None:
            conn = get_db_connection()
            close_conn = True
        
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
    finally:
        if close_conn and conn:
            conn.close()

def calculate_properties_for_consolidated(molecule_id: str, properties: List[str] = None) -> Dict[str, Any]:
    """
    Calculate properties for a molecule, handling consolidated molecules.
    
    If the molecule is a duplicate (consolidated) molecule, get the primary molecule
    and calculate properties for it instead.
    
    Args:
        molecule_id: The ID of the molecule
        properties: Optional list of properties to calculate
        
    Returns:
        Dictionary of calculated properties
    """
    try:
        # Get the primary molecule ID
        primary_id = get_primary_molecule_id(molecule_id)
        
        # If the molecule is a duplicate, use the primary
        if primary_id != molecule_id:
            logger.info(f"Molecule {molecule_id} is a duplicate, using primary molecule {primary_id}")
            
            # Calculate properties for the primary molecule
            return calculate_properties(primary_id, properties)
        
        # Otherwise, calculate properties for the original molecule
        return calculate_properties(molecule_id, properties)
        
    except Exception as e:
        logger.error(f"Error calculating properties for consolidated molecule: {str(e)}")
        return {}

def migrate_properties(source_id: str, target_id: str, property_types: List[str] = None) -> Dict[str, Any]:
    """
    Migrate properties from one molecule to another.
    
    This function is used when consolidating molecules to ensure that properties
    are preserved and migrated to the primary molecule.
    
    Args:
        source_id: The source molecule ID
        target_id: The target molecule ID
        property_types: Optional list of property types to migrate
        
    Returns:
        Dictionary with migration results
    """
    conn = None
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # If no property types specified, get all property types
        if not property_types:
            query = """
            SELECT DISTINCT pt.name
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE mp.molecule_id = %s
            """
            
            cursor.execute(query, (source_id,))
            results = cursor.fetchall()
            property_types = [row[0] for row in results]
        
        migrated_properties = []
        skipped_properties = []
        
        # For each property type, check if it exists in the target molecule
        for prop_type in property_types:
            # Get property type ID
            cursor.execute("SELECT id FROM property_types WHERE name = %s", (prop_type,))
            type_result = cursor.fetchone()
            if not type_result:
                skipped_properties.append({
                    'property_type': prop_type,
                    'reason': 'Property type not found'
                })
                continue
            
            property_type_id = type_result[0]
            
            # Check if property exists in source molecule
            cursor.execute(
                "SELECT id, property_value, calculation_method FROM molecular_properties WHERE molecule_id = %s AND property_type_id = %s",
                (source_id, property_type_id)
            )
            source_result = cursor.fetchone()
            if not source_result:
                skipped_properties.append({
                    'property_type': prop_type,
                    'reason': 'Property not found in source molecule'
                })
                continue
            
            source_property_id, property_value, calculation_method = source_result
            
            # Check if property exists in target molecule
            cursor.execute(
                "SELECT id FROM molecular_properties WHERE molecule_id = %s AND property_type_id = %s",
                (target_id, property_type_id)
            )
            target_result = cursor.fetchone()
            
            if target_result:
                # Property exists in target, update it
                target_property_id = target_result[0]
                cursor.execute(
                    """
                    UPDATE molecular_properties 
                    SET property_value = %s, calculation_method = %s, updated_at = NOW()
                    WHERE id = %s
                    """
                    ,
                    (property_value, calculation_method, target_property_id)
                )
            else:
                # Property doesn't exist in target, insert it
                cursor.execute(
                    """
                    INSERT INTO molecular_properties (molecule_id, property_type_id, property_value, calculation_method, created_at, updated_at)
                    VALUES (%s, %s, %s, %s, NOW(), NOW())
                    """
                    ,
                    (target_id, property_type_id, property_value, calculation_method)
                )
            
            migrated_properties.append({
                'property_type': prop_type,
                'value': property_value
            })
        
        # Commit the transaction
        conn.commit()
        
        return {
            'source_id': source_id,
            'target_id': target_id,
            'migrated_properties': migrated_properties,
            'skipped_properties': skipped_properties,
            'property_count': len(migrated_properties)
        }
        
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error migrating properties: {str(e)}")
        return {
            'source_id': source_id,
            'target_id': target_id,
            'error': str(e),
            'migrated_properties': [],
            'skipped_properties': [],
            'property_count': 0
        }
    finally:
        if conn:
            conn.close()

def generate_visualization_for_consolidated(molecule_id: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Generate visualization for a molecule, handling consolidated molecules.
    
    If the molecule is a duplicate (consolidated) molecule, get the primary molecule
    and generate visualization for it instead.
    
    Args:
        molecule_id: The ID of the molecule
        options: Optional visualization options
        
    Returns:
        Dictionary with visualization data
    """
    try:
        # Get the primary molecule ID
        primary_id = get_primary_molecule_id(molecule_id)
        
        # Generate visualization
        result = generate_visualization(primary_id, options)
        
        # If the molecule is a duplicate, add information about the consolidated status
        if primary_id != molecule_id:
            result['is_consolidated'] = True
            result['primary_molecule_id'] = primary_id
            result['consolidated_note'] = f"Visualization generated for the primary molecule {primary_id}"
        
        return result
        
    except Exception as e:
        logger.error(f"Error generating visualization for consolidated molecule: {str(e)}")
        return {
            'molecule_id': molecule_id,
            'error': str(e)
        }

# Other functions that need to be consolidated-molecule aware
def search_substructure_with_consolidated(query: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Search for substructures, handling consolidated molecules.
    
    When returning search results, indicate whether molecules are consolidated
    and provide their primary molecule ID if applicable.
    
    Args:
        query: The substructure query
        options: Optional search options
        
    Returns:
        Dictionary with search results
    """
    try:
        # Perform the substructure search
        results = search_substructure(query, options)
        
        if 'matches' in results:
            # Process the matches to add consolidated molecule information
            conn = get_db_connection()
            cursor = conn.cursor()
            
            for match in results['matches']:
                molecule_id = match.get('molecule_id')
                if molecule_id:
                    # Check if the molecule is consolidated
                    query = """
                    SELECT 
                        molecule_status,
                        CASE 
                            WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL
                            THEN primary_molecule_id
                            ELSE NULL
                        END as primary_id
                    FROM consolidated_molecules
                    WHERE id = %s
                    """
                    
                    cursor.execute(query, (molecule_id,))
                    result = cursor.fetchone()
                    
                    if result:
                        status, primary_id = result
                        match['molecule_status'] = status
                        
                        if primary_id:
                            match['primary_molecule_id'] = primary_id
                            match['is_consolidated'] = True
                        else:
                            match['is_consolidated'] = False
            
            conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"Error in substructure search with consolidated molecules: {str(e)}")
        return {
            'query': query,
            'error': str(e),
            'matches': []
        }

def search_similar_molecules_with_consolidated(query: str, options: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Search for similar molecules, handling consolidated molecules.
    
    When returning search results, indicate whether molecules are consolidated
    and provide their primary molecule ID if applicable.
    
    Args:
        query: The similarity query
        options: Optional search options
        
    Returns:
        Dictionary with search results
    """
    try:
        # Perform the similarity search
        results = search_similar_molecules(query, options)
        
        if 'matches' in results:
            # Process the matches to add consolidated molecule information
            conn = get_db_connection()
            cursor = conn.cursor()
            
            for match in results['matches']:
                molecule_id = match.get('molecule_id')
                if molecule_id:
                    # Check if the molecule is consolidated
                    query = """
                    SELECT 
                        molecule_status,
                        CASE 
                            WHEN molecule_status = 'duplicate' AND primary_molecule_id IS NOT NULL
                            THEN primary_molecule_id
                            ELSE NULL
                        END as primary_id
                    FROM consolidated_molecules
                    WHERE id = %s
                    """
                    
                    cursor.execute(query, (molecule_id,))
                    result = cursor.fetchone()
                    
                    if result:
                        status, primary_id = result
                        match['molecule_status'] = status
                        
                        if primary_id:
                            match['primary_molecule_id'] = primary_id
                            match['is_consolidated'] = True
                        else:
                            match['is_consolidated'] = False
            
            conn.close()
        
        return results
        
    except Exception as e:
        logger.error(f"Error in similarity search with consolidated molecules: {str(e)}")
        return {
            'query': query,
            'error': str(e),
            'matches': []
        }
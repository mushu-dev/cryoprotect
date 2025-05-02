"""
Module for populating the molecules table.

This module provides functions for populating the molecules table
from various data sources.
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)

def _load_molecule_data(
    environment: str,
    file_path: Optional[str] = None
) -> List[Dict]:
    """
    Load molecule data from appropriate source.
    
    Args:
        environment: Target environment
        file_path: Optional path to data file
        
    Returns:
        List of molecule records
    """
    # If a specific file is provided, use it
    if file_path and os.path.exists(file_path):
        with open(file_path, 'r') as f:
            data = json.load(f)
            return data
            
    # Otherwise, use environment-specific data
    base_dir = os.path.abspath(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.dirname(__file__)))))
        
    if environment == 'production':
        data_path = os.path.join(base_dir, 'data', 'production', 'molecules.json')
    elif environment == 'staging':
        data_path = os.path.join(base_dir, 'data', 'staging', 'molecules.json')
    else:  # development
        data_path = os.path.join(base_dir, 'data', 'development', 'molecules.json')
        
    if not os.path.exists(data_path):
        logger.warning(f"Data file not found: {data_path}")
        return []
        
    with open(data_path, 'r') as f:
        data = json.load(f)
        return data

def populate_molecules(
    conn: Any,
    environment: str = 'development',
    file_path: Optional[str] = None
) -> int:
    """
    Populate the molecules table.
    
    Args:
        conn: Database connection
        environment: Target environment
        file_path: Optional path to data file
        
    Returns:
        Number of records inserted
    """
    logger.info(f"Populating molecules table in {environment} environment")
    
    # Load data
    molecule_data = _load_molecule_data(environment, file_path)
    if not molecule_data:
        logger.warning("No molecule data found")
        return 0
        
    # Insert data
    count = 0
    for molecule in molecule_data:
        try:
            # Example using Supabase client
            result = conn.table('molecules').insert(molecule).execute()
            count += 1
        except Exception as e:
            logger.error(f"Error inserting molecule {molecule.get('name', 'unknown')}: {str(e)}")
            
    logger.info(f"Inserted {count} molecules")
    return count
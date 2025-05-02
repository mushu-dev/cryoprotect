"""
Module for populating the mixtures table.

This module provides functions for populating the mixtures table
from various data sources.
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)

def _load_mixture_data(
    environment: str,
    file_path: Optional[str] = None
) -> List[Dict]:
    """
    Load mixture data from appropriate source.
    
    Args:
        environment: Target environment
        file_path: Optional path to data file
        
    Returns:
        List of mixture records
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
        data_path = os.path.join(base_dir, 'data', 'production', 'mixtures.json')
    elif environment == 'staging':
        data_path = os.path.join(base_dir, 'data', 'staging', 'mixtures.json')
    else:  # development
        data_path = os.path.join(base_dir, 'data', 'development', 'mixtures.json')
        
    if not os.path.exists(data_path):
        logger.warning(f"Data file not found: {data_path}")
        return []
        
    with open(data_path, 'r') as f:
        data = json.load(f)
        return data

def populate_mixtures(
    conn: Any,
    environment: str = 'development',
    file_path: Optional[str] = None
) -> int:
    """
    Populate the mixtures table.
    
    Args:
        conn: Database connection
        environment: Target environment
        file_path: Optional path to data file
        
    Returns:
        Number of records inserted
    """
    logger.info(f"Populating mixtures table in {environment} environment")
    
    # Load data
    mixture_data = _load_mixture_data(environment, file_path)
    if not mixture_data:
        logger.warning("No mixture data found")
        return 0
        
    # Insert data
    count = 0
    for mixture in mixture_data:
        try:
            # Example using Supabase client
            result = conn.table('mixtures').insert(mixture).execute()
            
            # Handle mixture components if present
            if 'components' in mixture and mixture['components']:
                for component in mixture['components']:
                    component['mixture_id'] = result.data[0]['id']
                    conn.table('mixture_components').insert(component).execute()
                    
            count += 1
        except Exception as e:
            logger.error(f"Error inserting mixture {mixture.get('name', 'unknown')}: {str(e)}")
            
    logger.info(f"Inserted {count} mixtures")
    return count
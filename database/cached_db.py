#!/usr/bin/env python3
"""
Cached database module for CryoProtect.

This module provides a wrapped version of the database module with caching
for improved performance.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple, Union

# Import the base database module
from database.db import (
    get_connection,
    release_connection,
    init_connection_pool,
    close_all_connections,
    transaction
)

# Import the caching layer
from database.cache import (
    cache_query,
    cache_molecule,
    cache_property,
    invalidate_molecule,
    invalidate_property,
    clear_cache,
    get_cache_stats
)

# Import cache invalidation
from database.cache_invalidation import (
    setup_cache_invalidation,
    process_invalidation_events,
    invalidation_processor
)

# Configure logging
logger = logging.getLogger(__name__)

# Default TTL values (in seconds)
CACHE_TTL_SHORT = 60      # 1 minute
CACHE_TTL_MEDIUM = 600    # 10 minutes
CACHE_TTL_LONG = 3600     # 1 hour
CACHE_TTL_VERY_LONG = 86400  # 24 hours

# Cache configuration
_cache_initialized = False

def initialize_cache() -> bool:
    """
    Initialize the cache and invalidation system.
    
    Returns:
        True if successful, False otherwise
    """
    global _cache_initialized
    
    if _cache_initialized:
        return True
    
    try:
        # Setup cache invalidation infrastructure
        result = setup_cache_invalidation()
        
        if result:
            logger.info("Cache system initialized successfully")
            _cache_initialized = True
            return True
        else:
            logger.warning("Failed to initialize cache system")
            return False
    except Exception as e:
        logger.error(f"Error initializing cache: {e}")
        return False

# Process any pending invalidation events (asynchronously)
def check_invalidation_events() -> None:
    """
    Check for pending invalidation events and process them.
    """
    try:
        # Process any pending events
        process_invalidation_events()
    except Exception as e:
        logger.warning(f"Error processing invalidation events: {e}")

# Cached version of execute_query with different TTLs based on query type
@cache_query(ttl=CACHE_TTL_MEDIUM)
def execute_query(query: str, params: Any = None) -> List[Dict[str, Any]]:
    """
    Execute a SQL query and return the results (with caching).
    
    Args:
        query: SQL query string
        params: Query parameters
        
    Returns:
        Query results or None on error
    """
    from database.db import execute_query as base_execute_query
    
    # Initialize cache if needed
    if not _cache_initialized:
        initialize_cache()
    
    # Check if this is a read-only query that we can cache
    is_read_only = (
        query.strip().upper().startswith('SELECT') and
        'FOR UPDATE' not in query.upper() and
        'SET LOCAL' not in query.upper()
    )
    
    # If not read-only, bypass cache
    if not is_read_only:
        # Execute without caching
        return base_execute_query(query, params)
    
    # Periodically check for invalidation events (1% chance)
    import random
    if random.random() < 0.01:
        check_invalidation_events()
    
    # Execute with caching (the cache decorator handles the rest)
    return base_execute_query(query, params)

# Cached version of get_molecule_by_id
@cache_molecule(ttl=CACHE_TTL_LONG)
def get_molecule_by_id(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule by ID (with caching).
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        Molecule data or None if not found
    """
    query = """
        SELECT * FROM molecules 
        WHERE id = %s
    """
    result = execute_query(query, (molecule_id,))
    return result[0] if result else None

# Cached version of get_molecule_by_name
@cache_molecule(ttl=CACHE_TTL_LONG)
def get_molecule_by_name(name: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule by name (with caching).
    
    Args:
        name: Molecule name
        
    Returns:
        Molecule data or None if not found
    """
    query = """
        SELECT * FROM molecules 
        WHERE name = %s
    """
    result = execute_query(query, (name,))
    return result[0] if result else None

# Cached version of get_molecule_by_smiles
@cache_molecule(ttl=CACHE_TTL_LONG)
def get_molecule_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    """
    Get a molecule by SMILES (with caching).
    
    Args:
        smiles: SMILES notation
        
    Returns:
        Molecule data or None if not found
    """
    query = """
        SELECT * FROM molecules 
        WHERE smiles = %s
    """
    result = execute_query(query, (smiles,))
    return result[0] if result else None

# Cached version of get_molecular_properties
@cache_property(ttl=CACHE_TTL_MEDIUM)
def get_molecular_properties(molecule_id: str) -> List[Dict[str, Any]]:
    """
    Get molecular properties for a molecule (with caching).
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        List of properties
    """
    query = """
        SELECT mp.*, pt.name as property_name, pt.data_type 
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s
    """
    return execute_query(query, (molecule_id,))

# Cached version of get_properties_by_type
@cache_property(ttl=CACHE_TTL_MEDIUM)
def get_properties_by_type(property_type_id: str) -> List[Dict[str, Any]]:
    """
    Get properties by type (with caching).
    
    Args:
        property_type_id: Property type ID
        
    Returns:
        List of properties
    """
    query = """
        SELECT mp.*, m.name as molecule_name 
        FROM molecular_properties mp
        JOIN molecules m ON mp.molecule_id = m.id
        WHERE mp.property_type_id = %s
    """
    return execute_query(query, (property_type_id,))

# Cached version of get_property_types
@cache_query(ttl=CACHE_TTL_VERY_LONG)
def get_property_types() -> List[Dict[str, Any]]:
    """
    Get all property types (with caching).
    
    Returns:
        List of property types
    """
    query = """
        SELECT * FROM property_types
        ORDER BY name
    """
    return execute_query(query)

# Cached version of get_tables
@cache_query(ttl=CACHE_TTL_VERY_LONG)
def get_tables() -> List[str]:
    """
    Get list of tables in the database (with caching).
    
    Returns:
        List of table names or None on error
    """
    from database.db import get_tables as base_get_tables
    return base_get_tables()

# Cached version of get_table_info
@cache_query(ttl=CACHE_TTL_VERY_LONG)
def get_table_info(table_name: str) -> List[Dict[str, Any]]:
    """
    Get column information for a table (with caching).
    
    Args:
        table_name: Name of the table
        
    Returns:
        List of column details or None on error
    """
    from database.db import get_table_info as base_get_table_info
    return base_get_table_info(table_name)

# Batch operations that update the database (using the invalidation system)
def batch_update_molecules(molecules: List[Dict[str, Any]]) -> List[str]:
    """
    Batch update molecules.
    
    Args:
        molecules: List of molecule dictionaries with 'id' field
        
    Returns:
        List of updated molecule IDs
    """
    from database.db import execute_batch
    
    # Prepare the batch update
    updates = []
    for molecule in molecules:
        if 'id' not in molecule:
            continue
            
        # Prepare query and parameters
        set_clauses = []
        params = []
        for key, value in molecule.items():
            if key != 'id':
                set_clauses.append(f"{key} = %s")
                params.append(value)
                
        if not set_clauses:
            continue
            
        # Add molecule ID as the last parameter
        params.append(molecule['id'])
        
        # Create the complete query
        query = f"""
            UPDATE molecules 
            SET {', '.join(set_clauses)}, 
                updated_at = NOW() 
            WHERE id = %s
            RETURNING id
        """
        
        updates.append((query, tuple(params)))
    
    # Execute the batch update
    if not updates:
        return []
        
    results = execute_batch(updates)
    updated_ids = [result[0]['id'] for result in results if result]
        
    return updated_ids

def batch_update_properties(properties: List[Dict[str, Any]]) -> List[str]:
    """
    Batch update properties.
    
    Args:
        properties: List of property dictionaries with 'id' field
        
    Returns:
        List of updated property IDs
    """
    from database.db import execute_batch
    
    # Prepare the batch update
    updates = []
    
    for prop in properties:
        if 'id' not in prop:
            continue
            
        # Prepare query and parameters
        set_clauses = []
        params = []
        for key, value in prop.items():
            if key != 'id':
                set_clauses.append(f"{key} = %s")
                params.append(value)
                
        if not set_clauses:
            continue
            
        # Add property ID as the last parameter
        params.append(prop['id'])
        
        # Create the complete query
        query = f"""
            UPDATE molecular_properties 
            SET {', '.join(set_clauses)}, 
                updated_at = NOW() 
            WHERE id = %s
            RETURNING id
        """
        
        updates.append((query, tuple(params)))
    
    # Execute the batch update
    if not updates:
        return []
        
    results = execute_batch(updates)
    updated_ids = [result[0]['id'] for result in results if result]
    
    return updated_ids

# Cache statistics and management
def get_database_cache_stats() -> Dict[str, Any]:
    """
    Get statistics about the database cache.
    
    Returns:
        Dictionary with cache statistics
    """
    # Process any pending invalidation events
    check_invalidation_events()
    
    return get_cache_stats()

def clear_database_cache() -> bool:
    """
    Clear the entire database cache.
    
    Returns:
        True if successful, False otherwise
    """
    return clear_cache()

def process_cache_invalidation_events() -> int:
    """
    Process pending cache invalidation events.
    
    Returns:
        Number of events processed
    """
    return process_invalidation_events()

# Helper method to decide whether to use cache based on table name and operation
def should_cache_operation(table_name: str, operation: str) -> bool:
    """
    Determine if an operation on a table should be cached.
    
    Args:
        table_name: The table name
        operation: The operation (SELECT, INSERT, UPDATE, DELETE)
        
    Returns:
        True if the operation should be cached, False otherwise
    """
    # Only cache SELECT operations
    if operation != 'SELECT':
        return False
        
    # Cache most frequently accessed tables
    frequently_accessed = [
        'molecules',
        'molecular_properties',
        'property_types',
        'cryoprotection_scores',
        'mixtures'
    ]
    
    return table_name in frequently_accessed

# Initialize the cache on module import
initialize_cache()
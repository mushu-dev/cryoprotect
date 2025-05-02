#!/usr/bin/env python3
"""
Database utility functions for CryoProtect v2.
"""

import os
import logging
from typing import Any, Dict, List, Optional, Union, Tuple, Callable
import time
import functools

from .connection_manager import ConnectionManager

logger = logging.getLogger(__name__)

def get_db():
    """
    Get database connection manager.
    
    Returns:
        ConnectionManager: Database connection manager
    """
    return ConnectionManager.get_instance()

def with_connection(f: Callable) -> Callable:
    """
    Decorator to ensure database connection is established before function execution.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        db = get_db()
        
        # Ensure connection
        if not db.get_active_adapter():
            if not db.connect():
                raise ConnectionError("Failed to connect to database")
                
        return f(*args, **kwargs)
    return wrapper

def with_retry(max_retries: int = 3, backoff: float = 1.5) -> Callable:
    """
    Decorator to retry database operations with exponential backoff.
    
    Args:
        max_retries: Maximum number of retry attempts
        backoff: Backoff factor for retry delays
        
    Returns:
        Decorator function
    """
    def decorator(f: Callable) -> Callable:
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            last_exception = None
            retry_count = 0
            
            while retry_count <= max_retries:
                try:
                    return f(*args, **kwargs)
                except Exception as e:
                    retry_count += 1
                    last_exception = e
                    
                    # Check if we should retry
                    if retry_count > max_retries:
                        break
                        
                    # Calculate backoff time
                    delay = backoff ** (retry_count - 1)
                    
                    logger.warning(
                        f"Retrying database operation after error: {str(e)}. "
                        f"Retry {retry_count}/{max_retries} in {delay:.2f}s..."
                    )
                    
                    time.sleep(delay)
                    
            # Re-raise the last exception
            raise last_exception
        return wrapper
    return decorator

def with_transaction(f: Callable) -> Callable:
    """
    Decorator to execute function in a database transaction.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        db = get_db()
        
        # Ensure connection
        if not db.get_active_adapter():
            if not db.connect():
                raise ConnectionError("Failed to connect to database")
                
        # Begin transaction
        transaction = db.begin_transaction()
        
        try:
            # Execute function
            result = f(*args, transaction=transaction, **kwargs)
            
            # Commit transaction
            db.commit_transaction(transaction)
            
            return result
        except Exception as e:
            # Rollback transaction
            db.rollback_transaction(transaction)
            raise
    return wrapper

@with_connection
def execute_query(query: str, params: Optional[Union[Tuple, Dict]] = None) -> Any:
    """
    Execute SQL query.
    
    Args:
        query: SQL query to execute
        params: Query parameters
        
    Returns:
        Query results
    """
    db = get_db()
    return db.execute_query(query, params)

@with_connection
def execute_batch(queries: List[str]) -> List[Any]:
    """
    Execute multiple SQL queries.
    
    Args:
        queries: List of SQL queries to execute
        
    Returns:
        List of query results
    """
    db = get_db()
    return db.execute_batch(queries)

@with_connection
def get_molecule_by_id(molecule_id: str) -> Optional[Dict[str, Any]]:
    """
    Get molecule by ID.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        Molecule data or None if not found
    """
    db = get_db()
    results = db.execute_query(
        "SELECT * FROM molecules WHERE id = %s",
        (molecule_id,)
    )
    
    if results and len(results) > 0:
        return results[0]
    return None

@with_connection
def get_molecule_properties(molecule_id: str) -> List[Dict[str, Any]]:
    """
    Get properties for a molecule.
    
    Args:
        molecule_id: Molecule ID
        
    Returns:
        List of property data
    """
    db = get_db()
    return db.execute_query(
        """
        SELECT mp.*, pt.name as property_name, pt.unit 
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE mp.molecule_id = %s
        """,
        (molecule_id,)
    )

@with_connection
def get_molecules_by_inchikey(inchi_key: str) -> List[Dict[str, Any]]:
    """
    Get molecules by InChI Key.
    
    Args:
        inchi_key: InChI Key
        
    Returns:
        List of molecule data
    """
    db = get_db()
    return db.execute_query(
        "SELECT * FROM molecules WHERE inchi_key = %s",
        (inchi_key,)
    )

@with_connection
def insert_molecule(
    name: str,
    formula: Optional[str] = None,
    molecular_weight: Optional[float] = None,
    smiles: Optional[str] = None,
    inchi: Optional[str] = None,
    inchi_key: Optional[str] = None,
    chembl_id: Optional[str] = None,
    pubchem_cid: Optional[str] = None,
    data_source: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Insert new molecule.
    
    Args:
        name: Molecule name
        formula: Molecular formula
        molecular_weight: Molecular weight
        smiles: SMILES string
        inchi: InChI string
        inchi_key: InChI Key
        chembl_id: ChEMBL ID
        pubchem_cid: PubChem CID
        data_source: Data source
        
    Returns:
        Inserted molecule data
    """
    db = get_db()
    result = db.execute_query(
        """
        INSERT INTO molecules (
            name, formula, molecular_weight, smiles, inchi, inchi_key,
            chembl_id, pubchem_cid, data_source, created_at, updated_at
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW()
        )
        RETURNING *
        """,
        (
            name, formula, molecular_weight, smiles, inchi, inchi_key,
            chembl_id, pubchem_cid, data_source
        )
    )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def update_molecule(
    molecule_id: str,
    name: Optional[str] = None,
    formula: Optional[str] = None,
    molecular_weight: Optional[float] = None,
    smiles: Optional[str] = None,
    inchi: Optional[str] = None,
    inchi_key: Optional[str] = None,
    chembl_id: Optional[str] = None,
    pubchem_cid: Optional[str] = None,
    data_source: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Update molecule.
    
    Args:
        molecule_id: Molecule ID
        name: Molecule name
        formula: Molecular formula
        molecular_weight: Molecular weight
        smiles: SMILES string
        inchi: InChI string
        inchi_key: InChI Key
        chembl_id: ChEMBL ID
        pubchem_cid: PubChem CID
        data_source: Data source
        
    Returns:
        Updated molecule data
    """
    # Build update fields
    update_fields = []
    params = []
    
    if name is not None:
        update_fields.append("name = %s")
        params.append(name)
        
    if formula is not None:
        update_fields.append("formula = %s")
        params.append(formula)
        
    if molecular_weight is not None:
        update_fields.append("molecular_weight = %s")
        params.append(molecular_weight)
        
    if smiles is not None:
        update_fields.append("smiles = %s")
        params.append(smiles)
        
    if inchi is not None:
        update_fields.append("inchi = %s")
        params.append(inchi)
        
    if inchi_key is not None:
        update_fields.append("inchi_key = %s")
        params.append(inchi_key)
        
    if chembl_id is not None:
        update_fields.append("chembl_id = %s")
        params.append(chembl_id)
        
    if pubchem_cid is not None:
        update_fields.append("pubchem_cid = %s")
        params.append(pubchem_cid)
        
    if data_source is not None:
        update_fields.append("data_source = %s")
        params.append(data_source)
        
    # Add updated_at field
    update_fields.append("updated_at = NOW()")
    
    # If no fields to update, return current molecule
    if not update_fields:
        return get_molecule_by_id(molecule_id)
        
    # Build query
    query = f"""
        UPDATE molecules
        SET {", ".join(update_fields)}
        WHERE id = %s
        RETURNING *
    """
    
    # Add molecule_id to params
    params.append(molecule_id)
    
    # Execute query
    db = get_db()
    result = db.execute_query(query, tuple(params))
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def set_molecule_property(
    molecule_id: str,
    property_type_id: str,
    value: Union[float, str, bool],
    source: Optional[str] = None,
    confidence: Optional[float] = None
) -> Optional[Dict[str, Any]]:
    """
    Set molecule property.
    
    Args:
        molecule_id: Molecule ID
        property_type_id: Property type ID
        value: Property value
        source: Data source
        confidence: Confidence value
        
    Returns:
        Property data
    """
    db = get_db()
    
    # Check if property exists
    existing = db.execute_query(
        """
        SELECT * FROM molecular_properties
        WHERE molecule_id = %s AND property_type_id = %s
        """,
        (molecule_id, property_type_id)
    )
    
    if existing and len(existing) > 0:
        # Update existing property
        result = db.execute_query(
            """
            UPDATE molecular_properties
            SET value = %s, source = %s, confidence = %s, updated_at = NOW()
            WHERE molecule_id = %s AND property_type_id = %s
            RETURNING *
            """,
            (value, source, confidence, molecule_id, property_type_id)
        )
    else:
        # Insert new property
        result = db.execute_query(
            """
            INSERT INTO molecular_properties
            (molecule_id, property_type_id, value, source, confidence, created_at, updated_at)
            VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
            RETURNING *
            """,
            (molecule_id, property_type_id, value, source, confidence)
        )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def get_or_create_property_type(
    name: str,
    description: Optional[str] = None,
    data_type: str = 'numeric',
    unit: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Get or create property type.
    
    Args:
        name: Property type name
        description: Property type description
        data_type: Data type (numeric, text, boolean)
        unit: Unit of measurement
        
    Returns:
        Property type data
    """
    db = get_db()
    
    # Check if property type exists
    existing = db.execute_query(
        "SELECT * FROM property_types WHERE name = %s",
        (name,)
    )
    
    if existing and len(existing) > 0:
        return existing[0]
        
    # Create new property type
    result = db.execute_query(
        """
        INSERT INTO property_types
        (name, description, data_type, unit, created_at, updated_at)
        VALUES (%s, %s, %s, %s, NOW(), NOW())
        RETURNING *
        """,
        (name, description, data_type, unit)
    )
    
    if result and len(result) > 0:
        return result[0]
    return None

@with_connection
def test_database_connection() -> Dict[str, Any]:
    """
    Test database connection and return detailed status.
    
    Returns:
        Dict with connection status and test results
    """
    db = get_db()
    adapter_info = db.get_connection_info()
    test_results = db.test_all_connections()
    
    # Try simple query
    query_result = None
    query_error = None
    
    try:
        result = db.execute_query("SELECT 1 as test")
        if result and len(result) > 0:
            query_result = "Success"
    except Exception as e:
        query_error = str(e)
    
    return {
        'connection_info': adapter_info,
        'test_results': test_results,
        'query_test': {
            'result': query_result,
            'error': query_error
        }
    }
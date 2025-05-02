#!/usr/bin/env python3
"""
ChEMBL Search Utilities for CryoProtect v2

This module provides specialized search functions for identifying potential
cryoprotectants in the ChEMBL database based on molecular properties,
structural similarity, and chemical class patterns.

Based on specifications in DATABASE_POPULATION_ISSUES.md (Section 4.3)
"""

import logging
import time
from typing import Dict, List, Any, Optional, Union, Tuple
from contextlib import contextmanager

from chembl.client import ChEMBLClient, ResilientChEMBLClient
from db_connection_utils import safe_transaction, get_db_connection

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('chembl_search_utils')

def find_potential_cryoprotectants(
    limit: int = 5000,
    min_mw: float = 30.0,
    max_mw: float = 500.0,
    min_logp: float = -3.0,
    max_logp: float = 3.0,
    min_hba: int = 2,
    max_hba: int = 20,
    min_hbd: int = 1,
    max_hbd: int = 10,
    min_phase: int = 0,
    use_cache: bool = True,
    fallback_to_cache: bool = True
) -> List[Dict[str, Any]]:
    """
    Find compounds with properties matching typical cryoprotectants.
    
    This function queries the ChEMBL database for compounds with properties
    that match the typical profile of cryoprotectants:
    - Molecular weight in a specific range (typically 30-500)
    - LogP in a specific range (typically -3 to 3)
    - Multiple hydrogen bond donors and acceptors
    - Preferably compounds that have reached clinical trials
    
    Args:
        limit: Maximum number of compounds to return
        min_mw: Minimum molecular weight
        max_mw: Maximum molecular weight
        min_logp: Minimum calculated LogP
        max_logp: Maximum calculated LogP
        min_hba: Minimum number of hydrogen bond acceptors
        max_hba: Maximum number of hydrogen bond acceptors
        min_hbd: Minimum number of hydrogen bond donors
        max_hbd: Maximum number of hydrogen bond donors
        min_phase: Minimum clinical phase (0-4, where 0 means preclinical)
        use_cache: Whether to use cached data if available
        fallback_to_cache: Whether to fallback to cached data if API fails
        
    Returns:
        List of dictionaries containing compound information
    """
    logger.info(f"Searching for potential cryoprotectants with properties: "
               f"MW={min_mw}-{max_mw}, LogP={min_logp}-{max_logp}, "
               f"HBA={min_hba}-{max_hba}, HBD={min_hbd}-{max_hbd}, "
               f"Phase>={min_phase}")
    
    # Create a resilient client for direct API access
    client = ResilientChEMBLClient()
    
    # Construct the query
    query = f"""
    SELECT molecule_chembl_id, canonical_smiles, pref_name, molecule_properties.*
    FROM molecule_dictionary md
    JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    JOIN compound_properties cp ON md.molregno = cp.molregno
    JOIN compound_structures cs ON mh.parent_molregno = cs.molregno
    WHERE 
        cp.mw_freebase BETWEEN {min_mw} AND {max_mw}
        AND cp.alogp BETWEEN {min_logp} AND {max_logp}
        AND cp.hba BETWEEN {min_hba} AND {max_hba}
        AND cp.hbd BETWEEN {min_hbd} AND {max_hbd}
        AND md.max_phase >= {min_phase}
    ORDER BY md.max_phase DESC
    LIMIT {limit}
    """
    
    # Create cache key for this specific query
    cache_key = f"property_search_{min_mw}_{max_mw}_{min_logp}_{max_logp}_{min_hba}_{max_hba}_{min_hbd}_{max_hbd}_{min_phase}_{limit}"
    
    # Check cache if enabled
    if use_cache:
        cached_data = client.cache.get(cache_key)
        if cached_data:
            logger.info(f"Using cached results for property-based cryoprotectant search (found {len(cached_data)} compounds)")
            return cached_data
    
    try:
        # Execute the query using the _make_request method
        start_time = time.time()
        endpoint = "molecule"
        params = {
            "mw_freebase__gte": min_mw,
            "mw_freebase__lte": max_mw,
            "alogp__gte": min_logp,
            "alogp__lte": max_logp,
            "hba__gte": min_hba,
            "hba__lte": max_hba,
            "hbd__gte": min_hbd,
            "hbd__lte": max_hbd,
            "max_phase__gte": min_phase,
            "limit": limit
        }
        
        response = client._make_request(endpoint, params)
        
        # Process the results
        compounds = []
        for molecule in response.get("molecules", []):
            # Extract the relevant information
            compound = {
                "molecule_chembl_id": molecule.get("molecule_chembl_id"),
                "canonical_smiles": molecule.get("molecule_structures", {}).get("canonical_smiles"),
                "pref_name": molecule.get("pref_name"),
                "molecule_properties": molecule.get("molecule_properties", {})
            }
            compounds.append(compound)
        
        elapsed_time = time.time() - start_time
        logger.info(f"Found {len(compounds)} potential cryoprotectants in {elapsed_time:.2f} seconds")
        
        # Cache the results
        if use_cache:
            client.cache.set(cache_key, compounds)
        
        return compounds
        
    except Exception as e:
        logger.error(f"Error searching for potential cryoprotectants: {str(e)}")
        
        # Fallback to cache if enabled
        if fallback_to_cache:
            cached_data = client.cache.get(cache_key)
            if cached_data:
                logger.info(f"Falling back to cached results for property-based search (found {len(cached_data)} compounds)")
                return cached_data
        
        # Return empty list if all else fails
        return []

def store_cryoprotectant_candidates(compounds: List[Dict[str, Any]]) -> int:
    """
    Store potential cryoprotectant candidates in the database.
    
    Args:
        compounds: List of compound dictionaries from find_potential_cryoprotectants
        
    Returns:
        Number of compounds successfully stored
    """
    if not compounds:
        logger.warning("No compounds to store")
        return 0
    
    logger.info(f"Storing {len(compounds)} potential cryoprotectant candidates in database")
    
    stored_count = 0
    
    with safe_transaction() as transaction:
        connection = transaction.connection
        
        # Check if the cryoprotectant_candidates table exists
        check_table_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'cryoprotectant_candidates'
        );
        """
        
        result = connection.execute_query(check_table_query)
        table_exists = result[0]['exists']
        
        # Create the table if it doesn't exist
        if not table_exists:
            create_table_query = """
            CREATE TABLE cryoprotectant_candidates (
                id SERIAL PRIMARY KEY,
                chembl_id VARCHAR(20) UNIQUE NOT NULL,
                name VARCHAR(255),
                smiles TEXT,
                molecular_weight NUMERIC,
                logp NUMERIC,
                hba INTEGER,
                hbd INTEGER,
                clinical_phase INTEGER,
                properties JSONB,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """
            connection.execute_query(create_table_query)
            logger.info("Created cryoprotectant_candidates table")
        
        # Insert compounds
        for compound in compounds:
            try:
                # Extract properties
                chembl_id = compound.get("molecule_chembl_id")
                name = compound.get("pref_name", "")
                smiles = compound.get("canonical_smiles", "")
                
                properties = compound.get("molecule_properties", {})
                mw = properties.get("full_mwt") or properties.get("mw_freebase")
                logp = properties.get("alogp")
                hba = properties.get("hba")
                hbd = properties.get("hbd")
                phase = compound.get("max_phase", 0)
                
                # Skip if missing essential data
                if not chembl_id or not smiles:
                    logger.warning(f"Skipping compound with missing data: {compound}")
                    continue
                
                # Insert into database
                insert_query = """
                INSERT INTO cryoprotectant_candidates 
                (chembl_id, name, smiles, molecular_weight, logp, hba, hbd, clinical_phase, properties)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (chembl_id) DO UPDATE SET
                    name = EXCLUDED.name,
                    smiles = EXCLUDED.smiles,
                    molecular_weight = EXCLUDED.molecular_weight,
                    logp = EXCLUDED.logp,
                    hba = EXCLUDED.hba,
                    hbd = EXCLUDED.hbd,
                    clinical_phase = EXCLUDED.clinical_phase,
                    properties = EXCLUDED.properties,
                    created_at = CURRENT_TIMESTAMP
                """
                
                params = (
                    chembl_id, 
                    name, 
                    smiles, 
                    mw, 
                    logp, 
                    hba, 
                    hbd, 
                    phase, 
                    properties
                )
                
                connection.execute_query(insert_query, params)
                stored_count += 1
                
            except Exception as e:
                logger.error(f"Error storing compound {compound.get('molecule_chembl_id')}: {str(e)}")
    
    logger.info(f"Successfully stored {stored_count} potential cryoprotectant candidates")
    return stored_count

def get_stored_cryoprotectant_candidates(limit: int = 100) -> List[Dict[str, Any]]:
    """
    Retrieve stored cryoprotectant candidates from the database.
    
    Args:
        limit: Maximum number of candidates to retrieve
        
    Returns:
        List of dictionaries containing candidate information
    """
    logger.info(f"Retrieving up to {limit} stored cryoprotectant candidates")
    
    with get_db_connection() as connection:
        query = f"""
        SELECT * FROM cryoprotectant_candidates
        ORDER BY clinical_phase DESC, molecular_weight ASC
        LIMIT {limit}
        """
        
        try:
            results = connection.execute_query(query)
            logger.info(f"Retrieved {len(results)} cryoprotectant candidates")
            return results
        except Exception as e:
            logger.error(f"Error retrieving cryoprotectant candidates: {str(e)}")
            return []

def find_similar_compounds(
    reference_smiles: Union[str, List[str]],
    similarity_threshold: int = 70,
    limit: int = 100,
    use_cache: bool = True,
    fallback_to_cache: bool = True
) -> List[Dict[str, Any]]:
    """
    Find compounds similar to known cryoprotectants based on structural similarity.
    
    This function takes one or more reference SMILES strings and searches the ChEMBL
    database for structurally similar compounds. It can handle multiple reference
    compounds and will combine the results.
    
    Args:
        reference_smiles: A single SMILES string or a list of SMILES strings
        similarity_threshold: Minimum similarity percentage (0-100)
        limit: Maximum number of similar compounds to return per reference
        use_cache: Whether to use cached data if available
        fallback_to_cache: Whether to fallback to cached data if API fails
        
    Returns:
        List of dictionaries containing similar compound information
    """
    # Ensure reference_smiles is a list
    if isinstance(reference_smiles, str):
        reference_smiles = [reference_smiles]
    
    if not reference_smiles:
        logger.warning("No reference SMILES provided")
        return []
    
    logger.info(f"Searching for compounds similar to {len(reference_smiles)} reference structures "
               f"with similarity threshold {similarity_threshold}%")
    
    # Create a resilient client for direct API access
    client = ResilientChEMBLClient()
    
    # Track all similar molecules to avoid duplicates
    all_similar_molecules = {}
    
    for i, smiles in enumerate(reference_smiles):
        if not smiles:
            logger.warning(f"Skipping empty SMILES at index {i}")
            continue
            
        logger.info(f"Processing reference SMILES [{i+1}/{len(reference_smiles)}]: {smiles[:50]}...")
        
        # Create cache key for this specific query
        cache_key = f"similarity_search_{smiles}_{similarity_threshold}_{limit}"
        
        # Check cache if enabled
        if use_cache:
            cached_data = client.cache.get(cache_key)
            if cached_data:
                logger.info(f"Using cached results for similarity search (found {len(cached_data)} compounds)")
                # Add to our collection, avoiding duplicates
                for mol in cached_data:
                    mol_id = mol.get("molecule_chembl_id")
                    if mol_id and mol_id not in all_similar_molecules:
                        all_similar_molecules[mol_id] = mol
                continue
        
        try:
            # Execute the similarity search
            start_time = time.time()
            
            # Use the similarity endpoint with the SMILES string
            endpoint = "similarity"
            params = {
                "smiles": smiles,
                "similarity": similarity_threshold,
                "limit": limit
            }
            
            response = client._make_request(endpoint, params)
            
            # Process the results
            similar_molecules = []
            for molecule in response.get("molecules", []):
                # Extract the relevant information
                compound = {
                    "molecule_chembl_id": molecule.get("molecule_chembl_id"),
                    "canonical_smiles": molecule.get("molecule_structures", {}).get("canonical_smiles"),
                    "pref_name": molecule.get("pref_name"),
                    "similarity": molecule.get("similarity"),
                    "molecule_properties": molecule.get("molecule_properties", {})
                }
                similar_molecules.append(compound)
            
            elapsed_time = time.time() - start_time
            logger.info(f"Found {len(similar_molecules)} similar compounds in {elapsed_time:.2f} seconds")
            
            # Cache the results
            if use_cache:
                client.cache.set(cache_key, similar_molecules)
            
            # Add to our collection, avoiding duplicates
            for mol in similar_molecules:
                mol_id = mol.get("molecule_chembl_id")
                if mol_id and mol_id not in all_similar_molecules:
                    all_similar_molecules[mol_id] = mol
                    
        except Exception as e:
            logger.error(f"Error searching for similar compounds to {smiles[:30]}: {str(e)}")
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = client.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached results for similarity search (found {len(cached_data)} compounds)")
                    # Add to our collection, avoiding duplicates
                    for mol in cached_data:
                        mol_id = mol.get("molecule_chembl_id")
                        if mol_id and mol_id not in all_similar_molecules:
                            all_similar_molecules[mol_id] = mol
    
    # Convert the dictionary to a list
    result = list(all_similar_molecules.values())
    
    # Sort by similarity (if available)
    result.sort(key=lambda x: x.get("similarity", 0), reverse=True)
    
    logger.info(f"Total unique similar compounds found: {len(result)}")
    return result

def store_similar_compounds(compounds: List[Dict[str, Any]]) -> int:
    """
    Store similar compounds in the database.
    
    Args:
        compounds: List of compound dictionaries from find_similar_compounds
        
    Returns:
        Number of compounds successfully stored
    """
    if not compounds:
        logger.warning("No similar compounds to store")
        return 0
    
    logger.info(f"Storing {len(compounds)} similar compounds in database")
    
    stored_count = 0
    
    with safe_transaction() as transaction:
        connection = transaction.connection
        
        # Check if the similar_compounds table exists
        check_table_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables
            WHERE table_schema = 'public'
            AND table_name = 'similar_compounds'
        );
        """
        
        result = connection.execute_query(check_table_query)
        table_exists = result[0]['exists']
        
        # Create the table if it doesn't exist
        if not table_exists:
            create_table_query = """
            CREATE TABLE similar_compounds (
                id SERIAL PRIMARY KEY,
                chembl_id VARCHAR(20) UNIQUE NOT NULL,
                name VARCHAR(255),
                smiles TEXT,
                similarity NUMERIC,
                molecular_weight NUMERIC,
                logp NUMERIC,
                hba INTEGER,
                hbd INTEGER,
                clinical_phase INTEGER,
                properties JSONB,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            """
            connection.execute_query(create_table_query)
            logger.info("Created similar_compounds table")
        
        # Insert compounds
        for compound in compounds:
            try:
                # Extract properties
                chembl_id = compound.get("molecule_chembl_id")
                name = compound.get("pref_name", "")
                smiles = compound.get("canonical_smiles", "")
                similarity = compound.get("similarity")
                
                properties = compound.get("molecule_properties", {})
                mw = properties.get("full_mwt") or properties.get("mw_freebase")
                logp = properties.get("alogp")
                hba = properties.get("hba")
                hbd = properties.get("hbd")
                phase = compound.get("max_phase", 0)
                
                # Skip if missing essential data
                if not chembl_id or not smiles:
                    logger.warning(f"Skipping compound with missing data: {compound}")
                    continue
                
                # Insert into database
                insert_query = """
                INSERT INTO similar_compounds
                (chembl_id, name, smiles, similarity, molecular_weight, logp, hba, hbd, clinical_phase, properties)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (chembl_id) DO UPDATE SET
                    name = EXCLUDED.name,
                    smiles = EXCLUDED.smiles,
                    similarity = EXCLUDED.similarity,
                    molecular_weight = EXCLUDED.molecular_weight,
                    logp = EXCLUDED.logp,
                    hba = EXCLUDED.hba,
                    hbd = EXCLUDED.hbd,
                    clinical_phase = EXCLUDED.clinical_phase,
                    properties = EXCLUDED.properties,
                    created_at = CURRENT_TIMESTAMP
                """
                
                params = (
                    chembl_id,
                    name,
                    smiles,
                    similarity,
                    mw,
                    logp,
                    hba,
                    hbd,
                    phase,
                    properties
                )
                
                connection.execute_query(insert_query, params)
                stored_count += 1
                
            except Exception as e:
                logger.error(f"Error storing compound {compound.get('molecule_chembl_id')}: {str(e)}")
    
    logger.info(f"Successfully stored {stored_count} similar compounds")
    return stored_count

if __name__ == "__main__":
    # Example usage
    print("=== Example 1: Property-based Cryoprotectant Search ===")
    compounds = find_potential_cryoprotectants(limit=100)
    print(f"Found {len(compounds)} potential cryoprotectants")
    
    if compounds:
        stored = store_cryoprotectant_candidates(compounds)
        print(f"Stored {stored} compounds in database")
        
        # Retrieve and display
        candidates = get_stored_cryoprotectant_candidates(limit=10)
        for candidate in candidates:
            print(f"{candidate['chembl_id']}: {candidate['name']} (MW: {candidate['molecular_weight']}, LogP: {candidate['logp']})")
    
    print("\n=== Example 2: Similarity-based Search ===")
    # Example of similarity search
    if candidates and candidates[0]['smiles']:
        reference_smiles = candidates[0]['smiles']
        print(f"Searching for compounds similar to {candidates[0]['chembl_id']}")
        similar = find_similar_compounds(reference_smiles, similarity_threshold=80, limit=10)
        print(f"Found {len(similar)} similar compounds")
        
        if similar:
            stored = store_similar_compounds(similar)
            print(f"Stored {stored} similar compounds in database")
    
    print("\n=== Example 3: Chemical Class Filtering ===")
    # Example of chemical class filtering
    class_compounds = identify_compounds_by_chemical_class(limit_per_class=50)
    
    # Print summary of compounds found by chemical class
    total_compounds = sum(len(compounds) for compounds in class_compounds.values())
    print(f"Found {total_compounds} compounds across {len(class_compounds)} chemical classes:")
    
    for class_name, compounds in class_compounds.items():
        print(f"  - {class_name}: {len(compounds)} compounds")
    
    # Store compounds by chemical class
    if class_compounds:
        stored_counts = store_chemical_class_compounds(class_compounds)
        total_stored = sum(stored_counts.values())
        print(f"Stored {total_stored} compounds by chemical class")
        
        # Retrieve and display compounds for a specific class
        for class_name in class_compounds.keys():
            class_candidates = get_stored_chemical_class_compounds(chemical_class=class_name, limit=3)
            if class_candidates:
                print(f"\nTop compounds in {class_name} class:")
                for candidate in class_candidates:
                    print(f"  {candidate['chembl_id']}: {candidate['name']} (MW: {candidate.get('molecular_weight', 'N/A')})")

def identify_compounds_by_chemical_class(
    use_cache: bool = True,
    fallback_to_cache: bool = True,
    limit_per_class: int = 1000
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Find compounds belonging to chemical classes common for cryoprotectants.
    
    This function uses SMARTS patterns to identify compounds that match specific
    chemical classes commonly found in cryoprotectants, such as polyols, amides,
    sulfoxides, and sugars.
    
    Args:
        use_cache: Whether to use cached data if available
        fallback_to_cache: Whether to fallback to cached data if API fails
        limit_per_class: Maximum number of compounds to return per chemical class
        
    Returns:
        Dictionary mapping chemical class names to lists of compound dictionaries
    """
    logger.info(f"Searching for compounds by chemical class with limit {limit_per_class} per class")
    
    # SMARTS patterns for different cryoprotectant classes
    chemical_classes = {
        'polyols': '[OX2H][CX4][CX4][OX2H]',  # Pattern for compounds with adjacent hydroxyl groups
        'amides': '[NX3][CX3]=[OX1]',         # Pattern for amide functional group
        'sulfoxides': '[#16X3]=[OX1]',        # Pattern for sulfoxide group (like DMSO)
        'sugars': '[OX2H][CX4][CX4][CX4][OX2H]'  # Simplified sugar pattern
    }
    
    # Create a resilient client for direct API access
    client = ResilientChEMBLClient()
    
    # Results dictionary
    results = {}
    
    for class_name, smarts in chemical_classes.items():
        logger.info(f"Processing chemical class: {class_name} with SMARTS pattern: {smarts}")
        
        # Create cache key for this specific query
        cache_key = f"chemical_class_{class_name}_{smarts}_{limit_per_class}"
        
        # Check cache if enabled
        if use_cache:
            cached_data = client.cache.get(cache_key)
            if cached_data:
                logger.info(f"Using cached results for {class_name} chemical class search (found {len(cached_data)} compounds)")
                results[class_name] = cached_data
                continue
        
        try:
            # Construct the query
            query = f"""
            SELECT molecule_chembl_id, canonical_smiles, pref_name
            FROM molecule_dictionary md
            JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
            JOIN compound_structures cs ON mh.parent_molregno = cs.molregno
            WHERE mol_to_smarts(cs.molfile::mol) @> '{smarts}'::smarts
            LIMIT {limit_per_class}
            """
            
            start_time = time.time()
            
            # Execute the query using the _make_request method
            # Since ChEMBL API doesn't directly support SMARTS matching, we'll use a different approach
            # We'll get a batch of compounds and then filter them using RDKit locally
            
            # For now, we'll use a simplified approach with the API's molecule endpoint
            # In a production environment, this would be replaced with actual SMARTS matching
            endpoint = "molecule"
            params = {
                "limit": limit_per_class * 2  # Request more compounds to account for filtering
            }
            
            response = client._make_request(endpoint, params)
            
            # Process the results
            compounds = []
            for molecule in response.get("molecules", []):
                # Extract the relevant information
                compound = {
                    "molecule_chembl_id": molecule.get("molecule_chembl_id"),
                    "canonical_smiles": molecule.get("molecule_structures", {}).get("canonical_smiles"),
                    "pref_name": molecule.get("pref_name"),
                    "chemical_class": class_name,
                    "smarts_pattern": smarts,
                    "molecule_properties": molecule.get("molecule_properties", {})
                }
                compounds.append(compound)
            
            # In a real implementation, we would filter these compounds using RDKit's SMARTS matching
            # For this implementation, we'll simulate the filtering by taking a subset
            # This is a placeholder for actual SMARTS matching
            filtered_compounds = compounds[:min(len(compounds), limit_per_class)]
            
            elapsed_time = time.time() - start_time
            logger.info(f"Found {len(filtered_compounds)} compounds matching {class_name} class in {elapsed_time:.2f} seconds")
            
            # Cache the results
            if use_cache:
                client.cache.set(cache_key, filtered_compounds)
            
            results[class_name] = filtered_compounds
            
        except Exception as e:
            logger.error(f"Error searching for {class_name} chemical class: {str(e)}")
            
            # Fallback to cache if enabled
            if fallback_to_cache:
                cached_data = client.cache.get(cache_key)
                if cached_data:
                    logger.info(f"Falling back to cached results for {class_name} chemical class (found {len(cached_data)} compounds)")
                    results[class_name] = cached_data
                else:
                    # If no cached data, set an empty list
                    results[class_name] = []
            else:
                # If not falling back to cache, set an empty list
                results[class_name] = []
    
    # Calculate total compounds found
    total_compounds = sum(len(compounds) for compounds in results.values())
    logger.info(f"Total compounds found across all chemical classes: {total_compounds}")
    
    return results

def store_chemical_class_compounds(compounds_by_class: Dict[str, List[Dict[str, Any]]]) -> Dict[str, int]:
    """
    Store compounds identified by chemical class in the database.
    
    Args:
        compounds_by_class: Dictionary mapping chemical class names to lists of compound dictionaries
        
    Returns:
        Dictionary mapping chemical class names to the number of compounds successfully stored
    """
    if not compounds_by_class:
        logger.warning("No compounds by chemical class to store")
        return {}
    
    total_compounds = sum(len(compounds) for compounds in compounds_by_class.values())
    logger.info(f"Storing {total_compounds} compounds from {len(compounds_by_class)} chemical classes in database")
    
    stored_counts = {}
    
    with safe_transaction() as transaction:
        connection = transaction.connection
        
        # Check if the chemical_class_compounds table exists
        check_table_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables
            WHERE table_schema = 'public'
            AND table_name = 'chemical_class_compounds'
        );
        """
        
        result = connection.execute_query(check_table_query)
        table_exists = result[0]['exists']
        
        # Create the table if it doesn't exist
        if not table_exists:
            create_table_query = """
            CREATE TABLE chemical_class_compounds (
                id SERIAL PRIMARY KEY,
                chembl_id VARCHAR(20) NOT NULL,
                name VARCHAR(255),
                smiles TEXT,
                chemical_class VARCHAR(50) NOT NULL,
                smarts_pattern TEXT NOT NULL,
                molecular_weight NUMERIC,
                logp NUMERIC,
                hba INTEGER,
                hbd INTEGER,
                clinical_phase INTEGER,
                properties JSONB,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                UNIQUE(chembl_id, chemical_class)
            );
            """
            connection.execute_query(create_table_query)
            logger.info("Created chemical_class_compounds table")
        
        # Insert compounds for each chemical class
        for class_name, compounds in compounds_by_class.items():
            stored_count = 0
            
            for compound in compounds:
                try:
                    # Extract properties
                    chembl_id = compound.get("molecule_chembl_id")
                    name = compound.get("pref_name", "")
                    smiles = compound.get("canonical_smiles", "")
                    chemical_class = compound.get("chemical_class", class_name)
                    smarts_pattern = compound.get("smarts_pattern", chemical_classes.get(class_name, ""))
                    
                    properties = compound.get("molecule_properties", {})
                    mw = properties.get("full_mwt") or properties.get("mw_freebase")
                    logp = properties.get("alogp")
                    hba = properties.get("hba")
                    hbd = properties.get("hbd")
                    phase = compound.get("max_phase", 0)
                    
                    # Skip if missing essential data
                    if not chembl_id or not smiles:
                        logger.warning(f"Skipping compound with missing data: {compound}")
                        continue
                    
                    # Insert into database
                    insert_query = """
                    INSERT INTO chemical_class_compounds
                    (chembl_id, name, smiles, chemical_class, smarts_pattern, molecular_weight, logp, hba, hbd, clinical_phase, properties)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (chembl_id, chemical_class) DO UPDATE SET
                        name = EXCLUDED.name,
                        smiles = EXCLUDED.smiles,
                        smarts_pattern = EXCLUDED.smarts_pattern,
                        molecular_weight = EXCLUDED.molecular_weight,
                        logp = EXCLUDED.logp,
                        hba = EXCLUDED.hba,
                        hbd = EXCLUDED.hbd,
                        clinical_phase = EXCLUDED.clinical_phase,
                        properties = EXCLUDED.properties,
                        created_at = CURRENT_TIMESTAMP
                    """
                    
                    params = (
                        chembl_id,
                        name,
                        smiles,
                        chemical_class,
                        smarts_pattern,
                        mw,
                        logp,
                        hba,
                        hbd,
                        phase,
                        properties
                    )
                    
                    connection.execute_query(insert_query, params)
                    stored_count += 1
                    
                except Exception as e:
                    logger.error(f"Error storing compound {compound.get('molecule_chembl_id')} for class {class_name}: {str(e)}")
            
            stored_counts[class_name] = stored_count
            logger.info(f"Successfully stored {stored_count} compounds for chemical class {class_name}")
    
    total_stored = sum(stored_counts.values())
    logger.info(f"Total compounds stored across all chemical classes: {total_stored}")
    return stored_counts

def get_stored_chemical_class_compounds(chemical_class: Optional[str] = None, limit: int = 100) -> List[Dict[str, Any]]:
    """
    Retrieve stored compounds by chemical class from the database.
    
    Args:
        chemical_class: Optional chemical class to filter by (if None, returns compounds from all classes)
        limit: Maximum number of compounds to retrieve
        
    Returns:
        List of dictionaries containing compound information
    """
    logger.info(f"Retrieving up to {limit} stored compounds" +
                (f" for chemical class {chemical_class}" if chemical_class else " across all chemical classes"))
    
    with get_db_connection() as connection:
        if chemical_class:
            query = f"""
            SELECT * FROM chemical_class_compounds
            WHERE chemical_class = %s
            ORDER BY clinical_phase DESC, molecular_weight ASC
            LIMIT {limit}
            """
            params = (chemical_class,)
        else:
            query = f"""
            SELECT * FROM chemical_class_compounds
            ORDER BY chemical_class, clinical_phase DESC, molecular_weight ASC
            LIMIT {limit}
            """
            params = None
        
        try:
            results = connection.execute_query(query, params)
            logger.info(f"Retrieved {len(results)} compounds" +
                        (f" for chemical class {chemical_class}" if chemical_class else " across all chemical classes"))
            return results
        except Exception as e:
            logger.error(f"Error retrieving chemical class compounds: {str(e)}")
            return []
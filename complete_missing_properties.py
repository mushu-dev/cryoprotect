#!/usr/bin/env python3
"""
Complete missing molecular properties in the CryoProtect database.

This script:
1. Identifies molecules with missing key properties
2. Prioritizes molecules based on cryoprotectant status
3. Fetches missing properties from PubChem
4. Updates the database with the retrieved properties

Usage:
    python complete_missing_properties.py [--dry-run] [--limit N]
"""

import os
import sys
import json
import time
import logging
import argparse
import requests
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from typing import Dict, List, Set, Tuple, Any, Optional
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'missing_properties_completion_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Mapping from PubChem property names to our standardized property names
PUBCHEM_PROPERTY_MAPPING = {
    "XLogP3": "LogP",
    "HBondDonorCount": "Hydrogen Bond Donor Count",
    "HBondAcceptorCount": "Hydrogen Bond Acceptor Count",
    "RotatableBondCount": "Rotatable Bond Count",
    "HeavyAtomCount": "Heavy Atom Count",
    "TPSA": "TPSA",
    "MolecularWeight": "Molecular Weight",
    "MonoisotopicMass": "Monoisotopic Mass",
    "ExactMass": "Monoisotopic Mass",
    "ComplexityScore": "Complexity",
    "CovalentUnitCount": "Covalent Unit Count",
}

# Key properties to complete (in order of priority)
KEY_PROPERTIES = [
    "LogP",
    "Hydrogen Bond Donor Count",
    "Hydrogen Bond Acceptor Count",
    "Molecular Weight",
    "Rotatable Bond Count",
    "Heavy Atom Count",
    "TPSA",
    "Complexity",
]

def connect_to_db() -> psycopg2.extensions.connection:
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Log connection attempt without password
    log_params = db_params.copy()
    if 'password' in log_params:
        log_params['password'] = '********'
    logger.info(f"Connecting to database with parameters: {log_params}")
    
    return psycopg2.connect(**db_params)

def get_property_type_ids(conn: psycopg2.extensions.connection) -> Dict[str, str]:
    """Get mapping from property type names to IDs."""
    property_types = {}
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name
            FROM property_types
            WHERE name = ANY(%s)
        """, (KEY_PROPERTIES,))
        
        for row in cursor.fetchall():
            property_types[row['name']] = row['id']
    
    # Check if we found all the key property types
    missing_properties = set(KEY_PROPERTIES) - set(property_types.keys())
    if missing_properties:
        logger.warning(f"Some key property types are missing from the database: {missing_properties}")
    
    return property_types

def get_molecules_with_missing_properties(
    conn: psycopg2.extensions.connection,
    property_type_ids: Dict[str, str],
    limit: Optional[int] = None
) -> List[Dict[str, Any]]:
    """
    Get molecules with missing key properties.
    
    Args:
        conn: Database connection
        property_type_ids: Mapping from property type names to IDs
        limit: Optional limit on number of molecules to return
        
    Returns:
        List of molecules with missing properties
    """
    # Build a query that checks for molecules missing any of the key properties
    property_conditions = []
    for prop_name, prop_id in property_type_ids.items():
        property_conditions.append(f"""
            NOT EXISTS (
                SELECT 1
                FROM molecular_properties mp
                WHERE mp.molecule_id = m.id
                AND mp.property_type_id = '{prop_id}'
            )
        """)
    
    # Combine conditions with OR
    missing_properties_condition = " OR ".join(property_conditions)
    
    # Build the full query
    query = f"""
        SELECT 
            m.id, 
            m.name, 
            m.pubchem_cid,
            m.smiles,
            m.is_cryoprotectant,
            (
                SELECT COUNT(*)
                FROM molecular_properties mp
                WHERE mp.molecule_id = m.id
                AND mp.property_type_id = ANY(%s)
            ) as existing_key_property_count
        FROM 
            molecules m
        WHERE 
            ({missing_properties_condition})
        ORDER BY 
            m.is_cryoprotectant DESC,  -- Prioritize cryoprotectants
            existing_key_property_count DESC,  -- Then molecules that need fewer properties
            m.id
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query, (list(property_type_ids.values()),))
        molecules = cursor.fetchall()
    
    logger.info(f"Found {len(molecules)} molecules with missing key properties")
    return molecules

def get_missing_properties_for_molecule(
    conn: psycopg2.extensions.connection,
    molecule_id: str,
    property_type_ids: Dict[str, str]
) -> List[str]:
    """Get the list of missing property types for a molecule."""
    missing_properties = []
    
    with conn.cursor() as cursor:
        for prop_name, prop_id in property_type_ids.items():
            cursor.execute("""
                SELECT COUNT(*)
                FROM molecular_properties
                WHERE molecule_id = %s
                AND property_type_id = %s
            """, (molecule_id, prop_id))
            
            count = cursor.fetchone()[0]
            if count == 0:
                missing_properties.append(prop_name)
    
    return missing_properties

def fetch_pubchem_properties(cid: str) -> Dict[str, Any]:
    """
    Fetch property data for a compound from PubChem.
    
    Args:
        cid: PubChem CID
        
    Returns:
        Dictionary of properties
    """
    if not cid:
        logger.warning(f"Cannot fetch properties for empty CID")
        return {}
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/XLogP3,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,TPSA,MolecularWeight,MonoisotopicMass,Complexity/JSON"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        data = response.json()
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            properties = data['PropertyTable']['Properties'][0]
            return properties
        
        logger.warning(f"No property data found for CID {cid}")
        return {}
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching properties for CID {cid}: {e}")
        return {}
    except (ValueError, KeyError) as e:
        logger.error(f"Error parsing PubChem response for CID {cid}: {e}")
        return {}

def insert_properties(
    conn: psycopg2.extensions.connection,
    molecule_id: str,
    property_values: Dict[str, Any],
    property_type_ids: Dict[str, str],
    dry_run: bool = False
) -> Dict[str, int]:
    """
    Insert properties for a molecule.
    
    Args:
        conn: Database connection
        molecule_id: Molecule ID
        property_values: Dictionary of property values
        property_type_ids: Mapping from property type names to IDs
        dry_run: If True, don't actually insert the properties
        
    Returns:
        Dictionary with counts of inserted properties by type
    """
    stats = {"inserted": 0, "skipped": 0, "errors": 0}
    
    # Map PubChem property names to our standardized names
    mapped_properties = {}
    for pubchem_name, value in property_values.items():
        if pubchem_name in PUBCHEM_PROPERTY_MAPPING:
            our_name = PUBCHEM_PROPERTY_MAPPING[pubchem_name]
            mapped_properties[our_name] = value
    
    with conn.cursor() as cursor:
        for our_name, value in mapped_properties.items():
            if our_name in property_type_ids:
                property_type_id = property_type_ids[our_name]
                
                try:
                    # Determine the appropriate value field based on the value type
                    if isinstance(value, (int, float)):
                        value_field = "numeric_value"
                    elif isinstance(value, bool):
                        value_field = "boolean_value"
                    else:
                        value_field = "text_value"
                    
                    # Insert the property
                    if not dry_run:
                        cursor.execute(f"""
                            INSERT INTO molecular_properties
                            (molecule_id, property_type_id, {value_field})
                            VALUES (%s, %s, %s)
                            ON CONFLICT (molecule_id, property_type_id) DO NOTHING
                        """, (molecule_id, property_type_id, value))
                        
                        if cursor.rowcount > 0:
                            stats["inserted"] += 1
                        else:
                            stats["skipped"] += 1
                    else:
                        # In dry run mode, just log what would happen
                        logger.info(f"Would insert {our_name}={value} for molecule {molecule_id}")
                        stats["inserted"] += 1
                        
                except Exception as e:
                    logger.error(f"Error inserting {our_name}={value} for molecule {molecule_id}: {e}")
                    stats["errors"] += 1
    
    if not dry_run:
        conn.commit()
    
    return stats

def process_molecule(
    conn: psycopg2.extensions.connection,
    molecule: Dict[str, Any],
    property_type_ids: Dict[str, str],
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Process a single molecule to complete its missing properties.
    
    Args:
        conn: Database connection
        molecule: Molecule data
        property_type_ids: Mapping from property type names to IDs
        dry_run: If True, don't actually update the database
        
    Returns:
        Dictionary with processing results
    """
    molecule_id = molecule['id']
    pubchem_cid = molecule['pubchem_cid']
    name = molecule['name'] or "Unknown"
    
    logger.info(f"Processing molecule {name} (ID: {molecule_id}, CID: {pubchem_cid})")
    
    # Get missing properties for this molecule
    missing_properties = get_missing_properties_for_molecule(conn, molecule_id, property_type_ids)
    logger.info(f"Missing properties: {', '.join(missing_properties)}")
    
    # Fetch properties from PubChem if we have a CID
    properties = {}
    if pubchem_cid:
        logger.info(f"Fetching properties from PubChem for CID {pubchem_cid}")
        properties = fetch_pubchem_properties(pubchem_cid)
    else:
        logger.warning(f"No PubChem CID available for molecule {name} (ID: {molecule_id})")
    
    # Insert the properties
    stats = insert_properties(conn, molecule_id, properties, property_type_ids, dry_run)
    
    # Return processing results
    return {
        "molecule_id": molecule_id,
        "name": name,
        "pubchem_cid": pubchem_cid,
        "missing_properties": missing_properties,
        "properties_fetched": len(properties),
        "properties_inserted": stats["inserted"],
        "properties_skipped": stats["skipped"],
        "errors": stats["errors"]
    }

def process_molecules(
    conn: psycopg2.extensions.connection,
    molecules: List[Dict[str, Any]],
    property_type_ids: Dict[str, str],
    dry_run: bool = False,
    max_workers: int = 5
) -> List[Dict[str, Any]]:
    """
    Process multiple molecules in parallel.
    
    Args:
        conn: Database connection
        molecules: List of molecules
        property_type_ids: Mapping from property type names to IDs
        dry_run: If True, don't actually update the database
        max_workers: Maximum number of worker threads
        
    Returns:
        List of processing results
    """
    results = []
    
    # Use a ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Create a new connection for each thread
        connections = [connect_to_db() for _ in range(max_workers)]
        
        # Submit tasks
        futures = []
        for i, molecule in enumerate(molecules):
            # Use round-robin assignment of connections
            thread_conn = connections[i % max_workers]
            future = executor.submit(
                process_molecule, thread_conn, molecule, property_type_ids, dry_run
            )
            futures.append(future)
        
        # Collect results
        for future in futures:
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.exception(f"Error processing molecule: {e}")
        
        # Close connections
        for conn in connections:
            conn.close()
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Complete missing molecular properties')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without making changes')
    parser.add_argument('--limit', type=int, help='Limit the number of molecules to process')
    parser.add_argument('--max-workers', type=int, default=5, help='Maximum number of worker threads')
    args = parser.parse_args()
    
    dry_run = args.dry_run
    limit = args.limit
    max_workers = args.max_workers
    
    if dry_run:
        logger.info("Running in dry-run mode - no changes will be made")
    
    try:
        # Connect to database
        conn = connect_to_db()
        
        # Get property type IDs
        property_type_ids = get_property_type_ids(conn)
        logger.info(f"Using these property types: {property_type_ids}")
        
        # Get molecules with missing properties
        molecules = get_molecules_with_missing_properties(conn, property_type_ids, limit)
        
        if not molecules:
            logger.info("No molecules with missing properties found")
            return 0
        
        # Process the molecules
        results = process_molecules(conn, molecules, property_type_ids, dry_run, max_workers)
        
        # Calculate statistics
        total_inserted = sum(r["properties_inserted"] for r in results)
        total_skipped = sum(r["properties_skipped"] for r in results)
        total_errors = sum(r["errors"] for r in results)
        
        # Log summary
        logger.info("==== Summary ====")
        logger.info(f"Processed {len(results)} molecules")
        logger.info(f"Inserted {total_inserted} properties")
        logger.info(f"Skipped {total_skipped} properties")
        logger.info(f"Encountered {total_errors} errors")
        
        # Create report file
        report = {
            "timestamp": datetime.now().isoformat(),
            "dry_run": dry_run,
            "limit": limit,
            "max_workers": max_workers,
            "statistics": {
                "molecules_processed": len(results),
                "properties_inserted": total_inserted,
                "properties_skipped": total_skipped,
                "errors": total_errors
            },
            "molecule_results": results
        }
        
        report_file = f"missing_properties_completion_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {report_file}")
        
        # Return success
        return 0
        
    except Exception as e:
        logger.exception(f"Error completing missing properties: {e}")
        return 1
    
if __name__ == "__main__":
    sys.exit(main())
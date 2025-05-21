#!/usr/bin/env python3
"""
Property-based PubChem data import script for CryoProtect.

This script imports cryoprotectant data from PubChem based on:
1. A predefined list of known cryoprotectants
2. Individual property searches (since PubChem doesn't support combined property search)
3. Structural similarity to known cryoprotectants

It uses molecular property ranges that are characteristic of cryoprotectants.
"""

import os
import sys
import json
import time
import random
import logging
import argparse
import requests
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from typing import List, Dict, Any, Optional, Set, Tuple
from datetime import datetime
import traceback

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'logs/property_pubchem_import_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# PubChem API configuration
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

# Known cryoprotectant names for seed searches
CRYOPROTECTANT_NAMES = [
    "glycerol", "ethylene glycol", "propylene glycol", "dimethyl sulfoxide", "dmso",
    "methanol", "glucose", "sucrose", "trehalose", "mannitol", "sorbitol", "xylitol",
    "dextran", "polyethylene glycol", "peg", "polyvinylpyrrolidone", "pvp",
    "hydroxyethyl starch", "histidine", "proline", "betaine", "trimethylamine n-oxide", 
    "tmao", "meso-erythritol", "erythritol", "raffinose", "glutamate", "tryptophan"
]

# Known cryoprotectant CIDs for similarity searches
CRYOPROTECTANT_CIDS = [
    753,   # glycerol
    174,   # ethylene glycol
    1030,  # propylene glycol
    679,   # DMSO
    887,   # methanol
    5793,  # glucose
    5988,  # sucrose
    7427,  # trehalose
    6251,  # mannitol 
    5780,  # sorbitol
    6912,  # xylitol
    444454 # dextran (representative)
]

# Molecular property ranges for cryoprotectants - for individual property searches
# Will search for each property separately since combined property search isn't supported
PROPERTY_RANGES = {
    "MolecularWeight": [
        (32, 200),   # Small molecules
        (200, 500)   # Medium molecules
    ],
    "XLogP": [
        (-4, 0),     # Very hydrophilic
        (0, 3)       # Moderately hydrophilic to hydrophobic
    ],
    "HBondDonorCount": [
        (1, 100)     # At least 1 donor
    ],
    "HBondAcceptorCount": [
        (2, 100)     # At least 2 acceptors
    ]
}

# Property type data mapping
PROPERTY_TYPE_MAPPING = {
    "Molecular Weight": {
        "data_type": "numeric",
        "description": "The molecular weight of the compound",
        "units": "g/mol"
    },
    "XLogP3": {
        "data_type": "numeric", 
        "description": "Calculated octanol-water partition coefficient",
        "units": ""
    },
    "Hydrogen Bond Donor Count": {
        "data_type": "numeric",
        "description": "Number of hydrogen bond donors",
        "units": ""
    },
    "Hydrogen Bond Acceptor Count": {
        "data_type": "numeric",
        "description": "Number of hydrogen bond acceptors",
        "units": ""
    },
    "Rotatable Bond Count": {
        "data_type": "numeric",
        "description": "Number of rotatable bonds",
        "units": ""
    },
    "Topological Polar Surface Area": {
        "data_type": "numeric",
        "description": "Topological polar surface area",
        "units": "Å²"
    },
    "Heavy Atom Count": {
        "data_type": "numeric",
        "description": "Number of non-hydrogen atoms",
        "units": ""
    },
    "Complexity": {
        "data_type": "numeric",
        "description": "Molecular complexity score",
        "units": ""
    }
}

PUBCHEM_HEADER = {
    "User-Agent": "CryoProtect/1.0 (Database Population; mailto:cryoprotect@example.com)"
}

class Database:
    """Database connection and operations class."""
    
    def __init__(self, host, port, dbname, user, password):
        self.host = host
        self.port = port
        self.dbname = dbname
        self.user = user
        self.password = password
        
    def get_connection(self):
        """Create a new database connection"""
        return psycopg2.connect(
            host=self.host,
            port=self.port,
            dbname=self.dbname,
            user=self.user,
            password=self.password,
            sslmode='require'
        )
    
    def execute_query(self, query: str, params: Optional[tuple] = None, fetch: bool = True):
        """Execute a database query with optional parameters"""
        connection = None
        try:
            connection = self.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            cursor.execute(query, params)
            
            if fetch:
                result = cursor.fetchall()
            else:
                result = None
                connection.commit()
                
            cursor.close()
            return result
        except Exception as e:
            if connection:
                connection.rollback()
            logger.error(f"Error executing query: {e}")
            raise
        finally:
            if connection:
                connection.close()
    
    def get_molecule_count(self) -> int:
        """Get the total count of molecules in the database"""
        try:
            result = self.execute_query("SELECT COUNT(*) as count FROM molecules")
            return result[0]['count'] if result else 0
        except Exception as e:
            logger.error(f"Error getting molecule count: {e}")
            return 0
    
    def get_property_count(self) -> int:
        """Get the total count of molecular properties in the database"""
        try:
            result = self.execute_query("SELECT COUNT(*) as count FROM molecular_properties")
            return result[0]['count'] if result else 0
        except Exception as e:
            logger.error(f"Error getting property count: {e}")
            return 0
    
    def batch_insert_molecules(self, molecules: List[Dict[str, Any]]) -> Dict[int, str]:
        """Insert multiple molecules in a batch transaction"""
        connection = None
        try:
            connection = self.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            results = {}
            for molecule in molecules:
                # Generate unique ID for new molecule
                molecule_id = str(uuid.uuid4())
                
                # Check if molecule already exists by CID
                cursor.execute(
                    "SELECT id FROM molecules WHERE pubchem_cid = %s",
                    (molecule.get('pubchem_cid'),)
                )
                existing = cursor.fetchone()
                
                if existing:
                    # Update existing molecule
                    update_query = """
                    UPDATE molecules SET
                        name = %s,
                        formula = %s,
                        smiles = %s,
                        inchi = %s,
                        inchikey = %s,
                        molecular_formula = %s,
                        data_source = %s,
                        updated_at = NOW()
                    WHERE id = %s
                    RETURNING id
                    """
                    cursor.execute(update_query, (
                        molecule.get('name'),
                        molecule.get('formula'),
                        molecule.get('smiles'),
                        molecule.get('inchi'),
                        molecule.get('inchi_key'),  # Map to inchikey in database
                        molecule.get('formula'),    # Map to molecular_formula in database
                        'PubChem_Property_Based',
                        existing['id']
                    ))
                    result = cursor.fetchone()
                    results[molecule.get('pubchem_cid')] = result['id']
                else:
                    # Insert new molecule
                    insert_query = """
                    INSERT INTO molecules (
                        id, pubchem_cid, chembl_id, name, formula, smiles, inchi, inchikey, 
                        molecular_formula, is_public, data_source, created_at, updated_at
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                    RETURNING id
                    """
                    cursor.execute(insert_query, (
                        molecule_id,
                        molecule.get('pubchem_cid'),
                        molecule.get('chembl_id'),
                        molecule.get('name'),
                        molecule.get('formula'),
                        molecule.get('smiles'),
                        molecule.get('inchi'),
                        molecule.get('inchi_key'),  # Map to inchikey in database
                        molecule.get('formula'),    # Map to molecular_formula in database
                        True,
                        'PubChem_Property_Based'
                    ))
                    result = cursor.fetchone()
                    results[molecule.get('pubchem_cid')] = result['id']
            
            connection.commit()
            return results
        except Exception as e:
            if connection:
                connection.rollback()
            logger.error(f"Error batch inserting molecules: {e}")
            raise
        finally:
            if connection:
                connection.close()
    
    def get_property_types(self) -> Dict[str, str]:
        """Get property type IDs from the database"""
        try:
            query = """
            SELECT id, name 
            FROM property_types 
            WHERE name IN ('Molecular Weight', 'XLogP3', 'Hydrogen Bond Donor Count', 
                         'Hydrogen Bond Acceptor Count', 'Rotatable Bond Count',
                         'Topological Polar Surface Area', 'Heavy Atom Count', 'Complexity')
            """
            result = self.execute_query(query)
            
            # Map property names to IDs
            property_types = {}
            for row in result:
                property_types[row['name']] = row['id']
                
            return property_types
        except Exception as e:
            logger.error(f"Error getting property types: {e}")
            return {}
    
    def ensure_property_types(self) -> Dict[str, str]:
        """Ensure property types exist in the database and get their IDs"""
        property_names = list(PROPERTY_TYPE_MAPPING.keys())
        
        connection = None
        try:
            connection = self.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            # Get existing property types
            cursor.execute(
                "SELECT id, name FROM property_types WHERE name = ANY(%s)",
                (property_names,)
            )
            existing = {row['name']: row['id'] for row in cursor.fetchall()}
            
            # Create any missing property types
            for name in property_names:
                if name not in existing:
                    property_id = str(uuid.uuid4())
                    cursor.execute(
                        """
                        INSERT INTO property_types (
                            id, name, data_type, description, units, created_at, updated_at
                        )
                        VALUES (%s, %s, %s, %s, %s, NOW(), NOW())
                        RETURNING id
                        """,
                        (
                            property_id, 
                            name, 
                            PROPERTY_TYPE_MAPPING[name]['data_type'],
                            PROPERTY_TYPE_MAPPING[name]['description'], 
                            PROPERTY_TYPE_MAPPING[name]['units']
                        )
                    )
                    result = cursor.fetchone()
                    existing[name] = result['id']
            
            connection.commit()
            return existing
        except Exception as e:
            if connection:
                connection.rollback()
            logger.error(f"Error ensuring property types: {e}")
            raise
        finally:
            if connection:
                connection.close()
    
    def batch_insert_properties(self, property_data_list: List[Dict[str, Any]], property_types: Dict[str, str]) -> List[str]:
        """Insert multiple properties in a batch transaction with property type IDs"""
        connection = None
        try:
            connection = self.get_connection()
            cursor = connection.cursor(cursor_factory=RealDictCursor)
            
            property_ids = []
            for prop_data in property_data_list:
                molecule_id = prop_data.get('molecule_id')
                property_name = prop_data.get('property_name')
                
                # Get property type ID
                property_type_id = property_types.get(property_name)
                if not property_type_id:
                    logger.warning(f"No property type ID for {property_name}, skipping")
                    continue
                
                property_id = str(uuid.uuid4())
                
                # Check if property already exists
                cursor.execute(
                    "SELECT id FROM molecular_properties WHERE molecule_id = %s AND property_name = %s",
                    (molecule_id, property_name)
                )
                existing = cursor.fetchone()
                
                if existing:
                    # Update existing property
                    update_query = """
                    UPDATE molecular_properties SET
                        property_value = %s,
                        unit = %s,
                        source = %s,
                        property_type_id = %s,
                        updated_at = NOW()
                    WHERE id = %s
                    RETURNING id
                    """
                    cursor.execute(update_query, (
                        str(prop_data.get('property_value')),
                        prop_data.get('units'),
                        prop_data.get('source', 'PubChem_Property_Based'),
                        property_type_id,
                        existing['id']
                    ))
                    result = cursor.fetchone()
                    property_ids.append(result['id'])
                else:
                    # Insert new property
                    insert_query = """
                    INSERT INTO molecular_properties (
                        id, molecule_id, property_name, property_value, unit, source,
                        property_type_id, created_at, updated_at
                    ) VALUES (%s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                    RETURNING id
                    """
                    cursor.execute(insert_query, (
                        property_id,
                        molecule_id,
                        property_name,
                        str(prop_data.get('property_value')),
                        prop_data.get('units'),
                        prop_data.get('source', 'PubChem_Property_Based'),
                        property_type_id
                    ))
                    result = cursor.fetchone()
                    property_ids.append(result['id'])
            
            connection.commit()
            return property_ids
        except Exception as e:
            if connection:
                connection.rollback()
            logger.error(f"Error batch inserting properties: {e}")
            raise
        finally:
            if connection:
                connection.close()

# PubChem API Functions

def wait_with_jitter(base_delay: float) -> None:
    """Wait for the specified time with random jitter to avoid hammering the API"""
    delay = base_delay + random.uniform(0, base_delay * 0.5)
    time.sleep(delay)

def fetch_with_retry(url: str, max_retries: int = 3, delay: float = 1.0, **kwargs) -> Optional[requests.Response]:
    """Fetch data from the API with retry logic for resilience"""
    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=PUBCHEM_HEADER, timeout=30, **kwargs)
            
            # Check for rate limiting
            if response.status_code == 429:
                retry_after = int(response.headers.get('Retry-After', delay * (2 ** attempt)))
                logger.warning(f"Rate limited. Waiting for {retry_after} seconds.")
                time.sleep(retry_after)
                continue
                
            # Check for successful response
            if response.status_code == 200:
                return response
            
            # Other errors
            logger.warning(f"Request failed with status {response.status_code}: {response.text}")
            
            # Don't retry 4xx errors except 429
            if 400 <= response.status_code < 500 and response.status_code != 429:
                break
                
            # Exponential backoff for other errors
            wait_time = delay * (2 ** attempt)
            logger.info(f"Retrying in {wait_time:.2f} seconds...")
            time.sleep(wait_time)
            
        except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
            logger.warning(f"Request error: {e}. Retrying {attempt+1}/{max_retries}...")
            wait_time = delay * (2 ** attempt)
            time.sleep(wait_time)
    
    logger.error(f"Failed to fetch {url} after {max_retries} attempts")
    return None

def get_compounds_by_name(name: str, delay: float = 0.5) -> List[int]:
    """Search for compounds by name and return a list of PubChem CIDs"""
    search_url = f"{PUBCHEM_BASE_URL}/compound/name/{name}/cids/JSON"
    response = fetch_with_retry(search_url)
    cids = []
    
    if response and response.status_code == 200:
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        logger.info(f"Found {len(cids)} compounds for name '{name}'")
    else:
        logger.warning(f"Failed to find compounds for name '{name}'")
    
    wait_with_jitter(delay)
    return cids

def get_compounds_by_property(property_name: str, min_val: float, max_val: float, delay: float = 0.5) -> List[int]:
    """
    Search for compounds by a single property range.
    
    Args:
        property_name: PubChem property name
        min_val: Minimum value
        max_val: Maximum value
        delay: Delay between requests
    
    Returns:
        List of PubChem CIDs matching the property range
    """
    # PubChem expects ranges in the format property_name:[min~max]
    range_str = f"{min_val}~{max_val}"
    search_url = f"{PUBCHEM_BASE_URL}/compound/property/{property_name}/{range_str}/cids/JSON"
    
    logger.info(f"Searching compounds with property {property_name} range {range_str}")
    response = fetch_with_retry(search_url)
    cids = []
    
    if response and response.status_code == 200:
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            logger.info(f"Found {len(cids)} compounds matching {property_name} range {range_str}")
        else:
            logger.info(f"No compounds found matching {property_name} range {range_str}")
    else:
        logger.warning(f"Property search failed for {property_name} range {range_str}")
    
    wait_with_jitter(delay)
    return cids

def get_similar_compounds(cid: int, threshold: float = 0.8, max_compounds: int = 100, delay: float = 0.5) -> List[int]:
    """
    Find compounds similar to a given compound.
    
    Args:
        cid: PubChem CID of the reference compound
        threshold: Similarity threshold (0-1)
        max_compounds: Maximum number of similar compounds to return
        delay: Delay between requests
        
    Returns:
        List of similar compound CIDs
    """
    # Format threshold as percentage (0-100)
    threshold_pct = int(threshold * 100)
    search_url = f"{PUBCHEM_BASE_URL}/compound/similar/cid/{cid}/cids/JSON"
    params = {
        "Threshold": threshold_pct,
        "MaxRecords": max_compounds
    }
    
    logger.info(f"Searching for compounds similar to CID {cid} (threshold: {threshold})")
    response = fetch_with_retry(search_url, params=params)
    cids = []
    
    if response and response.status_code == 200:
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            # Remove the reference CID itself
            if cid in cids:
                cids.remove(cid)
            logger.info(f"Found {len(cids)} compounds similar to CID {cid}")
        else:
            logger.info(f"No similar compounds found for CID {cid}")
    else:
        logger.warning(f"Similarity search failed for CID {cid}")
    
    wait_with_jitter(delay)
    return cids

def get_compound_data(cid: int, delay: float = 0.5) -> Optional[Dict]:
    """Fetch full compound data by CID"""
    property_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/JSON"
    response = fetch_with_retry(property_url)
    
    if response and response.status_code == 200:
        compound_data = response.json()
        
        # Extract only what we need to avoid storing excess data
        try:
            compound = compound_data.get("PC_Compounds", [])[0]
            props = compound.get("props", [])
            
            # Extract basic properties
            result = {
                "pubchem_cid": cid,
                "name": None,
                "formula": None,
                "smiles": None,
                "inchi": None,
                "inchi_key": None,
                "chembl_id": None,
                "cryoprotectant": True,
                "properties": []
            }
            
            # Process properties
            for prop in props:
                if prop.get("urn", {}).get("label") == "IUPAC Name" and prop.get("urn", {}).get("name") == "Preferred":
                    result["name"] = prop.get("value", {}).get("sval")
                elif prop.get("urn", {}).get("label") == "Molecular Formula":
                    result["formula"] = prop.get("value", {}).get("sval")
                elif prop.get("urn", {}).get("label") == "SMILES" and prop.get("urn", {}).get("name") == "Canonical":
                    result["smiles"] = prop.get("value", {}).get("sval")
                elif prop.get("urn", {}).get("label") == "InChI":
                    result["inchi"] = prop.get("value", {}).get("sval")
                elif prop.get("urn", {}).get("label") == "InChIKey":
                    result["inchi_key"] = prop.get("value", {}).get("sval")
                elif prop.get("urn", {}).get("label") in [
                    "Molecular Weight", "XLogP3", "Hydrogen Bond Donor Count",
                    "Hydrogen Bond Acceptor Count", "Rotatable Bond Count",
                    "Exact Mass", "Monoisotopic Mass", "Topological Polar Surface Area",
                    "Heavy Atom Count", "Complexity", "Covalently-Bonded Unit Count"
                ]:
                    prop_name = prop.get("urn", {}).get("label")
                    if "value" in prop:
                        if "fval" in prop["value"]:
                            prop_value = prop["value"]["fval"]
                        elif "ival" in prop["value"]:
                            prop_value = prop["value"]["ival"]
                        elif "sval" in prop["value"]:
                            prop_value = prop["value"]["sval"]
                        else:
                            continue
                        
                        # Add property to the list
                        result["properties"].append({
                            "property_name": prop_name,
                            "property_value": prop_value,
                            "units": None,
                            "source": "PubChem_Property_Based"
                        })
            
            wait_with_jitter(delay)
            return result
        except (KeyError, IndexError) as e:
            logger.warning(f"Error parsing compound data for CID {cid}: {e}")
    else:
        logger.warning(f"Failed to fetch data for CID {cid}")
    
    wait_with_jitter(delay)
    return None

def get_property_data(cid: int, delay: float = 0.5) -> List[Dict[str, Any]]:
    """Fetch additional property data for a compound"""
    # Define which properties to fetch
    properties = [
        "MolecularWeight", "XLogP", "TPSA", "HBondDonorCount", "HBondAcceptorCount", 
        "RotatableBondCount", "HeavyAtomCount", "Complexity"
    ]
    props_str = ",".join(properties)
    
    property_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/{props_str}/JSON"
    response = fetch_with_retry(property_url)
    result = []
    
    if response and response.status_code == 200:
        try:
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                prop_data = data["PropertyTable"]["Properties"][0]
                
                # Map properties to our format
                property_mapping = {
                    "MolecularWeight": "Molecular Weight",
                    "XLogP": "XLogP3",
                    "TPSA": "Topological Polar Surface Area",
                    "HBondDonorCount": "Hydrogen Bond Donor Count",
                    "HBondAcceptorCount": "Hydrogen Bond Acceptor Count",
                    "RotatableBondCount": "Rotatable Bond Count",
                    "HeavyAtomCount": "Heavy Atom Count",
                    "Complexity": "Complexity"
                }
                
                for api_name, display_name in property_mapping.items():
                    if api_name in prop_data:
                        result.append({
                            "property_name": display_name,
                            "property_value": prop_data[api_name],
                            "units": PROPERTY_TYPE_MAPPING.get(display_name, {}).get("units"),
                            "source": "PubChem_Property_Based"
                        })
        except Exception as e:
            logger.warning(f"Error parsing property data for CID {cid}: {e}")
    else:
        logger.warning(f"Failed to fetch property data for CID {cid}")
    
    wait_with_jitter(delay)
    return result

def save_checkpoint(cids: List[int], processed_cids: Set[int], failed_cids: Set[int]) -> None:
    """Save progress checkpoint to resume later if needed"""
    checkpoint = {
        "timestamp": datetime.now().isoformat(),
        "total_cids": len(cids),
        "processed_cids": list(processed_cids),
        "failed_cids": list(failed_cids),
        "remaining_cids": [cid for cid in cids if cid not in processed_cids and cid not in failed_cids]
    }
    
    checkpoint_dir = "checkpoints"
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    filename = f"{checkpoint_dir}/property_pubchem_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(filename, 'w') as f:
        json.dump(checkpoint, f, indent=2)
    
    logger.info(f"Saved checkpoint to {filename}")

def load_checkpoint(checkpoint_path: str) -> Dict[str, Any]:
    """Load a previously saved checkpoint"""
    try:
        with open(checkpoint_path, 'r') as f:
            checkpoint = json.load(f)
            
        logger.info(f"Loaded checkpoint from {checkpoint_path}")
        logger.info(f"Checkpoint status: {len(checkpoint['processed_cids'])} processed, " 
                  f"{len(checkpoint['remaining_cids'])} remaining")
        
        return checkpoint
    except Exception as e:
        logger.error(f"Failed to load checkpoint: {e}")
        return {
            "processed_cids": [],
            "failed_cids": [],
            "remaining_cids": []
        }

def main():
    parser = argparse.ArgumentParser(description="Import cryoprotectant data from PubChem with property-based search")
    parser.add_argument("--host", required=True, help="Database host")
    parser.add_argument("--port", default="5432", help="Database port")
    parser.add_argument("--user", required=True, help="Database user")
    parser.add_argument("--password", required=True, help="Database password")
    parser.add_argument("--dbname", default="postgres", help="Database name")
    parser.add_argument("--target", type=int, default=100, help="Target number of compounds to import")
    parser.add_argument("--api-delay", type=float, default=1.0, help="Base delay between API calls (seconds)")
    parser.add_argument("--checkpoint", type=str, help="Path to checkpoint file to resume from")
    parser.add_argument("--batch-size", type=int, default=10, help="Number of compounds to process in each batch")
    parser.add_argument("--similarity-threshold", type=float, default=0.8, help="Similarity threshold for similar compounds")
    args = parser.parse_args()
    
    logger.info(f"Starting property-based PubChem import with database connection to {args.host}:{args.port}")
    
    # Initialize database connection
    db = Database(args.host, args.port, args.dbname, args.user, args.password)
    
    start_time = time.time()
    all_cids = set()
    processed_cids = set()
    failed_cids = set()
    
    # Check if we're resuming from a checkpoint
    if args.checkpoint:
        checkpoint = load_checkpoint(args.checkpoint)
        processed_cids = set(checkpoint.get("processed_cids", []))
        failed_cids = set(checkpoint.get("failed_cids", []))
        remaining_cids = checkpoint.get("remaining_cids", [])
        
        if remaining_cids:
            logger.info(f"Resuming import with {len(remaining_cids)} remaining compounds")
            all_cids = set(remaining_cids + list(processed_cids) + list(failed_cids))
    else:
        # Step 1: Collect seed compounds from known cryoprotectant names
        logger.info("Step 1: Collecting seed compounds from known cryoprotectant names")
        seed_cids = set()
        for name in CRYOPROTECTANT_NAMES:
            name_cids = get_compounds_by_name(name, delay=args.api_delay)
            logger.info(f"Found {len(name_cids)} compounds for '{name}'")
            seed_cids.update(name_cids)
        
        logger.info(f"Collected {len(seed_cids)} seed compounds from known cryoprotectants")
        all_cids.update(seed_cids)
        
        # Step 2: Collect compounds based on property ranges (individual searches)
        logger.info("Step 2: Collecting compounds based on property ranges")
        
        # Use a sample of seed compounds to search for similar property ranges
        # We use a limited number to avoid too many API calls
        property_cids = set()
        for property_name, ranges in PROPERTY_RANGES.items():
            for min_val, max_val in ranges:
                # For each property range, search for matching compounds
                range_cids = get_compounds_by_property(property_name, min_val, max_val, delay=args.api_delay)
                logger.info(f"Found {len(range_cids)} compounds with {property_name} in range {min_val}-{max_val}")
                
                # Limit the number of compounds per property search to avoid overwhelming the database
                if len(range_cids) > 100:
                    logger.info(f"Limiting to 100 random compounds from {len(range_cids)} total")
                    range_cids = random.sample(range_cids, 100)
                
                property_cids.update(range_cids)
        
        logger.info(f"Collected {len(property_cids)} compounds from property searches")
        all_cids.update(property_cids)
        
        # Step 3: Find similar compounds to known cryoprotectants
        logger.info("Step 3: Finding similar compounds to known cryoprotectants")
        similar_cids = set()
        
        # Use a limited number of CIDs to avoid too many API calls
        for cid in CRYOPROTECTANT_CIDS[:5]:  # Limit to first 5 cryoprotectants
            similar = get_similar_compounds(cid, threshold=args.similarity_threshold, max_compounds=20, delay=args.api_delay)
            similar_cids.update(similar)
            logger.info(f"Found {len(similar)} compounds similar to CID {cid}")
        
        all_cids.update(similar_cids)
        logger.info(f"Added {len(similar_cids)} similar compounds")
    
    # Convert set to list and limit to target number
    all_cids_list = list(all_cids)
    random.shuffle(all_cids_list)  # Shuffle to get a diverse set if we don't complete all
    
    if len(all_cids_list) > args.target:
        all_cids_list = all_cids_list[:args.target]
    
    logger.info(f"Processing {len(all_cids_list)} unique compounds (target: {args.target})")
    
    # Initial database counts
    initial_molecules = db.get_molecule_count()
    initial_properties = db.get_property_count()
    logger.info(f"Initial database state: {initial_molecules} molecules, {initial_properties} properties")
    
    # Get or create property types
    logger.info("Ensuring property types exist in the database")
    property_types = db.ensure_property_types()
    logger.info(f"Found {len(property_types)} property types: {list(property_types.keys())}")
    
    # Process compounds in batches
    batch_size = args.batch_size
    total_new_molecules = 0
    total_new_properties = 0
    
    for batch_start in range(0, len(all_cids_list), batch_size):
        # Create a fresh batch excluding already processed compounds
        batch_cids = [cid for cid in all_cids_list[batch_start:batch_start+batch_size] 
                     if cid not in processed_cids and cid not in failed_cids]
        
        if not batch_cids:
            continue  # Skip if all compounds in this batch were already processed
        
        logger.info(f"Processing batch {batch_start//batch_size + 1}/{(len(all_cids_list) + batch_size - 1)//batch_size} "
                   f"with {len(batch_cids)} compounds")
        
        # Process each compound in batch
        batch_molecules = []
        all_properties = []
        
        for cid in batch_cids:
            try:
                compound_data = get_compound_data(cid, delay=args.api_delay)
                
                if compound_data:
                    # Extract properties from compound_data
                    properties = compound_data.pop("properties", [])
                    batch_molecules.append(compound_data)
                    
                    processed_cids.add(cid)
                    logger.info(f"Successfully processed CID {cid}: {compound_data.get('name')}")
                else:
                    failed_cids.add(cid)
                    logger.warning(f"Failed to process CID {cid}")
            except Exception as e:
                failed_cids.add(cid)
                logger.error(f"Error processing CID {cid}: {e}")
                logger.debug(traceback.format_exc())
        
        # Insert compounds into database
        if batch_molecules:
            try:
                # Insert molecules
                molecule_map = db.batch_insert_molecules(batch_molecules)
                new_molecules = len(molecule_map)
                total_new_molecules += new_molecules
                logger.info(f"Inserted/updated {new_molecules} molecules in database")
                
                # Prepare properties with molecule IDs
                for molecule in batch_molecules:
                    cid = molecule.get("pubchem_cid")
                    if cid in molecule_map:
                        molecule_id = molecule_map[cid]
                        
                        # Get additional properties
                        additional_properties = get_property_data(cid, delay=args.api_delay)
                        
                        # Add molecule_id to each property
                        for prop in additional_properties:
                            prop["molecule_id"] = molecule_id
                            all_properties.append(prop)
                
                # Insert properties in batches with property type IDs
                if all_properties:
                    property_ids = db.batch_insert_properties(all_properties, property_types)
                    new_properties = len(property_ids)
                    total_new_properties += new_properties
                    logger.info(f"Inserted/updated {new_properties} properties in database")
            except Exception as e:
                logger.error(f"Database error in batch {batch_start//batch_size + 1}: {e}")
                logger.debug(traceback.format_exc())
        
        # Save checkpoint after each batch
        save_checkpoint(all_cids_list, processed_cids, failed_cids)
        
        # Check if we've reached or exceeded our target
        if len(processed_cids) >= args.target:
            logger.info(f"Reached target of {args.target} compounds, stopping")
            break
    
    # Calculate final statistics
    end_time = time.time()
    duration = end_time - start_time
    final_molecules = db.get_molecule_count()
    final_properties = db.get_property_count()
    
    # Generate summary
    summary = {
        "start_time": datetime.fromtimestamp(start_time).isoformat(),
        "end_time": datetime.fromtimestamp(end_time).isoformat(),
        "duration_seconds": duration,
        "total_cids_found": len(all_cids),
        "cids_processed": len(processed_cids),
        "cids_failed": len(failed_cids),
        "molecules_before": initial_molecules,
        "molecules_after": final_molecules,
        "molecules_added": final_molecules - initial_molecules,
        "properties_before": initial_properties,
        "properties_after": final_properties,
        "properties_added": final_properties - initial_properties,
        "total_new_molecules": total_new_molecules,
        "total_new_properties": total_new_properties
    }
    
    # Save summary
    summary_path = f"reports/property_pubchem_import_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    os.makedirs(os.path.dirname(summary_path), exist_ok=True)
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Import completed in {duration:.2f} seconds")
    logger.info(f"Processed {len(processed_cids)} compounds successfully, {len(failed_cids)} failed")
    logger.info(f"Added {final_molecules - initial_molecules} new molecules with {final_properties - initial_properties} properties")
    logger.info(f"Total inserted/updated: {total_new_molecules} molecules, {total_new_properties} properties")
    logger.info(f"Summary saved to {summary_path}")

if __name__ == "__main__":
    main()
"""
Simplified PubChem data import script for CryoProtect.

This script imports cryoprotectant data from PubChem based on:
1. A predefined list of known cryoprotectants
2. PubChem's classification data
3. Structural similarity to known cryoprotectants

It's designed to run in container environments with minimal dependencies.
"""

import os
import sys
import json
import time
import random
import logging
import argparse
import requests
from typing import List, Dict, Any, Optional, Set
from datetime import datetime
import traceback

# Local imports
from database.simple_db import (
    insert_molecule, insert_property, 
    batch_insert_molecules, batch_insert_properties,
    get_molecule_count, get_property_count
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'logs/pubchem_import_simple_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# PubChem API configuration
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CRYOPROTECTANT_CLASS_ID = "26001083"  # PubChem Classification for cryoprotectants
CRYOPROTECTANT_NAMES = [
    "glycerol", "ethylene glycol", "propylene glycol", "dimethyl sulfoxide", "dmso",
    "methanol", "glucose", "sucrose", "trehalose", "mannitol", "sorbitol", "xylitol",
    "dextran", "polyethylene glycol", "peg", "polyvinylpyrrolidone", "pvp",
    "hydroxyethyl starch", "histidine", "proline", "betaine", "trimethylamine n-oxide", 
    "tmao", "meso-erythritol", "erythritol", "raffinose", "glutamate", "tryptophan"
]

PUBCHEM_HEADER = {
    "User-Agent": "CryoProtect/1.0 (Database Population; mailto:cryoprotect@example.com)"
}

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

def get_compounds_by_class(class_id: str, delay: float = 0.5) -> List[int]:
    """Get compounds by PubChem classification ID"""
    class_url = f"{PUBCHEM_BASE_URL}/classification/tree/children/JSON?classification_id={class_id}&response_type=ids"
    response = fetch_with_retry(class_url)
    cids = []
    
    if response and response.status_code == 200:
        data = response.json()
        # The result includes both the compounds directly in this class and child classes
        # Extract all CIDs from the classification tree
        for ancestor in data.get("Hierarchies", {}).get("Ancestry", []):
            for node in ancestor.get("Node", []):
                if node.get("Response_Type") == "compound" and "ID" in node:
                    if isinstance(node["ID"], list):
                        cids.extend(node["ID"])
                    else:
                        cids.append(node["ID"])
        
        logger.info(f"Found {len(cids)} compounds for class ID {class_id}")
    else:
        logger.warning(f"Failed to find compounds for class ID {class_id}")
    
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
                            "source": "PubChem"
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
                            "units": None,
                            "source": "PubChem"
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
    
    filename = f"{checkpoint_dir}/pubchem_import_simple_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(filename, 'w') as f:
        json.dump(checkpoint, f, indent=2)
    
    logger.info(f"Saved checkpoint to {filename}")

def load_checkpoint(checkpoint_path: str) -> Dict[str, Any]:
    """Load a previously saved checkpoint"""
    try:
        with open(checkpoint_path, 'r') as f:
            checkpoint = json.load(f)
            
        logger.info(f"Loaded checkpoint from {checkpoint_path}")
        logger.info(f"Checkpoint status: {checkpoint['processed_cids']} processed, " 
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
    parser = argparse.ArgumentParser(description="Import cryoprotectant data from PubChem")
    parser.add_argument("--target", type=int, default=100, help="Target number of compounds to import")
    parser.add_argument("--api-delay", type=float, default=1.0, help="Base delay between API calls (seconds)")
    parser.add_argument("--checkpoint", type=str, help="Path to checkpoint file to resume from")
    parser.add_argument("--batch-size", type=int, default=10, help="Number of compounds to process in each batch")
    parser.add_argument("--workers", type=int, default=1, help="Number of worker processes (for future multiprocessing)")
    args = parser.parse_args()
    
    logger.info(f"Starting PubChem import with target: {args.target}, delay: {args.api_delay}s")
    
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
        # Step 1: Collect compounds from known cryoprotectant names
        logger.info("Step 1: Collecting compounds from known cryoprotectant names")
        for name in CRYOPROTECTANT_NAMES:
            name_cids = get_compounds_by_name(name, delay=args.api_delay)
            logger.info(f"Found {len(name_cids)} compounds for '{name}'")
            all_cids.update(name_cids)
        
        # Step 2: Collect compounds from PubChem classification
        logger.info("Step 2: Collecting compounds from PubChem classification")
        class_cids = get_compounds_by_class(CRYOPROTECTANT_CLASS_ID, delay=args.api_delay)
        logger.info(f"Found {len(class_cids)} compounds from classification")
        all_cids.update(class_cids)
    
    # Convert set to list and limit to target number
    all_cids_list = list(all_cids)
    random.shuffle(all_cids_list)  # Shuffle to get a diverse set if we don't complete all
    
    if len(all_cids_list) > args.target:
        all_cids_list = all_cids_list[:args.target]
    
    logger.info(f"Processing {len(all_cids_list)} unique compounds (target: {args.target})")
    
    # Initial database counts
    initial_molecules = get_molecule_count()
    initial_properties = get_property_count()
    logger.info(f"Initial database state: {initial_molecules} molecules, {initial_properties} properties")
    
    # Process compounds in batches
    batch_size = args.batch_size
    
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
                molecule_map = batch_insert_molecules(batch_molecules)
                logger.info(f"Inserted {len(molecule_map)} molecules to database")
                
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
                
                # Insert properties in batches
                if all_properties:
                    property_ids = batch_insert_properties(all_properties)
                    logger.info(f"Inserted {len(property_ids)} properties to database")
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
    final_molecules = get_molecule_count()
    final_properties = get_property_count()
    
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
        "properties_added": final_properties - initial_properties
    }
    
    # Save summary
    summary_path = f"reports/pubchem_import_report_simple_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    os.makedirs(os.path.dirname(summary_path), exist_ok=True)
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Import completed in {duration:.2f} seconds")
    logger.info(f"Processed {len(processed_cids)} compounds successfully, {len(failed_cids)} failed")
    logger.info(f"Added {final_molecules - initial_molecules} new molecules with {final_properties - initial_properties} properties")
    logger.info(f"Summary saved to {summary_path}")

if __name__ == "__main__":
    main()
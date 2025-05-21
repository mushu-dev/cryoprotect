#!/usr/bin/env python3
"""
expand_cid_list_classification.py

Expands the CID-Synonym-curated file to include at least 5,000 PubChem CIDs for
cryoprotectants and related compounds. This script uses the PubChem Classification
Browser to get compounds in specific classes related to cryoprotectants.

Usage:
    python expand_cid_list_classification.py --output CID-Synonym-curated --target 5000
"""

import os
import sys
import requests
import time
import argparse
import json
import random
from datetime import datetime
import logging
from pathlib import Path

# Set up logging
Path("logs").mkdir(exist_ok=True)
LOG_FILE = "logs/expand_cid_list_classification.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# PubChem API base URL
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# PubChem Classification Browser base URL
CLASSIFICATION_BASE = "https://pubchem.ncbi.nlm.nih.gov/classification/cgi/classifications.fcgi"

# Classification IDs for cryoprotectant-related classes
# These are the classification IDs for the PubChem Classification Browser
CLASSIFICATION_IDS = {
    # Chemical classes
    "Alcohols": "LCALCOHOLS",
    "Polyols": "LCPOLYOLS",
    "Glycols": "LCGLYCOLS",
    "Sugars": "LCMONOSACCHARIDES",
    "Disaccharides": "LCDISACCHARIDES",
    "Amides": "LCAMIDES",
    "Amines": "LCAMINES",
    "Ethers": "LCETHERS",
    "Esters": "LCESTERS",
    "Sulfoxides": "LCSULFOXIDES",
    "Amino Acids": "LCAMINOACIDS",
    "Peptides": "LCPEPTIDES",
    "Glycerides": "LCGLYCERIDES",
    "Phospholipids": "LCPHOSPHOLIPIDS",
    "Glycolipids": "LCGLYCOLIPIDS",
    
    # Functional groups
    "Hydroxyl": "LCHYDROXYL",
    "Carbonyl": "LCCARBONYL",
    "Carboxyl": "LCCARBOXYL",
    "Amine": "LCAMINE",
    "Amide": "LCAMIDE",
    "Ether": "LCETHER",
    "Ester": "LCESTER",
    "Sulfoxide": "LCSULFOXIDE",
    "Phosphate": "LCPHOSPHATE",
    "Glycosidic": "LCGLYCOSIDIC",
    
    # Biological roles
    "Cryoprotectants": "LCBR_CRYOPROTECTANTS",
    "Antifreeze": "LCBR_ANTIFREEZE",
    "Cell Preservatives": "LCBR_CELLPRESERVATIVES",
    "Tissue Preservatives": "LCBR_TISSUEPRESERVATIVES",
    "Organ Preservatives": "LCBR_ORGANPRESERVATIVES",
}

# Substructure SMARTS patterns for cryoprotectant-related functional groups
SUBSTRUCTURE_SMARTS = [
    # Alcohols
    "CO",  # Methanol
    "CCO",  # Ethanol
    "CCCO",  # Propanol
    "CCCCO",  # Butanol
    
    # Diols
    "OCCO",  # Ethylene glycol
    "OCCCO",  # Propylene glycol
    "OCCCCO",  # Butanediol
    
    # Triols
    "OCC(O)CO",  # Glycerol
    
    # Polyols
    "OCC(O)C(O)CO",  # Erythritol
    "OCC(O)C(O)C(O)CO",  # Xylitol, Arabitol
    "OCC(O)C(O)C(O)C(O)CO",  # Mannitol, Sorbitol
    
    # Sugars
    "OCC1OC(O)C(O)C(O)C1O",  # Glucose (pyranose form)
    "OCC1OC(O)C(O)C1O",  # Ribose (furanose form)
    
    # Disaccharides
    "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O",  # Sucrose-like
    
    # Amides
    "NC=O",  # Formamide
    "CNC=O",  # N-methylformamide
    "CC(=O)N",  # Acetamide
    
    # Sulfoxides
    "CS(=O)C",  # Dimethyl sulfoxide
    
    # Ethers
    "COC",  # Dimethyl ether
    "CCOC",  # Ethyl methyl ether
    "CCOCC",  # Diethyl ether
    
    # Glycol ethers
    "COCCOC",  # Diethylene glycol dimethyl ether
    "COCCO",  # Ethylene glycol monomethyl ether
    
    # Amino acids
    "NCC(=O)O",  # Glycine
    "NC(C)C(=O)O",  # Alanine
    "NC(CC(=O)O)C(=O)O",  # Aspartic acid
]

def fetch_cids_by_classification(classification_id, delay=0.2, max_retries=3):
    """Fetch CIDs for compounds in a specific classification."""
    url = f"{CLASSIFICATION_BASE}?format=json&hid={classification_id}&classification=chemical&count=5000"
    
    for retry in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code == 200:
                data = r.json()
                if "Hierarchies" in data and len(data["Hierarchies"]) > 0:
                    hierarchy = data["Hierarchies"][0]
                    if "Information" in hierarchy and "CID" in hierarchy["Information"]:
                        cids = hierarchy["Information"]["CID"]
                        logger.info(f"Found {len(cids)} CIDs for classification '{classification_id}'")
                        time.sleep(delay)  # Be polite to PubChem API
                        return set(str(cid) for cid in cids)
                logger.warning(f"No CIDs found in response for classification '{classification_id}'")
                time.sleep(delay)  # Be polite to PubChem API
                return set()
            elif r.status_code == 429:  # Rate limit
                wait_time = delay * (2 ** retry)
                logger.warning(f"Rate limit hit for classification '{classification_id}'. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                logger.warning(f"Error fetching CIDs for classification '{classification_id}'. Status code: {r.status_code}")
                time.sleep(delay * (2 ** retry))
        except Exception as e:
            logger.warning(f"Exception fetching CIDs for classification '{classification_id}': {str(e)}")
            time.sleep(delay * (2 ** retry))
    
    logger.error(f"Failed to fetch CIDs for classification '{classification_id}' after {max_retries} retries")
    return set()

def fetch_cids_by_substructure(smarts, delay=0.2, max_retries=3):
    """Fetch CIDs for compounds containing a specific substructure."""
    url = f"{PUBCHEM_BASE}/compound/substructure/smarts/{requests.utils.quote(smarts)}/cids/TXT"
    
    for retry in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code == 200:
                cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
                logger.info(f"Found {len(cids)} CIDs for substructure '{smarts}'")
                time.sleep(delay)  # Be polite to PubChem API
                return set(cids)
            elif r.status_code == 404:
                logger.warning(f"No results for substructure '{smarts}'. Status code: 404")
                time.sleep(delay)  # Be polite to PubChem API
                return set()
            elif r.status_code == 429:  # Rate limit
                wait_time = delay * (2 ** retry)
                logger.warning(f"Rate limit hit for substructure '{smarts}'. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                logger.warning(f"Error fetching CIDs for substructure '{smarts}'. Status code: {r.status_code}")
                time.sleep(delay * (2 ** retry))
        except Exception as e:
            logger.warning(f"Exception fetching CIDs for substructure '{smarts}': {str(e)}")
            time.sleep(delay * (2 ** retry))
    
    logger.error(f"Failed to fetch CIDs for substructure '{smarts}' after {max_retries} retries")
    return set()

def fetch_names_for_cids(cids, batch_size=100, delay=0.2, max_retries=3):
    """Fetch preferred names for a list of CIDs using PubChem."""
    cid_list = list(cids)
    cid_to_name = {}
    
    # Simple progress tracking
    total_batches = (len(cid_list) + batch_size - 1) // batch_size
    logger.info(f"Fetching names for {len(cid_list)} CIDs in {total_batches} batches")
    
    for i in range(0, len(cid_list), batch_size):
        batch_num = i // batch_size + 1
        batch = cid_list[i:i+batch_size]
        logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} CIDs)")
        
        ids = ",".join(str(cid) for cid in batch)
        url = f"{PUBCHEM_BASE}/compound/cid/{ids}/property/Title/TXT"
        
        for retry in range(max_retries):
            try:
                r = requests.get(url)
                if r.status_code == 200:
                    lines = r.text.splitlines()
                    if len(lines) == len(batch):
                        for j, line in enumerate(lines):
                            if line.strip():
                                cid_to_name[str(batch[j])] = line.strip()
                    else:
                        logger.warning(f"Mismatch in response length for batch {batch_num}. Expected {len(batch)}, got {len(lines)}")
                    break
                elif r.status_code == 429:  # Rate limit
                    wait_time = delay * (2 ** retry)
                    logger.warning(f"Rate limit hit for batch {batch_num}. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                    time.sleep(wait_time)
                else:
                    logger.warning(f"Error fetching names for batch {batch_num}. Status code: {r.status_code}")
                    time.sleep(delay * (2 ** retry))
            except Exception as e:
                logger.warning(f"Exception fetching names for batch {batch_num}: {str(e)}")
                time.sleep(delay * (2 ** retry))
        
        time.sleep(delay)  # Be polite to PubChem API
    
    logger.info(f"Retrieved names for {len(cid_to_name)} CIDs out of {len(cids)} requested")
    return cid_to_name

def load_existing_cids(file_path):
    """Load existing CIDs from a file."""
    if not os.path.exists(file_path):
        logger.warning(f"File {file_path} does not exist.")
        return set()
        
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            cids = [line.strip().split("\t")[0] for line in f if line.strip() and "\t" in line]
            return set(cids)
    except Exception as e:
        logger.warning(f"Error loading existing CIDs from {file_path}: {str(e)}")
        return set()

def save_cids_to_file(cid_to_name, output_path, mode="w"):
    """Save CIDs and names to a file."""
    try:
        with open(output_path, mode, encoding="utf-8") as f:
            for cid, name in cid_to_name.items():
                if name:
                    f.write(f"{cid}\t{name}\n")
        logger.info(f"Saved {len(cid_to_name)} CIDs to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error saving CIDs to {output_path}: {str(e)}")
        return False

def filter_cids_by_properties(cids, api_delay=0.2, max_retries=3, batch_size=100):
    """Filter CIDs based on molecular properties."""
    filtered_cids = set()
    cid_list = list(cids)
    
    # Core Filtering Criteria from import_pubchem_data_mcp.py
    CORE_CRITERIA = {
        "logP_range": (-5, 5),
        "mw_range": (0, 1000),
        "TPSA_range": (0, 200)
    }
    
    # Simple progress tracking
    total_batches = (len(cid_list) + batch_size - 1) // batch_size
    logger.info(f"Filtering {len(cid_list)} CIDs in {total_batches} batches")
    
    for i in range(0, len(cid_list), batch_size):
        batch_num = i // batch_size + 1
        batch = cid_list[i:i+batch_size]
        logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} CIDs)")
        
        ids = ",".join(str(cid) for cid in batch)
        url = f"{PUBCHEM_BASE}/compound/cid/{ids}/property/MolecularWeight,XLogP,TPSA/JSON"
        
        for retry in range(max_retries):
            try:
                r = requests.get(url)
                if r.status_code == 200:
                    data = r.json()
                    if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                        properties = data["PropertyTable"]["Properties"]
                        for prop in properties:
                            cid = str(prop["CID"])
                            mw = prop.get("MolecularWeight")
                            logp = prop.get("XLogP")
                            tpsa = prop.get("TPSA")
                            
                            # Apply filtering criteria
                            if (mw is not None and CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1] and
                                logp is not None and CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1] and
                                tpsa is not None and CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
                                filtered_cids.add(cid)
                    break
                elif r.status_code == 429:  # Rate limit
                    wait_time = delay * (2 ** retry)
                    logger.warning(f"Rate limit hit for batch {batch_num}. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                    time.sleep(wait_time)
                else:
                    logger.warning(f"Error fetching properties for batch {batch_num}. Status code: {r.status_code}")
                    time.sleep(delay * (2 ** retry))
            except Exception as e:
                logger.warning(f"Exception fetching properties for batch {batch_num}: {str(e)}")
                time.sleep(delay * (2 ** retry))
        
        time.sleep(api_delay)  # Be polite to PubChem API
    
    logger.info(f"Filtered {len(filtered_cids)} CIDs out of {len(cids)} based on properties")
    return filtered_cids

def main(output_path, target_count=5000, api_delay=0.2):
    """Main function to expand the CID list."""
    start_time = datetime.now()
    logger.info(f"Starting CID list expansion at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Target: {target_count} CIDs")
    
    # Load existing CIDs if the file exists
    existing_cids = load_existing_cids(output_path)
    logger.info(f"Loaded {len(existing_cids)} existing CIDs from {output_path}")
    
    all_cids = set(existing_cids)
    
    # Track the number of CIDs found by each method
    stats = {
        "existing": len(existing_cids),
        "classification": 0,
        "substructure": 0,
        "total": len(existing_cids)
    }
    
    # 1. Fetch CIDs by classification
    logger.info("Fetching CIDs by classification...")
    classification_cids = set()
    
    # Shuffle the classification IDs to get a more diverse set of compounds
    classification_items = list(CLASSIFICATION_IDS.items())
    random.shuffle(classification_items)
    
    for idx, (class_name, class_id) in enumerate(classification_items):
        logger.info(f"Processing classification {idx+1}/{len(classification_items)}: {class_name} ({class_id})")
        new_cids = fetch_cids_by_classification(class_id, api_delay)
        classification_cids.update(new_cids)
        
        # Check if we've reached the target
        if len(classification_cids) + len(all_cids) >= target_count:
            logger.info(f"Found enough CIDs via classification. Stopping search.")
            break
    
    # Add new CIDs to the set
    new_classification_cids = classification_cids - all_cids
    all_cids.update(new_classification_cids)
    stats["classification"] = len(new_classification_cids)
    stats["total"] = len(all_cids)
    logger.info(f"Found {len(new_classification_cids)} new CIDs via classification. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 2. Fetch CIDs by substructure
        logger.info("Fetching CIDs by substructure...")
        substructure_cids = set()
        
        # Shuffle the substructure SMARTS to get a more diverse set of compounds
        smarts_patterns = list(SUBSTRUCTURE_SMARTS)
        random.shuffle(smarts_patterns)
        
        for idx, smarts in enumerate(smarts_patterns):
            logger.info(f"Processing substructure {idx+1}/{len(smarts_patterns)}: {smarts}")
            new_cids = fetch_cids_by_substructure(smarts, api_delay)
            substructure_cids.update(new_cids)
            
            # Check if we've reached the target
            if len(substructure_cids) + len(all_cids) >= target_count:
                logger.info(f"Found enough CIDs via substructure. Stopping search.")
                break
        
        # Add new CIDs to the set
        new_substructure_cids = substructure_cids - all_cids
        all_cids.update(new_substructure_cids)
        stats["substructure"] = len(new_substructure_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_substructure_cids)} new CIDs via substructure. Total unique CIDs: {len(all_cids)}")
    
    # Filter CIDs based on properties
    logger.info(f"Filtering {len(all_cids)} CIDs based on properties...")
    filtered_cids = filter_cids_by_properties(all_cids, api_delay)
    logger.info(f"Filtered to {len(filtered_cids)} CIDs based on properties")
    
    # Fetch names for all CIDs
    logger.info(f"Fetching names for {len(filtered_cids)} CIDs...")
    cid_to_name = fetch_names_for_cids(filtered_cids, batch_size=100, delay=api_delay)
    logger.info(f"Retrieved names for {len(cid_to_name)} CIDs.")
    
    # Save to output file
    if save_cids_to_file(cid_to_name, output_path):
        logger.info(f"Successfully saved {len(cid_to_name)} CIDs to {output_path}")
    else:
        logger.error(f"Failed to save CIDs to {output_path}")
    
    # Generate a report
    end_time = datetime.now()
    duration = end_time - start_time
    
    report = {
        "timestamp": end_time.strftime("%Y-%m-%d %H:%M:%S"),
        "duration": str(duration),
        "target_count": target_count,
        "final_count": len(cid_to_name),
        "stats": stats,
        "output_file": output_path
    }
    
    report_file = f"reports/cid_expansion_report_{end_time.strftime('%Y%m%d_%H%M%S')}.json"
    Path("reports").mkdir(exist_ok=True)
    
    with open(report_file, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Report saved to {report_file}")
    logger.info(f"CID list expansion completed in {duration}")
    
    return len(cid_to_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Expand CID list for cryoprotectants using PubChem Classification Browser.")
    parser.add_argument("--output", type=str, default="CID-Synonym-curated", help="Output file path")
    parser.add_argument("--target", type=int, default=5000, help="Target number of CIDs")
    parser.add_argument("--delay", type=float, default=0.2, help="Delay between PubChem API calls in seconds")
    args = parser.parse_args()
    
    main(args.output, args.target, args.delay)
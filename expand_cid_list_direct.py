#!/usr/bin/env python3
"""
expand_cid_list_direct.py

A more direct approach to expand the CID-Synonym-curated file to include at least 5,000 
PubChem CIDs for cryoprotectants and related compounds. This script uses direct compound 
searches with the PubChem PUG REST API and focuses on specific properties relevant to 
cryoprotectants.

Usage:
    python expand_cid_list_direct.py --output CID-Synonym-curated --target 5000
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
LOG_FILE = "logs/expand_cid_list_direct.log"
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

# Direct search terms for cryoprotectants and related compounds
SEARCH_TERMS = [
    # Common cryoprotectants
    "dimethyl sulfoxide", "glycerol", "ethylene glycol", "propylene glycol", 
    "trehalose", "sucrose", "glucose", "fructose", "mannitol", "sorbitol",
    "methanol", "formamide", "acetamide", "propanediol", "butanediol",
    "erythritol", "xylitol", "inositol", "arabitol", "ribitol", "galactitol",
    
    # Chemical classes
    "polyol", "sugar alcohol", "disaccharide", "monosaccharide", "glycol",
    "diol", "triol", "tetrol", "pentol", "hexol", "amide", "sulfoxide",
    
    # Functional groups
    "hydroxyl", "alcohol", "ether", "glycoside", "glucoside", "pyranoside",
    "furanoside", "glycerol ether", "glycerol ester",
    
    # Properties
    "cryoprotective", "antifreeze", "freeze resistant", "ice inhibitor",
    "glass former", "vitrification agent",
    
    # Related compounds
    "polyethylene glycol", "hydroxyethyl starch", "dextran", "ficoll",
    "polyvinylpyrrolidone", "albumin", "proline", "betaine", "glycine",
    "alanine", "serine", "threonine", "glutamine", "asparagine",
    
    # Specific cryoprotectant derivatives
    "glycerol derivative", "ethylene glycol derivative", "propylene glycol derivative",
    "dimethyl sulfoxide derivative", "trehalose derivative", "sucrose derivative",
    
    # Expanded alcohol series
    "propanol", "butanol", "pentanol", "hexanol", "heptanol", "octanol",
    
    # Expanded glycol series
    "diethylene glycol", "triethylene glycol", "tetraethylene glycol",
    "dipropylene glycol", "tripropylene glycol", "tetrapropylene glycol",
    
    # Expanded sugar series
    "maltose", "lactose", "cellobiose", "maltotriose", "maltotetraose",
    "raffinose", "stachyose", "verbascose", "xylose", "arabinose", "ribose",
    "galactose", "mannose", "rhamnose", "fucose",
    
    # Expanded polyol series
    "glycerol", "erythritol", "threitol", "arabitol", "xylitol", "sorbitol",
    "mannitol", "dulcitol", "inositol", "maltitol", "lactitol", "isomalt",
    
    # Expanded amide series
    "formamide", "acetamide", "propionamide", "butyramide", "valeramide",
    "hexanamide", "heptanamide", "octanamide",
    
    # Expanded ether series
    "dimethyl ether", "diethyl ether", "dipropyl ether", "dibutyl ether",
    "methyl ethyl ether", "methyl propyl ether", "ethyl propyl ether",
    
    # Expanded ester series
    "methyl acetate", "ethyl acetate", "propyl acetate", "butyl acetate",
    "methyl propionate", "ethyl propionate", "propyl propionate",
    
    # Expanded amino acid series
    "glycine", "alanine", "valine", "leucine", "isoleucine", "proline",
    "phenylalanine", "tyrosine", "tryptophan", "serine", "threonine",
    "cysteine", "methionine", "asparagine", "glutamine", "aspartic acid",
    "glutamic acid", "lysine", "arginine", "histidine",
    
    # Expanded peptide series
    "dipeptide", "tripeptide", "tetrapeptide", "pentapeptide", "hexapeptide",
    "heptapeptide", "octapeptide", "nonapeptide", "decapeptide",
    
    # Expanded polymer series
    "polyethylene glycol", "polypropylene glycol", "polyvinyl alcohol",
    "polyvinylpyrrolidone", "polyacrylamide", "polyethylenimine",
    
    # Expanded carbohydrate series
    "glucose", "fructose", "galactose", "mannose", "xylose", "arabinose",
    "ribose", "sucrose", "lactose", "maltose", "trehalose", "raffinose",
    "stachyose", "verbascose", "cellulose", "starch", "glycogen", "inulin",
    
    # Expanded lipid series
    "phospholipid", "glycolipid", "sphingolipid", "ceramide", "ganglioside",
    "cerebroside", "sulfatide", "glycerophospholipid", "phosphatidylcholine",
    "phosphatidylethanolamine", "phosphatidylserine", "phosphatidylinositol",
    
    # Expanded nucleoside/nucleotide series
    "adenosine", "guanosine", "cytidine", "uridine", "thymidine",
    "inosine", "xanthosine", "adenosine monophosphate", "guanosine monophosphate",
    "cytidine monophosphate", "uridine monophosphate", "thymidine monophosphate",
]

# Property ranges for cryoprotectants
PROPERTY_RANGES = [
    # LogP range (-5 to 5)
    "xlogp:-5:5",
    
    # Molecular weight range (0 to 1000)
    "mw:0:1000",
    
    # Hydrogen bond donors (0 to 20)
    "hbonddonorcount:0:20",
    
    # Hydrogen bond acceptors (0 to 20)
    "hbondacceptorcount:0:20",
    
    # Rotatable bonds (0 to 20)
    "rotatablebondcount:0:20",
    
    # Topological polar surface area (0 to 200)
    "tpsa:0:200",
    
    # Heavy atom count (0 to 100)
    "heavyatomcount:0:100",
]

# Functional groups common in cryoprotectants
FUNCTIONAL_GROUPS = [
    "hydroxyl", "ether", "ester", "amide", "amine", "carboxyl", 
    "carbonyl", "sulfoxide", "sulfone", "phosphate", "phosphonate",
    "glycosidic", "acetal", "ketal", "hemiacetal", "hemiketal",
]

def fetch_cids_by_name(name, delay=0.2, max_retries=3):
    """Fetch CIDs for compounds matching a name."""
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/TXT"
    
    for retry in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code == 200:
                cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
                logger.info(f"Found {len(cids)} CIDs for name '{name}'")
                time.sleep(delay)  # Be polite to PubChem API
                return set(cids)
            elif r.status_code == 404:
                logger.warning(f"No results for name '{name}'. Status code: 404")
                time.sleep(delay)  # Be polite to PubChem API
                return set()
            elif r.status_code == 429:  # Rate limit
                wait_time = delay * (2 ** retry)
                logger.warning(f"Rate limit hit for name '{name}'. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                logger.warning(f"Error fetching CIDs for name '{name}'. Status code: {r.status_code}")
                time.sleep(delay * (2 ** retry))
        except Exception as e:
            logger.warning(f"Exception fetching CIDs for name '{name}': {str(e)}")
            time.sleep(delay * (2 ** retry))
    
    logger.error(f"Failed to fetch CIDs for name '{name}' after {max_retries} retries")
    return set()

def fetch_cids_by_property(property_query, delay=0.2, max_retries=3):
    """Fetch CIDs for compounds matching a property query."""
    url = f"{PUBCHEM_BASE}/compound/fastsubstructure/{property_query}/cids/TXT"
    
    for retry in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code == 200:
                cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
                logger.info(f"Found {len(cids)} CIDs for property query '{property_query}'")
                time.sleep(delay)  # Be polite to PubChem API
                return set(cids)
            elif r.status_code == 404:
                logger.warning(f"No results for property query '{property_query}'. Status code: 404")
                time.sleep(delay)  # Be polite to PubChem API
                return set()
            elif r.status_code == 429:  # Rate limit
                wait_time = delay * (2 ** retry)
                logger.warning(f"Rate limit hit for property query '{property_query}'. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                logger.warning(f"Error fetching CIDs for property query '{property_query}'. Status code: {r.status_code}")
                time.sleep(delay * (2 ** retry))
        except Exception as e:
            logger.warning(f"Exception fetching CIDs for property query '{property_query}': {str(e)}")
            time.sleep(delay * (2 ** retry))
    
    logger.error(f"Failed to fetch CIDs for property query '{property_query}' after {max_retries} retries")
    return set()

def fetch_cids_by_smarts(smarts, delay=0.2, max_retries=3):
    """Fetch CIDs for compounds matching a SMARTS pattern."""
    url = f"{PUBCHEM_BASE}/compound/substructure/smarts/{requests.utils.quote(smarts)}/cids/TXT"
    
    for retry in range(max_retries):
        try:
            r = requests.get(url)
            if r.status_code == 200:
                cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
                logger.info(f"Found {len(cids)} CIDs for SMARTS pattern '{smarts}'")
                time.sleep(delay)  # Be polite to PubChem API
                return set(cids)
            elif r.status_code == 404:
                logger.warning(f"No results for SMARTS pattern '{smarts}'. Status code: 404")
                time.sleep(delay)  # Be polite to PubChem API
                return set()
            elif r.status_code == 429:  # Rate limit
                wait_time = delay * (2 ** retry)
                logger.warning(f"Rate limit hit for SMARTS pattern '{smarts}'. Waiting {wait_time:.2f}s before retry {retry+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                logger.warning(f"Error fetching CIDs for SMARTS pattern '{smarts}'. Status code: {r.status_code}")
                time.sleep(delay * (2 ** retry))
        except Exception as e:
            logger.warning(f"Exception fetching CIDs for SMARTS pattern '{smarts}': {str(e)}")
            time.sleep(delay * (2 ** retry))
    
    logger.error(f"Failed to fetch CIDs for SMARTS pattern '{smarts}' after {max_retries} retries")
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
        "search_terms": 0,
        "property_ranges": 0,
        "functional_groups": 0,
        "total": len(existing_cids)
    }
    
    # 1. Fetch CIDs by search terms
    logger.info("Fetching CIDs by search terms...")
    search_term_cids = set()
    
    # Shuffle the search terms to get a more diverse set of compounds
    random.shuffle(SEARCH_TERMS)
    
    for idx, term in enumerate(SEARCH_TERMS):
        logger.info(f"Processing search term {idx+1}/{len(SEARCH_TERMS)}: {term}")
        new_cids = fetch_cids_by_name(term, api_delay)
        search_term_cids.update(new_cids)
        
        # Check if we've reached the target
        if len(search_term_cids) + len(all_cids) >= target_count:
            logger.info(f"Found enough CIDs via search terms. Stopping search.")
            break
    
    # Add new CIDs to the set
    new_search_term_cids = search_term_cids - all_cids
    all_cids.update(new_search_term_cids)
    stats["search_terms"] = len(new_search_term_cids)
    stats["total"] = len(all_cids)
    logger.info(f"Found {len(new_search_term_cids)} new CIDs via search terms. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 2. Fetch CIDs by property ranges
        logger.info("Fetching CIDs by property ranges...")
        property_range_cids = set()
        
        for idx, prop_range in enumerate(PROPERTY_RANGES):
            logger.info(f"Processing property range {idx+1}/{len(PROPERTY_RANGES)}: {prop_range}")
            new_cids = fetch_cids_by_property(prop_range, api_delay)
            property_range_cids.update(new_cids)
            
            # Check if we've reached the target
            if len(property_range_cids) + len(all_cids) >= target_count:
                logger.info(f"Found enough CIDs via property ranges. Stopping search.")
                break
        
        # Add new CIDs to the set
        new_property_range_cids = property_range_cids - all_cids
        all_cids.update(new_property_range_cids)
        stats["property_ranges"] = len(new_property_range_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_property_range_cids)} new CIDs via property ranges. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 3. Fetch CIDs by functional groups
        logger.info("Fetching CIDs by functional groups...")
        functional_group_cids = set()
        
        for idx, group in enumerate(FUNCTIONAL_GROUPS):
            logger.info(f"Processing functional group {idx+1}/{len(FUNCTIONAL_GROUPS)}: {group}")
            new_cids = fetch_cids_by_name(group, api_delay)
            functional_group_cids.update(new_cids)
            
            # Check if we've reached the target
            if len(functional_group_cids) + len(all_cids) >= target_count:
                logger.info(f"Found enough CIDs via functional groups. Stopping search.")
                break
        
        # Add new CIDs to the set
        new_functional_group_cids = functional_group_cids - all_cids
        all_cids.update(new_functional_group_cids)
        stats["functional_groups"] = len(new_functional_group_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_functional_group_cids)} new CIDs via functional groups. Total unique CIDs: {len(all_cids)}")
    
    # Fetch names for all CIDs
    logger.info(f"Fetching names for {len(all_cids)} CIDs...")
    cid_to_name = fetch_names_for_cids(all_cids, batch_size=100, delay=api_delay)
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
    parser = argparse.ArgumentParser(description="Expand CID list for cryoprotectants using direct PubChem searches.")
    parser.add_argument("--output", type=str, default="CID-Synonym-curated", help="Output file path")
    parser.add_argument("--target", type=int, default=5000, help="Target number of CIDs")
    parser.add_argument("--delay", type=float, default=0.2, help="Delay between PubChem API calls in seconds")
    args = parser.parse_args()
    
    main(args.output, args.target, args.delay)
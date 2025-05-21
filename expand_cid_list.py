#!/usr/bin/env python3
"""
expand_cid_list.py

Expands the CID-Synonym-curated file to include at least 5,000 PubChem CIDs for
cryoprotectants and related compounds. This script uses multiple strategies to find
relevant compounds:

1. MeSH terms related to cryoprotection and cryobiology
2. Keywords related to cryoprotectants and their properties
3. Similar compounds to known cryoprotectants (based on structure)
4. Compounds with similar properties to known cryoprotectants
5. Compounds in relevant chemical classes

The script outputs a tab-separated file: CID<TAB>Name, one per line, no header,
which can be used with import_pubchem_data_mcp.py to import the data into the database.

Requirements:
- Python 3.6+
- requests

Usage:
    python expand_cid_list.py --output CID-Synonym-curated --target 5000
"""

import os
import sys
import requests
import time
import argparse
import json
from datetime import datetime
import logging
from pathlib import Path

# Set up logging
Path("logs").mkdir(exist_ok=True)
LOG_FILE = "logs/expand_cid_list.log"
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

# MeSH terms related to cryoprotection and cryobiology
MESH_TERMS = [
    "D003440",  # Cryoprotective Agents
    "D003452",  # Cryobiology
    "D003453",  # Cryopreservation
    "D003457",  # Cryosurgery
    "D003463",  # Cryotherapy
    "D050051",  # Vitrification
    "D004041",  # Dimethyl Sulfoxide
    "D005930",  # Glycerol
    "D004791",  # Ethylene Glycol
    "D011498",  # Propylene Glycol
    "D014292",  # Trehalose
    "D013395",  # Sucrose
    "D008670",  # Methanol
    "D005557",  # Formamide
    "D000105",  # Acetamide
    "D011433",  # Propanediols
    "D008358",  # Mannitol
    "D013006",  # Sorbitol
    "D005947",  # Glucose
    "D005632",  # Fructose
    "D011725",  # Raffinose
    "D014948",  # Xylose
    "D004912",  # Erythritol
    "D005561",  # Formaldehyde
    "D001091",  # Antifreeze Proteins
    "D000262",  # Acetonitrile
    "D000096",  # Acetone
    "D002245",  # Butanediols
    "D007328",  # Isopropanol
    "D000431",  # Alcohols
    "D004040",  # Dimethylformamide
    "D013845",  # Thiourea
    "D013955",  # Thiodiglycol
    "D013997",  # Thiocyanates
    "D001323",  # Amides
    "D001564",  # Antioxidants
    "D004785",  # Ethers
    "D005680",  # Freezing
    "D005260",  # Ficoll
    "D011108",  # Polyvinylpyrrolidone
    "D010743",  # Polyethylene Glycols
    "D006820",  # Hydroxyethyl Starch
    "D003847",  # Dextrans
    "D000419",  # Albumins
    "D012770",  # Serum
]

# Keywords related to cryoprotectants and their properties
KEYWORDS = [
    "cryoprotectant",
    "cryoprotection",
    "cryopreservation",
    "antifreeze",
    "vitrification",
    "freeze protection",
    "ice inhibitor",
    "ice nucleation",
    "ice recrystallization",
    "glass transition",
    "glass former",
    "cell preservation",
    "tissue preservation",
    "organ preservation",
    "freezing point depression",
    "supercooling",
    "cold protection",
    "frost protection",
    "freeze-thaw",
    "cold tolerance",
    "cold acclimation",
    "cold hardiness",
    "cold resistance",
    "cold stress",
    "cold shock",
    "cold adaptation",
    "cold injury",
    "cold damage",
    "cold survival",
    "cold storage",
    "cold preservation",
    "cold stabilization",
    "cold protection",
    "cold hardening",
    "cold acclimatization",
    "cold tolerance",
    "cold resistance",
    "cold adaptation",
    "cold survival",
    "cold protection",
    "cold hardening",
    "cold acclimatization",
]

# Chemical classes related to cryoprotectants
CHEMICAL_CLASSES = [
    "polyols",
    "sugars",
    "sugar alcohols",
    "amino acids",
    "amides",
    "sulfoxides",
    "glycols",
    "diols",
    "triols",
    "alcohols",
    "glycerol ethers",
    "glycerol derivatives",
    "ethylene glycol derivatives",
    "propylene glycol derivatives",
    "dimethyl sulfoxide derivatives",
    "formamide derivatives",
    "acetamide derivatives",
    "trehalose derivatives",
    "sucrose derivatives",
    "glucose derivatives",
    "fructose derivatives",
    "mannitol derivatives",
    "sorbitol derivatives",
    "glycerol esters",
    "glycerol ethers",
]

# Core cryoprotectant CIDs (from the current CID-Synonym-curated file)
CORE_CRYOPROTECTANT_CIDS = [
    679,    # DMSO
    753,    # glycerol
    174,    # ethylene glycol
    1030,   # 1,2-propanediol
    7427,   # trehalose
    5988,   # sucrose
    887,    # methanol
    713,    # formamide
    178,    # acetamide
    10442,  # 1,3-propanediol
    6251,   # mannitol
    5780,   # sorbitol
    5793,   # glucose
    2723872,# fructose
    439242, # raffinose
    135191, # xylose
    222285, # erythritol
    302428, # EG
    1031,   # PROH
    4125253,# dextran
    6342,   # acetonitrile
    180,    # acetone
    165788, # butanediol
    3776,   # isopropanol
    702,    # ethanol
    8117,   # diethylene glycol
    8172,   # triethylene glycol
    8200,   # tetraethylene glycol
]

def fetch_cids_by_mesh(mesh_id, delay=0.2):
    """Fetch CIDs for compounds annotated with a given MeSH term."""
    url = f"{PUBCHEM_BASE}/compound/mesh/{mesh_id}/cids/TXT"
    try:
        r = requests.get(url)
        r.raise_for_status()
        cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
        time.sleep(delay)  # Be polite to PubChem API
        return set(cids)
    except Exception as e:
        logger.warning(f"Error fetching CIDs for MeSH term {mesh_id}: {str(e)}")
        time.sleep(delay * 2)  # Wait longer after an error
        return set()

def fetch_cids_by_keyword(keyword, delay=0.2):
    """Fetch CIDs for compounds matching a keyword in synonyms/titles."""
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(keyword)}/cids/TXT"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
            time.sleep(delay)  # Be polite to PubChem API
            return set(cids)
        else:
            logger.warning(f"No results for keyword '{keyword}'. Status code: {r.status_code}")
            time.sleep(delay * 2)  # Wait longer after an error
    except Exception as e:
        logger.warning(f"Error fetching CIDs for keyword '{keyword}': {str(e)}")
        time.sleep(delay * 2)  # Wait longer after an error
    return set()

def fetch_cids_by_chemical_class(chemical_class, delay=0.2):
    """Fetch CIDs for compounds in a specific chemical class."""
    return fetch_cids_by_keyword(chemical_class, delay)

def fetch_similar_compounds(cid, similarity=90, delay=0.2):
    """Fetch CIDs for compounds similar to a given CID based on structure."""
    url = f"{PUBCHEM_BASE}/compound/similarity/cid/{cid}/cids/TXT?Threshold={similarity}"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
            time.sleep(delay)  # Be polite to PubChem API
            return set(cids)
        else:
            logger.warning(f"No similar compounds found for CID {cid}. Status code: {r.status_code}")
            time.sleep(delay * 2)  # Wait longer after an error
    except Exception as e:
        logger.warning(f"Error fetching similar compounds for CID {cid}: {str(e)}")
        time.sleep(delay * 2)  # Wait longer after an error
    return set()

def fetch_names_for_cids(cids, batch_size=100, delay=0.2):
    """Fetch preferred names for a list of CIDs using PubChem."""
    cid_list = list(cids)
    cid_to_name = {}
    
    # Simple progress tracking
    total_batches = (len(cid_list) + batch_size - 1) // batch_size
    logger.info(f"Fetching names for {len(cid_list)} CIDs in {total_batches} batches")
    
    for i in range(0, len(cid_list), batch_size):
        batch_num = i // batch_size + 1
        logger.info(f"Processing batch {batch_num}/{total_batches} ({len(cid_list[i:i+batch_size])} CIDs)")
        batch = cid_list[i:i+batch_size]
        ids = ",".join(batch)
        url = f"{PUBCHEM_BASE}/compound/cid/{ids}/property/IUPACName,Title,CanonicalSMILES/TSV"
        
        try:
            r = requests.get(url)
            if r.status_code != 200:
                logger.warning(f"Error fetching names for batch {i//batch_size + 1}. Status code: {r.status_code}")
                time.sleep(delay * 2)  # Wait longer after an error
                continue
                
            lines = r.text.splitlines()
            if not lines:
                logger.warning(f"Empty response for batch {i//batch_size + 1}")
                time.sleep(delay)
                continue
                
            header = lines[0].split('\t')
            idx_cid = header.index("CID")
            idx_title = header.index("Title") if "Title" in header else -1
            idx_iupac = header.index("IUPACName") if "IUPACName" in header else -1
            
            for line in lines[1:]:
                fields = line.split('\t')
                if len(fields) <= max(idx_cid, idx_title, idx_iupac):
                    continue  # Skip malformed lines
                    
                cid = fields[idx_cid]
                
                # Try to get the Title first, then IUPAC name
                name = ""
                if idx_title >= 0 and idx_title < len(fields) and fields[idx_title]:
                    name = fields[idx_title]
                elif idx_iupac >= 0 and idx_iupac < len(fields) and fields[idx_iupac]:
                    name = fields[idx_iupac]
                    
                if name:
                    cid_to_name[cid] = name
            
            time.sleep(delay)  # Be polite to PubChem API
            
        except Exception as e:
            logger.warning(f"Error processing batch {i//batch_size + 1}: {str(e)}")
            time.sleep(delay * 2)  # Wait longer after an error
            
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

def save_cids_to_file(cid_to_name, output_path):
    """Save CIDs and names to a file."""
    try:
        with open(output_path, "w", encoding="utf-8") as f:
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
        "mesh_terms": 0,
        "keywords": 0,
        "chemical_classes": 0,
        "similar_compounds": 0,
        "total": len(existing_cids)
    }
    
    # 1. Fetch CIDs by MeSH terms
    logger.info("Fetching CIDs by MeSH terms...")
    mesh_cids = set()
    logger.info(f"Processing {len(MESH_TERMS)} MeSH terms...")
    for idx, mesh_id in enumerate(MESH_TERMS):
        logger.info(f"Processing MeSH term {idx+1}/{len(MESH_TERMS)}: {mesh_id}")
        new_cids = fetch_cids_by_mesh(mesh_id, api_delay)
        mesh_cids.update(new_cids)
        logger.info(f"Found {len(new_cids)} CIDs for MeSH term {mesh_id}")
    
    # Add new CIDs to the set
    new_mesh_cids = mesh_cids - all_cids
    all_cids.update(new_mesh_cids)
    stats["mesh_terms"] = len(new_mesh_cids)
    stats["total"] = len(all_cids)
    logger.info(f"Found {len(new_mesh_cids)} new CIDs via MeSH terms. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 2. Fetch CIDs by keywords
        logger.info("Fetching CIDs by keywords...")
        keyword_cids = set()
        logger.info(f"Processing {len(KEYWORDS)} keywords...")
        for idx, keyword in enumerate(KEYWORDS):
            logger.info(f"Processing keyword {idx+1}/{len(KEYWORDS)}: {keyword}")
            new_cids = fetch_cids_by_keyword(keyword, api_delay)
            keyword_cids.update(new_cids)
            logger.info(f"Found {len(new_cids)} CIDs for keyword '{keyword}'")
        
        # Add new CIDs to the set
        new_keyword_cids = keyword_cids - all_cids
        all_cids.update(new_keyword_cids)
        stats["keywords"] = len(new_keyword_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_keyword_cids)} new CIDs via keywords. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 3. Fetch CIDs by chemical classes
        logger.info("Fetching CIDs by chemical classes...")
        class_cids = set()
        logger.info(f"Processing {len(CHEMICAL_CLASSES)} chemical classes...")
        for idx, chemical_class in enumerate(CHEMICAL_CLASSES):
            logger.info(f"Processing chemical class {idx+1}/{len(CHEMICAL_CLASSES)}: {chemical_class}")
            new_cids = fetch_cids_by_chemical_class(chemical_class, api_delay)
            class_cids.update(new_cids)
            logger.info(f"Found {len(new_cids)} CIDs for chemical class '{chemical_class}'")
        
        # Add new CIDs to the set
        new_class_cids = class_cids - all_cids
        all_cids.update(new_class_cids)
        stats["chemical_classes"] = len(new_class_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_class_cids)} new CIDs via chemical classes. Total unique CIDs: {len(all_cids)}")
    
    # Check if we've reached the target
    if len(all_cids) >= target_count:
        logger.info(f"Reached target of {target_count} CIDs. Stopping search.")
    else:
        # 4. Fetch similar compounds to core cryoprotectants
        logger.info("Fetching similar compounds to core cryoprotectants...")
        similar_cids = set()
        logger.info(f"Processing {len(CORE_CRYOPROTECTANT_CIDS)} core cryoprotectant CIDs...")
        for idx, cid in enumerate(CORE_CRYOPROTECTANT_CIDS):
            logger.info(f"Processing core cryoprotectant {idx+1}/{len(CORE_CRYOPROTECTANT_CIDS)}: CID {cid}")
            # Start with a high similarity threshold and gradually lower it if needed
            for similarity in [95, 90, 85, 80]:
                new_cids = fetch_similar_compounds(cid, similarity, api_delay)
                similar_cids.update(new_cids)
                logger.info(f"Found {len(new_cids)} CIDs similar to CID {cid} at {similarity}% similarity")
                
                # If we found a good number of similar compounds, stop lowering the threshold
                if len(new_cids) > 10:
                    break
        
        # Add new CIDs to the set
        new_similar_cids = similar_cids - all_cids
        all_cids.update(new_similar_cids)
        stats["similar_compounds"] = len(new_similar_cids)
        stats["total"] = len(all_cids)
        logger.info(f"Found {len(new_similar_cids)} new CIDs via similar compounds. Total unique CIDs: {len(all_cids)}")
    
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
    parser = argparse.ArgumentParser(description="Expand CID list for cryoprotectants.")
    parser.add_argument("--output", type=str, default="CID-Synonym-curated", help="Output file path")
    parser.add_argument("--target", type=int, default=5000, help="Target number of CIDs")
    parser.add_argument("--delay", type=float, default=0.2, help="Delay between PubChem API calls in seconds")
    args = parser.parse_args()
    
    main(args.output, args.target, args.delay)
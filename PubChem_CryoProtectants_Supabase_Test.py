#!/usr/bin/env python3
"""
CryoProtect Analyzer - PubChem Cryoprotectant Data Importer (Supabase Version, Batch/Checkpoint/CLI Refactor)
Modified for testing with the new schema

This script retrieves cryoprotectant data from PubChem, filters molecules based on
predefined criteria, scores them, and stores the results in a Supabase database.

Features:
- Batch processing with parameterizable batch size
- Robust checkpointing for resumable operation
- CLI arguments for batch size, checkpoint path, resume/reset
- Efficient Supabase bulk inserts for molecules and properties
- Progress and ETA logging
- Error/skipped CID logging for review
"""

import os
import time
import requests
import json
import argparse
from datetime import datetime, timedelta
import logging
from dotenv import load_dotenv
from supabase import create_client, Client

# Load environment variables from .env file
load_dotenv()

# Supabase connection
supabase_url = os.getenv("SUPABASE_URL")
supabase_key = os.getenv("SUPABASE_KEY")

if not supabase_url or not supabase_key:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

supabase: Client = create_client(supabase_url, supabase_key)

# Test project ID (created for testing)
TEST_PROJECT_ID = "15e7f0ae-b909-42ff-b980-73a10ebbfcca"

# Scoring Weights (Total = 200)
WEIGHTS = {
    "hydrogen_bonding": 50,
    "solubility_polarity": 40,
    "membrane_permeability": 40,
    "toxicity_biocompatibility": 30,
    "protein_stabilization": 20,
    "stability_reactivity": 10,
    "environmental_safety": 10
}

# Stage 1: Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-5, 5),  # Relaxed for testing
    "mw_range": (0, 1000),  # Relaxed for testing
    "TPSA_range": (0, 200), # Relaxed for testing
    "functional_groups": [] # No functional group requirement for testing
}

# PubChem API delay
PUBCHEM_API_DELAY = float(os.getenv("PUBCHEM_API_DELAY", "0.2"))

# CID File
CID_FILE = "CID-Synonym-test"  # Test dataset for validation

# Set up logging
LOG_FILE = "cryoprotectant_analysis.log"
SKIPPED_CID_LOG = "skipped_cids.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def get_cid_list():
    """Retrieve all CIDs from the downloaded PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.warning(f"WARNING: CID file '{CID_FILE}' not found. Download it first.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

def get_molecule_properties(cid):
    """Fetch molecular properties from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IsomericSMILES,InChI,InChIKey/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        properties = data["PropertyTable"]["Properties"][0]
        return {
            "CID": cid,
            "Molecular Formula": properties.get("MolecularFormula"),
            "Molecular Weight": properties.get("MolecularWeight"),
            "LogP": properties.get("XLogP"),
            "TPSA": properties.get("TPSA"),
            "H-Bond Donors": properties.get("HBondDonorCount"),
            "H-Bond Acceptors": properties.get("HBondAcceptorCount"),
            "SMILES": properties.get("IsomericSMILES"),
            "InChI": properties.get("InChI"),
            "InChIKey": properties.get("InChIKey"),
            "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }
    logger.warning(f"WARNING: No molecular properties found for CID {cid}.")
    return {"CID": cid, "Error": "No data found"}

def get_additional_properties(cid):
    """Fetch toxicity, stability, and environmental impact from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        props = data["PC_Compounds"][0]["props"]

        extra_data = {
            "Toxicity": None,
            "Stability": None,
            "Environmental Safety": None
        }

        for prop in props:
            label = prop.get("urn", {}).get("label", "")
            value = prop.get("value", {}).get("sval", "")

            if "Toxicity" in label:
                extra_data["Toxicity"] = value
            elif "Stability" in label:
                extra_data["Stability"] = value
            elif "Environmental" in label:
                extra_data["Environmental Safety"] = value

        return extra_data
    logger.warning(f"WARNING: No additional properties found for CID {cid}.")
    return {"Error": "No additional data found"}

def filter_molecule(molecule):
    """Initial filtering based on core cryoprotectant properties, ensuring numerical values."""
    if "Error" in molecule:
        return False

    # Check for required fields for Supabase schema
    if not molecule.get("SMILES"):
        logger.warning(f"SKIPPED: Skipping CID {molecule['CID']} due to missing SMILES.")
        return False
    if not molecule.get("InChI"):
        logger.warning(f"SKIPPED: Skipping CID {molecule['CID']} due to missing InChI.")
        return False
    if not molecule.get("InChIKey"):
        logger.warning(f"SKIPPED: Skipping CID {molecule['CID']} due to missing InChIKey.")
        return False
    if not molecule.get("Molecular Formula"):
        logger.warning(f"SKIPPED: Skipping CID {molecule['CID']} due to missing Molecular Formula.")
        return False

    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError):
        logger.warning(f"SKIPPED: Skipping CID {molecule['CID']} due to invalid numerical values.")
        return False

    smiles = molecule["SMILES"]

    if mw is None or not (CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1]):
        return False
    if logp is None or not (CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1]):
        return False
    if tpsa is None or not (CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
        return False
    if CORE_CRITERIA["functional_groups"]:
        if not any(group in smiles for group in CORE_CRITERIA["functional_groups"]):
            return False

    return True

def score_molecule(molecule, extra_properties):
    """Compute final score out of 200 based on all properties."""
    score = 0

    # Hydrogen Bonding
    score += WEIGHTS["hydrogen_bonding"]

    # Solubility & Permeability
    score += WEIGHTS["solubility_polarity"]
    score += WEIGHTS["membrane_permeability"]

    # Toxicity & Stability
    if extra_properties.get("Toxicity"):
        score += WEIGHTS["toxicity_biocompatibility"]
    if extra_properties.get("Stability"):
        score += WEIGHTS["stability_reactivity"]

    # Environmental Safety
    if extra_properties.get("Environmental Safety"):
        score += WEIGHTS["environmental_safety"]

    return score

def fetch_property_types():
    """Fetch property types from Supabase once per run."""
    response = supabase.table("property_types").select("id, name, data_type").execute()
    logger.info(f"DEBUG: fetch_property_types() response type: {type(response)}; content: {response}")
    # Attempt to access .data and .error for debugging
    data = getattr(response, "data", None)
    error = getattr(response, "error", None)
    logger.info(f"DEBUG: response.data = {data}, response.error = {error}")
    if error:
        logger.error(f"Error fetching property types: {error}")
        return None
    return data

def bulk_insert_molecules(molecule_batch):
    """Bulk insert molecules, return list of inserted molecule IDs (in order)."""
    if not molecule_batch:
        return []
    try:
        response = supabase.table("molecule").insert(molecule_batch).execute()
        # Check if response has data attribute
        if hasattr(response, 'data') and response.data:
            return [row["id"] for row in response.data]
        else:
            logger.error(f"Error bulk inserting molecules: No data in response")
            return []
    except Exception as e:
        logger.error(f"Error bulk inserting molecules: {str(e)}")
        return []

def bulk_insert_properties(property_batch):
    """Bulk insert molecular properties."""
    if not property_batch:
        return True
    try:
        response = supabase.table("molecular_property").insert(property_batch).execute()
        # Check if response has data attribute
        if hasattr(response, 'data'):
            return True
        else:
            logger.error(f"Error bulk inserting properties: No data in response")
            return False
    except Exception as e:
        logger.error(f"Error bulk inserting properties: {str(e)}")
        return False

def read_checkpoint(checkpoint_path):
    if not os.path.exists(checkpoint_path):
        return None
    with open(checkpoint_path, "r") as f:
        return json.load(f)

def write_checkpoint(checkpoint_path, checkpoint_data):
    with open(checkpoint_path, "w") as f:
        json.dump(checkpoint_data, f)

def log_skipped_cid(cid, reason, skipped_log_path):
    with open(skipped_log_path, "a") as f:
        f.write(f"{cid}\t{reason}\n")

def estimate_time_remaining(start_time, batches_done, total_batches):
    elapsed = time.time() - start_time
    avg_per_batch = elapsed / batches_done if batches_done else 0
    remaining_batches = total_batches - batches_done
    eta_seconds = avg_per_batch * remaining_batches
    return str(timedelta(seconds=int(eta_seconds)))

def process_batches(
    cids,
    batch_size,
    checkpoint_path,
    resume,
    reset,
    skipped_log_path
):
    logger.info("STARTED: Starting batch dataset processing...")
    total_cids = len(cids)
    total_batches = (total_cids + batch_size - 1) // batch_size

    # Checkpoint logic
    checkpoint = None
    if resume and os.path.exists(checkpoint_path):
        checkpoint = read_checkpoint(checkpoint_path)
        start_batch = checkpoint.get("last_completed_batch", 0) + 1
        logger.info(f"Resuming from batch {start_batch} (checkpoint: {checkpoint_path})")
    elif reset and os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)
        logger.info(f"Resetting checkpoint at {checkpoint_path}")
        start_batch = 0
    else:
        start_batch = 0

    # Fetch property types once
    property_types = fetch_property_types()
    if property_types is None:
        logger.error("Could not fetch property types. Exiting.")
        return

    total_processed = 0
    total_imported = 0
    start_time = time.time()

    for batch_num in range(start_batch, total_batches):
        batch_start = batch_num * batch_size
        batch_end = min(batch_start + batch_size, total_cids)
        batch_cids = cids[batch_start:batch_end]

        molecules_to_insert = []
        molecule_cid_map = []
        batch_properties = []
        skipped_in_batch = 0

        for cid in batch_cids:
            try:
                molecule = get_molecule_properties(cid)
                logger.info(f"DEBUG: CID {cid} molecule properties before filtering: {molecule}")
                if not filter_molecule(molecule):
                    logger.info(f"DEBUG: CID {cid} did not pass filter criteria: {molecule}")
                    log_skipped_cid(cid, "Did not pass filter", skipped_log_path)
                    skipped_in_batch += 1
                    continue

                extra_properties = get_additional_properties(cid)
                score = score_molecule(molecule, extra_properties)

                # Modified for new schema
                molecules_to_insert.append({
                    "name": f"PubChem CID: {molecule['CID']}",
                    "smiles": molecule.get("SMILES"),
                    "inchi": molecule.get("InChI"),
                    "inchikey": molecule.get("InChIKey"),
                    "formula": molecule.get("Molecular Formula"),
                    "molecular_weight": float(molecule.get("Molecular Weight")) if molecule.get("Molecular Weight") else None,
                    "data_source": "PubChem",
                    "version": 1,
                    "project_id": TEST_PROJECT_ID,  # Use the test project ID
                    "modification_history": [{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "source": "PubChem_CryoProtectants_Supabase_Test.py"
                    }]
                })
                molecule_cid_map.append((cid, molecule, extra_properties, score))
            except Exception as e:
                logger.error(f"ERROR: Error processing CID {cid}: {str(e)}")
                log_skipped_cid(cid, f"Exception: {str(e)}", skipped_log_path)
                skipped_in_batch += 1

            time.sleep(PUBCHEM_API_DELAY)

        # Bulk insert molecules
        inserted_ids = bulk_insert_molecules(molecules_to_insert)
        if not inserted_ids or len(inserted_ids) != len(molecule_cid_map):
            logger.error(f"Bulk insert mismatch: {len(inserted_ids)} IDs for {len(molecule_cid_map)} molecules.")
            # Log all CIDs in this batch as skipped if bulk insert failed
            for cid, _, _, _ in molecule_cid_map:
                log_skipped_cid(cid, "Bulk insert failed", skipped_log_path)
            continue

        # Prepare property inserts for all molecules in batch
        for idx, (cid, molecule, extra_properties, score) in enumerate(molecule_cid_map):
            molecule_id = inserted_ids[idx]
            
            # Modified for new schema - using property_type TEXT and value columns
            properties_to_insert = [
                {
                    "molecule_id": molecule_id,
                    "property_type": "LogP",
                    "value": float(molecule.get("LogP")) if molecule.get("LogP") is not None else None,
                    "unit": "",
                    "data_source": "PubChem"
                },
                {
                    "molecule_id": molecule_id,
                    "property_type": "TPSA",
                    "value": float(molecule.get("TPSA")) if molecule.get("TPSA") is not None else None,
                    "unit": "Å²",
                    "data_source": "PubChem"
                },
                {
                    "molecule_id": molecule_id,
                    "property_type": "H-Bond Donors",
                    "value": float(molecule.get("H-Bond Donors")) if molecule.get("H-Bond Donors") is not None else None,
                    "unit": "",
                    "data_source": "PubChem"
                },
                {
                    "molecule_id": molecule_id,
                    "property_type": "H-Bond Acceptors",
                    "value": float(molecule.get("H-Bond Acceptors")) if molecule.get("H-Bond Acceptors") is not None else None,
                    "unit": "",
                    "data_source": "PubChem"
                },
                {
                    "molecule_id": molecule_id,
                    "property_type": "Total Score",
                    "value": float(score),
                    "unit": "",
                    "data_source": "PubChem_CryoProtectants_Supabase_Test.py"
                },
                {
                    "molecule_id": molecule_id,
                    "property_type": "PubChem CID",
                    "value": float(molecule.get("CID")),
                    "unit": "",
                    "data_source": "PubChem"
                }
            ]
            
            # Add text properties if they exist
            if extra_properties.get("Toxicity"):
                properties_to_insert.append({
                    "molecule_id": molecule_id,
                    "property_type": "Toxicity",
                    "value": 0,  # Placeholder numeric value
                    "unit": "",
                    "data_source": "PubChem",
                    "provenance": extra_properties.get("Toxicity")  # Store text in provenance field
                })
                
            if extra_properties.get("Stability"):
                properties_to_insert.append({
                    "molecule_id": molecule_id,
                    "property_type": "Stability",
                    "value": 0,  # Placeholder numeric value
                    "unit": "",
                    "data_source": "PubChem",
                    "provenance": extra_properties.get("Stability")  # Store text in provenance field
                })
                
            if extra_properties.get("Environmental Safety"):
                properties_to_insert.append({
                    "molecule_id": molecule_id,
                    "property_type": "Environmental Safety",
                    "value": 0,  # Placeholder numeric value
                    "unit": "",
                    "data_source": "PubChem",
                    "provenance": extra_properties.get("Environmental Safety")  # Store text in provenance field
                })
                
            batch_properties.extend(properties_to_insert)

        # Bulk insert properties
        if not bulk_insert_properties(batch_properties):
            logger.error(f"Failed to bulk insert properties for batch {batch_num}")
            # Log all CIDs in this batch as skipped for properties
            for cid, _, _, _ in molecule_cid_map:
                log_skipped_cid(cid, "Bulk property insert failed", skipped_log_path)
            continue

        total_processed += len(batch_cids)
        total_imported += len(inserted_ids)
        batches_done = batch_num + 1
        eta = estimate_time_remaining(start_time, batches_done, total_batches)
        logger.info(
            f"Batch {batch_num+1}/{total_batches} complete: {total_processed}/{total_cids} CIDs processed, "
            f"{total_imported} imported, {skipped_in_batch} skipped in this batch. ETA: {eta}"
        )

        # Save checkpoint
        checkpoint_data = {
            "last_completed_batch": batch_num,
            "total_processed": total_processed,
            "total_imported": total_imported,
            "timestamp": datetime.now().isoformat()
        }
        write_checkpoint(checkpoint_path, checkpoint_data)

    logger.info(f"SUCCESS: All batches complete! Processed {total_processed} molecules, imported {total_imported} to database.")

def main():
    parser = argparse.ArgumentParser(description="PubChem CryoProtectant Data Importer (Supabase, Batch/Checkpoint)")
    parser.add_argument("--batch-size", type=int, default=2, help="Batch size for processing (default: 2)")
    parser.add_argument("--checkpoint", type=str, default="checkpoint.json", help="Path to checkpoint file")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    parser.add_argument("--reset", action="store_true", help="Reset checkpoint and start from scratch")
    parser.add_argument("--skipped-log", type=str, default=SKIPPED_CID_LOG, help="Path to skipped CIDs log file")
    args = parser.parse_args()

    # Load CIDs
    cids = get_cid_list()
    if not cids:
        logger.error("No CIDs loaded. Exiting.")
        return

    # Remove skipped log if resetting
    if args.reset and os.path.exists(args.skipped_log):
        os.remove(args.skipped_log)

    process_batches(
        cids=cids,
        batch_size=args.batch_size,
        checkpoint_path=args.checkpoint,
        resume=args.resume,
        reset=args.reset,
        skipped_log_path=args.skipped_log
    )

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"ERROR: Fatal error: {str(e)}")
#!/usr/bin/env python3
"""
CryoProtect v2 - PubChem Incremental Data Updater

This script fetches new or updated PubChem entries and inserts them into the database,
supporting incremental updates, validation, and robust logging. Designed for scheduled
execution via cron (Linux/macOS) or Task Scheduler (Windows).

Configuration:
- Environment variables or .env file:
    SUPABASE_URL, SUPABASE_KEY, SUPABASE_USER, SUPABASE_PASSWORD
    PUBCHEM_CID_FILE (default: CID-Synonym-filtered)
    LOG_LEVEL, LOG_FILE, LOG_TO_FILE
    PUBCHEM_BATCH_SIZE (default: 50)
    PUBCHEM_REQUEST_DELAY (default: 0.2)
    UPDATE_FREQUENCY (for documentation only; not enforced by script)
    PUBCHEM_QUERY_PARAMS (for future extension)

Exit Codes:
- 0: Success
- 1: General error
- 2: Validation error
- 3: Database error

Example cron entry (run daily at 2:30am):
30 2 * * * /usr/bin/env python3 /path/to/update_pubchem_data.py

Example Task Scheduler:
- Action: Start a program
- Program/script: python
- Add arguments: C:\path\to\update_pubchem_data.py

"""

import os
import sys
import time
import logging
from dotenv import load_dotenv
from supabase import create_client, Client
import requests
from datetime import datetime
import traceback

# Import centralized logging setup
import logging_config

# Validation parameters (can be imported or duplicated from validate_cryoprotectants_data.py)
VALIDATION_PARAMS = {
    "molecular_weight_range": (0, 2000),
    "logp_range": (-10, 10),
    "tpsa_range": (0, 500),
    "hbond_donors_range": (0, 20),
    "hbond_acceptors_range": (0, 20),
}

REQUIRED_FIELDS = [
    "MolecularFormula", "MolecularWeight", "XLogP", "TPSA",
    "HBondDonorCount", "HBondAcceptorCount", "IsomericSMILES"
]

def validate_pubchem_data(data):
    """Validate a PubChem record for required fields and scientific plausibility."""
    for field in REQUIRED_FIELDS:
        if field not in data or data[field] is None:
            return False, f"Missing required field: {field}"

    mw = data["MolecularWeight"]
    logp = data["XLogP"]
    tpsa = data["TPSA"]
    hbd = data["HBondDonorCount"]
    hba = data["HBondAcceptorCount"]

    if not (VALIDATION_PARAMS["molecular_weight_range"][0] <= mw <= VALIDATION_PARAMS["molecular_weight_range"][1]):
        return False, f"Molecular weight {mw} out of range"
    if not (VALIDATION_PARAMS["logp_range"][0] <= logp <= VALIDATION_PARAMS["logp_range"][1]):
        return False, f"LogP {logp} out of range"
    if not (VALIDATION_PARAMS["tpsa_range"][0] <= tpsa <= VALIDATION_PARAMS["tpsa_range"][1]):
        return False, f"TPSA {tpsa} out of range"
    if not (VALIDATION_PARAMS["hbond_donors_range"][0] <= hbd <= VALIDATION_PARAMS["hbond_donors_range"][1]):
        return False, f"HBondDonorCount {hbd} out of range"
    if not (VALIDATION_PARAMS["hbond_acceptors_range"][0] <= hba <= VALIDATION_PARAMS["hbond_acceptors_range"][1]):
        return False, f"HBondAcceptorCount {hba} out of range"
    return True, "Valid"

def fetch_pubchem_properties(cid):
    """Fetch molecular properties from PubChem for a given CID."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IsomericSMILES/JSON"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        props = response.json()["PropertyTable"]["Properties"][0]
        return props
    except Exception as e:
        logging.error(f"Failed to fetch properties for CID {cid}: {e}")
        return None

def load_existing_cids(supabase: Client):
    """Fetch all CIDs currently in the molecules table."""
    try:
        result = supabase.table("molecules").select("cid").execute()
        if result.data is None:
            return set()
        return set(row["cid"] for row in result.data if "cid" in row)
    except Exception as e:
        logging.error(f"Error loading existing CIDs: {e}")
        return set()

def load_candidate_cids(cid_file):
    """Load candidate CIDs from the specified file."""
    if not os.path.exists(cid_file):
        logging.error(f"CID file '{cid_file}' not found.")
        return []
    with open(cid_file, "r") as f:
        cids = [int(line.strip().split("\t")[0]) for line in f if line.strip().split("\t")[0].isdigit()]
    return cids

def insert_molecule(supabase: Client, cid, props):
    """Insert a new molecule into the database."""
    try:
        molecule_data = {
            "cid": cid,
            "name": None,  # Name can be fetched/added if available
            "smiles": props.get("IsomericSMILES"),
            "molecular_weight": props.get("MolecularWeight"),
            "formula": props.get("MolecularFormula"),
            "logp": props.get("XLogP"),
            "tpsa": props.get("TPSA"),
            "hbond_donors": props.get("HBondDonorCount"),
            "hbond_acceptors": props.get("HBondAcceptorCount"),
            "created_at": datetime.utcnow().isoformat() + "Z",
        }
        result = supabase.table("molecules").insert(molecule_data).execute()
        if result.status_code not in (200, 201):
            raise Exception(f"Insert failed: {result.status_code} {result.data}")
        return True
    except Exception as e:
        logging.error(f"Failed to insert CID {cid}: {e}")
        return False

def main():
    # Load config and set up logging
    load_dotenv()
    logging_config.setup_logging()
    logger = logging.getLogger("update_pubchem_data")

    # Read config from env
    SUPABASE_URL = os.getenv("SUPABASE_URL")
    SUPABASE_KEY = os.getenv("SUPABASE_KEY")
    CID_FILE = os.getenv("PUBCHEM_CID_FILE", "CID-Synonym-filtered")
    BATCH_SIZE = int(os.getenv("PUBCHEM_BATCH_SIZE", "50"))
    REQUEST_DELAY = float(os.getenv("PUBCHEM_REQUEST_DELAY", "0.2"))

    if not SUPABASE_URL or not SUPABASE_KEY:
        logger.error("SUPABASE_URL and SUPABASE_KEY must be set in environment or .env file.")
        sys.exit(1)

    try:
        supabase: Client = create_client(SUPABASE_URL, SUPABASE_KEY)
    except Exception as e:
        logger.error(f"Failed to connect to Supabase: {e}")
        sys.exit(1)

    logger.info("Loading existing CIDs from database...")
    existing_cids = load_existing_cids(supabase)
    logger.info(f"Found {len(existing_cids)} CIDs in database.")

    logger.info(f"Loading candidate CIDs from {CID_FILE}...")
    candidate_cids = load_candidate_cids(CID_FILE)
    logger.info(f"Loaded {len(candidate_cids)} candidate CIDs.")

    new_cids = [cid for cid in candidate_cids if cid not in existing_cids]
    logger.info(f"{len(new_cids)} new CIDs to process.")

    success_count = 0
    validation_failures = 0
    db_failures = 0

    for i, cid in enumerate(new_cids):
        logger.info(f"[{i+1}/{len(new_cids)}] Processing CID {cid}...")
        props = fetch_pubchem_properties(cid)
        if not props:
            logger.error(f"Skipping CID {cid} due to fetch error.")
            continue

        valid, reason = validate_pubchem_data(props)
        if not valid:
            logger.warning(f"Validation failed for CID {cid}: {reason}")
            validation_failures += 1
            continue

        inserted = insert_molecule(supabase, cid, props)
        if inserted:
            logger.info(f"Inserted CID {cid} successfully.")
            success_count += 1
        else:
            logger.error(f"Database insert failed for CID {cid}.")
            db_failures += 1

        time.sleep(REQUEST_DELAY)

    logger.info(f"Update complete. Inserted: {success_count}, Validation failures: {validation_failures}, DB failures: {db_failures}")

    if db_failures > 0:
        sys.exit(3)
    if validation_failures > 0:
        sys.exit(2)
    sys.exit(0)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"Fatal error: {e}\n{traceback.format_exc()}")
        sys.exit(1)
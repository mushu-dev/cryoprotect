import requests
import time
import pandas as pd
import json
import os
import sys
from datetime import datetime

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
    "logP_range": (-1.5, 0),
    "mw_range": (50, 150),
    "TPSA_range": (40, 100),
    "functional_groups": ["OH", "CONH2", "S=O"]
}

# File Paths
CSV_FILE = "cryoprotectants.csv"
JSON_FILE = "cryoprotectants.json"
LOG_FILE = "cryoprotectant_analysis.log"
CID_FILE = "CID-Synonym-filtered"  # Updated CID list location

# Logging Function
def log_message(message):
    """Append message to log file with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as log:
        log.write(f"[{timestamp}] {message}\n")
    print(message)

# Fetch CIDs from the local NIH file instead of the PubChem API
def get_cid_list():
    """Retrieve all CIDs from the downloaded PubChem CID list."""
    if not os.path.exists(CID_FILE):
        log_message(f"⚠️ CID file '{CID_FILE}' not found. Download it first.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    log_message(f"✅ Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

# Fetch Molecular Properties
def get_molecule_properties(cid):
    """Fetch molecular properties from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IsomericSMILES/JSON"
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
            "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        }
    log_message(f"⚠️ Warning: No molecular properties found for CID {cid}.")
    return {"CID": cid, "Error": "No data found"}

# Fetch Additional Properties (Toxicity, Stability, Environmental Impact)
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
    log_message(f"⚠️ Warning: No additional properties found for CID {cid}.")
    return {"Error": "No additional data found"}

# Stage 1 Filtering with Float Conversion
def filter_molecule(molecule):
    """Initial filtering based on core cryoprotectant properties, ensuring numerical values."""
    if "Error" in molecule:
        return False

    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError):
        log_message(f"⚠️ Skipping CID {molecule['CID']} due to invalid numerical values.")
        return False

    smiles = molecule["SMILES"]

    if mw is None or not (50 <= mw <= 150):
        return False
    if logp is None or not (-1.5 <= logp <= 0):
        return False
    if tpsa is None or not (40 <= tpsa <= 100):
        return False
    if not any(group in smiles for group in ["OH", "CONH2", "S=O"]):
        return False

    return True

# Stage 2 Scoring
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

if __name__ == "__main__":
    log_message("🚀 Starting full dataset processing...")

    # Fetch all CIDs from the file
    cids = get_cid_list()

    for cid in cids:
        molecule = get_molecule_properties(cid)
        if filter_molecule(molecule):
            extra_properties = get_additional_properties(cid)
            molecule.update(extra_properties)
            molecule["Score"] = score_molecule(molecule, extra_properties)

            # Save progress after each molecule
            pd.DataFrame([molecule]).to_csv(CSV_FILE, mode="a", header=not os.path.exists(CSV_FILE), index=False)
            json.dump(molecule, open(JSON_FILE, "a"), indent=4)
            log_message(f"✅ Processed CID {cid}: Score = {molecule['Score']}")

        time.sleep(0.2)

    log_message(f"✅ Processing complete! Data saved as '{CSV_FILE}' and '{JSON_FILE}'.")


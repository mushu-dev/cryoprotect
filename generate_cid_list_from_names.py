#!/usr/bin/env python3
"""
generate_cid_list_from_names.py

Queries PubChem for a list of known cryoprotectant names and outputs a tab-separated
file of CID and name for production-scale import.

Usage:
    python generate_cid_list_from_names.py --output CID-Synonym-curated

Dependencies:
    - requests
"""

import requests
import argparse
import time

CRYOPROTECTANT_NAMES = [
    "dimethyl sulfoxide",
    "glycerol",
    "ethylene glycol",
    "propylene glycol",
    "trehalose",
    "sucrose",
    "methanol",
    "formamide",
    "acetamide",
    "1,2-propanediol",
    "ethylene glycol",
    "1,3-propanediol",
    "mannitol",
    "sorbitol",
    "glucose",
    "fructose",
    "raffinose",
    "xylose",
    "erythritol",
    "DMSO",
    "EG",
    "PG",
    "PROH",
    "Ficoll",
    "polyvinylpyrrolidone",
    "polyethylene glycol",
    "hydroxyethyl starch",
    "dextran",
    "albumin",
    "serum",
    "acetonitrile",
    "acetone",
    "butanediol",
    "isopropanol",
    "ethanol",
    "diethylene glycol",
    "triethylene glycol",
    "tetraethylene glycol"
]

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def fetch_cid_for_name(name):
    """Fetch the CID for a given compound name from PubChem."""
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/TXT"
    r = requests.get(url)
    if r.status_code == 200:
        lines = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
        if lines:
            return lines[0]  # Take the first CID if multiple
    return None

def main(output_path):
    cid_to_name = {}
    for name in CRYOPROTECTANT_NAMES:
        print(f"Querying PubChem for: {name}")
        cid = fetch_cid_for_name(name)
        if cid:
            cid_to_name[cid] = name
            print(f"  Found CID: {cid}")
        else:
            print(f"  No CID found for: {name}")
        time.sleep(0.2)  # Be polite to PubChem

    with open(output_path, "w", encoding="utf-8") as f:
        for cid, name in cid_to_name.items():
            f.write(f"{cid}\t{name}\n")
    print(f"Output written to {output_path} ({len(cid_to_name)} entries)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CID list from known cryoprotectant names.")
    parser.add_argument("--output", type=str, default="CID-Synonym-curated", help="Output file path")
    args = parser.parse_args()
    main(args.output)
#!/usr/bin/env python3
"""
generate_production_cid_list.py

Generates a comprehensive, production-grade list of PubChem CIDs and compound names
for all known cryoprotectants, suitable for use with PubChem_CryoProtectants_Supabase_Enhanced.py.

Approach:
- Uses PubChem PUG-REST API to identify cryoprotectants via:
    1. MeSH term "Cryoprotective Agents" (D003440)
    2. Keyword search for "cryoprotectant" in compound synonyms/titles
- Retrieves the preferred compound name for each CID
- Outputs a tab-separated file: CID<TAB>Name, one per line, no header
- Documents the process for reproducibility

Requirements:
- Python 3.6+
- requests

Usage:
    python generate_production_cid_list.py --output CID-Synonym-production

"""

import requests
import time
import argparse

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def fetch_cids_by_mesh(mesh_id):
    """Fetch CIDs for compounds annotated with a given MeSH term."""
    url = f"{PUBCHEM_BASE}/compound/mesh/{mesh_id}/cids/TXT"
    r = requests.get(url)
    r.raise_for_status()
    cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
    return set(cids)

def fetch_cids_by_keyword(keyword):
    """Fetch CIDs for compounds matching a keyword in synonyms/titles."""
    url = f"{PUBCHEM_BASE}/compound/name/{keyword}/cids/TXT"
    r = requests.get(url)
    if r.status_code == 200:
        cids = [line.strip() for line in r.text.splitlines() if line.strip().isdigit()]
        return set(cids)
    return set()

def fetch_names_for_cids(cids, batch_size=100):
    """Fetch preferred names for a list of CIDs using PubChem."""
    cid_list = list(cids)
    cid_to_name = {}
    for i in range(0, len(cid_list), batch_size):
        batch = cid_list[i:i+batch_size]
        ids = ",".join(batch)
        url = f"{PUBCHEM_BASE}/compound/cid/{ids}/property/IUPACName,Title,CanonicalSMILES/TSV"
        r = requests.get(url)
        if r.status_code != 200:
            time.sleep(1)
            continue
        lines = r.text.splitlines()
        header = lines[0].split('\t')
        idx_cid = header.index("CID")
        idx_title = header.index("Title")
        idx_iupac = header.index("IUPACName")
        for line in lines[1:]:
            fields = line.split('\t')
            cid = fields[idx_cid]
            name = fields[idx_title] if fields[idx_title] else fields[idx_iupac]
            cid_to_name[cid] = name
        time.sleep(0.2)
    return cid_to_name

def main(output_path):
    # 1. Fetch CIDs by MeSH term "Cryoprotective Agents" (D003440)
    mesh_cids = fetch_cids_by_mesh("D003440")
    print(f"Found {len(mesh_cids)} CIDs via MeSH term.")

    # 2. Fetch CIDs by keyword "cryoprotectant"
    keyword_cids = fetch_cids_by_keyword("cryoprotectant")
    print(f"Found {len(keyword_cids)} CIDs via keyword search.")

    # 3. Combine and deduplicate
    all_cids = mesh_cids.union(keyword_cids)
    print(f"Total unique CIDs: {len(all_cids)}")

    # 4. Fetch names for all CIDs
    cid_to_name = fetch_names_for_cids(all_cids)
    print(f"Retrieved names for {len(cid_to_name)} CIDs.")

    # 5. Write to output file
    with open(output_path, "w", encoding="utf-8") as f:
        for cid, name in cid_to_name.items():
            if name:
                f.write(f"{cid}\t{name}\n")
    print(f"Output written to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate production-grade CID list for cryoprotectants.")
    parser.add_argument("--output", type=str, default="CID-Synonym-production", help="Output file path")
    args = parser.parse_args()
    main(args.output)
#!/usr/bin/env python3
"""
Analyze duplicate molecule groups from the duplicate_molecules_report.json file.
"""

import json
import os

# Load the duplicate molecules report
with open('duplicate_molecules_report.json', 'r') as f:
    data = json.load(f)

# Analyze name duplicates
print("=== NAME DUPLICATE GROUPS ===")
for group in data['name_duplicates']['groups_data']:
    if group['name'] != 'none':  # Skip the 'None' group
        print(f"\nName Group: {group['name']}")
        print(f"Count: {group['count']}")
        for mol in group['molecules']:
            print(f"  - {mol['name']} (Formula: {mol['molecular_formula'] or 'None'}, PubChem CID: {mol['pubchem_cid'] or 'None'})")

# Analyze formula duplicates
print("\n\n=== FORMULA DUPLICATE GROUPS ===")
for group in data['formula_duplicates']['groups_data']:
    print(f"\nFormula Group: {group['formula']}")
    print(f"Count: {group['count']}")
    for mol in group['molecules']:
        print(f"  - {mol['name']} (Formula: {mol['molecular_formula'] or 'None'}, PubChem CID: {mol['pubchem_cid'] or 'None'})")

# Print summary statistics
print("\n\n=== SUMMARY STATISTICS ===")
print(f"Name duplicate groups: {data['name_duplicates']['groups']}")
print(f"Name duplicate molecules: {data['name_duplicates']['total_molecules']}")
print(f"Formula duplicate groups: {data['formula_duplicates']['groups']}")
print(f"Formula duplicate molecules: {data['formula_duplicates']['total_molecules']}")
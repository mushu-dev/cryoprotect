#!/usr/bin/env python3
"""
Analyze a specific complex merge group to understand the consolidation challenge.
"""

import json
import sys
import argparse

def load_analysis():
    """Load the duplicate groups analysis."""
    with open("duplicate_groups_analysis.json", "r") as f:
        return json.load(f)

def analyze_complex_group(analysis, group_id):
    """Analyze a specific complex merge group."""
    if group_id not in analysis:
        print(f"Group ID {group_id} not found in analysis")
        return
    
    group = analysis[group_id]
    
    if group.get("consolidation_type") != "COMPLEX_MERGE":
        print(f"Group {group_id} is not a complex merge group (type: {group.get('consolidation_type')})")
        return
    
    print(f"=== Complex Merge Group: {group_id} ===")
    print(f"Duplicate type: {group.get('duplicate_type')}")
    print(f"Duplicate value: {group.get('duplicate_value')}")
    print(f"Molecule count: {group.get('molecule_count')}")
    print(f"PubChem CIDs: {group.get('pubchem_cids')}")
    
    print(f"\nChemical Structure Differences:")
    print(f"  Has formula differences: {group.get('has_formula_differences')}")
    print(f"  Has SMILES differences: {group.get('has_smiles_differences')}")
    print(f"  Has property differences: {group.get('has_property_differences')}")
    
    print(f"\nRelationship Information:")
    print(f"  Has relationship conflicts: {group.get('has_relationship_conflicts')}")
    print(f"  Molecules with relationships: {group.get('molecules_with_relationships')}")
    print("  Relationship counts:")
    for rel_type, count in group.get("relationship_counts", {}).items():
        print(f"    {rel_type}: {count}")
    
    print(f"\nPrimary candidate: {group.get('primary_candidate')}")
    
    print(f"\nMolecule Details:")
    for i, mol in enumerate(group.get("molecules_detail", []), 1):
        mol_id = mol.get("id")
        primary_text = " (Primary Candidate)" if mol_id == group.get("primary_candidate") else ""
        
        print(f"\n  Molecule {i}{primary_text}:")
        print(f"    ID: {mol_id}")
        print(f"    Name: {mol.get('name')}")
        print(f"    Formula: {mol.get('molecular_formula')}")
        print(f"    PubChem CID: {mol.get('pubchem_cid')}")
        print(f"    Has relationships: {mol.get('has_relationships')}")
        print(f"    Relationship count: {mol.get('relationship_count')}")
        
        if mol.get("has_relationships"):
            print("    Relationship Details:")
            for rel_type, rels in mol.get("relationships", {}).items():
                if rels:
                    print(f"      {rel_type}: {len(rels)}")
                    if rel_type == "properties" and len(rels) > 0:
                        # Print property details
                        print("        Property Details:")
                        for prop in rels[:3]:  # Show first 3 properties
                            print(f"          {prop.get('property_type')}: {prop.get('property_value')} {prop.get('unit')}")
                        if len(rels) > 3:
                            print(f"          ... and {len(rels) - 3} more properties")
                    
                    if rel_type == "mixture_components" and len(rels) > 0:
                        # Print mixture details
                        print("        Mixture Details:")
                        for mix in rels[:3]:  # Show first 3 mixtures
                            print(f"          Mixture ID: {mix.get('mixture_id')}")
                        if len(rels) > 3:
                            print(f"          ... and {len(rels) - 3} more mixtures")

def main():
    parser = argparse.ArgumentParser(description="Analyze a specific complex merge group")
    parser.add_argument("group_id", help="Group ID to analyze")
    args = parser.parse_args()
    
    analysis = load_analysis()
    analyze_complex_group(analysis, args.group_id)

if __name__ == "__main__":
    main()
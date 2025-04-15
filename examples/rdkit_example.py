#!/usr/bin/env python3
"""
CryoProtect Analyzer - RDKit Integration Example

This script demonstrates the use of the RDKit integration module to calculate
molecular properties, generate visualizations, and perform structure searches.
"""

import os
import sys
import json
from pathlib import Path

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from api.rdkit_utils import (
    calculate_all_properties, generate_molecule_svg,
    perform_substructure_search, calculate_similarity
)

def print_section(title):
    """Print a section title."""
    print("\n" + "=" * 80)
    print(f" {title} ".center(80, "="))
    print("=" * 80)

def print_json(data):
    """Print data as formatted JSON."""
    print(json.dumps(data, indent=2))

def example_calculate_properties():
    """Example of calculating molecular properties."""
    print_section("Calculating Molecular Properties")
    
    # Define some example molecules
    molecules = {
        "Ethanol": "CCO",
        "Glycerol": "C(C(CO)O)O",
        "DMSO": "CS(=O)C",
        "Propylene Glycol": "CC(O)CO",
        "Trehalose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O"
    }
    
    # Calculate properties for each molecule
    for name, smiles in molecules.items():
        print(f"\nProperties for {name} (SMILES: {smiles}):")
        properties = calculate_all_properties(smiles)
        
        # Print selected properties
        print(f"  Molecular Weight: {properties['molecular_properties']['molecular_weight']:.2f} g/mol")
        print(f"  LogP: {properties['logp']:.2f}")
        print(f"  TPSA: {properties['tpsa']:.2f} Å²")
        print(f"  H-Bond Donors: {properties['hydrogen_bonding']['donors']}")
        print(f"  H-Bond Acceptors: {properties['hydrogen_bonding']['acceptors']}")
        
        # Print functional groups if any
        if properties['functional_groups']:
            print("  Functional Groups:")
            for group, count in properties['functional_groups'].items():
                print(f"    {group}: {count}")
        
        # Print permeability properties
        print("  Permeability:")
        print(f"    Rule of 5 Violations: {properties['permeability']['rule_of_5_violations']}")
        print(f"    BBB Permeant: {properties['permeability']['bbb_permeant']}")
        print(f"    Intestinal Absorption: {properties['permeability']['intestinal_absorption']}")

def example_generate_visualization():
    """Example of generating molecular visualizations."""
    print_section("Generating Molecular Visualizations")
    
    # Define some example molecules
    molecules = {
        "Ethanol": "CCO",
        "Glycerol": "C(C(CO)O)O",
        "DMSO": "CS(=O)C",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }
    
    # Generate SVG for each molecule
    for name, smiles in molecules.items():
        print(f"\nGenerating visualization for {name} (SMILES: {smiles}):")
        svg = generate_molecule_svg(smiles)
        
        # Save SVG to file
        output_dir = Path("examples/output")
        output_dir.mkdir(exist_ok=True)
        
        svg_path = output_dir / f"{name.lower().replace(' ', '_')}.svg"
        with open(svg_path, "w") as f:
            f.write(svg)
        
        print(f"  SVG saved to {svg_path}")

def example_substructure_search():
    """Example of performing substructure searches."""
    print_section("Performing Substructure Searches")
    
    # Define some example molecules and substructures
    molecules = {
        "Ethanol": "CCO",
        "Glycerol": "C(C(CO)O)O",
        "DMSO": "CS(=O)C",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }
    
    substructures = {
        "Hydroxyl Group": "[OH]",
        "Methyl Group": "[CH3]",
        "Carbonyl Group": "[CX3]=[OX1]",
        "Sulfoxide Group": "[#16X3](=[OX1])"
    }
    
    # Perform searches
    for substructure_name, smarts in substructures.items():
        print(f"\nSearching for {substructure_name} (SMARTS: {smarts}):")
        
        for molecule_name, smiles in molecules.items():
            result = perform_substructure_search(smarts, smiles)
            
            if result["match"]:
                print(f"  {molecule_name}: Found {result['match_count']} matches at positions {result['matches']}")
            else:
                print(f"  {molecule_name}: No matches found")

def example_similarity_search():
    """Example of calculating molecular similarity."""
    print_section("Calculating Molecular Similarity")
    
    # Define reference molecule
    reference = {
        "name": "Ethanol",
        "smiles": "CCO"
    }
    
    # Define molecules to compare against
    molecules = {
        "Methanol": "CO",
        "Propanol": "CCCO",
        "Isopropanol": "CC(C)O",
        "Butanol": "CCCCO",
        "Glycerol": "C(C(CO)O)O",
        "DMSO": "CS(=O)C",
        "Acetone": "CC(=O)C"
    }
    
    print(f"Reference molecule: {reference['name']} (SMILES: {reference['smiles']})")
    print("\nSimilarity scores (sorted by Tanimoto similarity):")
    
    # Calculate similarity for each molecule
    similarity_results = []
    for name, smiles in molecules.items():
        result = calculate_similarity(reference["smiles"], smiles)
        similarity_results.append({
            "name": name,
            "smiles": smiles,
            "tanimoto": result["tanimoto"],
            "dice": result["dice"]
        })
    
    # Sort by Tanimoto similarity (descending)
    similarity_results.sort(key=lambda x: x["tanimoto"], reverse=True)
    
    # Print results
    for result in similarity_results:
        print(f"  {result['name']} (SMILES: {result['smiles']}):")
        print(f"    Tanimoto: {result['tanimoto']:.4f}")
        print(f"    Dice: {result['dice']:.4f}")

def main():
    """Run all examples."""
    print("CryoProtect Analyzer - RDKit Integration Example")
    
    # Run examples
    example_calculate_properties()
    example_generate_visualization()
    example_substructure_search()
    example_similarity_search()
    
    print("\nAll examples completed successfully!")

if __name__ == "__main__":
    main()
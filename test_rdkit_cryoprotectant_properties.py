#!/usr/bin/env python3
"""
Test script for RDKit cryoprotectant property calculations.

This script tests the calculation of cryoprotectant-specific properties
using the rdkit_wrapper, without connecting to the database. It helps
verify that the scientific formulas are working correctly.
"""

import logging
import json
import sys
import os
from pprint import pprint

# Import our rdkit wrapper
from rdkit_wrapper import (
    create_molecule_from_smiles,
    calculate_properties,
    get_rdkit_status,
    RDKIT_AVAILABLE
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Cryoprotectant-specific property calculations
# These should be kept in sync with the calculations in rdkit_property_calculator.py
CRYOPROTECTANT_PROPERTIES = [
    {
        "name": "h_bond_donor_acceptor_ratio",
        "description": "Ratio of H-bond donors to acceptors",
        "calculation": lambda props: (
            props.get("h_donors", 0) / props.get("h_acceptors", 1)
            if props.get("h_acceptors", 0) > 0 else props.get("h_donors", 0)
        )
    },
    {
        "name": "total_h_bonding_capacity",
        "description": "Total hydrogen bonding capacity (donors + acceptors)",
        "calculation": lambda props: props.get("h_donors", 0) + props.get("h_acceptors", 0)
    },
    {
        "name": "polarity_index",
        "description": "Ratio of TPSA to molecular surface area",
        "calculation": lambda props: props.get("tpsa", 0) / (props.get("molecular_weight", 100) ** (2/3) * 10)
    },
    {
        "name": "membrane_interaction_score",
        "description": "Estimated ability to interact with cell membranes",
        "calculation": lambda props: (
            (5.0 - abs(props.get("logp", 0) - 1.5)) / 5.0 * 
            (1.0 if props.get("molecular_weight", 0) < 400 else 0.7) *
            min(1.0, props.get("h_donors", 0) / 3.0 + props.get("h_acceptors", 0) / 5.0)
        )
    },
    {
        "name": "ice_interaction_potential",
        "description": "Estimated ability to disrupt ice crystal formation",
        "calculation": lambda props: (
            min(1.5, props.get("h_donors", 0) / 2.0) * 
            (0.5 + min(0.5, props.get("tpsa", 0) / 150.0)) *
            (1.0 if props.get("molecular_weight", 0) < 500 else 
             0.7 if props.get("molecular_weight", 0) < 1000 else 0.3)
        )
    },
    {
        "name": "vitrification_potential",
        "description": "Estimated ability to promote vitrification",
        "calculation": lambda props: (
            min(1.0, props.get("h_donors", 0) / 3.0) *
            min(1.0, props.get("h_acceptors", 0) / 6.0) *
            (1.0 if -2.0 <= props.get("logp", 0) <= 2.0 else 0.5) *
            (1.2 if props.get("ring_count", 0) > 0 else 0.8)
        )
    },
    {
        "name": "estimated_toxicity",
        "description": "Simplified toxicity estimation",
        "calculation": lambda props: (
            0.5 + 
            (0.25 if props.get("logp", 0) > 4.0 else 0.0) +
            (0.25 if props.get("molecular_weight", 0) > 500 else 0.0) -
            min(0.25, (props.get("h_donors", 0) + props.get("h_acceptors", 0)) / 20.0)
        )
    },
]

def calculate_cryoprotectant_score(all_props):
    """Calculate the overall cryoprotectant score based on all properties."""
    return (
        # Weighted combination of key factors
        all_props.get("membrane_interaction_score", 0) * 0.25 +
        all_props.get("ice_interaction_potential", 0) * 0.25 +
        all_props.get("vitrification_potential", 0) * 0.25 + 
        # Inverse of toxicity (1 - toxicity) with lower weight
        (1.0 - all_props.get("estimated_toxicity", 0.5)) * 0.15 +
        # Bonus for balanced hydrogen bonding
        (0.1 if 0.5 <= all_props.get("h_bond_donor_acceptor_ratio", 0) <= 2.0 else 0.0)
    ) * 10  # Scale to 0-10 range

def calculate_cryoprotectant_properties(smiles):
    """Calculate all properties including cryoprotectant-specific ones for a molecule."""
    # Create molecule from SMILES
    mol = create_molecule_from_smiles(smiles)
    if not mol:
        logger.warning(f"Failed to create molecule from SMILES: {smiles}")
        return None
    
    # Get basic properties
    properties = calculate_properties(mol)
    
    # Calculate cryoprotectant-specific properties
    for prop in CRYOPROTECTANT_PROPERTIES:
        prop_name = prop["name"]
        try:
            prop_value = prop["calculation"](properties)
            properties[prop_name] = round(prop_value, 3) if prop_value is not None else None
        except Exception as e:
            logger.error(f"Error calculating {prop_name}: {e}")
            properties[prop_name] = None
    
    # Calculate overall score
    properties["cryoprotectant_score"] = round(calculate_cryoprotectant_score(properties), 2)
    
    return properties

def test_known_cryoprotectants():
    """Test property calculation on known cryoprotectants."""
    known_compounds = [
        {"name": "Glycerol", "smiles": "C(C(CO)O)O"},
        {"name": "DMSO", "smiles": "CS(=O)C"},
        {"name": "Ethylene glycol", "smiles": "C(CO)O"},
        {"name": "Propylene glycol", "smiles": "CC(CO)O"},
        {"name": "Formamide", "smiles": "C(=O)N"},
        {"name": "Methanol", "smiles": "CO"},
        {"name": "Trehalose", "smiles": "C(O1)C(OC2C(C(C(C(O2)CO)O)O)O)C(C(C(C1CO)O)O)O"},
        {"name": "Sucrose", "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O"},
        {"name": "Glucose", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
    ]
    
    # Ensure output directory exists
    os.makedirs("test_results", exist_ok=True)
    
    # Print RDKit status
    rdkit_status = get_rdkit_status()
    print(f"RDKit Status: {'Available' if rdkit_status['rdkit_available'] else 'Not Available'}")
    if rdkit_status['rdkit_available']:
        print(f"RDKit Version: {rdkit_status['rdkit_version']}")
    else:
        print("Using mock implementation - results will be less accurate")
    print("\nTesting property calculations on known cryoprotectants...\n")
    
    results = []
    
    for compound in known_compounds:
        name = compound["name"]
        smiles = compound["smiles"]
        print(f"Processing {name}...")
        
        properties = calculate_cryoprotectant_properties(smiles)
        if properties:
            # Print key properties
            print(f"  Molecular Weight: {properties.get('molecular_weight', 'N/A')}")
            print(f"  LogP: {properties.get('logp', 'N/A')}")
            print(f"  TPSA: {properties.get('tpsa', 'N/A')}")
            print(f"  H-Donors: {properties.get('h_donors', 'N/A')}")
            print(f"  H-Acceptors: {properties.get('h_acceptors', 'N/A')}")
            print(f"  H-Bond Donor/Acceptor Ratio: {properties.get('h_bond_donor_acceptor_ratio', 'N/A')}")
            print(f"  Total H-Bonding Capacity: {properties.get('total_h_bonding_capacity', 'N/A')}")
            print(f"  Membrane Interaction Score: {properties.get('membrane_interaction_score', 'N/A')}")
            print(f"  Ice Interaction Potential: {properties.get('ice_interaction_potential', 'N/A')}")
            print(f"  Vitrification Potential: {properties.get('vitrification_potential', 'N/A')}")
            print(f"  Estimated Toxicity: {properties.get('estimated_toxicity', 'N/A')}")
            print(f"  Cryoprotectant Score: {properties.get('cryoprotectant_score', 'N/A')}")
            print("")
            
            # Add to results
            results.append({
                "name": name,
                "smiles": smiles,
                "properties": properties
            })
    
    # Save results to file
    result_file = "test_results/cryoprotectant_property_calculations.json"
    with open(result_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to {result_file}")
    
    # Print ordered by cryoprotectant score
    print("\nCryoprotectants ranked by score:")
    for i, compound in enumerate(sorted(results, key=lambda x: x["properties"]["cryoprotectant_score"], reverse=True)):
        print(f"{i+1}. {compound['name']}: {compound['properties']['cryoprotectant_score']}")

def main():
    test_known_cryoprotectants()

if __name__ == "__main__":
    main()
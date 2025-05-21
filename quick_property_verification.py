#!/usr/bin/env python3
"""
Quick verification of property calculations by comparing mock_rdkit_formula results
with reference values for well-known molecules.
"""

import mock_rdkit_formula

# Define reference molecules with known properties
reference_molecules = [
    {
        "name": "Ethanol",
        "smiles": "CCO",
        "reference": {
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.07,
            "logp": -0.31,  # Reference from PubChem
            "tpsa": 20.23,  # Reference from PubChem
            "h_donors": 1,
            "h_acceptors": 1,
            "rotatable_bonds": 1,
            "ring_count": 0,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 3
        }
    },
    {
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "reference": {
            "molecular_formula": "C6H6",
            "molecular_weight": 78.11,
            "logp": 2.13,  # Reference from PubChem
            "tpsa": 0.00,  # Reference from PubChem
            "h_donors": 0,
            "h_acceptors": 0,
            "rotatable_bonds": 0,
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 6
        }
    },
    {
        "name": "Aspirin",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "reference": {
            "molecular_formula": "C9H8O4",
            "molecular_weight": 180.16,
            "logp": 1.19,  # Reference from PubChem
            "tpsa": 63.60,  # Reference from PubChem
            "h_donors": 1,
            "h_acceptors": 4,
            "rotatable_bonds": 3,
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 13
        }
    },
    {
        "name": "Glucose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
        "reference": {
            "molecular_formula": "C6H12O6",
            "molecular_weight": 180.16,
            "logp": -3.24,  # Reference from PubChem
            "tpsa": 110.38,  # Reference from PubChem
            "h_donors": 5,
            "h_acceptors": 6,
            "rotatable_bonds": 1,
            "ring_count": 1,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 12
        }
    },
    {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "reference": {
            "molecular_formula": "C8H10N4O2",
            "molecular_weight": 194.19,
            "logp": -0.07,  # Reference from PubChem
            "tpsa": 58.44,  # Reference from PubChem
            "h_donors": 0,
            "h_acceptors": 6,
            "rotatable_bonds": 0,
            "ring_count": 2,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 14
        }
    }
]

def calculate_mock_properties(smiles):
    """Calculate properties using mock_rdkit_formula module."""
    return {
        "molecular_formula": mock_rdkit_formula.calculate_molecular_formula(smiles),
        "molecular_weight": mock_rdkit_formula.calculate_molecular_weight(smiles),
        "logp": mock_rdkit_formula.calculate_logp(smiles),
        "tpsa": mock_rdkit_formula.calculate_tpsa(smiles),
        "h_donors": mock_rdkit_formula.calculate_h_donors(smiles),
        "h_acceptors": mock_rdkit_formula.calculate_h_acceptors(smiles),
        "rotatable_bonds": mock_rdkit_formula.calculate_rotatable_bonds(smiles),
        "ring_count": mock_rdkit_formula.calculate_ring_count(smiles),
        "aromatic_ring_count": mock_rdkit_formula.calculate_aromatic_ring_count(smiles),
        "heavy_atom_count": mock_rdkit_formula.calculate_heavy_atom_count(smiles)
    }

def compare_properties(reference, calculated):
    """Compare reference values with calculated values and compute accuracy."""
    results = {}
    
    for prop in reference:
        if reference[prop] is not None and calculated[prop] is not None:
            ref_val = reference[prop]
            calc_val = calculated[prop]
            
            # For numeric values, calculate percent error
            if isinstance(ref_val, (int, float)) and isinstance(calc_val, (int, float)):
                if ref_val == 0:
                    # Avoid division by zero
                    accuracy = 100.0 if calc_val == 0 else 0.0
                else:
                    percent_error = abs((calc_val - ref_val) / ref_val) * 100
                    accuracy = max(0, 100 - min(percent_error, 100))
            
            # For string values (like formula), exact match only
            elif isinstance(ref_val, str) and isinstance(calc_val, str):
                accuracy = 100.0 if ref_val == calc_val else 0.0
            
            # For other types, just check equality
            else:
                accuracy = 100.0 if ref_val == calc_val else 0.0
            
            results[prop] = {
                "reference": ref_val,
                "calculated": calc_val,
                "accuracy": accuracy
            }
    
    # Calculate overall accuracy
    if results:
        overall_accuracy = sum(item["accuracy"] for item in results.values()) / len(results)
    else:
        overall_accuracy = 0.0
    
    return {
        "properties": results,
        "overall_accuracy": overall_accuracy
    }

def main():
    """Main function to run verification of property calculations."""
    print("\nQUICK PROPERTY CALCULATION VERIFICATION")
    print("=====================================\n")
    
    total_accuracy = 0.0
    property_accuracies = {}
    
    for molecule in reference_molecules:
        name = molecule["name"]
        smiles = molecule["smiles"]
        reference = molecule["reference"]
        
        print(f"\nMolecule: {name}")
        print(f"SMILES: {smiles}")
        
        # Calculate properties
        calculated = calculate_mock_properties(smiles)
        
        # Compare with reference values
        comparison = compare_properties(reference, calculated)
        
        # Print results
        print(f"\nOverall Accuracy: {comparison['overall_accuracy']:.2f}%")
        print("\nProperty Comparison (Reference vs Calculated):")
        
        # Update running totals for each property type
        for prop, result in comparison["properties"].items():
            if prop not in property_accuracies:
                property_accuracies[prop] = []
            property_accuracies[prop].append(result["accuracy"])
            
            # Print property comparison
            print(f"  {prop}: {result['reference']} vs {result['calculated']} " +
                  f"(Accuracy: {result['accuracy']:.2f}%)")
        
        total_accuracy += comparison["overall_accuracy"]
    
    # Calculate and print overall averages
    average_accuracy = total_accuracy / len(reference_molecules)
    print("\n\nSUMMARY")
    print("=======")
    print(f"Average Overall Accuracy: {average_accuracy:.2f}%")
    
    print("\nAverage Property Accuracies:")
    for prop, accuracies in property_accuracies.items():
        avg = sum(accuracies) / len(accuracies)
        print(f"  {prop}: {avg:.2f}%")

if __name__ == "__main__":
    main()
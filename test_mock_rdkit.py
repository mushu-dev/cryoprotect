#!/usr/bin/env python3
"""
Script to test the mock RDKit functionality.
"""
import os
import sys
import json

print("Python version:", sys.version)

# Try to import RDKit - use mock if not available
try:
    import rdkit
    print("Using real RDKit version:", rdkit.__version__)
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    IS_MOCK = False
except ImportError:
    print("RDKit not available, using mock...")
    import mock_rdkit
    mock_rdkit.create_mock_rdkit()
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    print("Using mock RDKit")
    IS_MOCK = True

# Test basic RDKit functionality
def test_rdkit():
    """Test basic RDKit functionality with a simple molecule."""
    results = {}
    
    # Create a test molecule
    smiles = "CCO"  # Ethanol
    mol = Chem.MolFromSmiles(smiles)
    
    if mol:
        print(f"Successfully created molecule from SMILES: {smiles}")
        results["molecule_created"] = True
        
        # Calculate some basic properties
        mw = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hba = Lipinski.NumHAcceptors(mol)
        hbd = Lipinski.NumHDonors(mol)
        
        results["properties"] = {
            "molecular_weight": mw,
            "logP": logp,
            "h_bond_acceptors": hba,
            "h_bond_donors": hbd
        }
        
        print(f"Molecular weight: {mw}")
        print(f"LogP: {logp}")
        print(f"H-bond acceptors: {hba}")
        print(f"H-bond donors: {hbd}")
    else:
        print(f"Failed to create molecule from SMILES: {smiles}")
        results["molecule_created"] = False
    
    return results

# Main execution
if __name__ == "__main__":
    print("\n=== Testing RDKit Functionality ===\n")
    results = test_rdkit()
    
    print("\n=== Test Results ===\n")
    print(json.dumps(results, indent=2))
    
    print("\n=== Test Complete ===\n")
    
    if IS_MOCK:
        print("Note: Using mock RDKit - results are simulated")
    else:
        print("Note: Using real RDKit - results are accurate")
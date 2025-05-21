#!/usr/bin/env python3
"""
Simple test for RDKit wrapper functionality
"""
print("Testing RDKit Wrapper functionality...")

try:
    import rdkit_wrapper

    # Check if RDKit is available
    status = rdkit_wrapper.get_rdkit_status()
    print(f"RDKit Available: {status['rdkit_available']}")
    print(f"RDKit Version: {status['rdkit_version']}")
    print(f"Visualization Available: {status['visualization_available']}")

    # Test with a simple molecule
    smiles = "CCO"  # Ethanol
    print(f"\nTesting with Ethanol ({smiles}):")
    
    # Create molecule
    mol = rdkit_wrapper.create_molecule_from_smiles(smiles)
    print(f"Molecule created: {mol is not None}")
    
    # Calculate properties
    props = rdkit_wrapper.calculate_properties(smiles)
    print("Properties:")
    for key, value in props.items():
        if key not in ['rdkit_available', 'rdkit_version', 'calculation_method']:
            print(f"  {key}: {value}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
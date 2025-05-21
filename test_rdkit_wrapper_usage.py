import rdkit_wrapper

# Test with a common molecule
smiles = "CCO"  # Ethanol
properties = rdkit_wrapper.calculate_properties(smiles)

print(f"Molecule: Ethanol ({smiles})")
print(f"Properties:")
for prop, value in properties.items():
    print(f"  {prop}: {value}")

# Test if RDKit is available
if rdkit_wrapper.RDKIT_AVAILABLE:
    print("\nRDKit is available. Testing advanced features:")
    
    # Generate fingerprint
    fp = rdkit_wrapper.generate_fingerprint(smiles)
    print(f"  Fingerprint generated: {fp is not None}")
    
    # Generate SVG
    svg = rdkit_wrapper.generate_molecule_svg(smiles)
    print(f"  SVG generated: {svg is not None and len(svg) > 0}")
    
    # Check wrapper version info
    status = rdkit_wrapper.get_rdkit_status()
    print(f"  RDKit Version: {status[rdkit_version]}")
else:
    print("\nRDKit is not available. Using mock implementation.")

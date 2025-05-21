#!/usr/bin/env python3
"""
Test RDKit wrapper in the container environment.
This script demonstrates the use of the rdkit_wrapper module to perform
various cheminformatics operations.
"""
import sys
import json
import os

try:
    import rdkit_wrapper
    
    # Get RDKit status
    status = rdkit_wrapper.get_rdkit_status()
    print(f"RDKit Status:")
    print(f"  Available: {status['rdkit_available']}")
    print(f"  Version: {status['rdkit_version']}")
    print(f"  Visualization: {status['visualization_available']}")
    print(f"  Python Version: {status['python_version']}")
    
    # Test with complex molecules
    test_molecules = [
        ("Cryoprotectant-DMSO", "CS(=O)C"),  # Dimethyl sulfoxide 
        ("Cryoprotectant-Glycerol", "C(C(CO)O)O"),  # Glycerol
        ("Cryoprotectant-Sucrose", "C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O"),  # Sucrose
    ]
    
    results = {"molecules": {}}
    
    for name, smiles in test_molecules:
        print(f"\nTesting {name} ({smiles}):")
        
        # Create molecule
        mol = rdkit_wrapper.create_molecule_from_smiles(smiles)
        
        # Calculate properties
        props = rdkit_wrapper.calculate_properties(smiles)
        print(f"  Molecular Weight: {props.get('molecular_weight', 'N/A')}")
        print(f"  Molecular Formula: {props.get('molecular_formula', 'N/A')}")
        print(f"  LogP: {props.get('logp', 'N/A')}")
        print(f"  TPSA: {props.get('tpsa', 'N/A')}")
        print(f"  H-Donors: {props.get('h_donors', 'N/A')}")
        print(f"  H-Acceptors: {props.get('h_acceptors', 'N/A')}")
        
        # Generate fingerprint if available
        fp = None
        if rdkit_wrapper.RDKIT_AVAILABLE:
            fp = rdkit_wrapper.generate_fingerprint(smiles)
            print(f"  Fingerprint generated: {fp is not None}")
        
        # Generate SVG if available
        svg = None
        if rdkit_wrapper.VISUALIZATION_AVAILABLE:
            svg = rdkit_wrapper.generate_molecule_svg(smiles)
            print(f"  SVG generated: {svg is not None and len(svg) > 100}")
            
            # Save SVG to file for verification
            if svg:
                svg_name = name.lower().replace("-", "_")
                with open(f"{svg_name}.svg", "w") as f:
                    f.write(svg)
                print(f"  SVG saved to {svg_name}.svg")
        
        # Test 3D coordinates if available
        mol_3d = None
        if rdkit_wrapper.RDKIT_AVAILABLE:
            mol_3d = rdkit_wrapper.generate_molecule_3d_coordinates(smiles)
            print(f"  3D coordinates generated: {mol_3d is not None}")
        
        # Store results
        results["molecules"][name] = {
            "smiles": smiles,
            "properties": props,
            "svg_generated": svg is not None,
            "fingerprint_generated": fp is not None,
            "3d_generated": mol_3d is not None
        }
    
    # If molecules have SVGs, try highlighting substructures
    if rdkit_wrapper.VISUALIZATION_AVAILABLE:
        print("\nTesting substructure highlighting:")
        # Highlight hydroxyl groups in glycerol
        glycerol = rdkit_wrapper.create_molecule_from_smiles("C(C(CO)O)O")
        if glycerol:
            # Search for hydroxyl groups
            search_result = rdkit_wrapper.perform_substructure_search("[OH]", glycerol)
            if search_result["match"]:
                print(f"  Found {len(search_result['match_atoms'])} hydroxyl groups in glycerol")
                # Generate SVG with highlighting
                svg = rdkit_wrapper.generate_molecule_svg(
                    glycerol, 
                    highlight_atoms=search_result['match_atoms']
                )
                if svg:
                    with open("glycerol_highlighted.svg", "w") as f:
                        f.write(svg)
                    print("  Saved highlighted glycerol to glycerol_highlighted.svg")
    
    # Test similarity calculation
    if rdkit_wrapper.RDKIT_AVAILABLE:
        print("\nTesting molecular similarity:")
        similarity = rdkit_wrapper.calculate_similarity("CS(=O)C", "CS(=O)CC")  # DMSO vs. similar molecule
        print(f"  Tanimoto similarity: {similarity['tanimoto']:.4f}")
    
    # Save results
    with open("container_test_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to container_test_results.json")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
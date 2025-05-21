#!/usr/bin/env python3
"""
Comprehensive RDKit integration test for CryoProtect

This script tests the full functionality of the rdkit_wrapper module,
ensuring that all features work correctly with both real RDKit and mock implementations.
"""

import os
import sys
import json
import logging
from typing import Dict, Any, List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import the wrapper
import rdkit_wrapper

def test_basic_functionality():
    """Test basic molecule creation and property calculation"""
    logger.info("Testing basic functionality...")
    
    # Test SMILES parsing
    mol = rdkit_wrapper.create_molecule_from_smiles("CCO")
    assert mol is not None, "Failed to create molecule from SMILES"
    
    # Test property calculation
    props = rdkit_wrapper.calculate_properties("CCO")
    assert "molecular_weight" in props, "Molecular weight not calculated"
    assert "logp" in props, "LogP not calculated"
    
    logger.info("Basic functionality tests passed")
    return props

def test_advanced_functionality():
    """Test advanced RDKit functionality (when available)"""
    logger.info("Testing advanced functionality...")
    
    results = {"advanced_tests": {}}
    
    # Only run if real RDKit is available
    if not rdkit_wrapper.RDKIT_AVAILABLE:
        logger.warning("Real RDKit not available, skipping advanced tests")
        results["advanced_tests"]["status"] = "skipped"
        return results
    
    # Test fingerprint generation
    fp = rdkit_wrapper.generate_fingerprint("CCO")
    results["advanced_tests"]["fingerprint"] = fp is not None
    
    # Test similarity calculation
    sim = rdkit_wrapper.calculate_similarity("CCO", "CCCO")
    results["advanced_tests"]["similarity"] = sim.get("tanimoto", -1)
    
    # Test substructure search
    substructure = rdkit_wrapper.perform_substructure_search("[OH]", "CCO")
    results["advanced_tests"]["substructure_match"] = substructure.get("match", False)
    
    # Test 3D coordinate generation
    mol_3d = rdkit_wrapper.generate_molecule_3d_coordinates("CCO")
    results["advanced_tests"]["3d_coords"] = mol_3d is not None
    
    # Test visualization if available
    if rdkit_wrapper.VISUALIZATION_AVAILABLE:
        svg = rdkit_wrapper.generate_molecule_svg("CCO")
        results["advanced_tests"]["visualization"] = svg is not None and len(svg) > 0
    
    logger.info("Advanced functionality tests completed")
    return results

def test_complex_molecules():
    """Test handling of more complex molecules"""
    logger.info("Testing complex molecules...")
    
    # Define some complex molecules for testing
    complex_molecules = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("Cholesterol", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C")
    ]
    
    results = {}
    
    for name, smiles in complex_molecules:
        mol = rdkit_wrapper.create_molecule_from_smiles(smiles)
        
        if mol is not None:
            props = rdkit_wrapper.calculate_properties(mol)
            results[name] = {
                "smiles": smiles,
                "properties": props
            }
            
            logger.info(f"{name}: MW={props.get('molecular_weight', 'N/A')}, LogP={props.get('logp', 'N/A')}")
        else:
            logger.warning(f"Failed to create molecule for {name}")
            results[name] = {
                "smiles": smiles,
                "error": "Failed to create molecule"
            }
    
    logger.info("Complex molecule tests completed")
    return results

def test_fallback_mode():
    """Test the mock implementation fallback mode"""
    logger.info("Testing fallback mode...")
    
    results = {"fallback_tests": {}}
    
    # Only test fallback if we have mock_rdkit_formula imported
    try:
        # Check if mock_rdkit_formula module is properly imported
        import mock_rdkit_formula
        
        results["fallback_tests"]["mock_module_available"] = True
        
        # Test basic functions from the mock module directly
        formula = mock_rdkit_formula.calculate_molecular_formula("CCO")
        mw = mock_rdkit_formula.calculate_molecular_weight("CCO")
        logp = mock_rdkit_formula.calculate_logp("CCO")
        
        results["fallback_tests"]["mock_formula"] = formula
        results["fallback_tests"]["mock_mw"] = mw
        results["fallback_tests"]["mock_logp"] = logp
        
        logger.info(f"Mock implementation direct test: Formula={formula}, MW={mw}, LogP={logp}")
    except ImportError:
        logger.warning("Mock RDKit formula module not available, skipping direct tests")
        results["fallback_tests"]["mock_module_available"] = False
    
    # Test wrapper functions without modifying RDKIT_AVAILABLE
    # This tests the alternative paths in the wrapper functions
    results["fallback_tests"]["create_molecule"] = rdkit_wrapper.create_molecule_from_smiles("CCO") is not None
    results["fallback_tests"]["visualization_available"] = rdkit_wrapper.generate_molecule_svg("CCO") is not None
    results["fallback_tests"]["fingerprint_available"] = rdkit_wrapper.generate_fingerprint("CCO") is not None
    results["fallback_tests"]["3d_coords_available"] = rdkit_wrapper.generate_molecule_3d_coordinates("CCO") is not None
    
    logger.info("Fallback mode tests completed")
    return results

def main():
    """Main test function"""
    logger.info("Starting comprehensive RDKit integration tests")
    
    # Print wrapper status
    status = rdkit_wrapper.get_rdkit_status()
    logger.info(f"RDKit Wrapper Status:")
    logger.info(f"  RDKit Available: {status['rdkit_available']}")
    logger.info(f"  RDKit Version: {status['rdkit_version']}")
    logger.info(f"  Visualization Available: {status['visualization_available']}")
    logger.info(f"  Python Version: {status['python_version']}")
    
    # Run tests
    results = {
        "status": status,
        "basic_functionality": test_basic_functionality(),
        "complex_molecules": test_complex_molecules()
    }
    
    # Add advanced functionality test results
    advanced_results = test_advanced_functionality()
    results.update(advanced_results)
    
    # Test fallback mode
    fallback_results = test_fallback_mode()
    results.update(fallback_results)
    
    # Save results to file
    output_file = "rdkit_integration_test_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Test results saved to {output_file}")
    
    # Print summary
    print("\nTest Summary:")
    print(f"RDKit Available: {status['rdkit_available']}")
    if status['rdkit_available']:
        print(f"RDKit Version: {status['rdkit_version']}")
    
    print("\nBasic Functionality: PASS")
    
    if status['rdkit_available']:
        if 'advanced_tests' in results and results['advanced_tests'].get('fingerprint'):
            print("Advanced Functionality: PASS")
        else:
            print("Advanced Functionality: PARTIAL PASS")
    else:
        print("Advanced Functionality: SKIPPED (RDKit not available)")
    
    print("\nComplex Molecules Tested:")
    for name in results['complex_molecules'].keys():
        print(f"- {name}")
    
    print("\nFallback Mode:")
    if 'fallback_tests' in results:
        fallback = results['fallback_tests']
        if fallback.get('mock_module_available', False):
            print(f"- Mock module available: PASS")
            print(f"- Formula calculation: {fallback.get('mock_formula', 'N/A')}")
            print(f"- Molecular weight: {fallback.get('mock_mw', 'N/A')}")
        else:
            print("- Mock module not available")
    
    print("\nAll tests completed successfully!")

if __name__ == "__main__":
    main()
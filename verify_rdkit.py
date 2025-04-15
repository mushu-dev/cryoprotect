#!/usr/bin/env python3
"""
CryoProtect Analyzer - RDKit Verification Script

This script verifies that the RDKit integration works correctly by testing
key functionality and ensuring that the Docker container works properly.
"""

import os
import sys
import logging
import json
from typing import Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_rdkit_import() -> bool:
    """Check if RDKit can be imported."""
    logger.info("Checking RDKit import...")
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, Fragments
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        # Try to import IPythonConsole, but don't fail if it's not available
        try:
            from rdkit.Chem.Draw import IPythonConsole
            logger.info("IPython is available for RDKit visualization")
        except ImportError:
            logger.warning("IPython is not available. Some visualization features may be limited.")
        
        logger.info("RDKit import successful")
        return True
    except ImportError as e:
        logger.error(f"RDKit import failed: {str(e)}")
        return False

def test_molecule_parsing() -> bool:
    """Test molecule parsing functionality."""
    logger.info("Testing molecule parsing...")
    try:
        from rdkit import Chem
        
        # Test valid SMILES
        ethanol_smiles = "CCO"
        mol = Chem.MolFromSmiles(ethanol_smiles)
        if mol is None:
            logger.error("Failed to parse valid SMILES")
            return False
        
        logger.info("Molecule parsing successful")
        return True
    except Exception as e:
        logger.error(f"Error testing molecule parsing: {str(e)}")
        return False

def test_property_calculation() -> bool:
    """Test property calculation functionality."""
    logger.info("Testing property calculation...")
    try:
        from api.rdkit_utils import calculate_all_properties
        
        # Test with ethanol
        ethanol_smiles = "CCO"
        properties = calculate_all_properties(ethanol_smiles)
        
        # Check that we have all the property categories
        required_properties = [
            "hydrogen_bonding", "logp", "tpsa", "molecular_properties",
            "functional_groups", "permeability", "smiles", "inchi", "inchi_key"
        ]
        
        for prop in required_properties:
            if prop not in properties:
                logger.error(f"Missing property: {prop}")
                return False
        
        logger.info("Property calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing property calculation: {str(e)}")
        return False

def test_visualization() -> bool:
    """Test visualization functionality."""
    logger.info("Testing visualization...")
    try:
        from api.rdkit_utils import generate_molecule_svg
        
        # Test with ethanol
        ethanol_smiles = "CCO"
        svg = generate_molecule_svg(ethanol_smiles)
        
        # Check that we got an SVG string
        if not isinstance(svg, str) or not svg:
            logger.error("Failed to generate SVG")
            return False
        
        if "<svg" not in svg or "</svg>" not in svg:
            logger.error("Generated string is not a valid SVG")
            return False
        
        logger.info("Visualization successful")
        return True
    except Exception as e:
        logger.error(f"Error testing visualization: {str(e)}")
        return False

def test_substructure_search() -> bool:
    """Test substructure search functionality."""
    logger.info("Testing substructure search...")
    try:
        from api.rdkit_utils import perform_substructure_search
        
        # Search for alcohol group in ethanol
        ethanol_smiles = "CCO"
        result = perform_substructure_search("[OH]", ethanol_smiles, "smarts", "smiles")
        
        # Should find a match
        if not result.get("match", False):
            logger.error("Failed to find substructure match")
            return False
        
        logger.info("Substructure search successful")
        return True
    except Exception as e:
        logger.error(f"Error testing substructure search: {str(e)}")
        return False

def test_similarity_calculation() -> bool:
    """Test similarity calculation functionality."""
    logger.info("Testing similarity calculation...")
    try:
        from api.rdkit_utils import calculate_similarity
        
        # Compare ethanol to itself (should be identical)
        ethanol_smiles = "CCO"
        result = calculate_similarity(ethanol_smiles, ethanol_smiles)
        
        # Tanimoto similarity should be 1.0 for identical molecules
        if result.get("tanimoto", 0.0) != 1.0:
            logger.error(f"Unexpected similarity value: {result.get('tanimoto', 0.0)}")
            return False
        
        logger.info("Similarity calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing similarity calculation: {str(e)}")
        return False

def run_all_tests() -> Dict[str, bool]:
    """Run all tests and return results."""
    results = {
        "rdkit_import": check_rdkit_import(),
        "molecule_parsing": test_molecule_parsing(),
        "property_calculation": test_property_calculation(),
        "visualization": test_visualization(),
        "substructure_search": test_substructure_search(),
        "similarity_calculation": test_similarity_calculation()
    }
    
    return results

def main():
    """Main function."""
    logger.info("Starting RDKit verification...")
    
    # Add the parent directory to the path so we can import the api package
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    
    # Run all tests
    results = run_all_tests()
    
    # Print results
    logger.info("Verification results:")
    all_passed = True
    for test_name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        logger.info(f"  {test_name}: {status}")
        if not passed:
            all_passed = False
    
    # Print overall result
    if all_passed:
        logger.info("All tests passed! RDKit integration is working correctly.")
    else:
        logger.error("Some tests failed. Please check the logs for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Integration test for RDKit fallback mechanism in the CryoProtect application.
This test verifies that the application can use either real RDKit or the mock
implementation when needed.
"""

import os
import sys
import logging
import importlib
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_rdkit_availability():
    """Check if real RDKit is available"""
    try:
        import rdkit
        logger.info(f"Real RDKit is available (version: {rdkit.__version__})")
        return True
    except ImportError:
        logger.info("Real RDKit is not available")
        return False

def remove_rdkit_from_path():
    """Remove real RDKit from Python path"""
    # Clear RDKit modules if they're already imported
    for name in list(sys.modules.keys()):
        if name.startswith('rdkit'):
            logger.info(f"Removing {name} from sys.modules")
            del sys.modules[name]
    
    # Try to find and remove RDKit from Python path
    import site
    site_packages = site.getsitepackages()
    
    for path in list(sys.path):
        for site_path in site_packages:
            if path.startswith(site_path) and 'rdkit' in path.lower():
                logger.info(f"Removing {path} from sys.path")
                sys.path.remove(path)
    
    logger.info("Removed RDKit from Python path")

def setup_mock_rdkit():
    """Set up mock RDKit"""
    import enhanced_mock_rdkit
    mock_dir = enhanced_mock_rdkit.create_mock_rdkit()
    
    # Ensure it's at the beginning of Python path
    if mock_dir in sys.path:
        sys.path.remove(mock_dir)
    sys.path.insert(0, mock_dir)
    
    logger.info(f"Mock RDKit set up at {mock_dir}")
    return mock_dir

def test_app_behaviors():
    """Test application behavior with property calculations"""
    try:
        # Import application modules
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from api import rdkit_utils
        
        # Test basic property calculation
        smiles = "CCO"  # Ethanol
        logger.info(f"Testing property calculation for {smiles}...")
        
        properties = rdkit_utils.calculate_all_properties(smiles)
        logger.info("Property calculation result:")
        
        # Log a few selected properties
        logger.info(f"Molecular Weight: {properties.get('molecular_properties', {}).get('molecular_weight', 'N/A')}")
        logger.info(f"LogP: {properties.get('logp', 'N/A')}")
        logger.info(f"TPSA: {properties.get('tpsa', 'N/A')}")
        logger.info(f"H-Bond Donors: {properties.get('hydrogen_bonding', {}).get('donors', 'N/A')}")
        logger.info(f"H-Bond Acceptors: {properties.get('hydrogen_bonding', {}).get('acceptors', 'N/A')}")
        
        # Check that the key properties are present
        required_properties = [
            "hydrogen_bonding", "logp", "tpsa", "molecular_properties",
            "functional_groups", "permeability", "smiles", "inchi", "inchi_key"
        ]
        
        missing_properties = []
        for prop in required_properties:
            if prop not in properties:
                missing_properties.append(prop)
        
        if missing_properties:
            logger.error(f"Missing required properties: {', '.join(missing_properties)}")
            return False
        
        logger.info("All required properties are present")
        
        # Test visualization
        svg = rdkit_utils.generate_molecule_svg(smiles)
        if svg and "<svg" in svg:
            logger.info("Molecule visualization successful")
        else:
            logger.error("Molecule visualization failed")
            return False
        
        return True
    except Exception as e:
        logger.error(f"Error testing app behaviors: {str(e)}")
        return False

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Test RDKit fallback mechanisms")
    parser.add_argument("--force-mock", action="store_true", 
                        help="Force use of mock RDKit even if real RDKit is available")
    args = parser.parse_args()
    
    logger.info("Starting RDKit fallback integration test...")
    
    # Check if real RDKit is available
    has_real_rdkit = check_rdkit_availability()
    
    # Determine if we should force mock RDKit
    if args.force_mock:
        logger.info("Forcing use of mock RDKit as requested")
        if has_real_rdkit:
            remove_rdkit_from_path()
        mock_dir = setup_mock_rdkit()
        rdkit_type = "mock"
    else:
        if has_real_rdkit:
            rdkit_type = "real"
        else:
            logger.info("Using mock RDKit since real RDKit is not available")
            mock_dir = setup_mock_rdkit()
            rdkit_type = "mock"
    
    # Run application tests
    logger.info(f"Testing application with {rdkit_type} RDKit...")
    test_result = test_app_behaviors()
    
    if test_result:
        logger.info(f"Application test with {rdkit_type} RDKit: PASSED")
    else:
        logger.error(f"Application test with {rdkit_type} RDKit: FAILED")
    
    return test_result

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("RDKit fallback integration test passed!")
        sys.exit(0)
    else:
        logger.error("RDKit fallback integration test failed")
        sys.exit(1)
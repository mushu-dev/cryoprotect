#!/usr/bin/env python3
"""
Test script to verify RDKit fallback mechanisms when RDKit is not available
"""

import sys
import os
import logging
import importlib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_mock_rdkit():
    """Test the mock_rdkit implementation"""
    try:
        import mock_rdkit
        logger.info("Successfully imported mock_rdkit module")
        
        # Setup mock RDKit structure
        mock_dir = mock_rdkit.create_mock_rdkit()
        logger.info(f"Created mock RDKit structure at {mock_dir}")
        
        # Now try to import the mock RDKit
        sys.path.insert(0, mock_dir)
        logger.info(f"Added {mock_dir} to Python path")
        
        # Clear any cached imports
        if 'rdkit' in sys.modules:
            del sys.modules['rdkit']
        if 'rdkit.Chem' in sys.modules:
            del sys.modules['rdkit.Chem']
        
        # Try importing RDKit (should get the mock version)
        import rdkit
        from rdkit import Chem
        logger.info("Successfully imported mock rdkit module")
        
        # Test basic functionality
        smiles = "CCO"  # Ethanol
        mol = Chem.MolFromSmiles(smiles)
        logger.info(f"Parsed SMILES: {smiles}")
        
        # Test descriptors
        mw = Chem.Descriptors.MolWt(mol)
        logp = Chem.Descriptors.MolLogP(mol)
        h_donors = Chem.Lipinski.NumHDonors(mol)
        h_acceptors = Chem.Lipinski.NumHAcceptors(mol)
        
        logger.info("Mock property calculation:")
        logger.info(f"Molecular Weight: {mw}")  # Should be 100.0 (mock value)
        logger.info(f"LogP: {logp}")  # Should be 1.0 (mock value)
        logger.info(f"H-Bond Donors: {h_donors}")  # Should be 2 (mock value)
        logger.info(f"H-Bond Acceptors: {h_acceptors}")  # Should be 2 (mock value)
        
        return True
    except ImportError as e:
        logger.error(f"Failed to import mock_rdkit: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error testing mock_rdkit: {str(e)}")
        return False

def test_app_with_mock():
    """Test the application code with mock RDKit"""
    try:
        # First, we import mock_rdkit to create the structure
        import mock_rdkit
        mock_dir = mock_rdkit.create_mock_rdkit()
        sys.path.insert(0, mock_dir)
        
        # Clear any cached imports
        if 'rdkit' in sys.modules:
            del sys.modules['rdkit']
        if 'rdkit.Chem' in sys.modules:
            del sys.modules['rdkit.Chem']
        
        # Now try importing the application's rdkit_utils module
        sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
        
        # We need to import the module fresh
        spec = importlib.util.find_spec('api.rdkit_utils')
        if spec is None:
            logger.error("Could not find api.rdkit_utils module")
            return False
        
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        logger.info("Successfully imported api.rdkit_utils with mock RDKit")
        
        # Try calculating properties for a molecule
        try:
            properties = module.calculate_all_properties("CCO")
            logger.info("Successfully calculated properties using mock RDKit")
            logger.info(f"Properties: {properties}")
            return True
        except Exception as e:
            logger.error(f"Error calculating properties: {str(e)}")
            return False
            
    except ImportError as e:
        logger.error(f"Import error: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error testing application with mock RDKit: {str(e)}")
        return False

def main():
    logger.info("Testing RDKit fallback mechanisms...")
    
    # First, check if real RDKit is available
    try:
        import rdkit
        logger.warning(f"Real RDKit is available (version: {rdkit.__version__})")
        logger.warning("This test is meant to verify fallback when RDKit is not available.")
        logger.warning("For proper testing, run this script without RDKit in the environment.")
    except ImportError:
        logger.info("RDKit is not available, which is good for this test.")
    
    # Test mock_rdkit module
    mock_test_result = test_mock_rdkit()
    if mock_test_result:
        logger.info("Mock RDKit test: PASSED")
    else:
        logger.error("Mock RDKit test: FAILED")
    
    # Test application with mock
    app_test_result = test_app_with_mock()
    if app_test_result:
        logger.info("Application with mock RDKit test: PASSED")
    else:
        logger.error("Application with mock RDKit test: FAILED")
    
    return mock_test_result and app_test_result

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("All RDKit fallback tests passed!")
        sys.exit(0)
    else:
        logger.error("RDKit fallback test failed")
        sys.exit(1)
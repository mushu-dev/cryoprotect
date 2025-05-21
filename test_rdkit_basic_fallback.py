#!/usr/bin/env python3
"""
Basic test for RDKit fallback mechanism.
This test verifies that the code can handle switching between real RDKit and mock RDKit.
"""

import os
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Flag to track if we're using real or mock RDKit
RDKIT_AVAILABLE = False

def try_import_rdkit():
    """Try to import RDKit, falling back to mock if needed"""
    global RDKIT_AVAILABLE
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, AllChem
        logger.info("Successfully imported real RDKit")
        # Test that it's working by creating a molecule
        mol = Chem.MolFromSmiles("CCO")
        if mol is None:
            logger.error("Failed to create molecule with real RDKit")
            raise ImportError("RDKit available but not working properly")
        
        RDKIT_AVAILABLE = True
        logger.info("Real RDKit is working")
        return True
    except ImportError:
        logger.warning("RDKit not available, trying mock_rdkit")
        try:
            # Try to import and set up mock RDKit
            import enhanced_mock_rdkit
            mock_dir = enhanced_mock_rdkit.create_mock_rdkit()
            
            # Ensure mock dir is at the beginning of Python path
            if mock_dir in sys.path:
                sys.path.remove(mock_dir)
            sys.path.insert(0, mock_dir)
            
            # Clear any cached imports
            if 'rdkit' in sys.modules:
                del sys.modules['rdkit']
            for name in list(sys.modules.keys()):
                if name.startswith('rdkit.'):
                    del sys.modules[name]
            
            # Try to import the mock RDKit
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, AllChem
            
            # Test that it's working by creating a molecule
            mol = Chem.MolFromSmiles("CCO")
            if mol is None:
                logger.error("Failed to create molecule with mock RDKit")
                return False
            
            logger.info("Successfully set up and imported mock RDKit")
            RDKIT_AVAILABLE = False  # We're using mock, not real RDKit
            return True
        except Exception as e:
            logger.error(f"Failed to set up mock RDKit: {str(e)}")
            return False

def calculate_properties(smiles):
    """Calculate basic molecular properties using RDKit or mock"""
    if not try_import_rdkit():
        return {"error": "RDKit not available"}
    
    # Import RDKit modules (should be either real or mock by now)
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Failed to parse SMILES: {smiles}"}
    
    # Calculate basic properties
    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hydrogen_bond_donors": Lipinski.NumHDonors(mol),
        "hydrogen_bond_acceptors": Lipinski.NumHAcceptors(mol)
    }
    
    return properties

def main():
    logger.info("Testing RDKit fallback mechanism...")
    
    # Get initial RDKit status
    if try_import_rdkit():
        rdkit_status = "real" if RDKIT_AVAILABLE else "mock"
        logger.info(f"Using {rdkit_status} RDKit")
    else:
        logger.error("Failed to set up any RDKit implementation")
        return False
    
    # Test property calculation
    test_molecules = [
        ("Ethanol", "CCO"),
        ("Glycerol", "C(C(CO)O)O"),
        ("DMSO", "CS(=O)C")
    ]
    
    success = True
    for name, smiles in test_molecules:
        logger.info(f"Calculating properties for {name} ({smiles})...")
        properties = calculate_properties(smiles)
        
        if "error" in properties:
            logger.error(f"Error calculating properties: {properties['error']}")
            success = False
            continue
        
        logger.info(f"Properties for {name}:")
        for prop, value in properties.items():
            logger.info(f"  {prop}: {value}")
    
    return success

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("RDKit fallback test passed!")
        sys.exit(0)
    else:
        logger.error("RDKit fallback test failed")
        sys.exit(1)
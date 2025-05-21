#!/usr/bin/env python3
"""
Test script to verify the enhanced mock RDKit implementation
"""

import sys
import os
import logging
import importlib.util

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_mock_rdkit():
    """Set up the enhanced mock RDKit"""
    try:
        # First ensure real RDKit is not imported
        if 'rdkit' in sys.modules:
            logger.info("Removing existing rdkit from sys.modules")
            del sys.modules['rdkit']
            # Also remove any submodules
            for name in list(sys.modules.keys()):
                if name.startswith('rdkit.'):
                    del sys.modules[name]
        
        # Import enhanced mock
        import enhanced_mock_rdkit
        mock_dir = enhanced_mock_rdkit.create_mock_rdkit()
        
        # Ensure mock_dir is at the beginning of sys.path
        if mock_dir in sys.path:
            sys.path.remove(mock_dir)
        sys.path.insert(0, mock_dir)
        
        logger.info(f"Mock RDKit path: {mock_dir}")
        return True
    except Exception as e:
        logger.error(f"Error setting up mock RDKit: {str(e)}")
        return False

def test_basic_functionality():
    """Test basic functionality of the enhanced mock RDKit"""
    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, AllChem
        
        logger.info(f"Mock RDKit version: {rdkit.__version__}")
        
        # Create a molecule
        smiles = "CCO"
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            logger.error("Failed to create molecule from SMILES")
            return False
            
        logger.info(f"Successfully created molecule from SMILES: {smiles}")
        logger.info(f"Returned SMILES: {Chem.MolToSmiles(mol)}")
        
        # Test property calculations
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        
        logger.info("Property calculations:")
        logger.info(f"Molecular Weight: {mw}")
        logger.info(f"LogP: {logp}")
        logger.info(f"TPSA: {tpsa}")
        logger.info(f"H-Bond Donors: {h_donors}")
        logger.info(f"H-Bond Acceptors: {h_acceptors}")
        
        # Test 3D operations
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        
        logger.info(f"Number of conformers: {mol.GetNumConformers()}")
        
        # Test fingerprints
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
        sim = DataStructs.TanimotoSimilarity(fp, fp)
        
        logger.info(f"Self-similarity (Tanimoto): {sim}")
        
        # Test drawing
        from rdkit.Chem.Draw import rdMolDraw2D
        
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        if not svg or "<svg" not in svg:
            logger.error("Failed to generate SVG")
            return False
            
        logger.info("Successfully generated SVG visualization")
        
        return True
    except Exception as e:
        logger.error(f"Error testing basic functionality: {str(e)}")
        return False

def test_with_api_module():
    """Test importing and using the api.rdkit_utils module with mock RDKit"""
    try:
        # Make sure we can import the api module
        sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
        
        # Try to import directly first
        try:
            from api import rdkit_utils
            logger.info("Successfully imported api.rdkit_utils directly")
        except ImportError as e:
            logger.warning(f"Could not import api.rdkit_utils directly: {str(e)}")
            logger.info("Attempting spec-based import...")
            
            # Try spec-based import
            spec = importlib.util.find_spec('api.rdkit_utils')
            if spec is None:
                logger.error("Could not find api.rdkit_utils module spec")
                return False
                
            rdkit_utils = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(rdkit_utils)
            logger.info("Successfully imported api.rdkit_utils via spec")
        
        # Test calculating properties
        smiles = "CCO"  # Ethanol
        
        logger.info(f"Calculating properties for {smiles}...")
        properties = rdkit_utils.calculate_all_properties(smiles)
        
        if not properties or "error" in properties:
            logger.error(f"Error calculating properties: {properties.get('error', 'unknown error')}")
            return False
            
        logger.info("Successfully calculated properties:")
        logger.info(f"H-Bonding: {properties.get('hydrogen_bonding')}")
        logger.info(f"LogP: {properties.get('logp')}")
        logger.info(f"TPSA: {properties.get('tpsa')}")
        
        # Test generating SVG
        logger.info("Testing visualization...")
        svg = rdkit_utils.generate_molecule_svg(smiles)
        
        if not svg or "<svg" not in svg:
            logger.error("Failed to generate SVG visualization")
            return False
            
        logger.info("Successfully generated SVG visualization")
        
        # Test substructure search
        logger.info("Testing substructure search...")
        result = rdkit_utils.perform_substructure_search("[OH]", smiles)
        
        if not result.get("match", False):
            logger.error("Substructure search failed to find expected match")
            return False
            
        logger.info(f"Substructure search result: {result}")
        
        # Test similarity calculation
        logger.info("Testing similarity calculation...")
        sim_result = rdkit_utils.calculate_similarity(smiles, smiles)
        
        if sim_result.get("tanimoto", 0.0) < 0.9:
            logger.error(f"Unexpected similarity value: {sim_result.get('tanimoto', 0.0)}")
            return False
            
        logger.info(f"Similarity calculation result: {sim_result}")
        
        return True
    except Exception as e:
        logger.error(f"Error testing with api module: {str(e)}")
        return False

def main():
    logger.info("Testing enhanced mock RDKit implementation...")
    
    # Set up the mock RDKit
    setup_result = setup_mock_rdkit()
    if not setup_result:
        logger.error("Failed to set up mock RDKit")
        return False
    
    # Test basic functionality
    basic_result = test_basic_functionality()
    if basic_result:
        logger.info("Basic functionality test: PASSED")
    else:
        logger.error("Basic functionality test: FAILED")
    
    # Test with api module
    api_result = test_with_api_module()
    if api_result:
        logger.info("API module test: PASSED")
    else:
        logger.error("API module test: FAILED")
    
    return basic_result and api_result

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("All enhanced mock RDKit tests passed!")
        sys.exit(0)
    else:
        logger.error("Enhanced mock RDKit test failed")
        sys.exit(1)
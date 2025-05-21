#!/usr/bin/env python3
"""
Simple RDKit test script to verify basic property calculation
"""

import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    logger.info("Testing RDKit basic functionality...")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        
        logger.info("RDKit successfully imported")
        
        # Test molecule (Ethanol)
        smiles = "CCO"
        logger.info(f"Testing with molecule: Ethanol ({smiles})")
        
        # Parse molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error("Failed to parse SMILES")
            return False
        
        # Calculate basic properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        
        # Display results
        logger.info("Basic property calculation successful:")
        logger.info(f"Molecular Weight: {mw:.2f}")
        logger.info(f"LogP: {logp:.2f}")
        logger.info(f"H-Bond Donors: {h_donors}")
        logger.info(f"H-Bond Acceptors: {h_acceptors}")
        
        # Test generating SVG
        try:
            from rdkit.Chem import Draw
            from rdkit.Chem.Draw import rdMolDraw2D
            
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            if svg and "<svg" in svg:
                logger.info("Molecule visualization successful")
            else:
                logger.error("Failed to generate SVG visualization")
                
        except ImportError as e:
            logger.error(f"Visualization modules not available: {str(e)}")
        
        return True
        
    except ImportError as e:
        logger.error(f"RDKit import failed: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error during RDKit testing: {str(e)}")
        return False

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("All basic RDKit functionality tests passed!")
        sys.exit(0)
    else:
        logger.error("RDKit functionality test failed")
        sys.exit(1)
#!/usr/bin/env python3
"""
Simple test script to verify self-similarity in RDKit.
"""

import logging
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_self_similarity():
    """Test that a molecule compared to itself has a similarity of 1.0."""
    # Test with glycerol
    glycerol_smiles = "C(C(CO)O)O"
    
    # Parse molecule
    mol = Chem.MolFromSmiles(glycerol_smiles)
    if mol is None:
        logger.error("Failed to parse glycerol SMILES")
        return False
    
    # Get canonical SMILES
    canonical_smiles = Chem.MolToSmiles(mol)
    logger.info(f"Original SMILES: {glycerol_smiles}")
    logger.info(f"Canonical SMILES: {canonical_smiles}")
    
    # Generate fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    
    # Calculate self-similarity
    tanimoto = DataStructs.TanimotoSimilarity(fp, fp)
    logger.info(f"Self-similarity (Tanimoto): {tanimoto}")
    
    # Verify self-similarity is 1.0
    if tanimoto == 1.0:
        logger.info("PASS: Self-similarity is 1.0 as expected")
        return True
    else:
        logger.error(f"FAIL: Self-similarity is {tanimoto}, expected 1.0")
        return False

if __name__ == "__main__":
    logger.info("Testing self-similarity in RDKit...")
    result = test_self_similarity()
    if result:
        logger.info("All tests passed!")
    else:
        logger.error("Tests failed!")
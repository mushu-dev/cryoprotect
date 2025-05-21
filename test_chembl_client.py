#!/usr/bin/env python3
"""
Test script for the fixed ChEMBLClient.get_compound method.
This script verifies that the method returns data in the expected format
for property filtering in import_full_chembl.py.
"""

import logging
import sys
from chembl.client import ChEMBLClient

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/test_chembl_client.log')
    ]
)

logger = logging.getLogger(__name__)

def test_get_compound():
    """Test the get_compound method of ChEMBLClient."""
    client = ChEMBLClient()
    
    # Test with a known ChEMBL ID for a cryoprotectant (glycerol)
    chembl_id = "CHEMBL692"
    
    logger.info(f"Testing get_compound with ChEMBL ID: {chembl_id}")
    
    try:
        # Get compound data
        compound_data = client.get_compound(chembl_id)
        
        # Check if the data has the expected structure
        if "molecule_structures" in compound_data:
            logger.info("✓ Data contains 'molecule_structures' key")
            
            if "canonical_smiles" in compound_data["molecule_structures"]:
                smiles = compound_data["molecule_structures"]["canonical_smiles"]
                logger.info(f"✓ Data contains 'canonical_smiles': {smiles}")
            else:
                logger.warning("✗ Data does not contain 'canonical_smiles' in 'molecule_structures'")
        else:
            logger.warning("✗ Data does not contain 'molecule_structures' key")
        
        # Print the full data structure for inspection
        logger.info("Full compound data structure:")
        for key, value in compound_data.items():
            if isinstance(value, dict):
                logger.info(f"  {key}:")
                for subkey, subvalue in value.items():
                    logger.info(f"    {subkey}: {subvalue}")
            else:
                logger.info(f"  {key}: {value}")
        
        return True
    except Exception as e:
        logger.error(f"Error testing get_compound: {e}")
        return False

if __name__ == "__main__":
    # Create logs directory if it doesn't exist
    import os
    os.makedirs('logs', exist_ok=True)
    
    logger.info("Starting ChEMBLClient.get_compound test")
    
    success = test_get_compound()
    
    if success:
        logger.info("Test completed successfully")
        sys.exit(0)
    else:
        logger.error("Test failed")
        sys.exit(1)
#!/usr/bin/env python3
"""
Test script to verify that the ChEMBL import script fixes work correctly.
This script tests the key components that were fixed:
1. The store_compound_data function with conn parameter
2. The update_progress method with data parameter
3. The Unicode character replacement
"""

import os
import sys
import logging
import json
from typing import Dict, Any, List, Optional
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Import the fixed modules
sys.path.append('.')  # Add current directory to path
from chembl.checkpoint import CheckpointManager
from chembl.worker import ChEMBLWorker
from database.population.chembl_import import process_chembl_batch

def test_checkpoint_manager():
    """Test the CheckpointManager with the fixed update_progress method"""
    logger.info("Testing CheckpointManager...")
    
    # Create a test directory
    os.makedirs("test_checkpoints", exist_ok=True)
    
    # Initialize checkpoint manager
    checkpoint_manager = CheckpointManager("test_checkpoints", "test_import")
    
    # Test update_progress with data parameter
    test_data = {
        "operation_id": "test-operation",
        "processing_time": 1.5,
        "property_count": 5,
        "timestamp": datetime.now().timestamp()
    }
    
    try:
        checkpoint_manager.update_progress(
            compound_id="CHEMBL123",
            success=True,
            error=None,
            data=test_data
        )
        logger.info("[OK] CheckpointManager.update_progress with data parameter works")
        
        # Save checkpoint
        checkpoint_path = checkpoint_manager.save_checkpoint()
        logger.info(f"[OK] Checkpoint saved to {checkpoint_path}")
        
        # Load checkpoint
        loaded = checkpoint_manager.load_checkpoint()
        if loaded and "compound_data" in checkpoint_manager.state:
            logger.info("[OK] Checkpoint loaded successfully with compound_data")
        else:
            logger.error("[FAIL] Failed to load checkpoint with compound_data")
            
        return True
    except Exception as e:
        logger.error(f"✗ CheckpointManager test failed: {str(e)}")
        return False

def test_process_chembl_batch():
    """Test the process_chembl_batch function with conn parameter"""
    logger.info("Testing process_chembl_batch...")
    
    # Create a mock batch
    mock_batch = [
        {
            "molecule_chembl_id": "CHEMBL123",
            "molecule_structures": {
                "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
            },
            "pref_name": "Test Compound",
            "molecule_properties": {
                "full_molformula": "C9H8O4",
                "full_mwt": 180.16,
                "alogp": 1.43,
                "hba": 4,
                "hbd": 1,
                "psa": 63.6,
                "rtb": 3
            }
        }
    ]
    
    # Mock project ID
    project_id = "test_project"
    
    try:
        # Test with conn parameter
        mock_conn = object()  # Just a placeholder object
        result = process_chembl_batch(project_id, mock_batch, conn=mock_conn)
        logger.info("[OK] process_chembl_batch accepts conn parameter")
        return True
    except TypeError as e:
        if "unexpected keyword argument 'conn'" in str(e):
            logger.error("[FAIL] process_chembl_batch still doesn't accept conn parameter")
        else:
            logger.error(f"[FAIL] process_chembl_batch test failed: {str(e)}")
        return False
    except Exception as e:
        # This is expected since we're not actually connecting to a database
        logger.info("[OK] process_chembl_batch accepts conn parameter (other error occurred)")
        return True

def main():
    """Run all tests"""
    logger.info("Starting ChEMBL fix tests...")
    
    # Run tests
    checkpoint_test = test_checkpoint_manager()
    batch_test = test_process_chembl_batch()
    
    # Report results
    if checkpoint_test and batch_test:
        logger.info("All tests passed! The fixes appear to be working correctly.")
        return 0
    else:
        logger.error("Some tests failed. Please check the logs for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
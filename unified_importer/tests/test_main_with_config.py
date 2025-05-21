"""
Test script for verifying the main MolecularImporter with configuration.
"""

import os
import sys
import logging
from pathlib import Path

# Add the parent directory to the path to import the main module
parent_dir = str(Path(__file__).parent.parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from unified_importer.main import MolecularImporter

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('test_main')


def test_importer_initialization():
    """Test that the MolecularImporter can be properly initialized with configuration."""
    # Get the path to the test config file
    config_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
                              'config', 'test_config.json')
    
    logger.info(f"Loading configuration from {config_file}")
    
    try:
        # Initialize the importer
        importer = MolecularImporter(
            config_file=config_file,
            db_url="mock://localhost",
            db_key="mock_key",
            logger=logger
        )
        
        logger.info("Successfully initialized MolecularImporter")
        
        # Log some configuration details
        logger.info(f"Database pool min size: {importer.db.pool_min_size}")
        logger.info(f"Database pool max size: {importer.db.pool_max_size}")
        logger.info(f"Batch size: {importer.config.get('batch_size')}")
        logger.info(f"Transforms available: {importer.molecule_transformer is not None}")
        logger.info(f"Property transformer available: {importer.property_transformer is not None}")

        # Log the available data sources
        logger.info(f"Available data sources: {list(importer.sources.keys())}")
        
        return True
    except Exception as e:
        logger.error(f"Failed to initialize MolecularImporter: {str(e)}")
        return False


if __name__ == "__main__":
    if test_importer_initialization():
        print("Test passed: MolecularImporter initialized successfully")
        sys.exit(0)
    else:
        print("Test failed: Could not initialize MolecularImporter")
        sys.exit(1)
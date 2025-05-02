#!/usr/bin/env python3
"""
Configuration Verification Test Runner

This script demonstrates the configuration verification functionality implemented
in ChEMBL_CryoProtectants_Supabase.py. It verifies that:

1. All required configuration variables are validated
2. API connectivity tests for Supabase and ChEMBL are performed
3. Failures are properly logged and abort execution

Usage:
    python run_config_verification_test.py

The script will exit with code 0 if verification passes, or non-zero if it fails.
"""

import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import the module to test
import ChEMBL_CryoProtectants_Supabase as chembl_supabase
from config import active_config

def main():
    """Run the configuration verification test."""
    logger.info("=" * 80)
    logger.info("CONFIGURATION VERIFICATION TEST")
    logger.info("=" * 80)
    
    try:
        # Print active configuration environment
        env_name = active_config.__class__.__name__
        logger.info(f"Active configuration environment: {env_name}")
        
        # Run the verification
        logger.info("Running configuration verification...")
        result = chembl_supabase.verify_configuration()
        
        if result:
            logger.info("=" * 80)
            logger.info("VERIFICATION SUCCESSFUL!")
            logger.info("All required configuration variables are present and valid.")
            logger.info("Supabase and ChEMBL API connectivity tests passed.")
            logger.info("=" * 80)
            return 0
        else:
            logger.error("=" * 80)
            logger.error("VERIFICATION FAILED!")
            logger.error("See logs for details.")
            logger.error("=" * 80)
            return 1
    except Exception as e:
        logger.error("=" * 80)
        logger.error(f"VERIFICATION ERROR: {str(e)}")
        logger.error("=" * 80)
        return 1

if __name__ == "__main__":
    sys.exit(main())
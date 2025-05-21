#!/usr/bin/env python3
"""
Startup script for Heroku deployment.
This script is imported before the main app to set up the mock RDKit implementation.
"""

import os
import sys
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('heroku_startup')

logger.info("Heroku startup script running...")

# Try to import RDKit, if not available, use our mock implementation
try:
    import rdkit
    logger.info("Real RDKit found, will use actual implementation")
except ImportError:
    logger.warning("RDKit not available, loading mock implementation")
    # Import the mock rdkit module
    try:
        import mock_rdkit
        logger.info("Mock RDKit implementation loaded successfully")
    except Exception as e:
        logger.error(f"Error loading mock RDKit: {e}")
        sys.exit(1)

logger.info("Heroku startup complete, continuing with application startup")
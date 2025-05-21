#!/usr/bin/env python3
"""
CryoProtect Analyzer - Debug Script

This script helps identify issues with the Flask application startup.
It provides detailed error logging and runs the app in debug mode.
"""

import os
import logging
import traceback

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Check RDKit installation
logger.info("Checking RDKit installation...")
try:
    from rdkit import Chem
    logger.info("RDKit import successful")
except ImportError as e:
    logger.error(f"RDKit import failed: {str(e)}")
    logger.error("Please install RDKit using setup_environment script")
    exit(1)

# Import Flask app with error handling
logger.info("Importing Flask application...")
try:
    from app import create_app
    app = create_app()
    logger.info("Flask app created successfully!")
except Exception as e:
    logger.error(f"Error creating Flask app: {str(e)}")
    logger.error(traceback.format_exc())
    exit(1)

# Check environment variables
logger.info("Checking environment variables...")
required_vars = ['SUPABASE_URL', 'SUPABASE_KEY']
missing_vars = [var for var in required_vars if not os.environ.get(var)]
if missing_vars:
    logger.warning(f"Missing environment variables: {', '.join(missing_vars)}")
    logger.warning("Some functionality may not work correctly")

# Try to run the app
if __name__ == '__main__':
    try:
        logger.info("Starting Flask application in debug mode...")
        port = int(os.environ.get('PORT', 5000))
        app.run(host='0.0.0.0', port=port, debug=True)
    except Exception as e:
        logger.error(f"Error running Flask app: {str(e)}")
        logger.error(traceback.format_exc())
        exit(1)
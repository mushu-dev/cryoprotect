#!/usr/bin/env python3
"""
CryoProtect - Apply Toxicity Data Optimization

This script applies the toxicity data optimization improvements to the existing
API by updating the import statements and registering the optimized resources.

It automatically:
1. Backs up the original API configuration
2. Updates imports to use the optimized toxicity resources
3. Registers the optimized endpoints in the API

Usage:
    python apply_toxicity_optimization.py
"""

import os
import sys
import argparse
import shutil
import logging
import importlib
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(f'toxicity_optimization_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger()

def backup_file(file_path):
    """Create a backup of the specified file with timestamp."""
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        shutil.copy2(file_path, backup_path)
        logger.info(f"Created backup at {backup_path}")
        return backup_path
    except Exception as e:
        logger.error(f"Failed to create backup: {str(e)}")
        raise

def update_app_imports(app_file):
    """Update app.py to import the optimized toxicity resources."""
    try:
        # Create backup
        backup_file(app_file)
        
        # Read the file
        with open(app_file, 'r') as f:
            lines = f.readlines()
        
        # Find the toxicity resources import
        for i, line in enumerate(lines):
            if "from api.toxicity_resources import" in line:
                # Replace with optimized version
                logger.info(f"Updating import at line {i+1}")
                lines[i] = "from api.toxicity_resources_optimized import register_toxicity_resources\n"
                break
        
        # Write the file
        with open(app_file, 'w') as f:
            f.writelines(lines)
        
        logger.info(f"Updated imports in {app_file}")
    except Exception as e:
        logger.error(f"Failed to update app imports: {str(e)}")
        raise

def verify_optimized_resources():
    """Verify that the optimized toxicity resources module can be imported."""
    try:
        # Try to import the module
        module = importlib.import_module('api.toxicity_resources_optimized')
        
        # Check that the register_toxicity_resources function exists
        if not hasattr(module, 'register_toxicity_resources'):
            logger.error("Module does not contain register_toxicity_resources function")
            return False
        
        # Check that required resource classes exist
        required_classes = [
            'ToxicitySummaryResource', 
            'ToxicityDataResource',
            'ToxicityLD50Resource',
            'ToxicityTox21Resource',
            'ToxicityClassificationResource',
            'SimilarToxicityResource',
            'BulkToxicityResource'
        ]
        
        for cls_name in required_classes:
            if not hasattr(module, cls_name):
                logger.error(f"Module is missing required class: {cls_name}")
                return False
        
        logger.info("Optimized toxicity resources module verified successfully")
        return True
    except ImportError as e:
        logger.error(f"Failed to import optimized toxicity resources: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error during module verification: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Apply toxicity data optimization")
    parser.add_argument("--app-file", default="app.py", help="Path to the app.py file")
    args = parser.parse_args()
    
    app_file = args.app_file
    
    try:
        logger.info("Starting toxicity data optimization")
        
        # Check if the app file exists
        if not os.path.exists(app_file):
            logger.error(f"App file not found: {app_file}")
            return 1
        
        # Verify that the optimized resources module can be imported
        if not verify_optimized_resources():
            logger.error("Verification of optimized resources failed")
            return 1
        
        # Update app imports
        update_app_imports(app_file)
        
        logger.info("Toxicity data optimization completed successfully")
        return 0
    except Exception as e:
        logger.error(f"Optimization failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
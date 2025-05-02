#!/usr/bin/env python3
"""
CryoProtect v2 - Update Population Scripts

This script updates all population scripts to use the correct plural table names
(e.g., "molecules", "molecular_properties") instead of singular names.

Usage:
    python update_population_scripts.py
"""

import os
import re
import sys
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("update_population_scripts.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Table name mappings (singular to plural)
TABLE_MAPPINGS = {
    "molecule": "molecules",
    "molecular_property": "molecular_properties",
    "mixture": "mixtures",
    "mixture_component": "mixture_components",
    "prediction": "predictions",
    "experiment": "experiments",
    "experiment_property": "experiment_properties",
    "calculation_method": "calculation_methods",
    "property_type": "property_types",
    "project": "projects",
    "protocol": "protocols"
}

# Files to update
FILES_TO_UPDATE = [
    "populate_molecules.py",
    "populate_mixtures.py",
    "populate_predictions.py",
    "populate_experiments.py"
]

def update_file(file_path):
    """Update a file to use plural table names."""
    try:
        # Read the file
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Make a backup of the original file
        backup_path = f"{file_path}.bak"
        with open(backup_path, 'w') as f:
            f.write(content)
        logger.info(f"Created backup of {file_path} at {backup_path}")
        
        # Replace table names
        modified_content = content
        changes_made = False
        
        for singular, plural in TABLE_MAPPINGS.items():
            # Only replace when it's a table reference (preceded by 'table(')
            pattern = r'table\(["\']' + singular + r'["\']\)'
            replacement = r'table("' + plural + r'")'
            
            # Count occurrences before replacement
            count_before = len(re.findall(pattern, modified_content))
            
            # Perform the replacement
            modified_content = re.sub(pattern, replacement, modified_content)
            
            # Count occurrences after replacement
            count_after = len(re.findall(pattern, modified_content))
            
            # Check if any replacements were made
            if count_before != count_after:
                changes_made = True
                logger.info(f"Replaced {count_before - count_after} occurrences of '{singular}' with '{plural}' in {file_path}")
        
        # Write the modified content back to the file if changes were made
        if changes_made:
            with open(file_path, 'w') as f:
                f.write(modified_content)
            logger.info(f"Updated {file_path} with plural table names")
            return True
        else:
            logger.info(f"No changes needed for {file_path}")
            return False
    
    except Exception as e:
        logger.error(f"Error updating {file_path}: {str(e)}")
        return False

def main():
    """Update all population scripts to use plural table names."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Update Population Scripts")
    print("=" * 80)
    
    # Check if files exist
    missing_files = [f for f in FILES_TO_UPDATE if not os.path.exists(f)]
    if missing_files:
        logger.error(f"The following files are missing: {', '.join(missing_files)}")
        print(f"Error: The following files are missing: {', '.join(missing_files)}")
        return 1
    
    # Update each file
    updated_files = []
    for file_path in FILES_TO_UPDATE:
        print(f"\nUpdating {file_path}...")
        if update_file(file_path):
            updated_files.append(file_path)
    
    # Summary
    if updated_files:
        print("\n" + "=" * 60)
        print("Update Complete")
        print("=" * 60)
        print(f"\nSuccessfully updated {len(updated_files)} files:")
        for file in updated_files:
            print(f"- {file}")
        
        print("\nYou can now run the population scripts to populate the database with the correct table names.")
        return 0
    else:
        print("\nNo files were updated. All files may already be using the correct table names.")
        return 0

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""
CryoProtect v2 - Comprehensive Table Name Update

This script thoroughly updates all population scripts to use the correct plural table names
(e.g., "molecules", "molecular_properties") instead of singular names.
It checks for all possible ways a table name might be referenced.

Usage:
    python update_table_names_comprehensive.py [--dry-run]

Author: Claude
Date: April 18, 2025
"""

import os
import re
import sys
import argparse
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("update_table_names_comprehensive.log"),
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
    "protocol": "protocols",
    "user_profile": "user_profiles"
}

# Files to update
FILES_TO_UPDATE = [
    "populate_molecules.py",
    "populate_mixtures.py",
    "populate_predictions.py",
    "populate_experiments.py",
    "populate_molecules_production.py",
    "populate_mixtures_production.py",
    "populate_predictions_production.py",
    "populate_experiments_production.py",
    "populate_database.py",
    "populate_database_main.py",
    "populate_calculation_methods_production.py"
]

def update_file(file_path, dry_run=False):
    """Update a file to use plural table names with comprehensive pattern matching."""
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            logger.warning(f"File not found: {file_path}")
            return False
            
        # Read the file
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Make a backup of the original file
        backup_path = f"{file_path}.bak-table-names"
        if not dry_run:
            with open(backup_path, 'w', encoding='utf-8') as f:
                f.write(content)
            logger.info(f"Created backup of {file_path} at {backup_path}")
        else:
            logger.info(f"DRY RUN: Would create backup of {file_path} at {backup_path}")
        
        # Replace table names
        modified_content = content
        changes_made = False
        replacements = []
        
        for singular, plural in TABLE_MAPPINGS.items():
            # Patterns to match different ways tables might be referenced
            patterns = [
                # table("molecule") or table('molecule')
                (rf'table\(["\']({singular})["\']', f'table("{plural}"'),
                
                # supabase.from("molecule")
                (rf'from\(["\']({singular})["\']', f'from("{plural}"'),
                
                # "table_name": "molecule"
                (rf'"table_name":\s*["\']({singular})["\']', f'"table_name": "{plural}"'),
                
                # populate_molecule_table()
                (rf'populate_{singular}_table', f'populate_{plural}_table'),
                
                # response.molecule or response.data.molecule
                (rf'response\.{singular}', f'response.{plural}'),
                (rf'response\.data\.{singular}', f'response.data.{plural}'),
                
                # class MoleculeResource
                (rf'class {singular.title()}Resource', f'class {plural.title()}Resource'),
                
                # self.molecule = ...
                (rf'self\.{singular}\s*=', f'self.{plural} ='),
                
                # API paths like '/api/v1/molecule'
                (rf'[\'"/]api/v\d+/{singular}[\'"/]', f'"/api/v1/{plural}"'),
                
                # f"Creating {len(molecules)} molecule records"
                (rf'(["\'])(\{{\w+\}}|\d+) {singular}', f'\\1\\2 {plural}'),
                
                # Comment references: "Check if molecule already exists"
                (rf'# .* {singular} ', f'# \\g<0>'.replace(singular, plural)),
                
                # For explicit lookups: table="molecule"
                (rf'table\s*=\s*["\']({singular})["\']', f'table="{plural}"'),
                
                # Function names: def get_molecule_by_id
                (rf'def get_{singular}_by', f'def get_{plural}_by'),
                
                # Variable names that might include table references
                (rf'\b{singular}_id\b', f'{plural}_id'),
                (rf'\b{singular}_data\b', f'{plural}_data'),
                (rf'\b{singular}_response\b', f'{plural}_response')
            ]
            
            for pattern, replacement in patterns:
                # Count occurrences before replacement
                matches = re.findall(pattern, modified_content)
                
                if matches:
                    # Perform the replacement
                    modified_content = re.sub(pattern, replacement, modified_content)
                    changes_made = True
                    replacements.append((singular, plural, len(matches), pattern))
        
        # Write the modified content back to the file if changes were made
        if changes_made and not dry_run:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(modified_content)
            logger.info(f"Updated {file_path} with plural table names")
            
            # Log detailed replacements
            for singular, plural, count, pattern in replacements:
                logger.info(f"  - Replaced {count} occurrences of '{singular}' with '{plural}' (pattern: {pattern})")
                
            return True
        elif changes_made and dry_run:
            logger.info(f"DRY RUN: Would update {file_path} with plural table names")
            
            # Log detailed replacements
            for singular, plural, count, pattern in replacements:
                logger.info(f"  - Would replace {count} occurrences of '{singular}' with '{plural}' (pattern: {pattern})")
                
            return True
        else:
            logger.info(f"No changes needed for {file_path}")
            return False
    
    except Exception as e:
        logger.error(f"Error updating {file_path}: {str(e)}")
        return False

def find_additional_population_scripts():
    """Find any additional population scripts that might not be in the predefined list."""
    additional_files = []
    
    for f in os.listdir('.'):
        if f.endswith('.py') and f.startswith('populate_') and f not in FILES_TO_UPDATE:
            additional_files.append(f)
    
    return additional_files

def main():
    """Update all population scripts to use plural table names."""
    parser = argparse.ArgumentParser(description="Update population scripts to use plural table names.")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be changed without making actual changes")
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Comprehensive Table Name Update")
    print("=" * 80)
    
    # Find additional population scripts
    additional_files = find_additional_population_scripts()
    if additional_files:
        FILES_TO_UPDATE.extend(additional_files)
        print(f"\nFound {len(additional_files)} additional population scripts to check:")
        for f in additional_files:
            print(f"  - {f}")
    
    # Check if files exist
    existing_files = [f for f in FILES_TO_UPDATE if os.path.exists(f)]
    missing_files = [f for f in FILES_TO_UPDATE if not os.path.exists(f)]
    
    if missing_files:
        logger.warning(f"The following files were not found: {', '.join(missing_files)}")
        print(f"Warning: The following files were not found: {', '.join(missing_files)}")
        print(f"Continuing with the {len(existing_files)} files that were found.")
    
    # Update each file
    updated_files = []
    for file_path in existing_files:
        print(f"\nUpdating {file_path}...")
        if update_file(file_path, args.dry_run):
            updated_files.append(file_path)
    
    # Summary
    print("\n" + "=" * 60)
    if args.dry_run:
        print("Dry Run Complete")
    else:
        print("Update Complete")
    print("=" * 60)
    
    if updated_files:
        if args.dry_run:
            print(f"\nWould update {len(updated_files)} files:")
        else:
            print(f"\nSuccessfully updated {len(updated_files)} files:")
        for file in updated_files:
            print(f"- {file}")
        
        if not args.dry_run:
            print("\nYou can now run the population scripts to populate the database with the correct table names.")
        else:
            print("\nRun the script without the --dry-run flag to apply these changes.")
    else:
        print("\nNo files were updated. All files may already be using the correct table names.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

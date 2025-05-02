#!/usr/bin/env python3
"""
Runner script for ChEMBL data import.
This script adds the parent directory to the Python path and then runs the import_chembl_data function.
"""

import os
import sys

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath('../..'))

# Now we can import from the chembl module
from chembl.client import ResilientChEMBLClient
from chembl.checkpoint import CheckpointManager
from chembl.error_handler import ErrorCategory, classify_error, get_recovery_strategy
from chembl.worker import ChEMBLWorker

# Import the import_chembl_data function from chembl_import.py
from chembl_import import import_chembl_data, verify_imported_data

def main():
    """Main function to run the ChEMBL data import."""
    print("Starting ChEMBL data import...")
    
    # Get the Supabase project ID
    project_id = "tsdlmynydfuypiugmkev"
    
    # Import ChEMBL data with a target of 5000 molecules
    results = import_chembl_data(project_id, target_count=5000)
    
    print(f"Imported {results['molecules_imported']} molecules with {results['properties_imported']} properties")
    
    if results["errors"]:
        print(f"Encountered {len(results['errors'])} errors")
        for error in results["errors"][:10]:  # Log only the first 10 errors
            print(f"Error: {error}")
    
    # Verify imported data
    verification_results = verify_imported_data(project_id)
    
    print(f"Total molecules: {verification_results['total_molecules']}")
    print(f"Molecules with properties: {verification_results['molecules_with_properties']}")
    print(f"Property completeness: {verification_results['property_completeness_percentage']:.1f}%")
    
    # Check if success criteria are met
    success_criteria = verification_results["success_criteria"]
    for criterion, details in success_criteria.items():
        status = "✅ Met" if details["met"] else "❌ Not met"
        print(f"{criterion}: {details['actual']} / {details['target']} - {status}")
    
    print("ChEMBL data import completed successfully")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
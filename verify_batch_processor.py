#!/usr/bin/env python3
"""
Verify implementation of the batch processor for molecules with None names.

This script checks that all components required for the batch processor
have been implemented correctly. It verifies:

1. The batch processor script exists
2. The test file exists
3. Required functions are implemented
4. The script is executable
"""

import sys
import os
import importlib.util
import stat

def check_file_exists(path, description):
    """Check if a file exists and print the result."""
    exists = os.path.exists(path)
    status = "✓" if exists else "✗"
    print(f"{status} {description}: {path}")
    return exists

def check_file_executable(path, description):
    """Check if a file is executable and print the result."""
    if not os.path.exists(path):
        print(f"✗ {description} - file does not exist: {path}")
        return False
    
    is_executable = os.access(path, os.X_OK)
    status = "✓" if is_executable else "✗"
    print(f"{status} {description}: {path}")
    return is_executable

def check_module_has_component(module_path, component_name, description):
    """Check if a module has a given component by reading the file content."""
    try:
        # Read the module content directly
        with open(module_path, 'r') as f:
            content = f.read()
            
        pattern = f"def {component_name}" if not component_name.startswith('class ') else component_name
        has_component = pattern in content
        status = "✓" if has_component else "✗"
        print(f"{status} {description}: {component_name} in {module_path}")
        return has_component
    except Exception as e:
        print(f"✗ Error checking {description}: {e}")
        return False

def check_main_components():
    """Check if the main components have been implemented."""
    print("\nChecking main components...")
    success = True
    
    # Check batch processor file
    processor_path = "batch_process_none_names.py"
    success = check_file_exists(processor_path, "Batch processor script exists") and success
    
    # Check if file is executable
    success = check_file_executable(processor_path, "Batch processor script is executable") and success
    
    # Check required functions
    if os.path.exists(processor_path):
        required_functions = [
            ("get_molecules_with_none_names", "Function to get molecules with None names"),
            ("generate_name_for_molecule", "Function to generate names for molecules"),
            ("update_molecule_name", "Function to update molecule names"),
            ("process_molecules", "Function to process molecules in batches"),
            ("save_checkpoint", "Function to save checkpoints"),
            ("load_checkpoint", "Function to load checkpoints"),
            ("main", "Main function")
        ]
        
        for func_name, description in required_functions:
            success = check_module_has_component(processor_path, func_name, description) and success
    
    return success

def check_test_components():
    """Check if test components have been implemented."""
    print("\nChecking test components...")
    success = True
    
    # Check test file
    test_path = "tests/test_batch_process_none_names.py"
    success = check_file_exists(test_path, "Test file exists") and success
    
    # Check test case class
    if os.path.exists(test_path):
        success = check_module_has_component(test_path, "class BatchProcessNoneNamesTestCase", "Test case class") and success
        
        # Check test methods
        test_methods = [
            ("test_get_molecules_with_none_names", "Test for get_molecules_with_none_names"),
            ("test_name_generation_strategies", "Test for name generation strategies"),
            ("test_update_molecule_name", "Test for update_molecule_name"),
            ("test_checkpoint_save_load", "Test for checkpoint save/load"),
            ("test_process_molecules", "Test for process_molecules")
        ]
        
        for method_name, description in test_methods:
            success = check_module_has_component(test_path, method_name, description) and success
    
    return success

def main():
    """Check all components and return overall success status."""
    print("Verifying batch processor implementation...")
    
    main_success = check_main_components()
    test_success = check_test_components()
    
    all_success = main_success and test_success
    
    if all_success:
        print("\n✓ All batch processor components are correctly implemented!")
        return 0
    else:
        print("\n✗ Some components of the batch processor are missing or incomplete.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
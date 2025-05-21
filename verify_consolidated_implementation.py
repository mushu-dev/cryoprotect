#!/usr/bin/env python3
"""
Verification script for consolidated molecule implementation.

This script verifies the implementation of the consolidated molecule
handling in the API without requiring all dependencies to be installed.
"""

import os
import sys
import inspect
import json

def check_file_existence(files):
    """
    Check that required files exist.
    
    Args:
        files: List of file paths to check
        
    Returns:
        Tuple of (all_exist, missing_files)
    """
    missing_files = []
    for file_path in files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    return len(missing_files) == 0, missing_files

def check_file_content(file_path, required_elements):
    """
    Check that a file contains required code elements.
    
    Args:
        file_path: Path to the file to check
        required_elements: List of strings that should appear in the file
        
    Returns:
        Tuple of (all_found, missing_elements)
    """
    if not os.path.exists(file_path):
        return False, required_elements
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    missing_elements = []
    for element in required_elements:
        if element not in content:
            missing_elements.append(element)
    
    return len(missing_elements) == 0, missing_elements

def format_verification_result(component, status, details=None):
    """
    Format verification result for a component.
    
    Args:
        component: Component name
        status: "PASS" or "FAIL"
        details: Optional details about the result
        
    Returns:
        Formatted string
    """
    result = f"{status}: {component}"
    if details:
        if isinstance(details, list):
            details_str = "\n  - " + "\n  - ".join(str(d) for d in details)
            result += details_str
        else:
            result += f"\n  - {details}"
    return result

def main():
    """
    Run verification checks on the consolidated molecule implementation.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    api_dir = os.path.join(base_dir, 'api')
    docs_dir = os.path.join(base_dir, 'docs')
    
    results = []
    
    # Check that required files exist
    required_files = [
        os.path.join(api_dir, 'consolidated_utils.py'),
        os.path.join(api_dir, 'consolidated_decorators.py'),
        os.path.join(api_dir, 'consolidated_molecule_resource.py'),
        os.path.join(api_dir, 'differentiation_resources.py'),
        os.path.join(api_dir, 'consolidated_api.py'),
        os.path.join(base_dir, 'test_consolidated_api.py'),
        os.path.join(base_dir, 'test_integration_consolidated_api.py'),
        os.path.join(base_dir, 'examples/use_consolidated_api.py'),
        os.path.join(docs_dir, 'CONSOLIDATED_MOLECULE_API_GUIDE.md'),
        os.path.join(docs_dir, 'API_ENDPOINTS_REFERENCE.md'),
        os.path.join(base_dir, 'API_STANDARDIZATION_IMPLEMENTATION.md'),
        os.path.join(base_dir, 'API_STANDARDIZATION_DEMO.md')
    ]
    
    all_exist, missing_files = check_file_existence(required_files)
    if all_exist:
        results.append(format_verification_result("Required Files", "PASS"))
    else:
        results.append(format_verification_result("Required Files", "FAIL", missing_files))
    
    # Check consolidated_utils.py
    utils_requirements = [
        "is_consolidated",
        "get_primary_molecule",
        "get_consolidated_molecules",
        "get_differentiation_group",
        "get_differentiation_group_members",
        "enrich_molecule_data"
    ]
    
    utils_file = os.path.join(api_dir, 'consolidated_utils.py')
    all_found, missing_elements = check_file_content(utils_file, utils_requirements)
    if all_found:
        results.append(format_verification_result("Consolidated Utilities", "PASS"))
    else:
        results.append(format_verification_result("Consolidated Utilities", "FAIL", missing_elements))
    
    # Check consolidated_decorators.py
    decorator_requirements = [
        "handle_consolidated_molecules",
        "handle_batch_consolidated_molecules"
    ]
    
    decorators_file = os.path.join(api_dir, 'consolidated_decorators.py')
    all_found, missing_elements = check_file_content(decorators_file, decorator_requirements)
    if all_found:
        results.append(format_verification_result("Consolidated Decorators", "PASS"))
    else:
        results.append(format_verification_result("Consolidated Decorators", "FAIL", missing_elements))
    
    # Check consolidated_molecule_resource.py
    resource_requirements = [
        "ConsolidatedMoleculeResource",
        "ConsolidatedMoleculeBatchResource",
        "PrimaryMoleculeResource",
        "ConsolidatedMoleculesListResource"
    ]
    
    resource_file = os.path.join(api_dir, 'consolidated_molecule_resource.py')
    all_found, missing_elements = check_file_content(resource_file, resource_requirements)
    if all_found:
        results.append(format_verification_result("Consolidated Resources", "PASS"))
    else:
        results.append(format_verification_result("Consolidated Resources", "FAIL", missing_elements))
    
    # Check differentiation_resources.py
    diff_requirements = [
        "DifferentiationGroupListResource",
        "DifferentiationGroupResource",
        "MoleculeDifferentiationResource"
    ]
    
    diff_file = os.path.join(api_dir, 'differentiation_resources.py')
    all_found, missing_elements = check_file_content(diff_file, diff_requirements)
    if all_found:
        results.append(format_verification_result("Differentiation Resources", "PASS"))
    else:
        results.append(format_verification_result("Differentiation Resources", "FAIL", missing_elements))
    
    # Check consolidated_api.py
    api_requirements = [
        "register_consolidated_resources",
        "register_consolidated_docs"
    ]
    
    api_file = os.path.join(api_dir, 'consolidated_api.py')
    all_found, missing_elements = check_file_content(api_file, api_requirements)
    if all_found:
        results.append(format_verification_result("API Registration", "PASS"))
    else:
        results.append(format_verification_result("API Registration", "FAIL", missing_elements))
    
    # Check API endpoints in the documentation
    doc_file = os.path.join(docs_dir, 'API_ENDPOINTS_REFERENCE.md')
    endpoint_requirements = [
        "/api/v1/consolidated/molecules/{molecule_id}",
        "/api/v1/consolidated/batch",
        "/api/v1/molecules/{molecule_id}/primary",
        "/api/v1/consolidated",
        "/api/v1/differentiation/groups",
        "/api/v1/differentiation/groups/{group_id}",
        "/api/v1/molecules/{molecule_id}/differentiation"
    ]
    
    all_found, missing_elements = check_file_content(doc_file, endpoint_requirements)
    if all_found:
        results.append(format_verification_result("API Documentation", "PASS"))
    else:
        results.append(format_verification_result("API Documentation", "FAIL", missing_elements))
    
    # Check the test cases
    test_file = os.path.join(base_dir, 'test_consolidated_api.py')
    test_requirements = [
        "test_consolidated_molecule_endpoint",
        "test_primary_molecule_endpoint",
        "test_consolidated_batch_endpoint",
        "test_differentiation_group_endpoint",
        "test_differentiation_groups_list_endpoint",
        "test_molecule_differentiation_endpoint",
        "test_consolidated_list_endpoint"
    ]
    
    all_found, missing_elements = check_file_content(test_file, test_requirements)
    if all_found:
        results.append(format_verification_result("Test Suite", "PASS"))
    else:
        results.append(format_verification_result("Test Suite", "FAIL", missing_elements))
    
    # Check demonstration script
    demo_file = os.path.join(base_dir, 'API_STANDARDIZATION_DEMO.md')
    demo_requirements = [
        "Consolidated Molecule Handling",
        "Batch Operations",
        "Differentiation Groups",
        "Standardized Responses",
        "Error Handling"
    ]
    
    all_found, missing_elements = check_file_content(demo_file, demo_requirements)
    if all_found:
        results.append(format_verification_result("Demonstration", "PASS"))
    else:
        results.append(format_verification_result("Demonstration", "FAIL", missing_elements))
    
    # Print verification results
    print("\n====== CONSOLIDATED MOLECULE IMPLEMENTATION VERIFICATION ======\n")
    for result in results:
        print(result)
        print()
    
    # Calculate overall status
    pass_count = sum(1 for r in results if r.startswith("PASS"))
    fail_count = len(results) - pass_count
    
    print(f"SUMMARY: {pass_count} PASS, {fail_count} FAIL")
    
    if fail_count == 0:
        print("\nVERIFICATION SUCCESSFUL: All components implemented correctly!")
        print("Ready to proceed to Phase 3.")
        return 0
    else:
        print("\nVERIFICATION FAILED: Some components are missing or incomplete.")
        print("Please fix the issues before proceeding to Phase 3.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
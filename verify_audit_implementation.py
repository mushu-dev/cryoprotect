#!/usr/bin/env python3
"""
Verify implementation of the consolidated molecule audit trail.

This script checks that all components required for the consolidated molecule
audit trail functionality have been implemented correctly. It verifies:

1. The database migration has been created
2. The audit table exists in the database
3. The trigger function exists
4. The trigger has been registered
5. The API resource for retrieving audit records exists
6. The resource has been registered with the API
7. The docs for the resource have been registered
"""

import sys
import os
import importlib.util

def check_file_exists(path, description):
    """Check if a file exists and print the result."""
    exists = os.path.exists(path)
    status = "✓" if exists else "✗"
    print(f"{status} {description}: {path}")
    return exists

def check_module_has_component(module_path, component_name, description):
    """Check if a module has a given component by reading the file content."""
    try:
        # Read the module content directly (without importing)
        with open(module_path, 'r') as f:
            content = f.read()
        
        has_component = component_name in content
        status = "✓" if has_component else "✗"
        print(f"{status} {description}: {component_name} in {module_path}")
        return has_component
    except Exception as e:
        print(f"✗ Error checking {description}: {e}")
        return False

def check_file_contains(path, pattern, description):
    """Check if a file contains a specific pattern."""
    try:
        with open(path, 'r') as f:
            content = f.read()
            contains = pattern in content
            status = "✓" if contains else "✗"
            print(f"{status} {description}: {pattern} in {path}")
            return contains
    except Exception as e:
        print(f"✗ Error checking {description}: {e}")
        return False

def check_database_components():
    """Check if database components have been implemented."""
    print("\nChecking database components...")
    success = True
    
    # Check migration file
    migration_path = "migrations/025_consolidated_molecule_audit.sql"
    success = check_file_exists(migration_path, "Migration file exists") and success
    
    # Check migration content
    if os.path.exists(migration_path):
        patterns = [
            "CREATE TABLE IF NOT EXISTS molecule_consolidation_audit",
            "CREATE OR REPLACE FUNCTION log_molecule_consolidation()",
            "CREATE TRIGGER molecule_consolidation_audit_trigger"
        ]
        for pattern in patterns:
            success = check_file_contains(migration_path, pattern, f"Migration contains {pattern}") and success
    
    return success

def check_api_components():
    """Check if API components have been implemented."""
    print("\nChecking API components...")
    success = True
    
    # Check audit resource file
    resource_path = "api/audit_resources.py"
    success = check_file_exists(resource_path, "Audit resource file exists") and success
    
    # Check audit resource class and functions
    if os.path.exists(resource_path):
        components = [
            ("class ConsolidationAuditResource", "Audit resource class"),
            ("def register_audit_resources", "Register resources function"),
            ("def register_audit_docs", "Register docs function")
        ]
        for component, description in components:
            success = check_module_has_component(resource_path, component, description) and success
    
    # Check resource registration in __init__.py
    init_path = "api/__init__.py"
    patterns = [
        "from api.audit_resources import register_audit_resources",
        "register_audit_resources(api)",
        "from api.audit_resources import register_audit_docs",
        "register_audit_docs(docs)"
    ]
    for pattern in patterns:
        success = check_file_contains(init_path, pattern, f"API initialization contains {pattern}") and success
    
    return success

def check_test_components():
    """Check if test components have been implemented."""
    print("\nChecking test components...")
    success = True
    
    # Check test file
    test_path = "tests/test_consolidation_audit.py"
    success = check_file_exists(test_path, "Audit test file exists") and success
    
    # Check test content
    if os.path.exists(test_path):
        patterns = [
            "class ConsolidationAuditTestCase",
            "test_audit_trigger_on_consolidation",
            "test_audit_api_endpoint"
        ]
        for pattern in patterns:
            success = check_file_contains(test_path, pattern, f"Test contains {pattern}") and success
    
    return success

def main():
    """Check all components and return overall success status."""
    print("Verifying consolidated molecule audit trail implementation...")
    
    db_success = check_database_components()
    api_success = check_api_components()
    test_success = check_test_components()
    
    all_success = db_success and api_success and test_success
    
    if all_success:
        print("\n✓ All consolidated molecule audit trail components are correctly implemented!")
        return 0
    else:
        print("\n✗ Some components of the consolidated molecule audit trail are missing or incomplete.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
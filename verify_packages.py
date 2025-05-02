#!/usr/bin/env python
"""
Package Verification Script

This script verifies that all required packages are properly installed and accessible.
It attempts to import each package and reports success or failure.
"""

import sys
import importlib
import platform
import os
import datetime

# List of packages to verify
PACKAGES = [
    "scipy",
    "rdkit",
    "psutil",
    "xlsxwriter",
    "seaborn",
    "sklearn",               # scikit-learn is imported as sklearn
    "pythonjsonlogger.jsonlogger",  # Correct import path
    "ecs_logger",            # Use underscore instead of hyphen
    "prometheus_client",
    "yaml"                   # PyYAML is imported as yaml
]

def check_package(package_name):
    """Attempt to import a package and return success/failure status."""
    try:
        importlib.import_module(package_name)
        return True, None
    except ImportError as e:
        return False, str(e)

def get_python_info():
    """Get information about the Python environment."""
    info = {
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "python_path": sys.path,
        "current_directory": os.getcwd(),
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    return info

def main():
    """Main function to verify packages and print results."""
    print("\n===== PACKAGE VERIFICATION REPORT =====")
    print(f"Date/Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Python Version: {sys.version}")
    print(f"Python Executable: {sys.executable}")
    print(f"Platform: {platform.platform()}")
    print("\n")
    
    all_packages_available = True
    results = []
    
    for package in PACKAGES:
        success, error = check_package(package)
        results.append((package, success, error))
        if not success:
            all_packages_available = False
    
    # Print results in a table format
    print(f"{'Package':<20} {'Status':<10} {'Error Message'}")
    print("-" * 60)
    for package, success, error in results:
        status = "SUCCESS" if success else "FAILED"
        error_msg = error if error else ""
        print(f"{package:<20} {status:<10} {error_msg}")
    
    print("\n")
    if all_packages_available:
        print("[SUCCESS] All packages are successfully installed and accessible.")
    else:
        print("[FAILED] Some packages are not available. Please check the report above.")
    
    # Print Python path for debugging
    print("\nPython Path:")
    for path in sys.path:
        print(f"  - {path}")

if __name__ == "__main__":
    main()
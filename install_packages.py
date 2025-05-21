#!/usr/bin/env python
"""
Package Installation Script

This script installs all required packages for the CryoProtect application.
It ensures that packages are installed in the correct environment and persist between restarts.
"""

import sys
import subprocess
import os
import platform
import datetime

# List of packages to install
PACKAGES = [
    "scipy",
    "rdkit",
    "psutil",
    "xlsxwriter",
    "seaborn",
    "scikit-learn",         # Use latest available version
    "python_json_logger",   # Use underscore instead of hyphen
    "ecs_logger",           # Use underscore instead of hyphen
    "prometheus_client",
    "PyYAML"                # Use correct capitalization
]

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

def install_package(package_name):
    """Install a package using pip."""
    print(f"Installing {package_name}...")
    try:
        # Use the Python executable from the current environment to ensure packages are installed in the right place
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", package_name],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"Successfully installed {package_name}")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Failed to install {package_name}: {e.stderr}")
        return False, e.stderr

def main():
    """Main function to install packages."""
    print("\n===== PACKAGE INSTALLATION =====")
    print(f"Date/Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Print Python environment information
    info = get_python_info()
    print(f"Python Version: {info['python_version']}")
    print(f"Python Executable: {info['python_executable']}")
    print(f"Platform: {info['platform']}")
    print(f"Current Directory: {info['current_directory']}")
    print("\n")
    
    # Install each package
    results = []
    all_packages_installed = True
    
    for package in PACKAGES:
        success, output = install_package(package)
        results.append((package, success, output))
        if not success:
            all_packages_installed = False
    
    # Print results in a table format
    print("\n===== INSTALLATION RESULTS =====")
    print(f"{'Package':<20} {'Status':<10}")
    print("-" * 40)
    for package, success, _ in results:
        status = "SUCCESS" if success else "FAILED"
        print(f"{package:<20} {status:<10}")
    
    print("\n")
    if all_packages_installed:
        print("[SUCCESS] All packages were successfully installed.")
    else:
        print("[WARNING] Some packages failed to install. Check the output above for details.")
    
    # Print Python path for debugging
    print("\nPython Path:")
    for path in sys.path:
        print(f"  - {path}")

if __name__ == "__main__":
    main()
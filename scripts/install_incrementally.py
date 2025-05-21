#!/usr/bin/env python3
"""
Incremental Dependency Installer for CryoProtect

This script incrementally adds dependencies and tests if the app can import them.
"""

import sys
import importlib
import subprocess
import time

# ANSI colors
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[0;33m'
BLUE = '\033[0;34m'
RESET = '\033[0m'

def print_colored(color, message):
    """Print a colored message"""
    print(f"{color}{message}{RESET}")

def install_package(package_name):
    """Install a package using pip"""
    print_colored(YELLOW, f"Installing {package_name}...")
    result = subprocess.run([sys.executable, "-m", "pip", "install", package_name], 
                            capture_output=True, text=True)
    
    if result.returncode == 0:
        print_colored(GREEN, f"✓ Successfully installed {package_name}")
        return True
    else:
        print_colored(RED, f"✗ Failed to install {package_name}")
        print(result.stderr)
        return False

def test_import(module_name):
    """Test if a module can be imported"""
    try:
        importlib.import_module(module_name)
        print_colored(GREEN, f"✓ Successfully imported {module_name}")
        return True
    except ImportError as e:
        print_colored(RED, f"✗ Failed to import {module_name}: {str(e)}")
        return False

def main():
    print_colored(BLUE, "CryoProtect Incremental Dependency Installer")
    print("=" * 60)
    
    # Core Flask dependencies
    core_deps = [
        "flask==3.0.2",
        "flask-restful==0.3.10",
        "flask-cors==4.0.0",
        "python-dotenv==1.0.1",
        "requests==2.31.0"
    ]
    
    # Flask extensions
    flask_extensions = [
        "flask-apispec==0.11.4",
        "apispec==6.3.0",
        "marshmallow==3.21.0",
        "flask-mail==0.9.1",
    ]
    
    # Documentation packages
    doc_packages = [
        "sphinx==7.1.2",
        "sphinx-rtd-theme==1.3.0",
    ]
    
    # Database packages
    db_packages = [
        "psycopg2-binary==2.9.9",
        "sqlalchemy==2.0.27",
        "alembic==1.13.1",
        "supabase==2.1.0",
    ]
    
    # Monitoring/observability
    monitoring_packages = [
        "prometheus-client==0.19.0",
        "prometheus-flask-exporter==0.22.4",
    ]
    
    # Security packages
    security_packages = [
        "pyjwt==2.8.0",
        "cryptography==42.0.5",
    ]
    
    # Testing packages
    testing_packages = [
        "pytest==7.4.0",
        "pytest-cov==4.1.0",
    ]
    
    # All package groups
    all_groups = [
        ("Core Flask", core_deps),
        ("Flask Extensions", flask_extensions),
        ("Documentation", doc_packages),
        ("Database", db_packages),
        ("Monitoring", monitoring_packages),
        ("Security", security_packages),
        ("Testing", testing_packages),
    ]
    
    # Ask which groups to install
    print("Which package groups would you like to install?")
    for i, (name, _) in enumerate(all_groups):
        print(f"{i+1}. {name}")
    print("0. All groups")
    
    choice = input("Enter your choice (comma-separated for multiple): ")
    choices = [int(c.strip()) for c in choice.split(",") if c.strip().isdigit()]
    
    groups_to_install = []
    if 0 in choices:
        groups_to_install = all_groups
    else:
        groups_to_install = [all_groups[i-1] for i in choices if 1 <= i <= len(all_groups)]
    
    # Install selected package groups
    for group_name, packages in groups_to_install:
        print_colored(BLUE, f"\nInstalling {group_name} packages...")
        for package in packages:
            install_package(package)
            # Test import of the base package
            base_package = package.split('==')[0].split('[')[0]
            # Convert package name to module name (replace - with _)
            module_name = base_package.replace('-', '_')
            test_import(module_name)
            time.sleep(0.5)  # Small delay to avoid overwhelming output
    
    print_colored(GREEN, "\nDependency installation complete!")
    print("You can now try running the application again.")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
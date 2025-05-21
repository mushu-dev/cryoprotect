#!/usr/bin/env python3
"""
Dependency Analyzer for CryoProtect

This script analyzes app.py and attempts to identify all required dependencies.
It will create a minimal functional version of app.py for testing.
"""

import os
import re
import sys
import importlib
import subprocess
from pathlib import Path

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

def extract_imports(file_path):
    """Extract import statements from a Python file"""
    imports = []
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Find standard imports
    import_pattern = re.compile(r'^import\s+([\w.]+)', re.MULTILINE)
    for match in import_pattern.finditer(content):
        imports.append(match.group(1).split('.')[0])
    
    # Find from imports
    from_import_pattern = re.compile(r'^from\s+([\w.]+)\s+import', re.MULTILINE)
    for match in from_import_pattern.finditer(content):
        base_module = match.group(1).split('.')[0]
        if base_module != '':  # Skip relative imports
            imports.append(base_module)
    
    return list(set(imports))

def check_installed(module_name):
    """Check if a module is installed"""
    try:
        importlib.import_module(module_name)
        return True
    except ImportError:
        return False

def create_minimal_app(app_path, output_path):
    """Create a minimal version of app.py for testing"""
    print_colored(BLUE, f"Creating minimal version of {app_path}...")
    
    with open(app_path, 'r') as f:
        content = f.read()
    
    # Create a highly simplified version
    minimal_app = """#!/usr/bin/env python3
\"\"\"
Minimal CryoProtect App for Dependency Testing
\"\"\"
from flask import Flask, jsonify

app = Flask(__name__)

@app.route('/')
def index():
    return jsonify({
        'status': 'ok',
        'message': 'CryoProtect Minimal App is running'
    })

@app.route('/health')
def health():
    return jsonify({
        'status': 'healthy',
        'time': 'now'
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
"""
    
    with open(output_path, 'w') as f:
        f.write(minimal_app)
    
    print_colored(GREEN, f"Created minimal app at {output_path}")

def main():
    print_colored(BLUE, "CryoProtect Dependency Analyzer")
    print("=" * 60)
    
    repo_root = Path(__file__).parent.parent
    app_path = repo_root / "app.py"
    
    if not app_path.exists():
        print_colored(RED, f"Error: app.py not found at {app_path}")
        return 1
    
    # Extract imports from app.py
    print_colored(BLUE, "Analyzing app.py imports...")
    all_imports = extract_imports(app_path)
    
    # Filter out standard library imports
    std_lib_modules = sys.stdlib_module_names
    external_imports = [imp for imp in all_imports if imp not in std_lib_modules]
    
    print(f"Found {len(external_imports)} external modules:")
    for module in sorted(external_imports):
        if check_installed(module):
            print_colored(GREEN, f"✓ {module} (installed)")
        else:
            print_colored(RED, f"✗ {module} (not installed)")
    
    # Create requirements.txt
    requirements_path = repo_root / "requirements_minimal.txt"
    print_colored(BLUE, f"\nCreating {requirements_path}...")
    
    # Known pip package names that differ from import names
    package_mapping = {
        'flask': 'Flask',
        'sklearn': 'scikit-learn',
        'PIL': 'Pillow',
    }
    
    with open(requirements_path, 'w') as f:
        for module in sorted(external_imports):
            package_name = package_mapping.get(module, module)
            f.write(f"{package_name}\n")
    
    print_colored(GREEN, f"Created {requirements_path}")
    
    # Create a minimal app.py
    minimal_app_path = repo_root / "minimal_app.py"
    create_minimal_app(app_path, minimal_app_path)
    
    # Offer to install missing packages
    missing_modules = [m for m in external_imports if not check_installed(m)]
    if missing_modules:
        print_colored(YELLOW, f"\nFound {len(missing_modules)} missing modules.")
        install_all = input("Install all missing modules? (y/n): ").lower() == 'y'
        
        if install_all:
            for module in missing_modules:
                package_name = package_mapping.get(module, module)
                install_package(package_name)
        else:
            print_colored(YELLOW, "Skipping package installation.")
    
    print_colored(BLUE, "\nRecommended run command:")
    print(f"cd {repo_root} && python minimal_app.py")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
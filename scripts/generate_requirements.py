#!/usr/bin/env python3
"""
Generate comprehensive requirements.txt for CryoProtect

This script scans all Python files in the project and generates a requirements.txt
file with all external dependencies.
"""

import os
import re
import sys
import importlib
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

def find_python_files(root_dir):
    """Find all Python files in a directory recursively"""
    python_files = []
    for root, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    return python_files

def extract_imports(file_path):
    """Extract import statements from a Python file"""
    imports = []
    with open(file_path, 'r', errors='ignore') as f:
        try:
            content = f.read()
        except UnicodeDecodeError:
            print_colored(YELLOW, f"Warning: Could not read {file_path} due to encoding issues")
            return []
    
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

def is_standard_library(module_name):
    """Check if a module is part of the standard library"""
    if hasattr(sys, 'stdlib_module_names'):
        return module_name in sys.stdlib_module_names
    else:
        # For older Python versions that don't have stdlib_module_names
        try:
            spec = importlib.util.find_spec(module_name)
            if spec is None:
                return False
            return 'site-packages' not in spec.origin and 'dist-packages' not in spec.origin
        except (ImportError, AttributeError):
            return False

def classify_imports(imports):
    """Classify imports as standard library or external packages"""
    std_lib = []
    external = []
    
    for module in imports:
        if is_standard_library(module):
            std_lib.append(module)
        else:
            external.append(module)
    
    return std_lib, external

def main():
    print_colored(BLUE, "CryoProtect Requirements Generator")
    print("=" * 60)
    
    repo_root = Path(__file__).parent.parent
    
    # Find all Python files
    python_files = find_python_files(repo_root)
    print(f"Found {len(python_files)} Python files to analyze.")
    
    # Extract imports from all files
    all_imports = []
    for file in python_files:
        file_imports = extract_imports(file)
        all_imports.extend(file_imports)
    
    # Remove duplicates
    all_imports = list(set(all_imports))
    
    # Classify imports
    std_lib, external = classify_imports(all_imports)
    
    print(f"Found {len(std_lib)} standard library imports.")
    print(f"Found {len(external)} external package imports.")
    
    # Common package mappings (import name -> pip package name)
    package_mapping = {
        'flask': 'Flask',
        'sklearn': 'scikit-learn',
        'PIL': 'Pillow',
        'dotenv': 'python-dotenv',
        'psycopg2': 'psycopg2-binary',
        'yaml': 'PyYAML',
        'jwt': 'PyJWT',
        'bs4': 'beautifulsoup4',
        'jinja2': 'Jinja2',
        'sqlalchemy': 'SQLAlchemy',
        'werkzeug': 'Werkzeug',
    }
    
    # Additional known packages for this project
    known_packages = [
        'flask-cors',
        'flask-restful',
        'flask-apispec',
        'flask-mail',
        'apispec',
        'marshmallow',
        'supabase',
        'postgrest',
        'prometheus-client',
        'cryptography',
        'requests',
        'alembic',
    ]
    
    # Generate requirements.txt
    requirements_path = repo_root / "requirements_full.txt"
    print_colored(BLUE, f"\nGenerating {requirements_path}...")
    
    with open(requirements_path, 'w') as f:
        f.write("# Generated requirements for CryoProtect\n")
        f.write("# External packages detected from imports\n")
        
        # Write mapped packages
        for module in sorted(external):
            package_name = package_mapping.get(module, module)
            f.write(f"{package_name}\n")
        
        # Write additional known packages that might not be detected from imports
        f.write("\n# Additional known required packages\n")
        for package in sorted(known_packages):
            if package.lower() not in [p.lower() for p in external]:
                f.write(f"{package}\n")
    
    print_colored(GREEN, f"Generated {requirements_path}")
    print(f"You can now install all dependencies with: pip install -r {requirements_path}")
    
    # Generate install script
    install_script_path = repo_root / "setup_dependencies.sh"
    print_colored(BLUE, f"\nGenerating {install_script_path}...")
    
    with open(install_script_path, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Install all dependencies for CryoProtect\n\n")
        f.write("# ANSI colors\n")
        f.write("GREEN='\\033[0;32m'\n")
        f.write("YELLOW='\\033[0;33m'\n")
        f.write("BLUE='\\033[0;34m'\n")
        f.write("NC='\\033[0m'\n\n")
        f.write("echo -e \"${BLUE}Installing CryoProtect dependencies...${NC}\"\n\n")
        f.write("# System dependencies\n")
        f.write("echo -e \"${YELLOW}Installing system dependencies...${NC}\"\n")
        f.write("sudo dnf install -y postgresql-devel python3-devel gcc\n\n")
        f.write("# Create virtual environment\n")
        f.write("echo -e \"${BLUE}Creating virtual environment...${NC}\"\n")
        f.write("python -m venv venv\n")
        f.write("source venv/bin/activate\n\n")
        f.write("# Install Python dependencies\n")
        f.write("echo -e \"${BLUE}Installing Python dependencies...${NC}\"\n")
        f.write("pip install --upgrade pip\n")
        f.write("pip install -r requirements_full.txt\n\n")
        f.write("echo -e \"${GREEN}Dependencies installed successfully!${NC}\"\n")
        f.write("echo \"You can now run the application with: ./run_local.sh\"\n")
    
    os.chmod(install_script_path, 0o755)  # Make executable
    print_colored(GREEN, f"Generated {install_script_path}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
#!/bin/bash
# Check for and install missing dependencies

# ANSI colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Missing Dependency Installer${NC}"
echo "This script checks for missing dependencies and installs them"
echo "==========================================================================="

# Create a temporary Python script to check for and install missing dependencies
cat > check_deps.py << 'EOF'
import importlib
import subprocess
import sys

# List of required modules
required_modules = [
    # Core Flask
    ("flask", "flask==3.0.2"),
    ("flask_restful", "flask-restful==0.3.10"),
    ("flask_cors", "flask-cors==4.0.0"),
    ("flask_apispec", "flask-apispec==0.11.4"),
    ("flask_mail", "flask-mail==0.9.1"),
    
    # API and serialization
    ("apispec", "apispec==6.3.0"),
    ("marshmallow", "marshmallow==3.21.0"),
    ("webargs", "webargs==8.3.0"),
    
    # Data handling
    ("dotenv", "python-dotenv==1.0.1"),
    ("requests", "requests==2.31.0"),
    ("json", None),  # Standard library
    
    # Database
    ("postgrest", "postgrest==0.13.2"),
    ("supabase", "supabase==2.1.0"),
    ("psycopg2", "psycopg2-binary==2.9.9"),
    ("sqlalchemy", "sqlalchemy==2.0.27"),
    ("alembic", "alembic==1.13.1"),
    
    # Security
    ("jwt", "PyJWT==2.8.0"),
    ("cryptography", "cryptography==42.0.5"),
    
    # Monitoring
    ("prometheus_client", "prometheus-client==0.19.0"),
    
    # Scientific computing
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("pandas", "pandas"),
    ("matplotlib", "matplotlib"),
    ("sklearn", "scikit-learn"),
    ("rdkit", "rdkit-pypi"),  # This is optional
]

# Check for missing modules
missing_modules = []
for module_name, package_name in required_modules:
    try:
        importlib.import_module(module_name)
        print(f"\033[0;32m✓ {module_name} is installed\033[0m")
    except ImportError:
        if package_name:  # Skip standard library modules
            missing_modules.append((module_name, package_name))
            print(f"\033[0;31m✗ {module_name} is missing\033[0m")

# Install missing modules
if missing_modules:
    print(f"\n\033[0;34mInstalling {len(missing_modules)} missing modules...\033[0m")
    
    for module_name, package_name in missing_modules:
        print(f"\033[0;33mInstalling {package_name}...\033[0m")
        
        if package_name == "rdkit-pypi":
            # Ask before installing RDKit as it's large
            response = input("RDKit is a large package. Install it? (y/n): ")
            if response.lower() != 'y':
                print("\033[0;33mSkipping RDKit installation.\033[0m")
                continue
        
        try:
            subprocess.run([sys.executable, "-m", "pip", "install", package_name], check=True)
            print(f"\033[0;32m✓ Successfully installed {package_name}\033[0m")
        except subprocess.CalledProcessError:
            print(f"\033[0;31m✗ Failed to install {package_name}\033[0m")
            
    # Verify installation
    print("\n\033[0;34mVerifying installation...\033[0m")
    still_missing = []
    
    for module_name, package_name in missing_modules:
        try:
            importlib.import_module(module_name)
            print(f"\033[0;32m✓ {module_name} is now installed\033[0m")
        except ImportError:
            still_missing.append(module_name)
            print(f"\033[0;31m✗ {module_name} is still missing\033[0m")
    
    if still_missing:
        print(f"\n\033[0;31m{len(still_missing)} modules are still missing.\033[0m")
        sys.exit(1)
    else:
        print(f"\n\033[0;32mAll required modules are now installed.\033[0m")
else:
    print("\n\033[0;32mAll required modules are already installed.\033[0m")

sys.exit(0)
EOF

# Run the Python script
python check_deps.py

# Clean up
rm check_deps.py

echo -e "${GREEN}Dependency check complete!${NC}"
echo "You can now try running the application."
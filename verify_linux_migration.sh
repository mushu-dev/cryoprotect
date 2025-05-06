#!/bin/bash
# CryoProtect Linux Migration Verification Script

echo "=== CryoProtect Linux Migration Verification ==="
echo "Checking if all required Linux components are in place..."

# Check for shell script equivalents
echo -n "Checking run_app_with_fix.sh: "
if [ -f "run_app_with_fix.sh" ] && [ -x "run_app_with_fix.sh" ]; then
    echo "✅ Present and executable"
else
    echo "❌ Missing or not executable"
    echo "Fix with: chmod +x run_app_with_fix.sh"
fi

echo -n "Checking run_tests_conda.sh: "
if [ -f "run_tests_conda.sh" ] && [ -x "run_tests_conda.sh" ]; then
    echo "✅ Present and executable"
else
    echo "❌ Missing or not executable"
    echo "Fix with: chmod +x run_tests_conda.sh"
fi

echo -n "Checking batch_scripts/start_server.sh: "
if [ -f "batch_scripts/start_server.sh" ] && [ -x "batch_scripts/start_server.sh" ]; then
    echo "✅ Present and executable"
else
    echo "❌ Missing or not executable"
    echo "Fix with: chmod +x batch_scripts/start_server.sh"
fi

# Check environment and files
echo -n "Checking .env file: "
if [ -f ".env" ]; then
    echo "✅ Present"
else
    echo "❌ Missing"
    echo "Fix with: cp .env.template .env"
fi

echo -n "Checking conda: "
if command -v conda &> /dev/null; then
    echo "✅ Installed"
else
    echo "❌ Not installed"
    echo "Please install conda before continuing"
fi

echo -n "Checking PostgreSQL: "
if command -v psql &> /dev/null; then
    echo "✅ Installed"
else
    echo "❌ Not installed"
    echo "Fix with: sudo apt install postgresql postgresql-contrib"
fi

# Check if cryoprotect conda environment exists
echo -n "Checking cryoprotect conda environment: "
if conda env list | grep -q "cryoprotect"; then
    echo "✅ Created"
else
    echo "❌ Not created"
    echo "Fix with: ./setup_environment.sh"
fi

# Try to import rdkit to check if it's properly installed
echo -n "Checking RDKit installation: "
if source $(conda info --base)/etc/profile.d/conda.sh && conda activate cryoprotect && python -c "from rdkit import Chem; print('RDKit is properly installed')" &> /dev/null; then
    echo "✅ Properly installed"
else
    echo "❌ Not installed or not working properly"
    echo "Fix with: conda install -c conda-forge rdkit"
fi

# Check platform-specific code in ip_resolver.py
echo -n "Checking ip_resolver.py for Linux compatibility: "
if grep -q "_resolve_with_dig" ip_resolver.py; then
    echo "✅ Linux-compatible"
else
    echo "❌ May not be Linux-compatible"
    echo "Review ip_resolver.py for platform-specific code"
fi

echo -e "\nVerification complete! If all checks passed, your CryoProtect migration to Linux is successful."
echo "For detailed instructions, see LINUX_MIGRATION_GUIDE.md"
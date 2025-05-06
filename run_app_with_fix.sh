#!/bin/bash
echo "Starting CryoProtect Analyzer with Service Role Authentication Fix..."

# Check if auth_config.py exists
if [ ! -f auth_config.py ]; then
    echo "Service role authentication fix not found."
    echo "Running fix_auth_simple.py to apply the fix..."
    python fix_auth_simple.py
    if [ $? -ne 0 ]; then
        echo "Failed to apply service role authentication fix."
        echo "Please check the error messages above."
        exit 1
    fi
else
    echo "Service role authentication fix already applied."
fi

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Check if RDKit is installed
python -c "from rdkit import Chem; print('RDKit check passed')" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: RDKit is not installed or not working properly."
    echo "Please run setup_environment.sh first."
    exit 1
fi

# Test the authentication
echo "Testing Supabase authentication with service role approach..."
python test_service_role_auth.py
if [ $? -ne 0 ]; then
    echo "Authentication test failed with error code $?"
    echo "Please check the error messages above."
    exit 1
fi

# Run the Flask application with error handling
echo "Starting Flask application..."
python app.py
if [ $? -ne 0 ]; then
    echo "Application exited with error code $?"
    echo "Check the logs above for details."
fi
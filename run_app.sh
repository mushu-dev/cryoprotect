#!/bin/bash
echo "Starting CryoProtect Analyzer..."

# Activate the conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Verify and install required packages
echo "Verifying required packages..."
python verify_packages.py > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Some packages are missing or not accessible."
    echo "Running package installation script..."
    bash ./startup_packages.sh
fi

# Check if RDKit is installed
python -c "from rdkit import Chem; print('RDKit check passed')" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: RDKit is not installed or not working properly."
    echo "Please run setup_environment.sh first."
    exit 1
fi

# Run the Flask application with error handling
echo "Starting Flask application..."
python app.py
if [ $? -ne 0 ]; then
    echo "Application exited with error code $?"
    echo "Check the logs above for details."
fi
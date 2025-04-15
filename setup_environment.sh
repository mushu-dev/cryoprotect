#!/bin/bash
echo "Setting up CryoProtect Analyzer environment..."

# Try to create a new conda environment
conda create -n cryoprotect python=3.9 -y
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment. Trying with mamba..."
    conda install -c conda-forge mamba -y
    mamba create -n cryoprotect python=3.9 -y
    if [ $? -ne 0 ]; then
        echo "Failed to create environment with mamba. Exiting."
        exit 1
    fi
fi

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect

# Try to install RDKit with conda
echo "Installing RDKit..."
conda install -c conda-forge rdkit=2023.9.1 -y
if [ $? -ne 0 ]; then
    echo "Failed to install RDKit with conda. Trying with mamba..."
    mamba install -c conda-forge rdkit=2023.9.1 -y
    if [ $? -ne 0 ]; then
        echo "Failed to install RDKit with mamba. Trying with pip..."
        pip install rdkit
        if [ $? -ne 0 ]; then
            echo "Failed to install RDKit with pip. Please try manual installation."
            exit 1
        fi
    fi
fi

# Install other dependencies
echo "Installing other dependencies..."
pip install -r requirements.txt

# Verify RDKit installation
echo "Verifying RDKit installation..."
python -c "from rdkit import Chem; print('RDKit installation successful!')"
if [ $? -ne 0 ]; then
    echo "RDKit verification failed. Please check the installation."
    exit 1
fi

echo "Environment setup complete! Activate with: conda activate cryoprotect"
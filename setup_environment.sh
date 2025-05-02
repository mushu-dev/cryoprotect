#!/bin/bash
echo "Setting up CryoProtect v2 environment..."

# Check for .env file and create from template if it doesn't exist
if [ ! -f .env ]; then
    echo "Creating .env file from template..."
    cp .env.template .env
    echo "✅ Created .env file from template"
    echo "⚠️ IMPORTANT: You need to update the .env file with your credentials"
    echo "   Please edit the .env file and fill in the required values marked with [YOUR-*]"
else
    echo "✅ .env file already exists"
fi

# Try to create a new conda environment
echo "Creating conda environment..."
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

# Print next steps
echo ""
echo "🎉 Environment setup complete! 🎉"
echo ""
echo "Next steps:"
echo "1. Activate the environment with: conda activate cryoprotect"
echo "2. Edit the .env file to fill in your Supabase credentials and other required values"
echo "3. Run the application with: ./run_app.sh"
echo "4. For debugging, use: python debug_app.py"
echo ""
echo "For RDKit troubleshooting, see: README_RDKit_Troubleshooting.md"
echo "For API documentation, see: README_API.md"
echo "For authentication setup, see: README_Authentication.md"
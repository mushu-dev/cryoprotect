#!/bin/bash
# Script to set up a conda environment with RDKit for unified ChEMBL import

echo "Setting up a conda environment with RDKit for unified ChEMBL import..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "conda not found. Please install conda first."
    echo "You can install Miniconda from https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Create and activate a new conda environment
ENV_NAME="cryoprotect"
echo "Creating conda environment: $ENV_NAME"
conda create -y -n $ENV_NAME python=3.10

# Activate the environment
echo "Activating environment: $ENV_NAME"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install RDKit
echo "Installing RDKit..."
conda install -y -c conda-forge rdkit

# Install other dependencies
echo "Installing other dependencies..."
conda install -y -c conda-forge psycopg2
pip install chembl_webresource_client

# Create an activation script
ACTIVATE_SCRIPT="activate_rdkit_env.sh"
echo "Creating activation script: $ACTIVATE_SCRIPT"
cat > $ACTIVATE_SCRIPT << EOL
#!/bin/bash
# Source this script to activate the conda environment with RDKit
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME
echo "Activated $ENV_NAME environment with RDKit"
EOL
chmod +x $ACTIVATE_SCRIPT

echo "Setup complete!"
echo "To activate the environment, run: source $ACTIVATE_SCRIPT"
echo "Then you can run the unified ChEMBL import with RDKit enabled."
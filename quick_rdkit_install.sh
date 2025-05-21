#!/bin/bash
# Quick setup script for just RDKit in the existing container

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

CONTAINER_NAME="CryoProtect-RDKit-Conda"

echo -e "${BLUE}Setting up RDKit inside the container...${NC}"

# Check if container exists and is running
if ! podman container exists $CONTAINER_NAME || [ "$(podman container inspect -f '{{.State.Running}}' $CONTAINER_NAME 2>/dev/null)" != "true" ]; then
    echo -e "${YELLOW}Container $CONTAINER_NAME not found or not running.${NC}"
    echo "Please run ./quick_conda_container.sh first to create the container."
    exit 1
fi

# Create simple test script that will run inside the container
TEST_SCRIPT="/tmp/test_rdkit.py"
cat > $TEST_SCRIPT << 'EOF'
#!/usr/bin/env python3
# Simple script to verify RDKit installation

import sys
print(f"Python version: {sys.version}")

try:
    import rdkit
    print(f"RDKit version: {rdkit.__version__}")
    
    # Try to create a molecule
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    if mol:
        print("Successfully created molecule from SMILES")
        
    # Calculate a property
    from rdkit.Chem import Descriptors
    mw = Descriptors.MolWt(mol)
    print(f"Molecular weight of ethanol: {mw:.2f}")
    
    print("RDKit is working correctly!")
except ImportError as e:
    print(f"Error importing RDKit: {e}")
    sys.exit(1)
EOF

# Copy test script to container
podman cp $TEST_SCRIPT $CONTAINER_NAME:/app/test_rdkit.py

# Install RDKit and other dependencies directly with pip (faster than conda)
echo -e "${BLUE}Installing RDKit and other dependencies with pip (faster than waiting for conda)...${NC}"
podman exec -it $CONTAINER_NAME bash -c "
    # Downgrade numpy to 1.x for compatibility with rdkit-pypi
    pip uninstall -y numpy && \
    pip install numpy==1.24.4 && \
    pip install rdkit-pypi requests chembl_webresource_client psycopg2-binary && \
    python /app/test_rdkit.py
"

# Clean up temp file
rm $TEST_SCRIPT

echo -e "${GREEN}RDKit installation completed!${NC}"
echo 
echo -e "${BLUE}You can now run scripts with RDKit using:${NC}"
echo "  ./run_with_rdkit.sh your_script.py"
echo
echo -e "Or run the unified ChEMBL import with:${NC}"
echo "  ./run_unified_import_in_container.sh"
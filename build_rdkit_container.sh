#!/bin/bash
# Build and run the dedicated RDKit container for CryoProtect

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

IMAGE_NAME="cryoprotect-rdkit"
CONTAINER_NAME="CryoProtect-RDKit"

echo -e "${BLUE}Building dedicated RDKit container for CryoProtect${NC}"

# Make entrypoint script executable
chmod +x docker-entrypoint-rdkit.sh

# Build the container image
echo -e "${BLUE}Building container image...${NC}"
podman build -t $IMAGE_NAME -f Dockerfile.rdkit .

# Remove existing container if it exists
if podman container exists $CONTAINER_NAME; then
    echo -e "${YELLOW}Removing existing container...${NC}"
    podman rm -f $CONTAINER_NAME
fi

# Create app directory with required modules
TEMP_DIR=$(mktemp -d)
cp rdkit_wrapper.py mock_rdkit_formula.py $TEMP_DIR/

# Create test file
cat > $TEMP_DIR/test_rdkit.py << 'EOF'
#!/usr/bin/env python3
"""
Test RDKit functionality in the container
"""
import sys
print(f"Python version: {sys.version}")

# Import and test wrapper
import rdkit_wrapper

# Get status
status = rdkit_wrapper.get_rdkit_status()
print(f"RDKit Available: {status['rdkit_available']}")
print(f"RDKit Version: {status['rdkit_version']}")
print(f"Visualization Available: {status['visualization_available']}")

# Test with sample molecules
test_molecules = [
    ("Ethanol", "CCO"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
]

for name, smiles in test_molecules:
    print(f"\nTesting {name} ({smiles}):")
    props = rdkit_wrapper.calculate_properties(smiles)
    print(f"  Molecular Weight: {props['molecular_weight']}")
    print(f"  LogP: {props['logp']}")
    print(f"  TPSA: {props['tpsa']}")
    
    # Test fingerprint if available
    if rdkit_wrapper.RDKIT_AVAILABLE:
        fp = rdkit_wrapper.generate_fingerprint(smiles)
        print(f"  Fingerprint generated: {fp is not None}")
        
        # Test 3D coordinates
        mol_3d = rdkit_wrapper.generate_molecule_3d_coordinates(smiles)
        print(f"  3D coordinates generated: {mol_3d is not None}")
EOF

# Run the container interactively
echo -e "${BLUE}Running container...${NC}"
podman run -it --name $CONTAINER_NAME \
           -v $TEMP_DIR:/app \
           -p 5001:5000 \
           $IMAGE_NAME \
           python /app/test_rdkit.py

# Clean up
rm -rf $TEMP_DIR

echo -e "${GREEN}RDKit container test completed.${NC}"
echo
echo -e "${BLUE}Usage:${NC}"
echo "To run commands in the RDKit container:"
echo "  podman exec -it $CONTAINER_NAME python /app/your_script.py"
echo
echo "To use this container in your workflow:"
echo "  1. Mount your code directory: -v /path/to/your/code:/app"
echo "  2. Use the rdkit_wrapper.py module for all RDKit functionality"
echo "  3. For integration with existing code, use the wrapper pattern"
echo
echo -e "${BLUE}Next steps:${NC}"
echo "1. Update your scripts to use rdkit_wrapper.py"
echo "2. Set up CI/CD pipeline to build and test with this container"
echo "3. Use the container for production deployment"
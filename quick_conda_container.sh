#!/bin/bash
# Quick setup script for CryoProtect conda container
# This is a simplified version for testing without all the SELinux complexities

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Setting up quick CryoProtect conda container for testing...${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman"
    exit 1
fi

# Container name with clear indication it's the conda container with RDKit
CONTAINER_NAME="CryoProtect-RDKit-Conda"

# Remove existing container if it exists
if podman container exists $CONTAINER_NAME; then
    echo "Removing existing container..."
    podman rm -f $CONTAINER_NAME
fi

# Create a temp directory for testing
TEST_DIR="/tmp/cryoprotect-test"
mkdir -p "$TEST_DIR"

# Copy minimal files needed for testing
echo "Copying essential files to temp directory..."
cp "$(pwd)/environment.yml" "$TEST_DIR/"
cp "$(pwd)/mock_rdkit.py" "$TEST_DIR/"
cat > "$TEST_DIR/test_import.py" << EOF
# Test imports for CryoProtect dependencies
import sys
print("Python version:", sys.version)

try:
    import rdkit
    print("RDKit version:", rdkit.__version__)
except ImportError:
    print("RDKit not found, trying mock...")
    try:
        import mock_rdkit
        mock_rdkit.create_mock_rdkit()
        import rdkit
        print("Mock RDKit loaded successfully")
    except ImportError:
        print("Failed to load RDKit or mock")

try:
    import flask
    print("Flask version:", flask.__version__)
except ImportError:
    print("Flask not found")

try:
    import numpy
    print("NumPy version:", numpy.__version__)
except ImportError:
    print("NumPy not found")

print("Import test completed")
EOF

# Create and run the container with proper name and port forwarding
echo -e "${BLUE}Creating RDKit conda container...${NC}"
podman run -d \
    --name $CONTAINER_NAME \
    -p 5000:5000 \
    -v "$TEST_DIR:/app:Z" \
    continuumio/miniconda3:4.12.0 \
    tail -f /dev/null

echo -e "${GREEN}RDKit conda container created!${NC}"

# Create conda environment with RDKit
echo -e "${BLUE}Setting up conda environment with RDKit...${NC}"
podman exec -it $CONTAINER_NAME bash -c "
    cd /app && \
    echo 'Creating conda environment with RDKit...' && \
    conda env create -f environment.yml && \
    echo 'Conda environment created.'
"

echo -e "${GREEN}RDKit environment ready!${NC}"
echo
echo -e "${BLUE}Testing RDKit imports...${NC}"
podman exec -it $CONTAINER_NAME bash -c "
    cd /app && \
    conda init bash && \
    source ~/.bashrc && \
    conda activate cryoprotect && \
    python test_import.py || \
    /opt/conda/envs/cryoprotect/bin/python test_import.py
"

# Create a convenient label for easy identification
podman container label set $CONTAINER_NAME purpose="RDKit container for molecular property calculation"
podman container label set $CONTAINER_NAME description="Contains conda environment with RDKit for CryoProtect"

echo
echo -e "${BLUE}Commands for working with the RDKit conda container:${NC}"
echo "1. Run interactive shell with conda activated:"
echo "   podman exec -it $CONTAINER_NAME bash -c 'source ~/.bashrc && conda activate cryoprotect || echo \"Use /opt/conda/envs/cryoprotect/bin/python instead of conda activate\"'"
echo
echo "2. Run a Python script with conda environment:"
echo "   podman exec -it $CONTAINER_NAME bash -c 'cd /app && /opt/conda/envs/cryoprotect/bin/python your_script.py'"
echo
echo "3. Stop the container:"
echo "   podman stop $CONTAINER_NAME"
echo
echo "4. Remove the container when done:"
echo "   podman rm $CONTAINER_NAME"
echo
echo -e "${GREEN}RDKit conda container setup completed!${NC}"
echo -e "${YELLOW}Container name: ${GREEN}$CONTAINER_NAME${NC}"
#!/bin/bash
# Script to run the unified ChEMBL import inside the conda container with RDKit

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Preparing to run unified ChEMBL import in container...${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman"
    exit 1
fi

# Container name with clear indication it's the conda container with RDKit
CONTAINER_NAME="CryoProtect-RDKit-Conda"

# Check if RDKit container exists and is running
if ! podman container exists $CONTAINER_NAME || [ "$(podman container inspect -f '{{.State.Running}}' $CONTAINER_NAME 2>/dev/null)" != "true" ]; then
    echo -e "${YELLOW}RDKit conda container not found or not running.${NC}"
    echo -e "${BLUE}Setting up a new RDKit container...${NC}"

    # Run the quick setup script
    ./quick_conda_container.sh
else
    echo -e "${GREEN}RDKit container is already running.${NC}"
fi

# Copy the unified import script to the container
echo -e "${BLUE}Copying unified import script to RDKit container...${NC}"
podman cp unified_chembl_import.py $CONTAINER_NAME:/app/

# Run the script in the container with conda environment activated
echo -e "${BLUE}Running unified ChEMBL import in RDKit container...${NC}"
podman exec -it $CONTAINER_NAME bash -c "
    cd /app && \
    source ~/.bashrc && \
    conda activate cryoprotect || echo 'Using conda binary path directly' && \
    /opt/conda/envs/cryoprotect/bin/python unified_chembl_import.py \$@
"

echo -e "${GREEN}Import process completed.${NC}"
echo "You can check the logs and output files in the container at /app/"
echo "To access the RDKit container directly, run:"
echo "  podman exec -it $CONTAINER_NAME bash -c 'source ~/.bashrc && conda activate cryoprotect && bash'"
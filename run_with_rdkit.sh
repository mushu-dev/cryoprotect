#!/bin/bash
# Wrapper script to run Python scripts with RDKit support via the container

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

if [ $# -lt 1 ]; then
    echo -e "${RED}Error: No script specified.${NC}"
    echo "Usage: $0 script.py [arguments]"
    exit 1
fi

SCRIPT="$1"
shift
ARGS="$@"

echo -e "${BLUE}Running $SCRIPT with RDKit support in container...${NC}"

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

# Copy the script to the container
echo -e "${BLUE}Copying $SCRIPT to container...${NC}"
podman cp "$SCRIPT" $CONTAINER_NAME:/app/

# Run the script in the container with conda environment activated
echo -e "${BLUE}Running $SCRIPT in container with RDKit support...${NC}"
podman exec -it $CONTAINER_NAME bash -c "
    cd /app && \
    python3 -c 'import sys; print(\"Python path:\", sys.path)' && \
    python3 -c 'import rdkit; print(\"RDKit version:\", rdkit.__version__)' && \
    python3 $(basename "$SCRIPT") $ARGS
"

echo -e "${GREEN}Script execution completed.${NC}"
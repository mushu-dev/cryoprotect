#!/bin/bash
# Run a Python script with the dedicated RDKit container

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

CONTAINER_NAME="CryoProtect-RDKit"

# Check args
if [ $# -lt 1 ]; then
    echo -e "${RED}Error: No script specified${NC}"
    echo "Usage: $0 script.py [args...]"
    exit 1
fi

SCRIPT="$1"
shift
ARGS="$@"

# Check if script exists
if [ ! -f "$SCRIPT" ]; then
    echo -e "${RED}Error: Script $SCRIPT not found${NC}"
    exit 1
fi

# Check if container exists and is running
if ! podman container exists $CONTAINER_NAME; then
    echo -e "${YELLOW}RDKit container doesn't exist. Building it now...${NC}"
    ./build_rdkit_container.sh
elif [ "$(podman container inspect -f '{{.State.Running}}' $CONTAINER_NAME)" != "true" ]; then
    echo -e "${YELLOW}RDKit container exists but is not running. Starting it...${NC}"
    podman start $CONTAINER_NAME
fi

# Create a temporary directory for script execution
TEMP_DIR=$(mktemp -d)

# Copy the script, rdkit_wrapper.py, and mock_rdkit_formula.py
cp "$SCRIPT" $TEMP_DIR/
cp rdkit_wrapper.py mock_rdkit_formula.py $TEMP_DIR/ 2>/dev/null || true

# Copy any additional Python modules in the current directory that might be needed
SCRIPT_DIR=$(dirname "$SCRIPT")
if [ "$SCRIPT_DIR" != "." ]; then
    cp $SCRIPT_DIR/*.py $TEMP_DIR/ 2>/dev/null || true
fi

# Get absolute path and script name
SCRIPT_BASENAME=$(basename "$SCRIPT")

# Execute the script in the container
echo -e "${BLUE}Running $SCRIPT_BASENAME with RDKit...${NC}"
podman exec -it $CONTAINER_NAME bash -c "cd /app && rm -rf * && exit" || true
podman exec -it -w /app $CONTAINER_NAME rm -rf /app/* || true
podman cp $TEMP_DIR/. $CONTAINER_NAME:/app/

podman exec -it $CONTAINER_NAME python /app/$SCRIPT_BASENAME $ARGS

# Cleanup
rm -rf $TEMP_DIR

echo -e "${GREEN}Script execution completed.${NC}"
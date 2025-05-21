#!/bin/bash
# Cleanup script for CryoProtect containers

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Cleaning up CryoProtect containers...${NC}"

# Stop and remove containers
for container in "cryoprotect-rdkit" "cryoprotect-app" "CryoProtect-RDKit-Conda" "cryoprotect-minimal" "cryoprotect-dev" "CryoProtect"; do
    if podman container exists "$container"; then
        echo "Stopping $container..."
        podman stop "$container" >/dev/null 2>&1 || true
        
        echo "Removing $container..."
        podman rm "$container" >/dev/null 2>&1 || true
    fi
done

# Remove temporary directories
echo "Removing temporary directories..."
rm -rf /tmp/cryoprotect-app-data /tmp/cryoprotect-rdkit-data

echo -e "${GREEN}Cleanup complete!${NC}"

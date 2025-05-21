#!/bin/bash
#
# Container Cleanup Script for CryoProtect
# This script stops and removes the CryoProtect containers
#

echo "Cleaning up CryoProtect container environment..."

# Stop and remove containers
echo "Stopping containers..."
podman stop cryoprotect-app cryoprotect-rdkit

echo "Removing containers..."
podman rm cryoprotect-app cryoprotect-rdkit

# Optionally remove the network
if [ "$1" == "--all" ]; then
    echo "Removing cryoprotect-net network..."
    podman network rm cryoprotect-net
fi

echo "Cleanup completed."
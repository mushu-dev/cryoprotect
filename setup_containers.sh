#!/bin/bash
#
# Container Setup Script for CryoProtect
# This script sets up the complete container environment for CryoProtect
# with a Flask app container and a mock RDKit service container.
#

echo "Setting up CryoProtect container environment..."

# Create the network if it doesn't exist
if ! podman network inspect cryoprotect-net &>/dev/null; then
    echo "Creating cryoprotect-net network..."
    podman network create cryoprotect-net
fi

# Stop and remove existing containers if they exist
echo "Stopping existing containers if running..."
podman stop cryoprotect-app cryoprotect-rdkit &>/dev/null
podman rm cryoprotect-app cryoprotect-rdkit &>/dev/null

# Start the RDKit service container
echo "Starting RDKit service container..."
podman run -d --name=cryoprotect-rdkit \
    --network=cryoprotect-net \
    -p 5002:5000 \
    -v "/home/mushu/Projects/CryoProtect/mock_rdkit_service.py:/app/mock_rdkit_service.py:z" \
    python:3.10-slim \
    sh -c "cd /app && pip install flask && python mock_rdkit_service.py"

# Start the app container
echo "Starting app container..."
podman run -d --name=cryoprotect-app \
    --network=cryoprotect-net \
    -p 5001:5000 \
    -v "/home/mushu/Projects/CryoProtect:/app:z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt && python app.py"

echo "Waiting for containers to initialize..."
sleep 3

# Check if containers are running
echo "Checking container status..."
podman ps

echo "Setup completed. Use check_container_setup.sh to verify functionality."
echo "Access the app at: http://localhost:5001"
echo "Access the RDKit service at: http://localhost:5002"
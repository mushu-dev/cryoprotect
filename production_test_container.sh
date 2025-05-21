#!/bin/bash
#
# Production Test Container Setup Script for CryoProtect
# This script sets up a container environment for production testing with real data

echo "Setting up CryoProtect production test environment..."

# Create the network if it doesn't exist
if ! podman network inspect cryoprotect-net &>/dev/null; then
    echo "Creating cryoprotect-net network..."
    podman network create cryoprotect-net
fi

# Stop and remove existing containers if they exist
echo "Stopping existing containers if running..."
podman stop cryoprotect-test &>/dev/null
podman rm cryoprotect-test &>/dev/null
podman stop cryoprotect-rdkit &>/dev/null
podman rm cryoprotect-rdkit &>/dev/null
podman stop cryoprotect-app &>/dev/null
podman rm cryoprotect-app &>/dev/null

# Start the RDKit service container
echo "Starting RDKit service container..."
podman run -d --name=cryoprotect-rdkit --replace \
    --network=cryoprotect-net \
    -p 5002:5000 \
    -v "/home/mushu/Projects/CryoProtect/mock_rdkit_service.py:/app/mock_rdkit_service.py:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install flask && python mock_rdkit_service.py"

# Start the app container
echo "Starting app container..."
podman run -d --name=cryoprotect-app --replace \
    --network=cryoprotect-net \
    -p 5001:5000 \
    -e RDKIT_SERVICE_URL=http://cryoprotect-rdkit:5000 \
    -v "/home/mushu/Projects/CryoProtect:/app:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt && python app.py"

# Wait for the app and RDKit containers to start
echo "Waiting for containers to start..."
sleep 5

# Start the test container
echo "Starting test container..."
podman run -d --name=cryoprotect-test --replace \
    --network=cryoprotect-net \
    -v "/home/mushu/Projects/CryoProtect:/app:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt && pip install pytest pytest-html requests && sleep infinity"

echo "Container setup complete. Test container is now running."
echo "You can exec into the test container with:"
echo "  podman exec -it cryoprotect-test bash"
echo "Use the following commands inside the container to prepare data for testing:"
echo "  python tests/test_data/load_test_data.py --dataset all"
echo "Then run the test workflows with:"
echo "  python test_real_data_workflows.py"
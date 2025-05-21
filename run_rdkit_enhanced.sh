#!/bin/bash
# Run the enhanced RDKit service

# Build the container
echo "Building enhanced RDKit container..."
docker build -f Dockerfile.rdkit-enhanced -t cryoprotect-rdkit-enhanced .

# Run the container
echo "Starting enhanced RDKit service..."
docker run -d --name cryoprotect-rdkit-enhanced -p 5001:5000 cryoprotect-rdkit-enhanced

echo "Enhanced RDKit service started on port 5001"
echo "Health check: curl http://localhost:5001/health"
echo "Service info: curl http://localhost:5001/api/v1/rdkit/info"
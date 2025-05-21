#!/bin/bash
#
# Container Setup Verification Script for CryoProtect
# This script verifies that the container setup is working correctly
# and that the containers can communicate with each other.
#

echo "Checking container status..."
podman ps

echo -e "\nChecking if containers are on the same network..."
podman network inspect cryoprotect-net | grep -A 10 "containers"

echo -e "\nTesting RDKit service health endpoint..."
curl -s http://localhost:5002/health | jq

echo -e "\nTesting app's RDKit integration..."
curl -s http://localhost:5001/rdkit/check | jq

echo -e "\nTesting property calculation via RDKit service directly..."
curl -s http://localhost:5002/calculate/CCO | jq

echo -e "\nTesting property calculation via app..."
curl -s http://localhost:5001/molecule/CCO | jq

# Test a slightly more complex molecule (need to URL encode the SMILES)
echo -e "\nTesting property calculation for a more complex molecule (acetaminophen)..."
curl -s "http://localhost:5002/calculate/CC%28%3DO%29NC1%3DCC%3DC%28C%3DC1%29O" | jq
echo -e "\nTesting via app..."
curl -s "http://localhost:5001/molecule/CC%28%3DO%29NC1%3DCC%3DC%28C%3DC1%29O" | jq

echo -e "\nContainer setup verification completed."
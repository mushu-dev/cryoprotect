#!/bin/bash
#
# Simplified Test Environment Setup for CryoProtect
# This script sets up a more reliable test environment for production testing

set -e

echo "Setting up CryoProtect test environment..."

# Clean up any existing containers
echo "Cleaning up existing containers..."
podman stop cryoprotect-app cryoprotect-rdkit cryoprotect-test 2>/dev/null || true
podman rm cryoprotect-app cryoprotect-rdkit cryoprotect-test 2>/dev/null || true

# Create the network if it doesn't exist
if ! podman network inspect cryoprotect-net &>/dev/null; then
    echo "Creating cryoprotect-net network..."
    podman network create cryoprotect-net
fi

# Create test data directory if it doesn't exist
echo "Setting up test data directory..."
mkdir -p tests/test_data

# Create sample test data if it doesn't exist
if [ ! -f tests/test_data/core_cryoprotectants.json ]; then
    echo "Creating sample test data..."
    
    # Create core_cryoprotectants.json
    cat > tests/test_data/core_cryoprotectants.json << EOF
{
  "molecules": [
    {
      "name": "Dimethyl sulfoxide",
      "smiles": "CS(=O)C",
      "formula": "C2H6OS",
      "properties": {
        "molecular_weight": 78.13,
        "logp": -1.35,
        "tpsa": 17.07,
        "h_bond_donors": 0,
        "h_bond_acceptors": 1
      }
    },
    {
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "formula": "C3H8O3",
      "properties": {
        "molecular_weight": 92.09,
        "logp": -1.76,
        "tpsa": 60.69,
        "h_bond_donors": 3,
        "h_bond_acceptors": 3
      }
    }
  ]
}
EOF

    # Create edge_cases.json
    cat > tests/test_data/edge_cases.json << EOF
{
  "molecules": [
    {
      "name": "Sucrose",
      "smiles": "C(C1C(C(C(C(O1)O)O)O)O)OC2C(C(C(C(O2)CO)O)O)O",
      "formula": "C12H22O11",
      "properties": {
        "molecular_weight": 342.3,
        "logp": -3.76,
        "tpsa": 189.53,
        "h_bond_donors": 8,
        "h_bond_acceptors": 11
      }
    }
  ]
}
EOF

    # Create mixtures.json
    cat > tests/test_data/mixtures.json << EOF
{
  "mixtures": [
    {
      "name": "DMSO-Glycerol Mix",
      "description": "Common cryoprotectant mixture",
      "components": [
        {
          "molecule_name": "Dimethyl sulfoxide",
          "concentration": 70,
          "concentration_unit": "%v/v"
        },
        {
          "molecule_name": "Glycerol",
          "concentration": 30,
          "concentration_unit": "%v/v"
        }
      ]
    }
  ]
}
EOF

    echo "Test data files created."
fi

# Start the RDKit service container
echo "Starting RDKit service container..."
podman run -d --name=cryoprotect-rdkit --replace \
    --network=cryoprotect-net \
    -p 5002:5000 \
    -v "$(pwd)/mock_rdkit_service.py:/app/mock_rdkit_service.py:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install flask && python mock_rdkit_service.py"

echo "Waiting for RDKit service to start..."
sleep 5

# Start the app container
echo "Starting app container..."
podman run -d --name=cryoprotect-app --replace \
    --network=cryoprotect-net \
    -p 5001:5000 \
    -e RDKIT_SERVICE_URL=http://cryoprotect-rdkit:5000 \
    -v "$(pwd):/app:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt && python app.py"

echo "Waiting for app to start..."
sleep 5

# Start the test container
echo "Starting test container..."
podman run -d --name=cryoprotect-test --replace \
    --network=cryoprotect-net \
    -v "$(pwd):/app:Z" \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements_essential.txt requests pytest pytest-html && sleep infinity"

echo "Test environment setup complete."
echo 
echo "The following containers are now running:"
podman ps
echo
echo "You can test the environment with:"
echo "podman exec -it cryoprotect-test python /app/production_workflow_test.py --app http://cryoprotect-app:5000 --rdkit http://cryoprotect-rdkit:5000"
#!/bin/bash
# Simplified Container Setup for CryoProtect
# This script sets up the two essential containers with consistent naming,
# using a pre-built Python image with RDKit for faster setup.

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================="
echo "CryoProtect Simple Container Setup"
echo -e "==========================================${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman podman-compose"
    exit 1
fi

# Create the network if it doesn't exist
if ! podman network ls | grep -q "cryoprotect-net"; then
    echo -e "${YELLOW}Network 'cryoprotect-net' not found. Creating...${NC}"
    podman network create --ipv6=false cryoprotect-net
    echo -e "${GREEN}Created network 'cryoprotect-net'${NC}"
fi

# Clean up any existing containers
echo -e "${BLUE}Cleaning up any existing containers...${NC}"
for container in "cryoprotect-rdkit" "cryoprotect-app" "CryoProtect-RDKit-Conda" "cryoprotect-minimal" "cryoprotect-dev" "CryoProtect"; do
    if podman container exists "$container"; then
        echo "Stopping $container..."
        podman stop "$container" >/dev/null 2>&1 || true
        
        echo "Removing $container..."
        podman rm "$container" >/dev/null 2>&1 || true
    fi
done

# Step 1: Set up the App container (simple Python Flask)
echo -e "${BLUE}Setting up Flask App container (cryoprotect-app)...${NC}"

# Create the app directory
APP_DIR="/tmp/cryoprotect-app-data"
mkdir -p "$APP_DIR"

# Copy essential files to temp directory
echo "Copying essential files..."
mkdir -p "$APP_DIR/logs" "$APP_DIR/data" "$APP_DIR/cache"

# Create a test Flask app
echo "Creating test app..."
cat > "$APP_DIR/test_app.py" << EOF
#!/usr/bin/env python3
"""
Test CryoProtect Flask Application
"""

import os
import requests
from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'version': '1.0.0',
        'environment': os.environ.get('FLASK_ENV', 'development')
    })

@app.route('/rdkit/check', methods=['GET'])
def check_rdkit_service():
    """Check if RDKit service is available"""
    try:
        # Try to connect to RDKit service container
        rdkit_available = False
        
        try:
            response = requests.get("http://cryoprotect-rdkit:5000/health", timeout=2)
            if response.status_code == 200:
                rdkit_available = True
                rdkit_data = response.json()
            else:
                rdkit_data = {"error": "Service returned error code"}
        except requests.exceptions.RequestException as e:
            rdkit_data = {"error": str(e)}
            
        return jsonify({
            'rdkit_service': {
                'available': rdkit_available,
                'status': rdkit_data
            }
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/molecule/<smiles>', methods=['GET'])
def get_molecule(smiles):
    """Get molecule properties by calling RDKit service"""
    try:
        try:
            # Call RDKit service
            response = requests.get(f"http://cryoprotect-rdkit:5000/calculate/{smiles}", timeout=5)
            if response.status_code == 200:
                return jsonify(response.json())
            else:
                return jsonify({
                    'error': f"RDKit service error: {response.status_code}",
                    'smiles': smiles
                }), response.status_code
        except requests.exceptions.RequestException as e:
            return jsonify({
                'error': f"RDKit service unreachable: {str(e)}",
                'smiles': smiles
            }), 503
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF

# Create requirements file
cat > "$APP_DIR/requirements.txt" << EOF
flask==3.0.2
requests==2.32.3
EOF

# Create the container
echo "Creating Flask App container..."
podman run -d \
    --name "cryoprotect-app" \
    --network=cryoprotect-net \
    -v "$APP_DIR:/app:z" \
    -p 5001:5000 \
    -e FLASK_ENV=development \
    -e FLASK_DEBUG=1 \
    --security-opt label=disable \
    python:3.10-slim \
    sh -c "cd /app && pip install -r requirements.txt && python test_app.py"

echo -e "${GREEN}Flask App container created and started.${NC}"

# Step 2: Set up the RDKit container (using Python with RDKit pre-installed)
echo -e "${BLUE}Setting up RDKit service container (cryoprotect-rdkit)...${NC}"

# Create the RDKit directory
RDKIT_DIR="/tmp/cryoprotect-rdkit-data"
mkdir -p "$RDKIT_DIR"

# Create a simple RDKit service app
echo "Creating RDKit service app..."
cat > "$RDKIT_DIR/rdkit_service.py" << EOF
#!/usr/bin/env python3
"""
RDKit Service for CryoProtect
"""

import os
from flask import Flask, jsonify, request

# Import RDKit or use mock if not available
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
    RDKIT_VERSION = Chem.__version__
except ImportError:
    RDKIT_AVAILABLE = False
    RDKIT_VERSION = "Not available"

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy' if RDKIT_AVAILABLE else 'unhealthy',
        'rdkit_available': RDKIT_AVAILABLE,
        'rdkit_version': RDKIT_VERSION,
        'environment': os.environ.get('FLASK_ENV', 'development')
    })

@app.route('/calculate/<smiles>', methods=['GET'])
def calculate_properties(smiles):
    """Calculate molecular properties using RDKit"""
    if not RDKIT_AVAILABLE:
        return jsonify({'error': 'RDKit not available'}), 503
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
            
        # Calculate properties
        properties = {
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'logp': round(Descriptors.MolLogP(mol), 2),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'h_donors': Lipinski.NumHDonors(mol),
            'h_acceptors': Lipinski.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF

# Create requirements file
cat > "$RDKIT_DIR/requirements.txt" << EOF
flask==3.0.2
rdkit==2022.09.5
EOF

# Create a temporary Dockerfile for RDKit
cat > "$RDKIT_DIR/Dockerfile.rdkit" << EOF
FROM python:3.10-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        build-essential \
        libpq-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy service app
COPY rdkit_service.py .

# Expose port
EXPOSE 5000

# Run the service
CMD ["python", "rdkit_service.py"]
EOF

# Build and run the RDKit container
echo "Building and running RDKit container..."
podman build -t cryoprotect-rdkit-image -f "$RDKIT_DIR/Dockerfile.rdkit" "$RDKIT_DIR"
podman run -d \
    --name "cryoprotect-rdkit" \
    --network=cryoprotect-net \
    -p 5002:5000 \
    cryoprotect-rdkit-image

echo -e "${GREEN}RDKit container created and started.${NC}"

# Create a status check script
cat > "$(pwd)/check_containers.sh" << EOF
#!/bin/bash
# Status check script for CryoProtect containers

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "\${BLUE}CryoProtect Container Status\${NC}"
echo "================================="

# Check App container
echo -e "\${BLUE}App Container (cryoprotect-app):\${NC}"
if podman container exists "cryoprotect-app"; then
    if podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "  Status: \${GREEN}Running\${NC}"
        echo "  API URL: http://localhost:5001"
        echo "  Health Endpoint: http://localhost:5001/health"
    else
        echo -e "  Status: \${YELLOW}Stopped\${NC}"
    fi
else
    echo -e "  Status: \${RED}Not Created\${NC}"
fi

# Check RDKit container
echo -e "\${BLUE}RDKit Container (cryoprotect-rdkit):\${NC}"
if podman container exists "cryoprotect-rdkit"; then
    if podman container inspect "cryoprotect-rdkit" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "  Status: \${GREEN}Running\${NC}"
        echo "  API URL: http://localhost:5002"
        echo "  Health Endpoint: http://localhost:5002/health"
        
        # Check RDKit status
        echo "  RDKit Status:"
        curl -s http://localhost:5002/health | jq || echo "  Could not get RDKit status"
    else
        echo -e "  Status: \${YELLOW}Stopped\${NC}"
    fi
else
    echo -e "  Status: \${RED}Not Created\${NC}"
fi

# Check container communication
echo -e "\${BLUE}Container Communication:\${NC}"
# Check if app can reach RDKit
APP_CAN_REACH_RDKIT=false
if podman container exists "cryoprotect-app" && podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true"; then
    if curl -s http://localhost:5001/rdkit/check | jq -r '.rdkit_service.available' | grep -q "true"; then
        APP_CAN_REACH_RDKIT=true
        echo -e "  App → RDKit: \${GREEN}OK\${NC}"
    else
        echo -e "  App → RDKit: \${RED}Failed\${NC}"
    fi
else
    echo -e "  App → RDKit: \${YELLOW}App not running\${NC}"
fi

# Print overall status
echo
echo -e "\${BLUE}Overall Status:\${NC}"
if podman container exists "cryoprotect-app" && \
   podman container exists "cryoprotect-rdkit" && \
   podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true" && \
   podman container inspect "cryoprotect-rdkit" --format '{{.State.Running}}' | grep -q "true" && \
   [ "\$APP_CAN_REACH_RDKIT" = "true" ]; then
    echo -e "  \${GREEN}All containers running and communicating properly\${NC}"
else
    echo -e "  \${RED}Issues detected with container setup\${NC}"
fi
EOF
chmod +x "$(pwd)/check_containers.sh"

# Create a cleanup script
cat > "$(pwd)/clean_containers.sh" << EOF
#!/bin/bash
# Cleanup script for CryoProtect containers

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "\${BLUE}Cleaning up CryoProtect containers...\${NC}"

# Stop and remove containers
for container in "cryoprotect-rdkit" "cryoprotect-app" "CryoProtect-RDKit-Conda" "cryoprotect-minimal" "cryoprotect-dev" "CryoProtect"; do
    if podman container exists "\$container"; then
        echo "Stopping \$container..."
        podman stop "\$container" >/dev/null 2>&1 || true
        
        echo "Removing \$container..."
        podman rm "\$container" >/dev/null 2>&1 || true
    fi
done

# Remove temporary directories
echo "Removing temporary directories..."
rm -rf /tmp/cryoprotect-app-data /tmp/cryoprotect-rdkit-data

echo -e "\${GREEN}Cleanup complete!\${NC}"
EOF
chmod +x "$(pwd)/clean_containers.sh"

echo -e "${GREEN}Setup complete!${NC}"
echo
echo -e "${BLUE}Available services:${NC}"
echo "1. RDKit Service: http://localhost:5002"
echo "   - Health Check: http://localhost:5002/health"
echo "   - Calculate Properties: http://localhost:5002/calculate/CCO"
echo
echo "2. Flask App: http://localhost:5001"
echo "   - Health Check: http://localhost:5001/health"
echo "   - RDKit Check: http://localhost:5001/rdkit/check"
echo "   - Test Molecule: http://localhost:5001/molecule/CCO"
echo
echo "3. Container Communication:"
echo "   - App container can reach RDKit container via: http://cryoprotect-rdkit:5000"
echo "   - RDKit container can reach App container via: http://cryoprotect-app:5000"
echo
echo "4. Management Scripts:"
echo "   - Check containers: ./check_containers.sh"
echo "   - Clean up containers: ./clean_containers.sh"
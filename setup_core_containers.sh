#!/bin/bash
# CryoProtect Core Container Setup
# This script sets up the two essential containers with consistent naming:
# 1. cryoprotect-rdkit: For RDKit calculations
# 2. cryoprotect-app: For the Flask application

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================="
echo "CryoProtect Core Container Setup"
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

# Step 1: Set up the RDKit container
echo -e "${BLUE}Setting up RDKit container (cryoprotect-rdkit)...${NC}"

# Check for existing containers with old names
if podman container exists "CryoProtect-RDKit-Conda"; then
    echo "Found container with old name (CryoProtect-RDKit-Conda)"
    
    # Check if it's running
    if podman container inspect "CryoProtect-RDKit-Conda" --format '{{.State.Running}}' | grep -q "true"; then
        echo "Stopping container with old name..."
        podman stop "CryoProtect-RDKit-Conda"
    fi
    
    echo "Removing container with old name..."
    podman rm "CryoProtect-RDKit-Conda"
fi

# Check for our standardized container name
if podman container exists "cryoprotect-rdkit"; then
    echo -e "${GREEN}RDKit container already exists with standardized name.${NC}"
    
    # Check if it's running
    if podman container inspect "cryoprotect-rdkit" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "${GREEN}RDKit container is already running.${NC}"
    else
        echo "Starting RDKit container..."
        podman start "cryoprotect-rdkit"
    fi
else
    echo -e "${YELLOW}RDKit container not found. Creating...${NC}"
    
    # Create a temporary directory for RDKit data
    RDKIT_DIR="/tmp/cryoprotect-rdkit-data"
    mkdir -p "$RDKIT_DIR"
    
    # Create container with standardized name
    podman run -d \
        --name "cryoprotect-rdkit" \
        --network=cryoprotect-net \
        -v "$RDKIT_DIR:/data:Z" \
        -p 5002:5000 \
        --security-opt label=type:container_t \
        continuumio/miniconda3:4.12.0 \
        tail -f /dev/null
    
    echo -e "${GREEN}RDKit container created.${NC}"
    
    # Copy environment.yml to the temporary directory
    cp "$(pwd)/environment.yml" "$RDKIT_DIR/"
    
    # Install dependencies in the container
    echo "Installing conda environment (this may take several minutes)..."
    podman exec -it "cryoprotect-rdkit" bash -c "
        echo 'Creating conda environment from environment.yml...' && \
        conda env create -f /data/environment.yml && \
        echo 'Cleaning conda cache to reduce container size...' && \
        conda clean -afy && \
        echo 'Conda environment created successfully.'
    "
    
    echo -e "${GREEN}RDKit environment installed.${NC}"
fi

# Step 2: Set up the App container
echo -e "${BLUE}Setting up Flask App container (cryoprotect-app)...${NC}"

# Check for containers with inconsistent names
for OLD_NAME in "cryoprotect-minimal" "cryoprotect-dev" "CryoProtect"; do
    if podman container exists "$OLD_NAME"; then
        echo "Found container with inconsistent name ($OLD_NAME)"
        
        # Check if it's running
        if podman container inspect "$OLD_NAME" --format '{{.State.Running}}' | grep -q "true"; then
            echo "Stopping container..."
            podman stop "$OLD_NAME"
        fi
        
        echo "Removing container..."
        podman rm "$OLD_NAME"
    fi
done

# Create the app directory
APP_DIR="/tmp/cryoprotect-app-data"
mkdir -p "$APP_DIR"

# Create or start the app container
if podman container exists "cryoprotect-app"; then
    echo -e "${GREEN}App container already exists.${NC}"
    
    # Check if it's running
    if podman container inspect "cryoprotect-app" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "${GREEN}App container is already running.${NC}"
    else
        echo "Starting App container..."
        podman start "cryoprotect-app"
    fi
else
    echo "Creating Flask App container..."
    
    # Copy essential files to temp directory
    echo "Copying essential files..."
    cp -r $(pwd)/*.py "$APP_DIR/" 2>/dev/null || true
    cp -r $(pwd)/api "$APP_DIR/" 2>/dev/null || true
    mkdir -p "$APP_DIR/static" "$APP_DIR/templates" "$APP_DIR/logs" "$APP_DIR/data" "$APP_DIR/cache"
    
    # Create the container
    podman run -d \
        --name "cryoprotect-app" \
        --network=cryoprotect-net \
        -v "$APP_DIR:/app:z" \
        -p 5001:5000 \
        -e FLASK_APP=test_app.py \
        -e FLASK_ENV=development \
        -e FLASK_DEBUG=1 \
        -e LOG_LEVEL=DEBUG \
        --security-opt label=disable \
        python:3.10-slim \
        sh -c "cd /app && pip install -q flask flask-restful flask-cors python-dotenv requests && tail -f /dev/null"
    
    echo -e "${GREEN}Flask App container created.${NC}"
fi

# Create a test app that communicates with the RDKit container
echo "Creating test app that communicates with RDKit container..."
cat > "$APP_DIR/test_app.py" << EOF
#!/usr/bin/env python3
"""
Test CryoProtect Flask Application with RDKit Integration
"""

import os
import requests
from flask import Flask, jsonify, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

# RDKit service URL (container name as hostname in the podman network)
RDKIT_SERVICE_URL = "http://cryoprotect-rdkit:5000"

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
        # First try directly in the app container
        local_rdkit = False
        rdkit_version = "Not installed"
        
        try:
            from rdkit import Chem
            rdkit_version = Chem.__version__
            local_rdkit = True
        except ImportError:
            pass
            
        # Try to connect to RDKit service container
        service_available = False
        service_status = "Unknown"
        
        try:
            # This will only work once we have the RDKit service running
            # response = requests.get(f"{RDKIT_SERVICE_URL}/health", timeout=2)
            # service_available = response.status_code == 200
            # service_status = response.json() if service_available else "Service unreachable"
            
            # For now, just check if the container is reachable
            import socket
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.settimeout(2)
            s.connect(("cryoprotect-rdkit", 5000))
            s.close()
            service_available = True
            service_status = "Container reachable"
        except Exception as e:
            service_status = str(e)
        
        return jsonify({
            'local_rdkit': local_rdkit,
            'rdkit_version': rdkit_version,
            'rdkit_service': {
                'available': service_available,
                'status': service_status
            }
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/molecule/<smiles>', methods=['GET'])
def get_molecule(smiles):
    """Get molecule properties using local implementation"""
    try:
        properties = {
            'smiles': smiles,
            'properties': {
                'molecular_weight': 0.0,
                'logp': 0.0,
                'h_donors': 0,
                'h_acceptors': 0
            }
        }
        
        # Try to calculate using local RDKit if available
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                properties['properties'] = {
                    'molecular_weight': round(Descriptors.MolWt(mol), 2),
                    'logp': round(Descriptors.MolLogP(mol), 2),
                    'h_donors': Lipinski.NumHDonors(mol),
                    'h_acceptors': Lipinski.NumHAcceptors(mol)
                }
        except ImportError:
            properties['calculation_method'] = 'default'
        
        return jsonify(properties)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF

# Create a simple RDKit service app in the RDKit container
echo "Creating RDKit service app..."
mkdir -p "$RDKIT_DIR/app"
cat > "$RDKIT_DIR/app/rdkit_service.py" << EOF
#!/usr/bin/env python3
"""
RDKit Service for CryoProtect
"""

import os
from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    try:
        # Check if RDKit is available
        from rdkit import Chem
        return jsonify({
            'status': 'healthy',
            'rdkit_version': Chem.__version__,
            'environment': os.environ.get('FLASK_ENV', 'development')
        })
    except ImportError:
        return jsonify({
            'status': 'unhealthy',
            'error': 'RDKit not available'
        }), 500

@app.route('/calculate/<smiles>', methods=['GET'])
def calculate_properties(smiles):
    """Calculate molecular properties using RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
            
        # Calculate properties
        properties = {
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'exact_mass': round(Descriptors.ExactMolWt(mol), 4),
            'logp': round(Descriptors.MolLogP(mol), 2),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'h_donors': Lipinski.NumHDonors(mol),
            'h_acceptors': Lipinski.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Lipinski.NumAromaticRings(mol),
            'heavy_atoms': mol.GetNumHeavyAtoms()
        }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties
        })
    except ImportError:
        return jsonify({'error': 'RDKit not available'}), 500
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF

# Install Flask in the RDKit container and start the service
echo "Installing dependencies in RDKit container..."
podman exec -it "cryoprotect-rdkit" bash -c "
    export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH && \
    pip install flask flask-cors requests && \
    echo 'Dependencies installed'
"

# Start the Flask app in the App container
echo -e "${BLUE}Starting Flask app in app container...${NC}"
podman exec -d "cryoprotect-app" sh -c "cd /app && python test_app.py"

# Start the RDKit service in the RDKit container
echo -e "${BLUE}Starting RDKit service in RDKit container...${NC}"
podman exec -d "cryoprotect-rdkit" bash -c "
    export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH && \
    cd /data/app && python rdkit_service.py
"

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
rm -rf /tmp/cryoprotect-app-data /tmp/cryoprotect-rdkit-data /tmp/cryoprotect-test

echo -e "\${GREEN}Cleanup complete!\${NC}"
EOF
chmod +x "$(pwd)/clean_containers.sh"

echo -e "${GREEN}Setup complete!${NC}"
echo
echo -e "${BLUE}Available services:${NC}"
echo "1. RDKit Container: podman exec -it cryoprotect-rdkit bash"
echo "   - RDKit version: $(podman exec -it cryoprotect-rdkit bash -c "export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH && python -c 'import rdkit; print(rdkit.__version__)'") "
echo "   - RDKit API: http://localhost:5002/health"
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
echo "4. Cleanup:"
echo "   - To clean up all containers: ./clean_containers.sh"
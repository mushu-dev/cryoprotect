#!/bin/bash
# Minimal Container Setup for CryoProtect
# This script sets up only the necessary containers for development and testing

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}CryoProtect Minimal Container Setup${NC}"
echo "Setting up only the necessary containers for development and testing"
echo "=================================================================="

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

# Check for the RDKit container
RDKIT_CONTAINER="CryoProtect-RDKit-Conda"

if podman container exists $RDKIT_CONTAINER; then
    echo -e "${GREEN}RDKit container already exists.${NC}"
    
    # Check if it's running
    if podman container inspect $RDKIT_CONTAINER --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "${GREEN}RDKit container is already running.${NC}"
    else
        echo "Starting RDKit container..."
        podman start $RDKIT_CONTAINER
    fi
else
    echo -e "${YELLOW}RDKit container not found. Creating...${NC}"
    ./quick_conda_container.sh
fi

# Create the app directory
APP_DIR="/tmp/cryoprotect-app"
mkdir -p "$APP_DIR"

# Set up minimal container for the Flask app
echo -e "${BLUE}Setting up minimal Flask app container...${NC}"

# Create container if it doesn't exist
if ! podman container exists "cryoprotect-minimal"; then
    echo "Creating minimal Flask app container..."
    
    # Copy essential files to temp directory
    echo "Copying essential files..."
    cp -r $(pwd)/*.py "$APP_DIR/" 2>/dev/null || true
    cp -r $(pwd)/api "$APP_DIR/" 2>/dev/null || true
    mkdir -p "$APP_DIR/static" "$APP_DIR/templates" "$APP_DIR/logs" "$APP_DIR/data" "$APP_DIR/cache"
    
    # Create the container
    podman run -d \
        --name cryoprotect-minimal \
        --network=cryoprotect-net \
        -v "$APP_DIR:/app:z" \
        -p 5001:5000 \
        -e FLASK_APP=minimal_app.py \
        -e FLASK_ENV=development \
        -e FLASK_DEBUG=1 \
        -e LOG_LEVEL=DEBUG \
        --security-opt label=disable \
        python:3.10-slim \
        sh -c "cd /app && pip install -q flask flask-restful flask-cors python-dotenv && tail -f /dev/null"
    
    echo -e "${GREEN}Minimal Flask app container created.${NC}"
else
    echo -e "${GREEN}Minimal Flask app container already exists.${NC}"
    
    # Check if it's running
    if podman container inspect "cryoprotect-minimal" --format '{{.State.Running}}' | grep -q "true"; then
        echo -e "${GREEN}Minimal Flask app container is running.${NC}"
    else
        echo "Starting minimal Flask app container..."
        podman start "cryoprotect-minimal"
    fi
fi

# Create the minimal app file if it doesn't exist
if [ ! -f "$APP_DIR/minimal_app.py" ]; then
    echo "Creating minimal Flask app..."
    cat > "$APP_DIR/minimal_app.py" << EOF
#!/usr/bin/env python3
"""
Minimal CryoProtect Flask Application for Testing
"""

import os
from flask import Flask, jsonify, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'version': '1.0.0',
        'environment': os.environ.get('FLASK_ENV', 'development')
    })

@app.route('/api/v1/molecule/<smiles>', methods=['GET'])
def get_molecule(smiles):
    """Get molecule properties"""
    try:
        # Create a simple response
        return jsonify({
            'smiles': smiles,
            'properties': {
                'molecular_weight': 100.0,
                'logp': 1.0,
                'h_donors': 1,
                'h_acceptors': 1
            }
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/status', methods=['GET'])
def rdkit_status():
    """Check RDKit status"""
    try:
        # Try to import RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            rdkit_available = True
            rdkit_version = Chem.__version__
        except ImportError:
            rdkit_available = False
            rdkit_version = "Not available"
            
        return jsonify({
            'rdkit_available': rdkit_available,
            'rdkit_version': rdkit_version
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF
fi

# Start the minimal app
echo -e "${BLUE}Starting minimal Flask app...${NC}"
podman exec -it cryoprotect-minimal sh -c "cd /app && python minimal_app.py &"

echo -e "${GREEN}Setup complete!${NC}"
echo
echo -e "${BLUE}Available services:${NC}"
echo "1. RDKit Container: podman exec -it $RDKIT_CONTAINER bash"
echo "2. Flask App: http://localhost:5001"
echo "3. Health Check: http://localhost:5001/health"
echo
echo -e "${YELLOW}Note: For a full setup, use podman-compose with the appropriate compose file.${NC}"
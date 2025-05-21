#!/bin/bash
# Simplified container setup that uses mock RDKit 
# This allows for quick testing without the full conda environment

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Setting up simplified CryoProtect container...${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman"
    exit 1
fi

# Remove existing container if it exists
if podman container exists CryoProtect-Slim; then
    echo "Removing existing container..."
    podman rm -f CryoProtect-Slim
fi

# Create a temp directory for app files
TEMP_DIR="/tmp/cryoprotect-slim"
mkdir -p $TEMP_DIR

# Copy essential files
echo "Copying essential files..."
cp $(pwd)/mock_rdkit.py $TEMP_DIR/
cp $(pwd)/requirements_essential.txt $TEMP_DIR/

# Create a simple test script
cat > $TEMP_DIR/test_app.py << EOF
#!/usr/bin/env python3
"""
Simple test application for CryoProtect without RDKit.
"""
import os
import sys

# Try to import RDKit, use mock if not available
try:
    import rdkit
    print(f"RDKit version: {rdkit.__version__}")
except ImportError:
    print("RDKit not available, using mock...")
    import mock_rdkit
    mock_rdkit.create_mock_rdkit()
    import rdkit
    print("Mock RDKit loaded successfully")

# Import other dependencies
try:
    import flask
    print(f"Flask version: {flask.__version__}")
except ImportError:
    print("Flask not available")

try:
    import numpy as np
    print(f"NumPy version: {np.__version__}")
except ImportError:
    print("NumPy not available")

# Simple minimal Flask app
from flask import Flask, jsonify

app = Flask(__name__)

@app.route('/')
def home():
    return jsonify({
        "status": "ok",
        "message": "CryoProtect Minimal API Running",
        "python_version": sys.version,
        "working_directory": os.getcwd()
    })

@app.route('/health')
def health():
    return jsonify({"status": "ok"})

if __name__ == '__main__':
    print("Starting minimal Flask application...")
    app.run(host='0.0.0.0', port=5000, debug=True)
EOF

# Create a minimal Dockerfile
cat > $TEMP_DIR/Dockerfile << EOF
FROM python:3.10-slim

WORKDIR /app

COPY requirements_essential.txt .
RUN pip install --no-cache-dir -r requirements_essential.txt

COPY . .

ENV FLASK_APP=test_app.py
ENV PYTHONPATH=/app

EXPOSE 5000

CMD ["python", "test_app.py"]
EOF

# Build and run the container
echo -e "${BLUE}Building and running simplified container...${NC}"
cd $TEMP_DIR
podman build -t cryoprotect-slim .
podman run -d --name CryoProtect-Slim -p 5001:5000 cryoprotect-slim

# Create a helper script to interact with the container
cat > $(pwd)/run_slim_container.sh << EOF
#!/bin/bash
# Helper script to run commands in the simplified container

# Check if container exists
if ! podman container exists CryoProtect-Slim; then
    echo "Error: CryoProtect-Slim container does not exist."
    echo "Please run simplified_container_setup.sh first."
    exit 1
fi

# Check if container is running
if ! podman container inspect CryoProtect-Slim --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting CryoProtect-Slim container..."
    podman start CryoProtect-Slim
fi

# Execute the command in the container
if [ \$# -eq 0 ]; then
    # If no arguments, start an interactive shell
    echo "Starting interactive shell in container..."
    podman exec -it CryoProtect-Slim bash
else
    # Otherwise, run the specified command
    echo "Running command in container: \$@"
    podman exec -it CryoProtect-Slim bash -c "\$*"
fi
EOF
chmod +x $(pwd)/run_slim_container.sh

echo -e "${GREEN}Simplified container setup complete!${NC}"
echo 
echo -e "${BLUE}Container is running at:${NC} http://localhost:5000"
echo
echo "You can access the container using the run_slim_container.sh script:"
echo "  ./run_slim_container.sh"
echo
echo "To check the container logs:"
echo "  podman logs CryoProtect-Slim"
echo
echo "To stop the container:"
echo "  podman stop CryoProtect-Slim"
echo
echo "To start the container again:"
echo "  podman start CryoProtect-Slim"
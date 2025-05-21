#!/bin/bash
# Safely mount the project directory into the CryoProtect container
# This script uses a temporary directory with proper SELinux context

set -e

TEMP_APP_DIR="/tmp/cryoprotect-app"
PROJ_DIR="$(pwd)"

# Check if container exists
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_conda_container.sh first."
    exit 1
fi

# Prepare temporary directory with proper SELinux context
echo "Setting up temporary directory with proper SELinux context..."
mkdir -p "$TEMP_APP_DIR"
sudo chcon -Rt container_file_t "$TEMP_APP_DIR"

# Copy important files to temporary directory
echo "Copying project files to temporary directory..."
rsync -av --exclude="*/__pycache__" --exclude="*.pyc" \
    "$PROJ_DIR/"*.py \
    "$PROJ_DIR/api" \
    "$PROJ_DIR/database" \
    "$PROJ_DIR/migrations" \
    "$PROJ_DIR/chembl" \
    "$PROJ_DIR/pubchem" \
    "$TEMP_APP_DIR/"

mkdir -p "$TEMP_APP_DIR/reports"
mkdir -p "$TEMP_APP_DIR/data"

# Recreate the container with the app directory mounted
echo "Stopping existing container..."
podman stop CryoProtect
podman rm CryoProtect

echo "Creating new container with project mount..."
podman run -d \
    --name CryoProtect \
    -v "$TEMP_APP_DIR:/app:Z" \
    -p 5000:5000 \
    --security-opt label=type:container_t \
    -w /app \
    continuumio/miniconda3:4.12.0 \
    tail -f /dev/null

echo "Project directory mounted securely in container."
echo "Use './run_in_cryoprotect.sh' to run commands in the container."

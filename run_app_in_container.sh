#!/bin/bash
# Run the CryoProtect application in the container with conda environment

# Check if container exists and is running
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_conda_container.sh first."
    exit 1
fi

if ! podman container inspect CryoProtect --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting CryoProtect container..."
    podman start CryoProtect
fi

# Create a directory for mounted app code (with proper SELinux context)
TEMP_APP_DIR="/tmp/cryoprotect-app"
mkdir -p "$TEMP_APP_DIR"
sudo chcon -Rt container_file_t "$TEMP_APP_DIR"

# Copy essential application files to the temporary directory
echo "Copying application files to temporary directory..."
cp -r /home/mushu/Projects/CryoProtect/*.py "$TEMP_APP_DIR/"
cp -r /home/mushu/Projects/CryoProtect/api "$TEMP_APP_DIR/"
cp -r /home/mushu/Projects/CryoProtect/database "$TEMP_APP_DIR/"
mkdir -p "$TEMP_APP_DIR/reports"

# Stop the container and recreate it with the app mount
echo "Recreating container with application mount..."
podman rm -f CryoProtect

# Recreate the container with the app directory mounted
podman run -d \
    --name CryoProtect \
    -v "$TEMP_APP_DIR:/app:Z" \
    -p 5000:5000 \
    --security-opt label=type:container_t \
    -w /app \
    continuumio/miniconda3:4.12.0 \
    tail -f /dev/null

# Set environment variables and run the app
echo "Running CryoProtect application..."
podman exec -it CryoProtect bash -c "
    cd /app && \
    export PATH=/opt/conda/envs/cryoprotect/bin:$PATH && \
    export FLASK_APP=app.py && \
    export FLASK_ENV=development && \
    python app.py
"

#!/bin/bash
# Create a Podman container for CryoProtect with conda and all required dependencies
# This script addresses SELinux permission issues and sets up the container properly

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================"
echo "Creating CryoProtect Conda Podman Container"
echo -e "======================================${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman podman-compose"
    exit 1
fi

# Check if container already exists
if podman container exists CryoProtect; then
    echo -e "${YELLOW}Container 'CryoProtect' already exists.${NC}"
    read -p "Do you want to remove it and create a new one? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing container..."
        podman rm -f CryoProtect
    else
        echo "Exiting without changes."
        exit 0
    fi
fi

# SELinux preparation - Create a temporary directory with appropriate context
TEMP_DIR="/tmp/cryoprotect-data"
PROJ_DIR="$(pwd)"

echo -e "${BLUE}Setting up temporary directory with appropriate SELinux context...${NC}"
mkdir -p $TEMP_DIR

# Try setting SELinux context if possible, otherwise continue with warning
if command -v chcon &> /dev/null; then
    chcon -Rt container_file_t $TEMP_DIR 2>/dev/null || \
    echo -e "${YELLOW}Warning: Could not set SELinux context. If you encounter permission issues, run: sudo chcon -Rt container_file_t $TEMP_DIR${NC}"
else
    echo -e "${YELLOW}Warning: chcon command not found. SELinux context not set.${NC}"
fi

# Create container - using a specific continuumio/miniconda3 version that's known to work well
echo -e "${BLUE}Creating CryoProtect container from miniconda image...${NC}"
podman run -d \
    --name CryoProtect \
    -v "$TEMP_DIR:/data:Z" \
    -p 5000:5000 \
    --security-opt label=type:container_t \
    continuumio/miniconda3:4.12.0 \
    tail -f /dev/null

echo -e "${GREEN}Container created!${NC}"

# Copy environment.yml to the temporary directory
echo -e "${BLUE}Copying environment configuration...${NC}"
cp "$PROJ_DIR/environment.yml" "$TEMP_DIR/"

# Install dependencies in the container
echo -e "${BLUE}Installing conda environment (this may take several minutes)...${NC}"
podman exec -it CryoProtect bash -c "
    echo 'Creating conda environment from environment.yml...' && \
    conda env create -f /data/environment.yml && \
    echo 'Cleaning conda cache to reduce container size...' && \
    conda clean -afy && \
    echo 'Conda environment created successfully.'
"

# Create a convenience script to run commands in the container
echo -e "${BLUE}Creating convenience script...${NC}"
cat > "$PROJ_DIR/run_in_cryoprotect.sh" << EOL
#!/bin/bash
# Run commands in the CryoProtect container with conda environment activated

# Check if container exists
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_conda_container.sh first."
    exit 1
fi

# Check if container is running
if ! podman container inspect CryoProtect --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting CryoProtect container..."
    podman start CryoProtect
fi

# Execute the command in the container with conda environment activated
if [ \$# -eq 0 ]; then
    # If no arguments, start an interactive shell with conda environment activated
    echo "Starting interactive shell in container with conda environment..."
    podman exec -it CryoProtect bash -c "
        # Try to initialize conda if not already done
        conda init bash &>/dev/null || true
        source ~/.bashrc &>/dev/null || true
        # Try to activate with conda or use the direct path to the environment's Python
        echo 'Activating conda environment...'
        conda activate cryoprotect &>/dev/null
        if [ \$? -ne 0 ]; then
            echo 'Using direct path to conda environment instead'
            export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH
        fi
        bash
    "
else
    # Otherwise, run the specified command with conda environment activated
    echo "Running command in container with conda environment: \$@"
    podman exec -it CryoProtect bash -c "
        # Try using direct path to the environment's Python first
        export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH
        \$*
    "
fi
EOL
chmod +x "$PROJ_DIR/run_in_cryoprotect.sh"

# Create a script to run the application
echo -e "${BLUE}Creating application run script...${NC}"
cat > "$PROJ_DIR/run_app_in_container.sh" << EOL
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
mkdir -p "\$TEMP_APP_DIR"
sudo chcon -Rt container_file_t "\$TEMP_APP_DIR"

# Copy essential application files to the temporary directory
echo "Copying application files to temporary directory..."
cp -r $(pwd)/*.py "\$TEMP_APP_DIR/"
cp -r $(pwd)/api "\$TEMP_APP_DIR/"
cp -r $(pwd)/database "\$TEMP_APP_DIR/"
mkdir -p "\$TEMP_APP_DIR/reports"

# Stop the container and recreate it with the app mount
echo "Recreating container with application mount..."
podman rm -f CryoProtect

# Recreate the container with the app directory mounted
podman run -d \\
    --name CryoProtect \\
    -v "\$TEMP_APP_DIR:/app:Z" \\
    -p 5000:5000 \\
    --security-opt label=type:container_t \\
    -w /app \\
    continuumio/miniconda3:4.12.0 \\
    tail -f /dev/null

# Set environment variables and run the app
echo "Running CryoProtect application..."
podman exec -it CryoProtect bash -c "
    cd /app && \\
    export PATH=/opt/conda/envs/cryoprotect/bin:\$PATH && \\
    export FLASK_APP=app.py && \\
    export FLASK_ENV=development && \\
    python app.py
"
EOL
chmod +x "$PROJ_DIR/run_app_in_container.sh"

# Create a script to mount the project files securely
echo -e "${BLUE}Creating mount helper script...${NC}"
cat > "$PROJ_DIR/mount_project_in_container.sh" << EOL
#!/bin/bash
# Safely mount the project directory into the CryoProtect container
# This script uses a temporary directory with proper SELinux context

set -e

TEMP_APP_DIR="/tmp/cryoprotect-app"
PROJ_DIR="\$(pwd)"

# Check if container exists
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_conda_container.sh first."
    exit 1
fi

# Prepare temporary directory with proper SELinux context
echo "Setting up temporary directory with proper SELinux context..."
mkdir -p "\$TEMP_APP_DIR"
sudo chcon -Rt container_file_t "\$TEMP_APP_DIR"

# Copy important files to temporary directory
echo "Copying project files to temporary directory..."
rsync -av --exclude="*/__pycache__" --exclude="*.pyc" \\
    "\$PROJ_DIR/"*.py \\
    "\$PROJ_DIR/api" \\
    "\$PROJ_DIR/database" \\
    "\$PROJ_DIR/migrations" \\
    "\$PROJ_DIR/chembl" \\
    "\$PROJ_DIR/pubchem" \\
    "\$TEMP_APP_DIR/"

mkdir -p "\$TEMP_APP_DIR/reports"
mkdir -p "\$TEMP_APP_DIR/data"

# Recreate the container with the app directory mounted
echo "Stopping existing container..."
podman stop CryoProtect
podman rm CryoProtect

echo "Creating new container with project mount..."
podman run -d \\
    --name CryoProtect \\
    -v "\$TEMP_APP_DIR:/app:Z" \\
    -p 5000:5000 \\
    --security-opt label=type:container_t \\
    -w /app \\
    continuumio/miniconda3:4.12.0 \\
    tail -f /dev/null

echo "Project directory mounted securely in container."
echo "Use './run_in_cryoprotect.sh' to run commands in the container."
EOL
chmod +x "$PROJ_DIR/mount_project_in_container.sh"

# Create a systemd service file for user-level service
echo -e "${BLUE}Creating systemd service file...${NC}"
mkdir -p ~/.config/systemd/user/
cat > ~/.config/systemd/user/cryoprotect-container.service << EOL
[Unit]
Description=CryoProtect Application Container Service
After=network.target
Documentation=https://github.com/yourusername/cryoprotect

[Service]
Type=simple
WorkingDirectory=$PROJ_DIR
ExecStartPre=$PROJ_DIR/run_in_cryoprotect.sh "echo 'Starting CryoProtect service...'"
ExecStart=$PROJ_DIR/run_app_in_container.sh
ExecStop=podman stop CryoProtect
Restart=on-failure
RestartSec=10s

[Install]
WantedBy=default.target
EOL

# Create a mock_rdkit script for testing without RDKit
echo -e "${BLUE}Creating RDKit mock script for testing...${NC}"
cat > "$PROJ_DIR/setup_mock_rdkit.sh" << EOL
#!/bin/bash
# Setup mock RDKit for testing when RDKit is not available

set -e

echo "Creating mock RDKit module for testing..."
python $PROJ_DIR/mock_rdkit.py

# Get the mock module path
MOCK_PATH=\$(python -c "import sys; print(f'/tmp/mock_modules:\{\":\".\".join(sys.path)}')")

echo "Mock RDKit created successfully."
echo
echo "To use the mock RDKit in the container, run:"
echo "./run_in_cryoprotect.sh \"export PYTHONPATH=$MOCK_PATH && python your_script.py\""
EOL
chmod +x "$PROJ_DIR/setup_mock_rdkit.sh"

# Create a guide document
echo -e "${BLUE}Creating usage guide...${NC}"
cat > "$PROJ_DIR/CONDA_CONTAINER_GUIDE.md" << EOL
# CryoProtect Conda Container Guide

This guide explains how to use the CryoProtect conda container setup with SELinux compatibility.

## Setup Scripts

1. **create_cryoprotect_conda_container.sh**
   - Creates the main container with conda environment
   - Installs all dependencies including RDKit
   - Sets up helper scripts
   - Handles SELinux permissions properly

2. **run_in_cryoprotect.sh**
   - Runs commands in the container with conda environment activated
   - Can be used interactively or with specific commands

3. **run_app_in_container.sh**
   - Runs the CryoProtect application in the container
   - Handles mounting the application code securely

4. **mount_project_in_container.sh**
   - Safely mounts project files into the container
   - Uses temporary directory with proper SELinux context

5. **setup_mock_rdkit.sh**
   - Creates a mock RDKit module for testing without RDKit

## Common Tasks

### Running Interactive Shell

```bash
./run_in_cryoprotect.sh
```

### Running Specific Commands

```bash
./run_in_cryoprotect.sh "python -c 'import rdkit; print(rdkit.__version__)'"
```

### Running the Application

```bash
./run_app_in_container.sh
```

### Mounting Project Files

```bash
./mount_project_in_container.sh
```

### Using Mock RDKit

```bash
./setup_mock_rdkit.sh
./run_in_cryoprotect.sh "export PYTHONPATH=/tmp/mock_modules:\$PYTHONPATH && python verify_rdkit.py"
```

### Setting Up as User Service

```bash
systemctl --user daemon-reload
systemctl --user enable cryoprotect-container.service
systemctl --user start cryoprotect-container.service
```

## Troubleshooting SELinux Issues

If you encounter SELinux permission issues:

1. Check SELinux denials:
   ```bash
   sudo ausearch -m AVC -ts recent
   ```

2. Run the container with permissive SELinux type:
   ```bash
   podman run --security-opt label=type:container_t,label=level:s0:c1023,label=disable ...
   ```

3. Apply SELinux configurations:
   ```bash
   sudo ./setup_podman_selinux.sh --apply
   ```

## Container Management

- Stop container: `podman stop CryoProtect`
- Start container: `podman start CryoProtect`
- Remove container: `podman rm CryoProtect`
- View container logs: `podman logs CryoProtect`
- Inspect container: `podman inspect CryoProtect`
EOL

echo -e "${GREEN}Setup completed successfully!${NC}"
echo 
echo -e "${BLUE}Usage:${NC}"
echo "1. Run commands in the container with conda environment:"
echo "   ./run_in_cryoprotect.sh [command]"
echo 
echo "2. Start an interactive shell in the container with conda environment:"
echo "   ./run_in_cryoprotect.sh"
echo 
echo "3. Run the application in the container:"
echo "   ./run_app_in_container.sh"
echo 
echo "4. Mount project files securely:"
echo "   ./mount_project_in_container.sh"
echo
echo "5. Install as a systemd service (user mode):"
echo "   systemctl --user daemon-reload"
echo "   systemctl --user enable cryoprotect-container.service"
echo "   systemctl --user start cryoprotect-container.service"
echo 
echo "6. Check service status:"
echo "   systemctl --user status cryoprotect-container.service"
echo
echo "See CONDA_CONTAINER_GUIDE.md for more information and troubleshooting."
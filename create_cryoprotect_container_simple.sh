#!/bin/bash
# Create a Podman container for CryoProtect with all required dependencies

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================"
echo "Creating CryoProtect Podman Container"
echo -e "======================================${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman podman-compose"
    exit 1
fi

# Check if the container already exists
if podman container exists CryoProtect; then
    echo -e "${YELLOW}Container 'CryoProtect' already exists.${NC}"
    
    # Ask if we should remove it
    read -p "Do you want to remove it and create a new one? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing container..."
        podman stop CryoProtect || true
        podman rm CryoProtect
    else
        echo "Keeping existing container."
        exit 0
    fi
fi

# Create a container from miniconda
echo -e "${BLUE}Creating CryoProtect container from Python image...${NC}"
podman run -d \
    --name CryoProtect \
    -v "$(pwd):/app" \
    -p 5000:5000 \
    -w /app \
    docker.io/python:3.10-slim \
    tail -f /dev/null

echo -e "${GREEN}Container created!${NC}"

# Install dependencies in the container
echo -e "${BLUE}Installing basic dependencies in the container...${NC}"
podman exec -it CryoProtect bash -c "
    echo 'Installing dependencies...' && \
    pip install flask flask-restful psycopg2-binary requests python-dotenv
"

# Create a convenience script to run commands in the container
echo -e "${BLUE}Creating convenience script...${NC}"
cat > run_in_cryoprotect.sh << EOL
#!/bin/bash
# Run commands in the CryoProtect container

# Check if container exists
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_container.sh first."
    exit 1
fi

# Check if container is running
if ! podman container inspect CryoProtect --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting CryoProtect container..."
    podman start CryoProtect
fi

# Execute the command in the container
if [ \$# -eq 0 ]; then
    # If no arguments, start an interactive shell
    echo "Starting interactive shell in container..."
    podman exec -it CryoProtect bash
else
    # Otherwise, run the specified command
    echo "Running command in container: \$@"
    podman exec -it CryoProtect bash -c "\$*"
fi
EOL
chmod +x run_in_cryoprotect.sh

# Create a script to run the application
echo -e "${BLUE}Creating application run script...${NC}"
cat > run_app_in_container.sh << EOL
#!/bin/bash
# Run the CryoProtect application in the container

# Set environment variables
export SUPABASE_DB_HOST=host.containers.internal
export SUPABASE_DB_PORT=5432
export SUPABASE_DB_NAME=postgres
export SUPABASE_DB_USER=postgres
export SUPABASE_DB_PASSWORD=postgres

# Run the app
./run_in_cryoprotect.sh "cd /app && python app.py"
EOL
chmod +x run_app_in_container.sh

# Create a systemd service file
echo -e "${BLUE}Creating systemd service file...${NC}"
cat > cryoprotect-container.service << EOL
[Unit]
Description=CryoProtect Application Container Service
After=network.target
Documentation=https://github.com/yourusername/cryoprotect

[Service]
Type=simple
User=$(whoami)
WorkingDirectory=$(pwd)
ExecStartPre=$(pwd)/run_in_cryoprotect.sh "echo 'Starting CryoProtect service...'"
ExecStart=$(pwd)/run_app_in_container.sh
ExecStop=podman stop CryoProtect
Restart=on-failure
RestartSec=10s

[Install]
WantedBy=multi-user.target
EOL

echo -e "${GREEN}Setup completed successfully!${NC}"
echo 
echo -e "${BLUE}Usage:${NC}"
echo "1. Run commands in the container:"
echo "   ./run_in_cryoprotect.sh [command]"
echo 
echo "2. Start an interactive shell in the container:"
echo "   ./run_in_cryoprotect.sh"
echo 
echo "3. Run the application in the container:"
echo "   ./run_app_in_container.sh"
echo 
echo "4. Install as a systemd service (user mode):"
echo "   mkdir -p ~/.config/systemd/user/"
echo "   cp cryoprotect-container.service ~/.config/systemd/user/"
echo "   systemctl --user daemon-reload"
echo "   systemctl --user enable cryoprotect-container.service"
echo "   systemctl --user start cryoprotect-container.service"
echo 
echo "5. Check service status:"
echo "   systemctl --user status cryoprotect-container.service"
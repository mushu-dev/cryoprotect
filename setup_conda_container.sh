#!/bin/bash
# Setup script for CryoProtect conda container

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================"
echo "CryoProtect Conda Container Setup"
echo -e "======================================${NC}"

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo -e "${RED}Error: Podman is not installed.${NC}"
    echo "Please install Podman first using:"
    echo "  sudo dnf install podman podman-compose"
    exit 1
fi

# Create a Containerfile
echo -e "${BLUE}Creating Containerfile...${NC}"
cat > Containerfile << EOL
FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy environment file
COPY environment.prod.yml /app/environment.yml

# Create conda environment
RUN conda env create -f environment.yml && \\
    conda clean -afy

# Add conda init to .bashrc
RUN conda init bash

# Set default shell to bash
SHELL ["/bin/bash", "-c"]

# Create a startup script
RUN echo '#!/bin/bash' > /app/start.sh && \\
    echo 'eval "$(conda shell.bash hook)"' >> /app/start.sh && \\
    echo 'conda activate cryoprotect' >> /app/start.sh && \\
    echo 'exec bash "\$@"' >> /app/start.sh && \\
    chmod +x /app/start.sh

ENTRYPOINT ["/app/start.sh"]
CMD ["--login"]
EOL

# Build the container
echo -e "${BLUE}Building the container...${NC}"
podman build -t cryoprotect-conda -f Containerfile .

# Create a convenience script to run the container
echo -e "${BLUE}Creating run script...${NC}"
cat > run_conda_container.sh << EOL
#!/bin/bash
# Run the CryoProtect conda container with the current directory mounted

# Check if container exists, if not create it
if ! podman container exists cryoprotect-conda; then
    echo "Creating cryoprotect-conda container..."
    podman run -d \\
        --name cryoprotect-conda \\
        -v "\$(pwd):/app" \\
        -p 5000:5000 \\
        --entrypoint sleep \\
        cryoprotect-conda infinity
    
    echo "Container created and running in background"
fi

# Check if container is running
if ! podman container inspect cryoprotect-conda --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting cryoprotect-conda container..."
    podman start cryoprotect-conda
fi

# Execute the command in the container
if [ \$# -eq 0 ]; then
    # If no arguments, start an interactive shell
    echo "Starting interactive shell in container..."
    podman exec -it cryoprotect-conda /app/start.sh
else
    # Otherwise, run the specified command
    echo "Running command in container: \$@"
    podman exec -it cryoprotect-conda /app/start.sh -c "\$*"
fi
EOL

chmod +x run_conda_container.sh

# Create a script to run the application
echo -e "${BLUE}Creating application run script...${NC}"
cat > run_app_in_container.sh << EOL
#!/bin/bash
# Run the CryoProtect application in the conda container

# Set environment variables
export SUPABASE_DB_HOST=host.containers.internal
export SUPABASE_DB_PORT=5432
export SUPABASE_DB_NAME=postgres
export SUPABASE_DB_USER=postgres
export SUPABASE_DB_PASSWORD=postgres

# Run the app
./run_conda_container.sh "cd /app && python app.py"
EOL

chmod +x run_app_in_container.sh

# Create a systemd service file
echo -e "${BLUE}Creating systemd service file...${NC}"
cat > cryoprotect-app-container.service << EOL
[Unit]
Description=CryoProtect Application Container
After=network.target
Documentation=https://github.com/yourusername/cryoprotect

[Service]
Type=forking
User=$(whoami)
WorkingDirectory=$(pwd)
ExecStart=$(pwd)/run_app_in_container.sh
Restart=on-failure
RestartSec=10s
TimeoutStartSec=180
TimeoutStopSec=60

[Install]
WantedBy=multi-user.target
EOL

echo -e "${GREEN}Setup completed successfully!${NC}"
echo 
echo -e "${BLUE}Usage:${NC}"
echo "1. Start an interactive shell in the container:"
echo "   ./run_conda_container.sh"
echo 
echo "2. Run the application in the container:"
echo "   ./run_app_in_container.sh"
echo 
echo "3. Install as a systemd service (user mode):"
echo "   mkdir -p ~/.config/systemd/user/"
echo "   cp cryoprotect-app-container.service ~/.config/systemd/user/"
echo "   systemctl --user daemon-reload"
echo "   systemctl --user enable cryoprotect-app-container.service"
echo "   systemctl --user start cryoprotect-app-container.service"
echo 
echo "4. Check service status:"
echo "   systemctl --user status cryoprotect-app-container.service"
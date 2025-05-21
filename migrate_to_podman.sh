#!/bin/bash
#
# migrate_to_podman.sh - Docker to Podman migration script for CryoProtect on Fedora
#
# This script assists with migrating from Docker to Podman on Fedora systems.
# It checks for Podman installation, provides installation instructions,
# migrates Docker configuration files, and provides basic usage guidance.

set -e

# Text formatting
BOLD="\033[1m"
GREEN="\033[32m"
YELLOW="\033[33m"
RED="\033[31m"
RESET="\033[0m"

# Function to log messages with timestamp
log() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] $1"
}

# Function to log success messages
log_success() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${GREEN}${BOLD}$1${RESET}"
}

# Function to log warning messages
log_warning() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${YELLOW}${BOLD}$1${RESET}"
}

# Function to log error messages
log_error() {
    echo -e "[$(date +"%Y-%m-%d %H:%M:%S")] ${RED}${BOLD}$1${RESET}"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" &> /dev/null
}

# Function to check Fedora version
check_fedora() {
    if [ -f /etc/fedora-release ]; then
        FEDORA_VERSION=$(cat /etc/fedora-release | grep -oP '(?<=release )[0-9]+')
        log "Detected Fedora $FEDORA_VERSION"
        return 0
    else
        log_warning "This script is designed for Fedora Linux. Your system may not be fully supported."
        return 1
    fi
}

# Function to check and install Podman
check_podman() {
    if command_exists podman; then
        PODMAN_VERSION=$(podman --version | awk '{print $3}')
        log_success "Podman is already installed (version $PODMAN_VERSION)"
        return 0
    else
        log_warning "Podman is not installed on this system."
        read -p "Would you like to install Podman now? (y/n): " install_podman
        if [[ "$install_podman" =~ ^[Yy]$ ]]; then
            log "Installing Podman..."
            sudo dnf install -y podman podman-compose container-selinux
            if command_exists podman; then
                PODMAN_VERSION=$(podman --version | awk '{print $3}')
                log_success "Podman installed successfully (version $PODMAN_VERSION)"
                return 0
            else
                log_error "Failed to install Podman. Please install it manually."
                return 1
            fi
        else
            log_warning "Podman installation skipped. You will need to install it manually."
            return 1
        fi
    fi
}

# Function to check Docker installation
check_docker() {
    if command_exists docker; then
        DOCKER_VERSION=$(docker --version | awk '{print $3}' | tr -d ',')
        log "Docker is installed (version $DOCKER_VERSION)"
        log_warning "It's recommended to remove Docker after migrating to Podman."
        return 0
    else
        log "Docker is not installed. Proceeding with Podman-only setup."
        return 1
    fi
}

# Function to create Podman-compatible compose files
create_podman_compose_files() {
    log "Creating Podman-compatible compose files..."
    
    # Check if docker-compose.yml exists
    if [ -f docker-compose.yml ]; then
        log "Converting docker-compose.yml to podman-compose.yml..."
        
        # Create a backup
        cp docker-compose.yml docker-compose.yml.bak
        log "Created backup of docker-compose.yml to docker-compose.yml.bak"
        
        # Create podman-compose.yml with SELinux context labels
        cat docker-compose.yml | sed 's/\(volume:.*\)/\1:Z/g' > podman-compose.yml
        
        # Check if development compose file exists
        if [ -f docker-compose.dev.yml ]; then
            cp docker-compose.dev.yml docker-compose.dev.yml.bak
            log "Created backup of docker-compose.dev.yml to docker-compose.dev.yml.bak"
            
            cat docker-compose.dev.yml | sed 's/\(volume:.*\)/\1:Z/g' > podman-compose.dev.yml
            log_success "Created podman-compose.dev.yml"
        fi
        
        log_success "Created podman-compose.yml"
    else
        log_warning "docker-compose.yml not found. Cannot create podman-compose.yml."
        log_warning "Please run this script from the CryoProtect project root directory."
        return 1
    fi
    
    # Create minimal compose file for quick startup
    cat > podman-compose.minimal.yml << 'EOL'
version: '3.8'

services:
  cryoprotect:
    build: .
    image: cryoprotect:latest
    container_name: cryoprotect
    ports:
      - "5000:5000"
    volumes:
      - ./:/app:Z
    environment:
      - FLASK_APP=app.py
      - FLASK_ENV=development
      - PYTHONUNBUFFERED=1
    restart: unless-stopped
    command: python app.py
EOL
    log_success "Created podman-compose.minimal.yml for simplified deployment"
    
    return 0
}

# Function to create a quickstart script
create_quickstart_script() {
    log "Creating Podman quickstart script..."
    
    cat > quickstart_podman.sh << 'EOL'
#!/bin/bash
#
# quickstart_podman.sh - Quick start script for CryoProtect with Podman
#
# This script provides a simplified way to get started with CryoProtect using Podman.

set -e

# Check if Podman is installed
if ! command -v podman &> /dev/null; then
    echo "Error: Podman is not installed. Please install it first."
    echo "Run: sudo dnf install -y podman podman-compose"
    exit 1
fi

# Make sure we're in the project root
if [ ! -f "app.py" ] || [ ! -f "podman-compose.minimal.yml" ]; then
    echo "Error: This script must be run from the CryoProtect project root directory."
    exit 1
fi

# Create necessary directories if they don't exist
mkdir -p logs
mkdir -p backup/data
mkdir -p redis-data

# Set SELinux context if available
if command -v chcon &> /dev/null; then
    echo "Setting SELinux context for volume mount points..."
    chcon -Rt container_file_t ./logs
    chcon -Rt container_file_t ./backup/data
    chcon -Rt container_file_t ./redis-data
fi

# Build and start the application
echo "Building and starting CryoProtect..."
podman-compose -f podman-compose.minimal.yml up -d

echo ""
echo "CryoProtect is now running!"
echo "Access the application at: http://localhost:5000"
echo ""
echo "To view logs: podman logs -f cryoprotect"
echo "To stop: podman-compose -f podman-compose.minimal.yml down"
EOL
    
    chmod +x quickstart_podman.sh
    log_success "Created quickstart_podman.sh"
    
    return 0
}

# Function to provide Docker to Podman command equivalents
show_command_equivalents() {
    log "Docker to Podman Command Equivalents:"
    echo -e "\n${BOLD}Common Docker and Podman Command Equivalents${RESET}\n"
    echo -e "| ${BOLD}Task${RESET}                     | ${BOLD}Docker${RESET}                          | ${BOLD}Podman${RESET}                          |"
    echo -e "|--------------------------|----------------------------------|----------------------------------|"
    echo -e "| Pull image               | docker pull image:tag            | podman pull image:tag            |"
    echo -e "| List images              | docker images                    | podman images                    |"
    echo -e "| List containers          | docker ps                        | podman ps                        |"
    echo -e "| Build image              | docker build -t name:tag .       | podman build -t name:tag .       |"
    echo -e "| Run container            | docker run image                 | podman run image                 |"
    echo -e "| Stop container           | docker stop container            | podman stop container            |"
    echo -e "| Remove container         | docker rm container              | podman rm container              |"
    echo -e "| Execute in container     | docker exec -it container cmd    | podman exec -it container cmd    |"
    echo -e "| View logs                | docker logs container            | podman logs container            |"
    echo -e "| Compose up               | docker-compose up -d             | podman-compose up -d             |"
    echo -e "| Compose down             | docker-compose down              | podman-compose down              |"
    echo -e "| Login to registry        | docker login registry            | podman login registry            |"
    echo -e "| Prune unused resources   | docker system prune              | podman system prune              |"
    echo -e "\n${BOLD}Major Differences${RESET}\n"
    echo -e "1. Podman runs containers rootless by default (no daemon)"
    echo -e "2. Podman uses SELinux labels (:Z or :z) for volume mounts on Fedora"
    echo -e "3. Podman manages containers per-user rather than system-wide"
    echo -e "4. Podman supports Kubernetes YAML format natively"
}

# Function to test a basic Podman container
test_podman() {
    log "Testing basic Podman functionality..."
    
    if podman run --rm hello-world | grep -q "Hello from Docker"; then
        log_success "Podman test successful. Basic functionality is working."
        return 0
    else
        log_error "Podman test failed. Please check your installation."
        return 1
    fi
}

# Function to create full Podman deployment guide
create_deployment_guide() {
    log "Creating comprehensive Podman deployment guide..."
    
    cat > PODMAN_DEPLOYMENT_GUIDE.md << 'EOL'
# Podman Deployment Guide for CryoProtect on Fedora

This guide provides detailed instructions for deploying CryoProtect using Podman on Fedora Linux systems. Podman is a daemonless container engine that is compatible with Docker but provides additional security features, especially on SELinux-enabled systems like Fedora.

## Prerequisites

- Fedora Linux (version 35 or higher recommended)
- Podman and podman-compose installed
- Git (to clone the repository if needed)
- Basic familiarity with container concepts

## Installation

### Install Podman and Related Tools

```bash
# Install Podman and podman-compose
sudo dnf install -y podman podman-compose container-selinux

# Verify installation
podman --version
podman-compose --version
```

## Quick Start

The simplest way to get started is to use the included quickstart script:

```bash
./quickstart_podman.sh
```

This script will:
1. Check for Podman installation
2. Create necessary directories
3. Set appropriate SELinux contexts
4. Build and start CryoProtect using podman-compose.minimal.yml

## Deployment Options

### Option 1: Minimal Deployment

For a simple deployment with just the CryoProtect application:

```bash
# Create necessary directories
mkdir -p logs backup/data redis-data

# Set SELinux context for volume mounts
chcon -Rt container_file_t ./logs
chcon -Rt container_file_t ./backup/data
chcon -Rt container_file_t ./redis-data

# Start the application
podman-compose -f podman-compose.minimal.yml up -d
```

### Option 2: Full Deployment

For a complete deployment with all services:

```bash
# Start the full stack
podman-compose -f podman-compose.yml up -d
```

### Option 3: Development Environment

For development purposes with live reloading:

```bash
# Start the development environment
podman-compose -f podman-compose.dev.yml up
```

## SELinux Considerations

Podman on Fedora enforces SELinux policies by default. Volume mounts require proper context labeling:

- Use `:Z` suffix for dedicated volumes (used by a single container)
- Use `:z` suffix for shared volumes (used by multiple containers)

Example in compose file:
```yaml
volumes:
  - ./data:/app/data:Z
```

## Common Operations

### View Logs

```bash
# View logs from the CryoProtect container
podman logs -f cryoprotect
```

### Access Shell in Container

```bash
# Execute a shell in the running container
podman exec -it cryoprotect /bin/bash
```

### Stop and Remove Containers

```bash
# Stop all containers defined in compose file
podman-compose -f podman-compose.yml down

# Remove stopped containers
podman container prune
```

### Update the Application

```bash
# Pull latest code changes
git pull

# Rebuild and restart containers
podman-compose -f podman-compose.yml down
podman-compose -f podman-compose.yml up -d --build
```

## Troubleshooting

### Permission Denied on Volume Mounts

If you encounter permission issues with volume mounts:

```bash
# Check SELinux contexts
ls -lZ ./logs ./backup/data

# Reset context if needed
chcon -Rt container_file_t ./logs
chcon -Rt container_file_t ./backup/data
```

### Network Connectivity Issues

If containers can't connect to external services:

```bash
# Check if firewalld is blocking connections
sudo firewall-cmd --list-all

# Allow necessary ports
sudo firewall-cmd --permanent --add-port=5000/tcp
sudo firewall-cmd --reload
```

### Container Fails to Start

If the container exits immediately:

```bash
# Check container logs
podman logs cryoprotect

# Inspect container details
podman inspect cryoprotect

# Start container interactively for debugging
podman run -it --rm cryoprotect /bin/bash
```

## Podman vs Docker Differences

- **Rootless Containers**: Podman runs containers as your user by default
- **No Daemon**: Podman doesn't require a running daemon
- **SELinux Integration**: Podman has better integration with SELinux
- **Compatible CLI**: Most Docker commands work with Podman
- **Compose Differences**: podman-compose has some subtle differences from docker-compose

## Resources

- [Podman Documentation](https://podman.io/docs/)
- [Podman-Compose GitHub](https://github.com/containers/podman-compose)
- [Fedora Documentation on Containers](https://docs.fedoraproject.org/en-US/fedora-silverblue/docker-podman/)
EOL
    
    log_success "Created PODMAN_DEPLOYMENT_GUIDE.md"
    
    return 0
}

# Main script execution
log "Starting Docker to Podman migration script for CryoProtect..."
echo ""

# Check if running on Fedora
check_fedora

# Check Podman and Docker installation
check_podman
check_docker

# Create Podman-compatible compose files
create_podman_compose_files

# Create quickstart script
create_quickstart_script

# Test Podman functionality
test_podman

# Create comprehensive deployment guide
create_deployment_guide

# Show command equivalents
show_command_equivalents

# Summary and next steps
echo ""
log_success "Migration preparation completed!"
echo ""
echo -e "${BOLD}Next Steps:${RESET}"
echo "1. Review the generated files:"
echo "   - podman-compose.yml (main Podman compose file)"
echo "   - podman-compose.minimal.yml (simplified compose file)"
echo "   - quickstart_podman.sh (easy startup script)"
echo "   - PODMAN_DEPLOYMENT_GUIDE.md (comprehensive guide)"
echo ""
echo "2. Start CryoProtect with Podman using one of these methods:"
echo "   - Run the quickstart script: ./quickstart_podman.sh"
echo "   - Manual start: podman-compose -f podman-compose.yml up -d"
echo "   - Development mode: podman-compose -f podman-compose.dev.yml up"
echo ""
echo "3. Access the application at: http://localhost:5000"
echo ""
echo -e "${BOLD}For more information, see PODMAN_DEPLOYMENT_GUIDE.md${RESET}"
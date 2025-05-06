#!/bin/bash
# CryoProtect - Docker Setup Script
# This script sets up Docker for the CryoProtect project

set -e  # Exit on error

echo "================================================"
echo "CryoProtect Docker Setup"
echo "================================================"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Detect Linux distribution
if [ -f /etc/os-release ]; then
    . /etc/os-release
    DISTRO=$ID
    echo "Detected Linux distribution: $DISTRO"
else
    echo "Unable to determine Linux distribution. Assuming Debian/Ubuntu compatible."
    DISTRO="ubuntu"
fi

# Install Docker if not already installed
echo -e "\n[1/4] Installing Docker..."
if ! command_exists docker; then
    echo "Docker not found. Installing..."
    
    case $DISTRO in
        ubuntu|debian|pop|mint|kali|elementary)
            # Update package index
            sudo apt update
            
            # Install prerequisites
            sudo apt install -y apt-transport-https ca-certificates curl software-properties-common
            
            # Add Docker's official GPG key
            curl -fsSL https://download.docker.com/linux/$DISTRO/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
            
            # Add the Docker repository
            echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/$DISTRO $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
            
            # Install Docker
            sudo apt update
            sudo apt install -y docker-ce docker-ce-cli containerd.io
            ;;
            
        fedora|rhel|centos|rocky)
            # Install Docker
            sudo dnf -y install dnf-plugins-core
            sudo dnf config-manager --add-repo https://download.docker.com/linux/$DISTRO/docker-ce.repo
            sudo dnf install -y docker-ce docker-ce-cli containerd.io
            ;;
            
        arch|manjaro|endeavouros)
            # Install Docker
            sudo pacman -Sy --noconfirm docker
            ;;
            
        *)
            echo "Unsupported distribution for automatic Docker installation."
            echo "Please install Docker manually following the official documentation:"
            echo "https://docs.docker.com/engine/install/"
            exit 1
            ;;
    esac
    
    # Start and enable Docker service
    sudo systemctl start docker
    sudo systemctl enable docker
    
    echo "Docker installed successfully"
else
    echo "Docker is already installed"
fi

# Add current user to the docker group to run docker without sudo
echo -e "\n[2/4] Adding user to docker group..."
if ! groups | grep -q docker; then
    sudo usermod -aG docker $(whoami)
    echo "Added $(whoami) to the docker group"
    echo "NOTE: You may need to log out and log back in for this to take effect"
else
    echo "User is already in the docker group"
fi

# Install Docker Compose if not already installed
echo -e "\n[3/4] Installing Docker Compose..."
if ! command_exists docker-compose; then
    echo "Docker Compose not found. Installing..."
    
    case $DISTRO in
        ubuntu|debian|pop|mint|kali|elementary|fedora|rhel|centos|rocky)
            # Install Docker Compose
            DOCKER_COMPOSE_VERSION=$(curl -s https://api.github.com/repos/docker/compose/releases/latest | grep 'tag_name' | cut -d '"' -f 4)
            sudo curl -L "https://github.com/docker/compose/releases/download/${DOCKER_COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
            sudo chmod +x /usr/local/bin/docker-compose
            ;;
            
        arch|manjaro|endeavouros)
            # Install Docker Compose
            sudo pacman -Sy --noconfirm docker-compose
            ;;
    esac
    
    echo "Docker Compose installed successfully"
else
    echo "Docker Compose is already installed"
fi

# Set up environment for Docker
echo -e "\n[4/4] Setting up Docker environment..."

# Create .env file if it doesn't exist
if [ ! -f .env ] && [ -f .env.template ]; then
    cp .env.template .env
    echo "Created .env file from template"
    echo "⚠️ Please edit .env file to add your Supabase credentials"
elif [ ! -f .env ] && [ ! -f .env.template ]; then
    echo "⚠️ No .env.template found. Please create a .env file manually."
else
    echo ".env file already exists"
fi

echo "✅ Docker environment setup complete"

echo -e "\n================================================"
echo "Docker Setup Complete!"
echo "================================================"

echo -e "\nTo run CryoProtect with Docker:"
echo "1. Edit the .env file with your Supabase credentials"
echo "2. Build and start the containers:"
echo "   docker-compose up -d"
echo "3. Access the application at http://localhost:5000"
echo -e "\nTo run database migrations in Docker:"
echo "   docker-compose exec app node migrations/apply_migration.js"
echo -e "\nTo view logs:"
echo "   docker-compose logs -f"
echo -e "\nTo stop the containers:"
echo "   docker-compose down"
echo -e "\nIf you made changes to your environment, rebuild the containers:"
echo "   docker-compose up -d --build"
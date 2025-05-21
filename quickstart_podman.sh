#!/bin/bash
set -e

# CryoProtect Podman Quickstart for Fedora
# A simplified setup script focused on getting the core app running quickly

# Function to log messages
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Check for required dependencies
if ! command -v podman &> /dev/null; then
  log "Error: podman is not installed. Please install podman first:"
  log "  sudo dnf install -y podman podman-docker"
  exit 1
fi

# Create necessary directories
mkdir -p ./logs
mkdir -p ./cache

# Set SELinux context if enabled
if command -v getenforce &> /dev/null && [ "$(getenforce)" != "Disabled" ]; then
  log "Setting SELinux context for volume mount points..."
  chcon -Rt container_file_t ./logs || true
  chcon -Rt container_file_t ./cache || true
fi

# Create minimal .env file if it doesn't exist
if [ ! -f ./.env ]; then
  log "Creating .env file with placeholder values (you should customize these)..."
  cat > ./.env << EOF
# Supabase credentials - replace with your own
SUPABASE_URL=https://your-project.supabase.co
SUPABASE_KEY=your-supabase-anon-key
SECRET_KEY=dev-secret-key-please-change-in-production
EOF
  log "Created default .env file. Please edit it with your configuration."
  log "For now, we'll start with placeholder values..."
fi

# Source environment variables
set -a
source ./.env
set +a

# Check if podman-compose is installed
if ! command -v podman-compose &> /dev/null; then
  log "podman-compose not found, checking for docker-compose as an alternative..."
  
  if command -v docker-compose &> /dev/null; then
    log "Using docker-compose (it works with podman on Fedora)..."
    COMPOSE_CMD="docker-compose"
  else
    log "Installing podman-compose..."
    pip3 install --user podman-compose
    
    if ! command -v podman-compose &> /dev/null; then
      log "Error: Failed to install podman-compose. Please install it manually:"
      log "  pip3 install --user podman-compose"
      exit 1
    fi
    
    COMPOSE_CMD="podman-compose"
  fi
else
  COMPOSE_CMD="podman-compose"
fi

# Start the minimal application
log "Starting CryoProtect in development mode..."
log "Using the minimal configuration for simplicity..."
$COMPOSE_CMD -f podman-compose.minimal.yml up

# This script will terminate when the user presses Ctrl+C to stop the containers
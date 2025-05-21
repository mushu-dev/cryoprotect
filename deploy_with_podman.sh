#!/bin/bash
set -e

# CryoProtect Podman Deployment Script for Fedora
# This script sets up and launches CryoProtect using podman-compose

# Function to log messages
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Ensure required directories exist
setup_directories() {
  log "Setting up directories..."
  
  # Create logs directory with proper SELinux context
  mkdir -p ./logs
  mkdir -p ./backup/data
  mkdir -p ./redis-data
  
  # Set SELinux context for volume mount points
  # This is critical for Fedora with SELinux enforcing
  if command -v chcon &> /dev/null; then
    log "Setting SELinux context for volume mount points..."
    chcon -Rt container_file_t ./logs
    chcon -Rt container_file_t ./backup/data
    chcon -Rt container_file_t ./redis-data
  else
    log "Warning: chcon command not found. If SELinux is enabled, volume mounts may fail."
  fi
}

# Setup environment variables
setup_env() {
  log "Setting up environment variables..."
  
  # Create .env file if it doesn't exist
  if [ ! -f ./.env ]; then
    log "Creating .env file with default values (you should customize these)..."
    cat > ./.env << EOF
# CryoProtect environment variables
FLASK_ENV=development
HOST_PORT=5000
DEV_HOST_PORT=5000
REDIS_PORT=6379

# Supabase credentials - replace with your own
SUPABASE_URL=https://your-project.supabase.co
SUPABASE_KEY=your-supabase-anon-key
SECRET_KEY=your-secret-key

# Volume paths (absolute paths recommended for production)
LOGS_VOLUME_PATH=./logs
BACKUPS_VOLUME_PATH=./backup/data
REDIS_DATA_PATH=./redis-data

# Container registry (default to local)
PODMAN_REGISTRY=localhost
PODMAN_NAMESPACE=cryoprotect
IMAGE_TAG=latest
DEV_IMAGE_TAG=dev
EOF
    log "Created default .env file. Please edit it with your configuration."
    exit 1
  fi
  
  # Source environment variables
  set -a
  source ./.env
  set +a
}

# Check podman and podman-compose
check_dependencies() {
  log "Checking dependencies..."
  
  if ! command -v podman &> /dev/null; then
    log "Error: podman is not installed. Please install podman first."
    log "  Run: sudo dnf install -y podman podman-docker"
    exit 1
  fi
  
  if ! command -v podman-compose &> /dev/null; then
    log "Warning: podman-compose is not installed. Will try to use python-based podman-compose."
    if ! command -v pip3 &> /dev/null; then
      log "Error: pip3 is not installed. Please install python3-pip first."
      log "  Run: sudo dnf install -y python3-pip"
      exit 1
    fi
    
    log "Installing podman-compose using pip..."
    pip3 install --user podman-compose
    
    # Check if installation was successful
    if ! command -v podman-compose &> /dev/null; then
      log "Error: Failed to install podman-compose. Please install it manually."
      log "  Run: sudo pip3 install podman-compose"
      exit 1
    fi
  fi
  
  # Check if we can connect to the Podman socket
  if ! podman info &> /dev/null; then
    log "Error: Cannot connect to Podman. Make sure podman is running and you have permission to use it."
    log "  You may need to enable and start the podman.socket service:"
    log "  Run: systemctl --user enable --now podman.socket"
    exit 1
  }
}

# Start services using podman-compose
start_services() {
  local mode=$1
  
  case "$mode" in
    dev)
      log "Starting CryoProtect in development mode..."
      podman-compose -f podman-compose.yml up cryoprotect-dev redis
      ;;
    prod)
      log "Starting CryoProtect in production mode..."
      podman-compose -f podman-compose.yml up -d cryoprotect redis
      ;;
    *)
      log "Error: Invalid mode. Use 'dev' or 'prod'."
      exit 1
      ;;
  esac
}

# Stop services
stop_services() {
  log "Stopping CryoProtect services..."
  podman-compose -f podman-compose.yml down
}

# Create systemd service for autostart (user level)
create_systemd_service() {
  local service_file="$HOME/.config/systemd/user/cryoprotect.service"
  
  log "Creating systemd user service for automatic startup..."
  
  # Create directory if it doesn't exist
  mkdir -p "$HOME/.config/systemd/user"
  
  # Create service file
  cat > "$service_file" << EOF
[Unit]
Description=CryoProtect Application
After=network.target

[Service]
WorkingDirectory=$(pwd)
Environment=PODMAN_USERNS=keep-id
ExecStart=/usr/bin/podman-compose -f podman-compose.yml up
ExecStop=/usr/bin/podman-compose -f podman-compose.yml down
Restart=always
RestartSec=10s
Type=simple

[Install]
WantedBy=default.target
EOF
  
  log "Reloading systemd user daemon..."
  systemctl --user daemon-reload
  
  log "To enable autostart at login, run: systemctl --user enable cryoprotect.service"
  log "To start the service now, run: systemctl --user start cryoprotect.service"
}

# Show help
show_help() {
  echo "CryoProtect Podman Deployment Script for Fedora"
  echo ""
  echo "Usage: $0 [command]"
  echo ""
  echo "Commands:"
  echo "  setup       Setup directories, SELinux contexts, and check dependencies"
  echo "  dev         Start CryoProtect in development mode"
  echo "  prod        Start CryoProtect in production mode"
  echo "  stop        Stop all CryoProtect services"
  echo "  systemd     Create systemd user service for auto-starting CryoProtect"
  echo "  help        Show this help message"
  echo ""
  echo "Examples:"
  echo "  $0 setup    # Initial setup, run this first"
  echo "  $0 dev      # Start in development mode"
  echo "  $0 prod     # Start in production mode"
  echo ""
}

# Main logic
main() {
  local command="${1:-help}"
  
  case "$command" in
    setup)
      check_dependencies
      setup_directories
      setup_env
      log "Setup completed successfully."
      ;;
    dev)
      setup_env
      start_services "dev"
      ;;
    prod)
      setup_env
      start_services "prod"
      ;;
    stop)
      stop_services
      ;;
    systemd)
      create_systemd_service
      ;;
    help|*)
      show_help
      ;;
  esac
}

# Execute main function with command line arguments
main "$@"
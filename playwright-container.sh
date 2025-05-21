#!/bin/bash

# playwright-container.sh
# A robust solution for running Playwright in a container on any Linux distribution
# This script creates and manages a container for Playwright operations

# Configuration
CONTAINER_NAME="cryoprotect-playwright"
IMAGE_NAME="mcr.microsoft.com/playwright:v1.40.0-jammy"
HOST_PORT=8080
CONTAINER_PORT=8080
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="$SCRIPT_DIR/playwright-shared"
CONTAINER_SHARED_DIR="/app/shared"
CONTAINER_LOG_FILE="$SHARED_DIR/container.log"

# Create shared directory if it doesn't exist
mkdir -p "$SHARED_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Log message with timestamp
log() {
  local level=$1
  local message=$2
  local color=$NC
  
  case $level in
    "INFO") color=$GREEN ;;
    "WARN") color=$YELLOW ;;
    "ERROR") color=$RED ;;
    "DEBUG") color=$BLUE ;;
  esac
  
  echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message${NC}"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" >> "$SCRIPT_DIR/playwright-container.log"
}

# Check if podman is installed
check_podman() {
  if ! command -v podman &> /dev/null; then
    log "ERROR" "Podman is not installed. Please install podman first."
    exit 1
  fi
  log "INFO" "Podman is installed"
}

# Pull the container image if it doesn't exist
pull_image() {
  if ! podman image exists "$IMAGE_NAME"; then
    log "INFO" "Pulling image $IMAGE_NAME..."
    podman pull "$IMAGE_NAME"
    if [ $? -ne 0 ]; then
      log "ERROR" "Failed to pull image $IMAGE_NAME"
      exit 1
    fi
    log "INFO" "Image pulled successfully"
  else
    log "INFO" "Image $IMAGE_NAME already exists"
  fi
}

# Check if container exists
check_container() {
  if podman container exists "$CONTAINER_NAME"; then
    log "INFO" "Container $CONTAINER_NAME exists"
    return 0
  else
    log "INFO" "Container $CONTAINER_NAME does not exist"
    return 1
  fi
}

# Check if container is running
check_running() {
  if podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' 2>/dev/null | grep -q "true"; then
    log "INFO" "Container $CONTAINER_NAME is running"
    return 0
  else
    log "INFO" "Container $CONTAINER_NAME is not running"
    return 1
  fi
}

# Create and start container
start_container() {
  if check_container; then
    if ! check_running; then
      log "INFO" "Starting existing container $CONTAINER_NAME..."
      podman start "$CONTAINER_NAME"
      if [ $? -ne 0 ]; then
        log "ERROR" "Failed to start container $CONTAINER_NAME"
        exit 1
      fi
      log "INFO" "Container started successfully"
    fi
  else
    log "INFO" "Creating and starting container $CONTAINER_NAME..."
    podman run -d --name "$CONTAINER_NAME" \
      -p "$HOST_PORT:$CONTAINER_PORT" \
      --security-opt label=disable \
      -v "$SHARED_DIR:$CONTAINER_SHARED_DIR:Z" \
      "$IMAGE_NAME" \
      sleep infinity
    
    if [ $? -ne 0 ]; then
      log "ERROR" "Failed to create and start container $CONTAINER_NAME"
      exit 1
    fi
    log "INFO" "Container created and started successfully"
  fi
}

# Install necessary tools in container
setup_container() {
  log "INFO" "Setting up container environment..."
  
  # Update and install basic tools
  podman exec "$CONTAINER_NAME" bash -c "apt-get update && apt-get install -y jq curl socat"
  if [ $? -ne 0 ]; then
    log "WARN" "Could not install some tools in container. Some features may be limited."
  else
    log "INFO" "Basic tools installed successfully"
  fi
  
  # Install playwright if needed
  podman exec "$CONTAINER_NAME" bash -c "cd && npm list playwright || npm install -g playwright"
  if [ $? -ne 0 ]; then
    log "ERROR" "Failed to install Playwright in container"
    return 1
  fi
  
  # Verify Playwright installation
  podman exec "$CONTAINER_NAME" bash -c "node -e \"const { chromium } = require('playwright'); console.log('Playwright installed successfully');\""
  if [ $? -ne 0 ]; then
    log "ERROR" "Playwright installation verification failed"
    return 1
  fi
  
  log "INFO" "Container environment set up successfully"
  return 0
}

# Stop container
stop_container() {
  if check_container && check_running; then
    log "INFO" "Stopping container $CONTAINER_NAME..."
    podman stop "$CONTAINER_NAME"
    if [ $? -ne 0 ]; then
      log "ERROR" "Failed to stop container $CONTAINER_NAME"
      exit 1
    fi
    log "INFO" "Container stopped successfully"
  else
    log "INFO" "Container is not running"
  fi
}

# Remove container
remove_container() {
  if check_container; then
    stop_container
    log "INFO" "Removing container $CONTAINER_NAME..."
    podman rm "$CONTAINER_NAME"
    if [ $? -ne 0 ]; then
      log "ERROR" "Failed to remove container $CONTAINER_NAME"
      exit 1
    fi
    log "INFO" "Container removed successfully"
  else
    log "INFO" "Container does not exist"
  fi
}

# Run a Playwright script in the container
run_playwright_script() {
  local script_path=$1
  local script_name=$(basename "$script_path")
  local container_script_path="$CONTAINER_SHARED_DIR/$script_name"
  
  # Copy the script to the shared directory
  cp "$script_path" "$SHARED_DIR/"
  
  log "INFO" "Running Playwright script $script_name in container..."
  podman exec "$CONTAINER_NAME" bash -c "cd $CONTAINER_SHARED_DIR && node $script_name"
  local exit_code=$?
  
  if [ $exit_code -ne 0 ]; then
    log "ERROR" "Script execution failed with exit code $exit_code"
    return $exit_code
  else
    log "INFO" "Script executed successfully"
    return 0
  fi
}

# Run a Playwright command directly
run_playwright_command() {
  local command=$1
  
  log "INFO" "Running Playwright command in container: $command"
  podman exec "$CONTAINER_NAME" bash -c "cd $CONTAINER_SHARED_DIR && $command"
  local exit_code=$?
  
  if [ $exit_code -ne 0 ]; then
    log "ERROR" "Command execution failed with exit code $exit_code"
    return $exit_code
  else
    log "INFO" "Command executed successfully"
    return 0
  fi
}

# Parse command line arguments
if [ $# -eq 0 ]; then
  echo "Playwright Container Utility"
  echo "Usage: $0 <command> [arguments]"
  echo ""
  echo "Commands:"
  echo "  start             Start the Playwright container"
  echo "  stop              Stop the Playwright container"
  echo "  restart           Restart the Playwright container"
  echo "  status            Show container status"
  echo "  remove            Remove the container"
  echo "  run-script <path> Run a Node.js script with Playwright in the container"
  echo "  exec <command>    Execute a command in the container"
  echo "  logs              Show container logs"
  echo ""
  exit 0
fi

command=$1
shift

# Execute command
case $command in
  "start")
    check_podman
    pull_image
    start_container
    setup_container
    log "INFO" "Playwright container is ready"
    ;;
  "stop")
    stop_container
    ;;
  "restart")
    stop_container
    start_container
    setup_container
    log "INFO" "Playwright container restarted"
    ;;
  "status")
    if check_container; then
      if check_running; then
        log "INFO" "Container is running"
        podman exec "$CONTAINER_NAME" bash -c "node -e \"const { chromium } = require('playwright'); (async () => { try { const browser = await chromium.launch(); await browser.close(); console.log('Playwright is working properly'); } catch (e) { console.error('Playwright error:', e); process.exit(1); } })();\""
        if [ $? -eq 0 ]; then
          log "INFO" "Playwright is functioning correctly in container"
        else
          log "ERROR" "Playwright is not functioning correctly in container"
        fi
      else
        log "WARN" "Container exists but is not running"
      fi
    else
      log "WARN" "Container does not exist"
    fi
    ;;
  "remove")
    remove_container
    ;;
  "run-script")
    if [ $# -eq 0 ]; then
      log "ERROR" "No script path provided"
      echo "Usage: $0 run-script <path>"
      exit 1
    fi
    
    script_path=$1
    if [ ! -f "$script_path" ]; then
      log "ERROR" "Script file $script_path does not exist"
      exit 1
    fi
    
    check_podman
    pull_image
    start_container
    setup_container
    run_playwright_script "$script_path"
    ;;
  "exec")
    if [ $# -eq 0 ]; then
      log "ERROR" "No command provided"
      echo "Usage: $0 exec <command>"
      exit 1
    fi
    
    command="$*"
    check_podman
    pull_image
    start_container
    run_playwright_command "$command"
    ;;
  "logs")
    if check_container; then
      podman logs "$CONTAINER_NAME"
    else
      log "WARN" "Container does not exist"
    fi
    ;;
  *)
    log "ERROR" "Unknown command: $command"
    echo "Usage: $0 <command> [arguments]"
    exit 1
    ;;
esac

exit 0
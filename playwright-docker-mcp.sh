#!/bin/bash

# playwright-docker-mcp.sh
# A direct container mapping for MCP Playwright commands
# This script is designed to directly forward MCP commands to the Playwright Docker container

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_NAME="cryoprotect-playwright-mcp"
IMAGE_NAME="mcr.microsoft.com/playwright:v1.40.0-jammy"
SHARED_DIR="$SCRIPT_DIR/playwright-mcp-shared"
LOG_FILE="$SCRIPT_DIR/playwright-docker-mcp.log"
CONTAINER_SHARED_DIR="/app/shared"

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
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" >> "$LOG_FILE"
}

# Ensure container exists and is running
ensure_container() {
  # Check if container exists
  if podman container exists "$CONTAINER_NAME" 2>/dev/null; then
    log "INFO" "Container exists"
    
    # Check if container is running
    if podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' 2>/dev/null | grep -q "true"; then
      log "INFO" "Container is running"
    else
      log "INFO" "Starting container..."
      podman start "$CONTAINER_NAME"
    fi
  else
    log "INFO" "Creating container..."
    podman run -d --name "$CONTAINER_NAME" \
      --security-opt label=disable \
      -v "$SHARED_DIR:$CONTAINER_SHARED_DIR:Z" \
      "$IMAGE_NAME" \
      sleep infinity
    
    # Install dependencies
    log "INFO" "Installing dependencies..."
    podman exec "$CONTAINER_NAME" npm install -g playwright
  fi
}

# Execute a Playwright command
execute_browser_navigate() {
  ensure_container
  
  local url=$1
  local output_file="$SHARED_DIR/browser_navigate_result.json"
  
  log "INFO" "Navigating to URL: $url"
  
  # Create script file
  cat > "$SHARED_DIR/browser_navigate.js" << EOF
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    await page.goto('$url');
    const title = await page.title();
    const result = { success: true, data: { title, url: page.url() } };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/browser_navigate_result.json', JSON.stringify(result));
  } catch (error) {
    const result = { success: false, error: error.message };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/browser_navigate_result.json', JSON.stringify(result));
  } finally {
    await browser.close();
  }
})();
EOF

  # Execute script
  podman exec "$CONTAINER_NAME" node "$CONTAINER_SHARED_DIR/browser_navigate.js"
  
  # Return result
  if [ -f "$output_file" ]; then
    cat "$output_file"
  else
    echo '{"success":false,"error":"Failed to create result file"}'
  fi
}

execute_browser_take_screenshot() {
  ensure_container
  
  local url=$1
  local filename=$2
  local output_file="$SHARED_DIR/screenshot_result.json"
  
  log "INFO" "Taking screenshot of URL: $url"
  
  # Create script file
  cat > "$SHARED_DIR/take_screenshot.js" << EOF
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    await page.goto('$url');
    await page.screenshot({ path: '$CONTAINER_SHARED_DIR/screenshot.png' });
    const result = { success: true, data: { filename: 'screenshot.png' } };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/screenshot_result.json', JSON.stringify(result));
  } catch (error) {
    const result = { success: false, error: error.message };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/screenshot_result.json', JSON.stringify(result));
  } finally {
    await browser.close();
  }
})();
EOF

  # Execute script
  podman exec "$CONTAINER_NAME" node "$CONTAINER_SHARED_DIR/take_screenshot.js"
  
  # Copy screenshot to desired location if provided
  if [ -f "$SHARED_DIR/screenshot.png" ] && [ -n "$filename" ]; then
    cp "$SHARED_DIR/screenshot.png" "$filename"
  fi
  
  # Return result
  if [ -f "$output_file" ]; then
    cat "$output_file"
  else
    echo '{"success":false,"error":"Failed to create result file"}'
  fi
}

execute_browser_snapshot() {
  ensure_container
  
  local url=$1
  local output_file="$SHARED_DIR/snapshot_result.json"
  
  log "INFO" "Taking accessibility snapshot of URL: $url"
  
  # Create script file
  cat > "$SHARED_DIR/snapshot.js" << EOF
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    await page.goto('$url');
    const snapshot = await page.accessibility.snapshot();
    const result = { success: true, data: { snapshot } };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/snapshot_result.json', JSON.stringify(result));
  } catch (error) {
    const result = { success: false, error: error.message };
    fs.writeFileSync('$CONTAINER_SHARED_DIR/snapshot_result.json', JSON.stringify(result));
  } finally {
    await browser.close();
  }
})();
EOF

  # Execute script
  podman exec "$CONTAINER_NAME" node "$CONTAINER_SHARED_DIR/snapshot.js"
  
  # Return result
  if [ -f "$output_file" ]; then
    cat "$output_file"
  else
    echo '{"success":false,"error":"Failed to create result file"}'
  fi
}

# Parse command line arguments
if [ $# -eq 0 ]; then
  echo "Playwright Docker MCP"
  echo "Usage: $0 <command> [arguments]"
  echo ""
  echo "Commands:"
  echo "  browser_navigate <url>                       Navigate to a URL"
  echo "  browser_take_screenshot <url> [filename]     Take a screenshot"
  echo "  browser_snapshot <url>                       Get accessibility snapshot"
  echo "  status                                       Show container status"
  echo ""
  exit 0
fi

command=$1
shift

# Execute command
case $command in
  "browser_navigate")
    if [ $# -eq 0 ]; then
      log "ERROR" "URL is required"
      echo '{"success":false,"error":"URL is required"}'
      exit 1
    fi
    execute_browser_navigate "$1"
    ;;
  "browser_take_screenshot")
    if [ $# -eq 0 ]; then
      log "ERROR" "URL is required"
      echo '{"success":false,"error":"URL is required"}'
      exit 1
    fi
    execute_browser_take_screenshot "$1" "$2"
    ;;
  "browser_snapshot")
    if [ $# -eq 0 ]; then
      log "ERROR" "URL is required"
      echo '{"success":false,"error":"URL is required"}'
      exit 1
    fi
    execute_browser_snapshot "$1"
    ;;
  "status")
    # Check if container exists
    if podman container exists "$CONTAINER_NAME" 2>/dev/null; then
      # Check if container is running
      if podman container inspect "$CONTAINER_NAME" --format '{{.State.Running}}' 2>/dev/null | grep -q "true"; then
        log "INFO" "Container is running"
        
        # Test if Playwright works
        podman exec "$CONTAINER_NAME" bash -c "node -e \"const { chromium } = require('playwright'); console.log('Playwright is working');\""
        if [ $? -eq 0 ]; then
          log "INFO" "Playwright is functioning correctly"
          echo '{"success":true,"status":"running","playwright":true}'
        else
          log "ERROR" "Playwright is not functioning correctly"
          echo '{"success":true,"status":"running","playwright":false}'
        fi
      else
        log "WARN" "Container exists but is not running"
        echo '{"success":true,"status":"stopped"}'
      fi
    else
      log "WARN" "Container does not exist"
      echo '{"success":false,"status":"not_created"}'
    fi
    ;;
  *)
    log "ERROR" "Unknown command: $command"
    echo '{"success":false,"error":"Unknown command"}'
    exit 1
    ;;
esac

exit 0
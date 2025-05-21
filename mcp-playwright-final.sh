#!/bin/bash

# mcp-playwright-final.sh
# Final solution for MCP Playwright on Linux
# This script uses the official Microsoft Playwright Docker image with npx commands

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$SCRIPT_DIR/mcp-playwright-final.log"
TEMP_DIR="$SCRIPT_DIR/playwright-final-temp"
OUTPUT_DIR="$SCRIPT_DIR/playwright-final-output"
CONTAINER_TAG="v1.52.0-jammy"

# Create directories if they don't exist
mkdir -p "$TEMP_DIR" "$OUTPUT_DIR"

# Function to log messages
log() {
  local level=$1
  local message=$2
  
  echo "[$level] $message"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" >> "$LOG_FILE"
}

# Ensure container image is available
ensure_container() {
  if ! podman image exists "mcr.microsoft.com/playwright:$CONTAINER_TAG"; then
    log "INFO" "Pulling Playwright container image..."
    podman pull "mcr.microsoft.com/playwright:$CONTAINER_TAG"
    if [ $? -ne 0 ]; then
      log "ERROR" "Failed to pull container image"
      echo "{\"success\":false,\"error\":\"Failed to pull container image\"}"
      exit 1
    fi
  fi
}

# Execute a command in the Playwright container
run_in_container() {
  local command="$1"
  
  podman run --rm \
    --security-opt label=disable \
    -v "$TEMP_DIR:/input" \
    -v "$OUTPUT_DIR:/output" \
    --network host \
    "mcr.microsoft.com/playwright:$CONTAINER_TAG" \
    bash -c "$command"
  
  return $?
}

# Navigate to URL
browser_navigate() {
  local url="$1"
  log "INFO" "Navigating to URL: $url"
  
  # Create a script file to navigate to the URL
  cat > "$TEMP_DIR/navigate.js" << EOF
// navigate.js
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to URL: $url');
    await page.goto('$url');
    const title = await page.title();
    console.log(JSON.stringify({ 
      success: true, 
      data: { 
        title: title, 
        url: page.url() 
      } 
    }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ 
      success: false, 
      error: error.message 
    }));
  } finally {
    await browser.close();
  }
})();
EOF

  # Run the script in the container
  ensure_container
  run_in_container "cd /input && npm init -y && npm install playwright && node navigate.js"
}

# Take screenshot
browser_take_screenshot() {
  local url="$1"
  local filename="$2"
  
  # If no filename provided, use a default one
  if [ -z "$filename" ]; then
    filename="$OUTPUT_DIR/screenshot-$(date +%s).png"
  fi
  
  log "INFO" "Taking screenshot of URL: $url"
  
  # Create a script file to take a screenshot
  cat > "$TEMP_DIR/screenshot.js" << EOF
// screenshot.js
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to URL: $url');
    await page.goto('$url');
    console.log('Taking screenshot...');
    await page.screenshot({ path: '/output/screenshot.png', fullPage: true });
    console.log(JSON.stringify({ 
      success: true, 
      data: { 
        filename: '/output/screenshot.png'
      } 
    }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ 
      success: false, 
      error: error.message 
    }));
  } finally {
    await browser.close();
  }
})();
EOF

  # Run the script in the container
  ensure_container
  run_in_container "cd /input && npm init -y && npm install playwright && node screenshot.js"
  
  # Copy the screenshot to the desired location
  if [ -f "$OUTPUT_DIR/screenshot.png" ]; then
    cp "$OUTPUT_DIR/screenshot.png" "$filename"
    log "INFO" "Screenshot saved to $filename"
    echo "{\"success\":true,\"data\":{\"filename\":\"$filename\"}}"
  else
    log "ERROR" "Failed to create screenshot"
    echo "{\"success\":false,\"error\":\"Failed to create screenshot\"}"
  fi
}

# Get accessibility snapshot
browser_snapshot() {
  local url="$1"
  log "INFO" "Getting accessibility snapshot of URL: $url"
  
  # Create a script file to get an accessibility snapshot
  cat > "$TEMP_DIR/snapshot.js" << EOF
// snapshot.js
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to URL: $url');
    await page.goto('$url');
    console.log('Getting accessibility snapshot...');
    const snapshot = await page.accessibility.snapshot();
    console.log(JSON.stringify({ 
      success: true, 
      data: { 
        snapshot: snapshot
      } 
    }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ 
      success: false, 
      error: error.message 
    }));
  } finally {
    await browser.close();
  }
})();
EOF

  # Run the script in the container
  ensure_container
  run_in_container "cd /input && npm init -y && npm install playwright && node snapshot.js"
}

# Check status
check_status() {
  log "INFO" "Checking Playwright container status"
  
  # Create a simple script to test Playwright
  cat > "$TEMP_DIR/status.js" << EOF
// status.js
const { chromium } = require('playwright');

(async () => {
  try {
    const browser = await chromium.launch();
    await browser.close();
    console.log(JSON.stringify({ success: true, status: "Playwright is working correctly" }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ success: false, error: error.message }));
  }
})();
EOF

  # Run the script in the container
  ensure_container
  run_in_container "cd /input && npm init -y && npm install playwright && node status.js"
}

# Main function to handle commands
main() {
  case "$1" in
    browser_navigate)
      if [ -z "$2" ]; then
        log "ERROR" "URL required"
        echo "{\"success\":false,\"error\":\"URL required\"}"
        exit 1
      fi
      browser_navigate "$2"
      ;;
    browser_take_screenshot)
      if [ -z "$2" ]; then
        log "ERROR" "URL required"
        echo "{\"success\":false,\"error\":\"URL required\"}"
        exit 1
      fi
      browser_take_screenshot "$2" "$3"
      ;;
    browser_snapshot)
      if [ -z "$2" ]; then
        log "ERROR" "URL required"
        echo "{\"success\":false,\"error\":\"URL required\"}"
        exit 1
      fi
      browser_snapshot "$2"
      ;;
    status)
      check_status
      ;;
    help|--help|-h)
      echo "MCP Playwright Final Integration"
      echo "-------------------------------"
      echo "Usage: $0 <command> [arguments]"
      echo ""
      echo "Commands:"
      echo "  browser_navigate <url>                 Navigate to a URL"
      echo "  browser_take_screenshot <url> [file]   Take a screenshot"
      echo "  browser_snapshot <url>                 Get accessibility snapshot"
      echo "  status                                 Check Playwright container status"
      echo "  help                                   Show this help message"
      ;;
    *)
      echo "Unknown command: $1"
      echo "Run '$0 help' for usage information"
      exit 1
      ;;
  esac
}

# Execute main function with all arguments
main "$@"
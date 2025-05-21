#!/bin/bash

# mcp-playwright-direct.sh
# Direct MCP Playwright integration using Microsoft's official container
# This script provides a direct interface to Playwright functionality in a container

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$SCRIPT_DIR/mcp-playwright-direct.log"
TEMP_DIR="$SCRIPT_DIR/playwright-temp"
SCREENSHOT_DIR="$SCRIPT_DIR/playwright-screenshots"

# Create directories if they don't exist
mkdir -p "$TEMP_DIR" "$SCREENSHOT_DIR"

# Function to log messages
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"
  echo "$1"
}

# Function to execute a playwright command in the container
run_playwright() {
  local script="$1"
  log "Running Playwright script: $script"
  
  # Write script to temp file
  local script_file="$TEMP_DIR/script-$(date +%s).js"
  echo "$script" > "$script_file"
  
  # Run in container
  podman run --rm \
    --security-opt label=disable \
    -v "$TEMP_DIR:/app" \
    -v "$SCREENSHOT_DIR:/screenshots" \
    --network host \
    mcr.microsoft.com/playwright:v1.40.0-jammy \
    node "/app/$(basename "$script_file")"
    
  return $?
}

# Navigate to URL
browser_navigate() {
  local url="$1"
  log "Navigating to URL: $url"
  
  local script="
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to $url');
    await page.goto('$url');
    const title = await page.title();
    console.log(JSON.stringify({ success: true, data: { title, url: page.url() } }));
  } catch (error) {
    console.error(error);
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();
"
  run_playwright "$script"
}

# Take screenshot
browser_take_screenshot() {
  local url="$1"
  local filename="$2"
  
  if [ -z "$filename" ]; then
    filename="screenshot-$(date +%s).png"
  fi
  
  log "Taking screenshot of URL: $url, saving to: $filename"
  
  local script="
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to $url');
    await page.goto('$url');
    console.log('Taking screenshot');
    await page.screenshot({ path: '/screenshots/temp-screenshot.png', fullPage: true });
    console.log(JSON.stringify({ success: true, data: { filename: '$filename' } }));
  } catch (error) {
    console.error(error);
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();
"
  run_playwright "$script"
  
  # Copy screenshot to desired location
  if [ -f "$SCREENSHOT_DIR/temp-screenshot.png" ]; then
    cp "$SCREENSHOT_DIR/temp-screenshot.png" "$filename"
    log "Screenshot saved to: $filename"
  else
    log "Error: Screenshot not created"
  fi
}

# Get accessibility snapshot
browser_snapshot() {
  local url="$1"
  log "Getting accessibility snapshot of URL: $url"
  
  local script="
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to $url');
    await page.goto('$url');
    console.log('Getting accessibility snapshot');
    const snapshot = await page.accessibility.snapshot();
    console.log(JSON.stringify({ success: true, data: { snapshot } }));
  } catch (error) {
    console.error(error);
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();
"
  run_playwright "$script"
}

# Check if playwright container is available
check_playwright() {
  if ! podman image exists mcr.microsoft.com/playwright:v1.40.0-jammy; then
    log "Pulling Playwright image..."
    podman pull mcr.microsoft.com/playwright:v1.40.0-jammy
    if [ $? -ne 0 ]; then
      log "Failed to pull Playwright image"
      return 1
    fi
  fi
  
  # Test running a simple script
  local script="
const { chromium } = require('playwright');
console.log(JSON.stringify({ success: true, message: 'Playwright container working properly' }));
"
  run_playwright "$script"
  return $?
}

# Main function
main() {
  case "$1" in
    browser_navigate)
      if [ -z "$2" ]; then
        echo "Error: URL required"
        echo "Usage: $0 browser_navigate <url>"
        exit 1
      fi
      browser_navigate "$2"
      ;;
    browser_take_screenshot)
      if [ -z "$2" ]; then
        echo "Error: URL required"
        echo "Usage: $0 browser_take_screenshot <url> [filename]"
        exit 1
      fi
      browser_take_screenshot "$2" "$3"
      ;;
    browser_snapshot)
      if [ -z "$2" ]; then
        echo "Error: URL required"
        echo "Usage: $0 browser_snapshot <url>"
        exit 1
      fi
      browser_snapshot "$2"
      ;;
    status)
      echo "Checking Playwright container..."
      check_playwright
      ;;
    help|--help|-h)
      echo "MCP Playwright Direct Integration"
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

# Run main function with all arguments
main "$@"
#!/bin/bash

# mcp-playwright.sh
# Integration script for MCP Playwright commands
# This script serves as the primary entry point for MCP Playwright operations

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BRIDGE_SCRIPT="$SCRIPT_DIR/mcp-playwright-bridge.js"
LOG_FILE="$SCRIPT_DIR/mcp-playwright.log"

# Log message with timestamp
log() {
  local level=$1
  local message=$2
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" >> "$LOG_FILE"
}

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
  log "ERROR" "Node.js is not installed. Please install Node.js first."
  echo "Node.js is not installed. Please install Node.js first."
  exit 1
fi

# Check if bridge script exists
if [ ! -f "$BRIDGE_SCRIPT" ]; then
  log "ERROR" "Bridge script not found: $BRIDGE_SCRIPT"
  echo "Bridge script not found: $BRIDGE_SCRIPT"
  exit 1
fi

# Parse command line arguments
operation=$1
shift

# Check if operation is provided
if [ -z "$operation" ]; then
  echo "MCP Playwright Integration"
  echo "Usage: $0 <operation> [parameters]"
  echo ""
  echo "Operations:"
  echo "  browser_navigate         Navigate to a URL"
  echo "  browser_take_screenshot  Take a screenshot"
  echo "  browser_snapshot         Get accessibility snapshot"
  echo "  browser_click            Click on an element"
  echo "  browser_type             Type text into an element"
  echo ""
  exit 0
fi

# Handle parameters
params="{}"
if [ $# -gt 0 ]; then
  params="$*"
fi

# Log the operation
log "INFO" "Executing operation: $operation with params: $params"

# Execute the bridge script
node "$BRIDGE_SCRIPT" "$operation" "$params"

# Return the exit code from the bridge script
exit $?
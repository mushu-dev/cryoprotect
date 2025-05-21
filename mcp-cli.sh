#!/bin/bash

# mcp-cli.sh - CLI wrapper for MCP tools with library compatibility
# Usage: ./mcp-cli.sh [command]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LINK_DIR="$SCRIPT_DIR/lib_links"

# Ensure lib_links directory exists with necessary symlinks
if [ ! -d "$LINK_DIR" ] || [ ! -f "$LINK_DIR/libicudata.so.66" ]; then
    echo "Setting up library symlinks..."
    "$SCRIPT_DIR/setup_mcp_playwright.sh"
fi

# Set the library path
export LD_LIBRARY_PATH="$LINK_DIR:$LD_LIBRARY_PATH"

# If no arguments provided, show help
if [ $# -eq 0 ]; then
    echo "MCP CLI Tool Wrapper"
    echo "--------------------"
    echo "Usage: ./mcp-cli.sh [command]"
    echo ""
    echo "Examples:"
    echo "  ./mcp-cli.sh npx playwright install # Install Playwright browsers"
    echo "  ./mcp-cli.sh node test_playwright.js # Run a Playwright test script"
    echo "  ./mcp-cli.sh node your-app.js # Run Node.js app with MCP compatibility"
    echo ""
    echo "Environment:"
    echo "  LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
    exit 0
fi

# Handle special case for 'mcp' command
if [ "$1" = "mcp" ]; then
    echo "The 'mcp' command is meant to be run from Claude Code CLI interface."
    echo "This script provides the library compatibility layer needed for Playwright."
    echo "To test Playwright functionality with our compatibility layer, run:"
    echo "  ./mcp-cli.sh node test_mcp_playwright.js"
    exit 0
fi

# Execute the command with the library path
echo "Running with MCP compatibility: $@"
exec "$@"
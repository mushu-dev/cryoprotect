#!/bin/bash

# setup_mcp_playwright.sh
# This script sets up the environment for MCP Playwright on Fedora

# Create necessary symbolic links for library compatibility
LINK_DIR="/home/mushu/Projects/cryoprotect/lib_links"
mkdir -p "$LINK_DIR"

# Create symbolic links for required libraries
ln -sf /usr/lib64/libicudata.so.76 "$LINK_DIR/libicudata.so.66"
ln -sf /usr/lib64/libicui18n.so.76 "$LINK_DIR/libicui18n.so.66"
ln -sf /usr/lib64/libicuuc.so.76 "$LINK_DIR/libicuuc.so.66"
ln -sf /usr/lib64/libjpeg.so.62 "$LINK_DIR/libjpeg.so.8"
ln -sf /usr/lib64/libwebp.so.7 "$LINK_DIR/libwebp.so.6"
ln -sf /usr/lib64/libffi.so.8 "$LINK_DIR/libffi.so.7"

# Export the library path
export LD_LIBRARY_PATH="$LINK_DIR:$LD_LIBRARY_PATH"

echo "MCP Playwright environment configured!"
echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"
echo ""
echo "To use in your current shell session, run:"
echo "source setup_mcp_playwright.sh"
echo ""
echo "To test, run:"
echo "node test_mcp_playwright.js"

# Execute the command passed as arguments if any
if [ $# -gt 0 ]; then
    echo "Running command: $@"
    exec "$@"
fi
# Playwright Linux Compatibility Guide

This guide addresses compatibility issues with Playwright on modern Linux distributions, specifically Fedora 42 and other newer distros that have library versions that differ from what Playwright expects.

## Problem

Playwright expects specific library versions that aren't always available on modern Linux distributions. On Fedora 42, we encountered the following missing libraries:

- `libicudata.so.66` (Fedora has version 76)
- `libicui18n.so.66` (Fedora has version 76)
- `libicuuc.so.66` (Fedora has version 76)
- `libjpeg.so.8` (Fedora has version 62)
- `libwebp.so.6` (Fedora has version 7)
- `libffi.so.7` (Fedora has version 8)

## Solution

We've created symbolic links from the newer library versions to the older ones that Playwright expects. This allows Playwright to run without requiring system-wide changes or installing older library versions.

## Quick Start

1. Run the setup script once to create the necessary symbolic links:
   ```bash
   ./setup_mcp_playwright.sh
   ```

2. Use the MCP CLI wrapper to run commands with Playwright compatibility:
   ```bash
   ./mcp-cli.sh mcp                     # Run MCP with Playwright compatibility
   ./mcp-cli.sh node your-script.js     # Run a Node.js script with compatibility
   ```

## Files

- `setup_mcp_playwright.sh`: Sets up symbolic links and configures the environment
- `mcp-cli.sh`: CLI wrapper for running commands with the correct environment
- `lib_links/`: Directory containing symbolic links to required libraries
- `test_playwright.js`: Simple test script to verify Playwright functionality
- `test_mcp_playwright.js`: Test script to verify MCP integration

## Manual Environment Setup

If you prefer to set up the environment manually, you can:

1. Create the symbolic links:
   ```bash
   mkdir -p ~/Projects/cryoprotect/lib_links
   cd ~/Projects/cryoprotect/lib_links
   ln -s /usr/lib64/libicudata.so.76 libicudata.so.66
   ln -s /usr/lib64/libicui18n.so.76 libicui18n.so.66
   ln -s /usr/lib64/libicuuc.so.76 libicuuc.so.66
   ln -s /usr/lib64/libjpeg.so.62 libjpeg.so.8
   ln -s /usr/lib64/libwebp.so.7 libwebp.so.6
   ln -s /usr/lib64/libffi.so.8 libffi.so.7
   ```

2. Set the environment variable:
   ```bash
   export LD_LIBRARY_PATH=~/Projects/cryoprotect/lib_links:$LD_LIBRARY_PATH
   ```

3. Run your commands:
   ```bash
   mcp
   node your-script.js
   ```

## Troubleshooting

If you encounter issues:

1. Verify the symbolic links are correctly created:
   ```bash
   ls -la lib_links/
   ```

2. Check the environment variable is set:
   ```bash
   echo $LD_LIBRARY_PATH
   ```

3. Run the test script to verify Playwright functionality:
   ```bash
   ./mcp-cli.sh node test_playwright.js
   ```

## Alternative Approaches

If this solution doesn't work for your needs, consider:

1. Using a container with Playwright pre-installed
2. Installing older versions of the required libraries (not recommended)
3. Using a different browser automation tool

## Notes

- This solution doesn't modify system libraries
- The symbolic links are contained within the project directory
- Performance impact should be minimal
- This is a workaround - future Playwright versions may support newer library versions directly